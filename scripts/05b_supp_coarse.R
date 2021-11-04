
### Estimate LTLA prevalence with partial/coarsened REACT data ###

library(dplyr)
library(prevdebiasr)
library(parallel)
library(foreach)
source("scripts/SIR_utils.R")

trans_mats <- readRDS("transmats/poisson_SIR_epi_gamma_1.RDS")
vax_df <- readr::read_csv("data/vaccination.csv")
pcr_infectious_df <- readr::read_csv("data/moment_match_infectious.csv")

############################
### Remove REACT round 8 ###
############################

round8_start_midweek <- as.Date("2021-01-03")
round8_end_midweek <- as.Date("2021-01-24")
region_df <- readr::read_csv("data/region.csv") %>%
  left_join(pcr_infectious_df, by = c("phe_region", "mid_week")) %>%
  mutate(
    Nr = if_else(mid_week <= round8_end_midweek &
                   mid_week >= round8_start_midweek, 0, Nr),
    nr = if_else(mid_week <= round8_end_midweek &
                   mid_week >= round8_start_midweek, 0, nr)
  )

ltla_df <- readr::read_csv("data/ltla.csv") %>%
  left_join(vax_df, by = c("ltla", "mid_week")) %>%
  left_join(select(region_df, -c(Nt, nt, M)),
            by = c("phe_region", "mid_week")
  ) %>%
  select(ltla, phe_region, mid_week, Nt, nt, M, V, alpha, beta)

alpha_testing <- 3e-4
control_debias <- prevdebiasr::get_control_parameters(
  alpha_testing = alpha_testing
)

mid_week_unique <- sort(unique(ltla_df$mid_week))

n_cores <- 4
clust <- makeCluster(n_cores)
doParallel::registerDoParallel(clust)

out_dir <- paste0(
  "output/AR", control_debias$delta_AR_rho,
  "sd", control_debias$delta_AR_sd,
  "Rsd", control_SIR$R_AR_sd
)
control <- c(control_debias, control_SIR)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
#saveRDS(control, file.path(out_dir, "control.RDS"), version = 2)

clusterExport(clust, c("control_debias", "control_SIR"))

### PCR positivity, perfect testing ###
imperfect <- FALSE
type <- "PCR_positive"

# Calculate regional bias parameters
delta_df <- region_df %>%
  group_by(phe_region) %>%
  group_modify(~ cbind(
    mid_week = .x$mid_week,
    specify_delta_prior(.x, control_debias, imperfect)
  ))

# Estimate local prevalence for each LTLA
ltla_list <- ltla_df %>%
  left_join(delta_df, by = c("phe_region", "mid_week")) %>%
  group_by(ltla) %>%
  group_split()
ltla_names <- sapply(ltla_list, function(x) x$ltla[1])
ltla_prevalence <- parLapply(clust, ltla_list, local_prevalence,
                             control_debias, imperfect, type)
names(ltla_prevalence) <- ltla_names
# Save output

delta_out_file <- file.path(out_dir, "delta_no_round8.csv")
readr::write_csv(delta_df, delta_out_file)

ltla_out_file <- file.path(out_dir, "ltla_prevalence_no_round8.RDS")
saveRDS(ltla_prevalence, ltla_out_file, version = 2)

######################
### National delta ###
######################

region_df <- readr::read_csv("data/region.csv") %>%
  select(-phe_region) %>%
  group_by(mid_week) %>%
  summarise_all(sum)

ltla_df <- readr::read_csv("data/ltla.csv") %>%
  left_join(vax_df, by = c("ltla", "mid_week")) %>%
  left_join(select(region_df, -c(Nt, nt, M)),
            by = c("mid_week")
  ) %>%
  select(ltla, mid_week, Nt, nt, M, V)

delta_df <- cbind(
  region_df,
  specify_delta_prior(region_df, control_debias, imperfect)
) %>%
  select(mid_week, delta_prior_mean, delta_prior_sd)

# Estimate local prevalence for each LTLA
ltla_list <- ltla_df %>%
  left_join(delta_df, by = c("mid_week")) %>%
  group_by(ltla) %>%
  group_split()
ltla_names <- sapply(ltla_list, function(x) x$ltla[1])
ltla_prevalence <- parLapply(
  clust, ltla_list, local_prevalence,
  control_debias, imperfect, type
)
names(ltla_prevalence) <- ltla_names
# Save output

delta_out_file <- file.path(out_dir, "delta_national.csv")
readr::write_csv(delta_df, delta_out_file)

ltla_out_file <- file.path(out_dir, "ltla_prevalence_national.RDS")
saveRDS(ltla_prevalence, ltla_out_file, version = 2)

stopCluster(clust)
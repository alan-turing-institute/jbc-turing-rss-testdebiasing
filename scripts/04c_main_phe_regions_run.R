#!/apps/eb/skylake/software/R/3.6.2-foss-2019b/bin/Rscript
#$ -cwd
#$ -N main_run
#$ -q short.qe
#$ -pe shmem 12
#$ -P holmes.prjc
#$ -o logs/main.log
#$ -e logs/main.log

### Estimate LTLA prevalence ###

library(dplyr)
library(prevdebiasr)
library(parallel)
library(foreach)
source("scripts/SIR_utils.R")

trans_mats <- readRDS("transmats/poisson_SIR_epi_gamma_1.RDS")
pcr_infectious_df <- readr::read_csv("data/moment_match_infectious.csv")
ltla_df <- readr::read_csv("data/ltla.csv")

# %>%
#   left_join(vax_df, by = c("ltla", "mid_week")) %>%
#   left_join(select(region_df, -c(Nt, nt, M)),
#     by = c("phe_region", "mid_week")
#   ) %>%
#   select(ltla, phe_region, mid_week, Nt, nt, M, V, alpha, beta)
# 
vax_df <- readr::read_csv("data/vaccination.csv")
vax_df_to_aggreg <- as.data.frame(vax_df)
vax_df_to_aggreg$phe_region <- unlist(ltla_df[match(vax_df_to_aggreg$ltla, vax_df_to_aggreg$ltla), "phe_region"])
vax_df_to_aggreg$phe_region_mid_week <- paste0(vax_df_to_aggreg$phe_region, "_", vax_df_to_aggreg$mid_week)
phe_vax_df <- vax_df_to_aggreg[, c("phe_region", "mid_week", "phe_region_mid_week")]
phe_vax_df <- unique(phe_vax_df)
phe_vax_df$V <- sapply(phe_vax_df$phe_region_mid_week, function(x) sum(vax_df_to_aggreg[vax_df_to_aggreg$phe_region_mid_week == x, "V"]))


region_df <- readr::read_csv("data/region.csv") %>%
  left_join(pcr_infectious_df, by = c("phe_region", "mid_week")) %>%
  left_join(phe_vax_df, by = c("phe_region", "mid_week"))
region_df$ltla <- region_df$phe_region


str(region_df)
alpha_testing <- 3e-4
control_debias <- prevdebiasr::get_control_parameters(
    alpha_testing = alpha_testing
)

mid_week_unique <- sort(unique(region_df$mid_week))

n_cores <- 12
clust <- makeCluster(n_cores)
doParallel::registerDoParallel(clust)

out_dir <- paste0(
    "output/AR", control_debias$delta_AR_rho,
    "sd", control_debias$delta_AR_sd,
    "Rsd", control_SIR$R_AR_sd
)
control <- c(control_debias, control_SIR)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(control, file.path(out_dir, "control.RDS"), version = 2)

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

# Estimate local prevalence for each PHE region
region_list <- region_df %>%
left_join(delta_df, by = c("phe_region", "mid_week")) %>%
group_by(phe_region) %>%
group_split()
region_names <- sapply(region_list, function(x) x$phe_region[1])
region_prevalence <- parLapply(clust, region_list, local_prevalence,
                            control_debias, imperfect, type)
names(region_prevalence) <- region_names
# Save output

delta_out_file <- file.path(out_dir, "delta_pcr_perfect.csv")
# readr::write_csv(delta_df, delta_out_file)

region_out_file <- file.path(out_dir, "phe_region_prevalence_pcr_perfect.RDS")
saveRDS(region_prevalence, region_out_file, version = 2)


### Imperfect testing, PCR positivity and infectiousness ###
imperfect <- FALSE
type <- "PCR_positive" # for (type in c("Infectious", "PCR_positive")) {

    # Calculate regional bias parameters
    delta_df <- region_df %>%
        group_by(phe_region) %>%
        group_modify(~ cbind(
            mid_week = .x$mid_week,
            specify_delta_prior(.x, control_debias, imperfect)
        ))

    # # Estimate local prevalence for each LTLA
    # ltla_list <- ltla_df %>%
    #     left_join(delta_df, by = c("phe_region", "mid_week")) %>%
    #     # left_join(pcr_infectious_df, by = c("phe_region", "mid_week")) %>%
    #     group_by(ltla) %>%
    #     group_split()
    # ltla_names <- sapply(ltla_list, function(x) x$ltla[1])

    region_prevalence <- parLapply(
        clust, region_list, local_prevalence,
        control_debias, imperfect, type
    )
    names(region_prevalence) <- region_names

    # Save output
    dir.create(file.path(out_dir, type), showWarnings = FALSE, recursive = TRUE)
    # delta_out_file <- file.path(out_dir, type, "delta.csv")
    # readr::write_csv(delta_df, delta_out_file)

    region_out_file <- file.path(out_dir, type, "phe_region_prevalence.RDS")
    saveRDS(region_prevalence, region_out_file, version = 2)

    # Fit SIR to latest available date
    foreach(region_name = region_names, .packages = "dplyr") %dopar% {
        d_region <- region_df %>%
            filter(phe_region == region_name)

        I_log_lik <- region_prevalence[[region_name]]$log_post
        SIR_model_out_region <- sample_I_R_SIR_from_pre_calc_lik(
            d_region, I_log_lik,
            trans_mats, control_debias,
            control_SIR
        )

        out_file <- file.path(out_dir, type, "SIR", paste0(region_name, ".RDS"))
        dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
        saveRDS(SIR_model_out_region, out_file, version = 2)
    }
}

stopCluster(clust)

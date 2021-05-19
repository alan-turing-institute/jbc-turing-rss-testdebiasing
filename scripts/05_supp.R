#!/apps/eb/skylake/software/R/3.6.2-foss-2019b/bin/Rscript
#$ -cwd
#$ -N supp_run
#$ -q short.qe
#$ -P holmes.prjc
#$ -o logs/supp.log
#$ -e logs/supp.log

### Estimate LTLA prevalence ###

library(dplyr)
library(prevdeviasr)
library(parallel)
library(foreach)
source("scripts/SIR_utils.R")

trans_mats <- readRDS("transmats/poisson_SIR_epi_gamma_1.RDS")
vax_df <- readr::read_csv("data/vaccination.csv")
pcr_infectious_df <- readr::read_csv("data/moment_match_infectious.csv")
region_df <- readr::read_csv("data/region.csv") %>%
  left_join(pcr_infectious_df, by = c("phe_region", "mid_week"))
ltla_df <- readr::read_csv("data/ltla.csv") %>%
  left_join(vax_df, by = c("ltla", "mid_week")) %>%
  left_join(select(region_df, -c(Nt, nt, M)),
            by = c("phe_region", "mid_week")) %>%
  select(ltla, phe_region, mid_week, Nt, nt, M, V, alpha, beta)

control_debias <- prevdebiasr::get_control_parameters()

mid_week_unique <- sort(unique(ltla_df$mid_week))

param_df <- expand.grid(alpha_testing = c(0.001, 0.003),
                        beta_testing = c(0.05, 0.1, 0.3))

n.cores <- 12
clust <- makeCluster(n.cores)
doParallel::registerDoParallel(clust)

control_debias$delta_AR_rho <- 0.99
control_debias$delta_AR_sd <- 1
control_SIR$R_AR_sd <- 0.2

for (i in 1:nrow(param_df)) {
  
  control_debias$alpha_testing <- param_df$alpha_testing[i]
  control_debias$beta_testing <- param_df$beta_testing[i]

  out_dir <- paste0(
    "output/alpha", control_debias$alpha_testing,
    "beta", control_debias$beta_testing
  )
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  clusterExport(clust, c("control_debias", "control_SIR"))
  
  ### Imperfect testing, PCR positivity and infectiousness ###
  imperfect <- TRUE
  
  # Calculate regional bias parameters
  delta_df <- region_df %>%
    group_by(phe_region) %>%
    group_modify(~ cbind(mid_week = .x$mid_week, 
                         specify_delta_prior(.x, control_debias, imperfect)))

  for (type in c("PCR_positive", "Infectious")) {
    # Estimate local prevalence for each LTLA
    ltla_list <- ltla_df %>%
      left_join(delta_df, by = c("phe_region", "mid_week")) %>%
      left_join(pcr_infectious_df,  by = c("phe_region", "mid_week")) %>%
      group_by(ltla) %>%
      group_split()
    ltla_names <- sapply(ltla_list, function(x) x$ltla[1])
    
    ltla_prevalence <- parLapply(clust, ltla_list, local_prevalence, 
                                 control_debias, imperfect, type)
    names(ltla_prevalence) <- ltla_names

    # Save output
    dir.create(file.path(out_dir, type), showWarnings = FALSE, recursive = TRUE)
    delta_out_file <- file.path(out_dir, type, "delta.csv")
    readr::write_csv(delta_df, delta_out_file)

    ltla_out_file <- file.path(out_dir, type, "ltla_prevalence.RDS")
    saveRDS(ltla_prevalence, ltla_out_file)

    # Fit SIR to latest available date
    foreach(ltla_name = ltla_names, .packages = "dplyr") %dopar% {
      d_ltla <- ltla_df %>%
        filter(ltla == ltla_name)

      I_log_lik <- ltla_prevalence[[ltla_name]]$log_post
      SIR_model_out_ltla <- sample_I_R_SIR_from_pre_calc_lik(
        d_ltla, I_log_lik,
        trans_mats, control_debias,
        control_SIR
      )

      out_file <- file.path(out_dir, type, "SIR", paste0(ltla_name, ".RDS"))
      dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
      saveRDS(SIR_model_out_ltla, out_file)
    }
  }
}

stopCluster(clust)
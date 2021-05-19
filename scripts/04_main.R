#!/apps/eb/skylake/software/R/3.6.2-foss-2019b/bin/Rscript
#$ -cwd
#$ -N main_run
#$ -q short.qe
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
vax_df <- readr::read_csv("data/vaccination.csv")
ltla_df <- readr::read_csv("data/ltla.csv") %>%
  left_join(vax_df, by = c("ltla", "mid_week"))
pcr_infectious_df <- readr::read_csv("data/moment_match_infectious.csv")
region_df <- readr::read_csv("data/region.csv") %>%
  left_join(pcr_infectious_df, by = c("phe_region", "mid_week"))

control_debias <- prevdebiasr::get_control_parameters()

mid_week_unique <- sort(unique(ltla_df$mid_week))
cut_dates <- rev(mid_week_unique)[2:12]

param_df <- expand.grid(
  delta_AR_rho = c(0.975, 0.99, 0.999),
  delta_AR_sd = c(1, 2),
  R_AR_sd = 0.2
)

n.cores <- 2
clust <- makeCluster(n.cores)
doParallel::registerDoParallel(clust)

for (i in 1:nrow(param_df)) {
  control_debias$delta_AR_rho <- param_df$delta_AR_rho[i]
  control_debias$delta_AR_sd <- param_df$delta_AR_sd[i]
  control_SIR$R_AR_sd <- param_df$R_AR_sd[i]
  out_dir <- paste0(
    "output/AR", control_debias$delta_AR_rho,
    "sd", control_debias$delta_AR_sd,
    "Rsd", control_SIR$R_AR_sd
  )
  control <- c(control_debias, control_SIR)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  saveRDS(control, file.path(out_dir, "control.RDS"))

  clusterExport(clust, c("control_debias", "control_SIR"))

  ### PCR positivity, perfect testing ###
  imperfect <- FALSE
  type <- "PCR_positive"

  # Calculate regional bias parameters
  delta_df <- region_df %>%
    group_by(phe_region) %>%
    group_modify(~ cbind(mid_week = .x$mid_week, 
                         specify_delta_prior(.x, control_debias, imperfect)))
  
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

  delta_out_file <- file.path(out_dir, "delta_pcr_perfect.csv")
  readr::write_csv(delta_df, delta_out_file)

  ltla_out_file <- file.path(out_dir, "ltla_prevalence_pcr_perfect.RDS")
  saveRDS(ltla_prevalence, ltla_out_file)


  ### Imperfect testing, PCR positivity and infectiousness ###
  imperfect <- TRUE
  for (type in c("PCR_positive", "Infectious")) {
    
    # Calculate regional bias parameters
    delta_df <- region_df %>%
      group_by(phe_region) %>%
      group_modify(~ cbind(mid_week = .x$mid_week, 
                           specify_delta_prior(.x, control_debias, imperfect)))
    
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
       SIR_model_out_ltla <- sample_I_R_SIR_from_pre_calc_lik(d_ltla, I_log_lik, 
                                                              trans_mats, control_debias,
                                                              control_SIR)
      
       out_file <- file.path(out_dir, type, "SIR", paste0(ltla_name, ".RDS"))
       dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
       saveRDS(SIR_model_out_ltla, out_file)
    }

    # Fit SIR with cut dates (for n-week ahead predictions)
    for (mid_week_int in cut_dates) {
      mid_week_cut <- as.Date(mid_week_int, origin = "1970-01-01")
      print(c(mid_week_cut = mid_week_cut))
      
      results_file <- file.path(out_dir, type, paste0("cut_", mid_week_cut, ".RDS"))
      
      delta_df <- region_df %>%
        filter(mid_week <= mid_week_cut) %>%
        group_by(phe_region) %>%
        group_modify(~ cbind(mid_week = .x$mid_week, 
                             specify_delta_prior(.x, control_debias, imperfect)))
      
      ltla_list <- ltla_df %>%
         filter(mid_week <= mid_week_cut) %>%
         left_join(delta_df, by = c("phe_region", "mid_week")) %>%
         left_join(pcr_infectious_df,  by = c("phe_region", "mid_week")) %>%
         group_by(ltla) %>%
         group_split()
      ltla_names <- sapply(ltla_list, function(x) x$ltla[1]) 
      
      ltla_prevalence <- parLapply(clust, ltla_list, local_prevalence, 
                                   control_debias, imperfect, type)
      names(ltla_prevalence) <- ltla_names
      
      SIR_out <- foreach(ltla_curr = ltla_names, .verbose = TRUE,
                          .packages = "dplyr") %dopar% {
                           
                            d_ltla <- ltla_df %>%
                              filter(ltla == ltla_curr &
                                       mid_week <= mid_week_cut)
                           
                            I_log_lik <- ltla_prevalence[[ltla_curr]]$log_post
                            SIR_model_out_ltla <- sample_I_R_SIR_from_pre_calc_lik(d_ltla, I_log_lik, 
                                                                                   trans_mats, control_debias,
                                                                                   control_SIR)
                           
                            return(SIR_model_out_ltla)
                          }
      names(SIR_out) <- ltla_names
      saveRDS(SIR_out, results_file)
    }
  }
}
stopCluster(clust)


library(prevdebiasr)
source("scripts/SIR_utils.R")


trans_mats <- readRDS("transmats/poisson_SIR_epi_gamma_1.RDS")
pcr_infectious_df <- readr::read_csv("data/moment_match_infectious.csv")
ltla_df <- readr::read_csv("data/ltla.csv")

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


unique_regions <- unique(region_df$phe_region)
str(region_df)
pdf("/mnt/c/Temp/bsu_checking.pdf", 9, 9)
par(mfrow = c(3, 3), oma = c(2, 2, 5, 2))
alpha_testing <- 3e-4
imperfect <- FALSE#for (imperfect in c(FALSE, TRUE)) {
  type <- "PCR_positive"
  type_in_file_path <- paste0(type, "_", ifelse(imperfect, "Imperfect", "Perfect"))
  control_debias <- prevdebiasr::get_control_parameters(
    alpha_testing = alpha_testing
  )
  out_dir <- paste0(
    "output/AR", control_debias$delta_AR_rho,
    "sd", control_debias$delta_AR_sd,
    "Rsd", control_SIR$R_AR_sd
  )
  
  
  
    type_in_file_path <- paste0(type, "_", ifelse(imperfect, "Imperfect", "Perfect"))
    region_out_file <- file.path(out_dir, type_in_file_path, "phe_region_prevalence.RDS")
    region_prevalence <- readRDS(region_out_file)
  control_debias$quant_approx
    str(region_prevalence)
    region_prevalence[[1]]$log_post[1, ]
  # ###############
  # pdf("/mnt/c/Temp/bsu_checking_SIR_fitting_issue.pdf", 12, 9)
  # ind_curr <- 23
  # par(mfrow = c(2, 1), oma = c(2, 2, 6, 2))
  # for (imperfect in c(FALSE, TRUE)) {
  #   type_in_file_path <- paste0(type, "_", ifelse(imperfect, "Imperfect", "Perfect"))
  #   region_out_file <- file.path(out_dir, type_in_file_path, "phe_region_prevalence.RDS")
  #   region_prevalence <- readRDS(region_out_file)
  #   region_name <- "London"
  #   source("scripts/SIR_utils.R")
  #   d_region <- region_df %>%
  #     filter(phe_region == region_name)
  #     
  #   I_log_lik <- region_prevalence[[region_name]]$log_post
  #   plot(control_debias$I_seq,  I_log_lik[ind_curr, ], log = "x", ylim = c(-100, 0), ty = "l", ylab = "log posterior", xlab = "I", 
  #        main = paste0("imperfect = ", imperfect))
  #   mtext(outer = T, text = paste0(region_name, " ", d_region$mid_week[ind_curr]))
  # }  
  # dev.off()
  # 
  # plot(I_log_lik[1, ])
  # apply(I_log_lik, 1, max)
  # SIR_model_out_region <- sample_I_R_SIR_from_pre_calc_lik(
  #   d_region, I_log_lik,
  #   trans_mats, control_debias,
  #   control_SIR
  # )
  # 
  # I_seq = control_debias$I_seq
  # log_lik = I_log_lik
  # d = d_region; 
  # ##############
    
  # region_name <- unique_regions[1]#for (region_name in unique_regions) {
  #   out_file <- file.path(out_dir, type_in_file_path, "SIR", paste0(region_name, ".RDS"))
  #   SIR_model_out_region <- readRDS(out_file)
  #   
  # }
  
  out_files <- file.path(out_dir, type_in_file_path, "SIR", paste0(unique_regions, ".RDS"))
  SIR_model_results <- lapply(out_files, readRDS)
  names(SIR_model_results) <- unique_regions
  
  R_all <- prevent_all <- I_all <- data.frame()
  out_df <- data.frame()
  quant_plot <- c(0.025, 0.5, 0.975)
  for (region_curr in unique_regions) {#region_curr <- "London"#
    
    this_M <- region_df %>%
      filter(phe_region == region_curr) %>%
      pull(M)
    this_M <- this_M[1]
    saml_biased <- SIR_model_results[[region_curr]]
    its_keep <- (control_SIR$burn_in + 1):control_SIR$n_iters
    I_quant_curr <- t(apply(saml_biased$I[, its_keep], 1, function(v) quantile(v, quant_plot, na.rm = T))) / this_M * 100
    region_curr_df <- region_df[region_df$phe_region == region_curr, ]
    # region_curr_df <- region_curr_df[order(region_curr_df$mid_week), ]
    # prob_recover_in_week <- pexp(7, 1 / 7)
    # E_num_infections_prevented_mcmc_out <- (saml_biased$I[, its_keep] - region_curr_df$nt) / this_M * saml_biased$R[, its_keep] * prob_recover_in_week
    # SIR_model_results[[region_curr]]$prevent_quant <- t(apply(E_num_infections_prevented_mcmc_out, 1, 
    #                                                         function(v) quantile(v, quant_plot, na.rm = T)))
    SIR_model_results[[region_curr]]$I_quant <- I_quant_curr
    SIR_model_results[[region_curr]]$R_quant <- t(apply(saml_biased$R[, its_keep], 1, 
                                                      function(v) quantile(v, quant_plot, na.rm = T)))
    # SIR_model_results[[region_curr]]$R_quant[n_weeks, ]
    
    R_add <- as.data.frame((SIR_model_results[[region_curr]]$R_quant))
    I_add <- as.data.frame((SIR_model_results[[region_curr]]$I_quant))
    # prevent_add <- as.data.frame(t(SIR_model_results[[region_curr]]$prevent_quant[n_weeks, ]))
    names(I_add) <- paste0(c("l", "m", "u"), "_I")
    names(R_add) <- paste0(c("l", "m", "u"), "_R")
    region_curr_df <- cbind(region_curr_df, I_add, R_add)
    out_df <- rbind(out_df, region_curr_df)
    # names(prevent_add) <- paste0(c("l", "m", "u"), "_prevent")
    # R_add$region <- I_add$region <- region_curr
    # R_all <- rbind(R_all, R_add)
    # I_all <- rbind(I_all, I_add)
    # prevent_all <- rbind(prevent_all, prevent_add)
  }
  pdf("/mnt/c/Temp/phe_regions_bsu_plot.pdf", 9, 9)
  par(mfrow = c(3, 3))
  for (region_curr in unique_regions) {# region_curr <- "London"#
    region_curr_df <- out_df[out_df$phe_region == region_curr, ]
    plot(1:nrow(region_curr_df), region_curr_df$m_I, ty = "l", lwd = 2, ylim = c(0, 4), main = region_curr, xlab = "Week index", ylab = "Prevalence (%)")
    lines(1:nrow(region_curr_df), region_curr_df$l_I)
    lines(1:nrow(region_curr_df), region_curr_df$u_I)
  }
  dev.off()
  # mtext(outer = T, text = paste0("Imperfect testing ", imperfect))
# }
dev.off()






































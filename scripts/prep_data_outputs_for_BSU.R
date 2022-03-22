library(prevdebiasr)
library(readr)
library(readODS)
library(readxl)
library(dplyr)
library(tidyr)
library(zoo)
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

alpha_testing <- 3e-4
imperfect <- FALSE
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

out_files <- file.path(out_dir, type_in_file_path, "SIR", paste0(unique_regions, ".RDS"))
SIR_model_results <- lapply(out_files, readRDS)
names(SIR_model_results) <- unique_regions
saveRDS(SIR_model_results, file = "output/SIR_MCMC_outputs_for_Epi_ensemble_team.RDS")
out_df <- data.frame()
quant_plot <- c(0.025, 0.5, 0.975)
for (region_curr in unique_regions) {#region_curr <- "London"#
  this_M <- region_df %>%
    filter(phe_region == region_curr) %>%
    pull(M)
  this_M <- this_M[1]
  region_curr_df <- region_df[region_df$phe_region == region_curr, ]
  saml_biased <- SIR_model_results[[region_curr]]
  its_keep <- (control_SIR$burn_in + 1):control_SIR$n_iters
  # I SIR quantiles
  I_quant_curr <- t(apply(saml_biased$I[, its_keep], 1, function(v) quantile(v, quant_plot, na.rm = T))) / this_M * 100
  I_add <- as.data.frame(I_quant_curr)
  names(I_add) <- paste0(c("l", "m", "u"), "_I")
  
  # R SIR quantiles
  R_quant_curr <- t(apply(saml_biased$R[, its_keep], 1, function(v) quantile(v, quant_plot, na.rm = T)))
  R_add <- as.data.frame(R_quant_curr)
  names(R_add) <- paste0(c("l", "m", "u"), "_R")
  
  # I cross-sectional posterior quantiles
  norm_post_curr <- region_prevalence[[region_curr]]$norm_post
  cum_post_curr <- t(apply(norm_post_curr, 1, cumsum))
  Ix_quant_curr <- t(apply(cum_post_curr, 1, function(v) control_debias$I_seq[findInterval(quant_plot, v)])) / this_M * 100
  Ix_add <- as.data.frame(Ix_quant_curr)
  names(Ix_add) <- paste0(c("l", "m", "u"), "_Ix")

  # REACT cross-sectional mean and exact binomial CI
  react_quant_curr <- Hmisc::binconf(x = region_curr_df$nr, n = region_curr_df$Nr, method = "exact", alpha = .05)[, c(2, 1, 3)] * 100
  react_quant_curr[rowSums(is.na(react_quant_curr)) > 0, ] <- NA
  react_add <- as.data.frame(react_quant_curr)
  names(react_add) <- paste0(c("l", "m", "u"), "_react")
  
  region_curr_df <- cbind(region_curr_df, I_add, R_add, Ix_add, react_add)
  out_df <- rbind(out_df, region_curr_df)
}

pdf("plots/phe_regions_bsu_plot.pdf", 9, 9)
  par(mfrow = c(3, 3), oma = c(2, 2, 6, 2))
  for (region_curr in unique_regions) {# region_curr <- "London"#
    region_curr_df <- out_df[out_df$phe_region == region_curr, ]
    plot(1:nrow(region_curr_df), region_curr_df$m_I, ty = "n", lwd = 2, ylim = c(0, 4), main = region_curr, xlab = "Week index", ylab = "Prevalence (%)")
    lines(1:nrow(region_curr_df), region_curr_df$m_Ix, col = 2, lwd = 2)
    lines(1:nrow(region_curr_df), region_curr_df$l_Ix, col = 2, lwd = 1)
    lines(1:nrow(region_curr_df), region_curr_df$u_Ix, col = 2, lwd = 1)
    lines(1:nrow(region_curr_df), region_curr_df$m_I, lwd = 2)
    lines(1:nrow(region_curr_df), region_curr_df$l_I)
    lines(1:nrow(region_curr_df), region_curr_df$u_I)
    points(1:nrow(region_curr_df), region_curr_df$m_react, col = 4)
    for(j in 1:nrow(region_curr_df)) {
      lines(rep(j, 2), c(region_curr_df$l_react[j], region_curr_df$u_react[j]), col = 4, lwd = 1)
    }
    if (match(region_curr, unique_regions) == 1) {
      par(xpd = NA)
      legend(x = 0, y = 7.5, legend = c("REACT raw cross-sectional", "Debiasing cross-sectional", "Debiasing SIR"), title = "Prevalence outputs",
             col = c(4, 2, 1), lty = 1, lwd = 2, bg = "white")
      par(xpd = F)
    }
  }
  mtext(outer = T, text = "PCR+ prevalence")
  for (region_curr in unique_regions) {# region_curr <- "London"#
    region_curr_df <- out_df[out_df$phe_region == region_curr, ]
    plot(1:nrow(region_curr_df), region_curr_df$m_R, ty = "l", lwd = 2, ylim = c(0, 2.5), main = region_curr, xlab = "Week index", ylab = "Rt")
    lines(1:nrow(region_curr_df), region_curr_df$l_R)
    lines(1:nrow(region_curr_df), region_curr_df$u_R)
  }
  mtext(outer = T, text = "Effective R number")
dev.off()


dput(names(out_df))
names_df <- data.frame(keep = c("phe_region", "mid_week", "Nt", "nt", "Nr", "nr", "M", "l_I", "m_I", "u_I", "l_R", "m_R", "u_R", "l_Ix", 
  "m_Ix", "u_Ix", "l_react", "m_react", "u_react"),
  name_out = c("region", "mid_week", "P12_N", "P12_n", "react_N", "react_n", "pop_size", "l_I", "m_I", "u_I", "l_R", "m_R", "u_R", "l_Ix", 
               "m_Ix", "u_Ix", "l_react", "m_react", "u_react"))
out_df_write <- out_df[, names_df$keep]
names(out_df_write) <- names_df$name_out
write.csv(out_df_write, file = "data/outputs_to_BSU.csv", row.names = F)

d <- read.csv(file = "data/outputs_to_BSU.csv")
#






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



























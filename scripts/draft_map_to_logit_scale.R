
### Estimate LTLA prevalence ###
library(dplyr)
library(prevdebiasr)
library(parallel)
library(foreach)
library(truncnorm)
source("scripts/SIR_utils.R")

trans_mats <- readRDS("transmats/poisson_SIR_epi_gamma_1.RDS")
vax_df <- readr::read_csv("data/vaccination.csv")
pcr_infectious_df <- readr::read_csv("data/moment_match_infectious.csv")
region_df <- readr::read_csv("data/region.csv") %>%
  left_join(pcr_infectious_df, by = c("phe_region", "mid_week"))
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

reg_df_curr <- region_df[region_df$phe_region == "South West", ]
reg_df_curr[nrow(reg_df_curr) - 10:0, ]
mid_week_unique <- sort(unique(ltla_df$mid_week))

n_cores <- 12
run_type <- c("fast", "full")[2]
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

# Estimate local prevalence for each LTLA
ltla_list <- ltla_df %>%
  left_join(delta_df, by = c("phe_region", "mid_week")) %>%
  group_by(ltla) %>%
  group_split()
ltla_names <- sapply(ltla_list, function(x) x$ltla[1])
ltla_prevalence <- parLapply(clust, ltla_list, local_prevalence,
                             control_debias, imperfect, type)
names(ltla_prevalence) <- ltla_names

# str(ltla_prevalence[[ltla_curr]])

# Code from mapping to logit scale

ltla_curr <- "Adur"
ltla_unique <- unique(ltla_df$ltla)
out_all <- data.frame()
for (ltla_curr in ltla_unique) {
  M_curr <- unlist(ltla_df[match(ltla_curr, ltla_df$ltla), "M"])
  # NB: this next line will import log_lik
  output_to_transform <- c("log_post", "log_lik")[2]
  log_lik <- ltla_prevalence[[ltla_curr]][[output_to_transform]]
  dimnames(log_lik) <- dimnames(ltla_prevalence[[ltla_curr]]$log_post)
  log_lik <- log_lik - apply(log_lik, 1, max)
  norm_lik <- exp(log_lik) / rowSums(exp(log_lik))
  norm_lik_dezeroed <- norm_lik
  norm_lik_dezeroed[, "0"] <- 0
  norm_lik_dezeroed[, "1"] <- norm_lik[, "0"] + norm_lik[, "1"]
  logit_scale_prev_propn_seq <- boot::logit(control$I_seq / M_curr)
  mom_match_1 <- rowSums(sweep(norm_lik_dezeroed, 2, logit_scale_prev_propn_seq, '*'), na.rm = TRUE)
  mom_match_2 <- rowSums(sweep(norm_lik_dezeroed, 2, logit_scale_prev_propn_seq^2, '*'), na.rm = TRUE)
  mean_on_logit_scale <- mom_match_1
  sd_on_logit_scale <- sqrt(mom_match_2 - mom_match_1^2)
  add <- data.frame(mean = mean_on_logit_scale, sd = sd_on_logit_scale, ltla = ltla_curr, mid_week = mid_week_unique)
  out_all <- rbind(out_all, add)
}
out_all <- out_all[order(out_all$ltla, out_all$mid_week), ]
# Sanity check :)
# pdf("/mnt/c/Temp/logit_sanity_check.pdf", 9, 9)
par(mfrow = c(4, 4))
for (j in 1:16) {
  plot(logit_scale_prev_propn_seq, norm_lik_dezeroed[j, ], ty = "l")
  logit_scale_for_gaussian <- seq(-20, 0, len = 200)
  raw_gauss_dens <- dnorm(logit_scale_for_gaussian, mean = mean_on_logit_scale[j], sd = sd_on_logit_scale[j])
  gauss_cdf <- pnorm(logit_scale_prev_propn_seq, mean = mean_on_logit_scale[j], sd = sd_on_logit_scale[j])
  map_back_to_bins <- diff(c(0, gauss_cdf))
  map_back_to_bins[is.na(map_back_to_bins)] <- 0
  lines(logit_scale_prev_propn_seq, map_back_to_bins, col = 2)
}



# dev.off()
# 
logit_mom <- function(x, values) {
  logit_values <- boot::logit(values)
  norm_prob <- exp(x) / sum(exp(x))
  ind_zero <- which(values == 0)
  ind_one <- which(values == 1)
  norm_prob_dezeroed <- norm_prob
  norm_prob_dezeroed[ind_zero] <- 0
  norm_prob_dezeroed[ind_one] <- norm_prob[ind_zero] + norm_prob[ind_one]
  mom_match1 <- sum(norm_prob_dezeroed * logit_values, na.rm = TRUE)
  mom_match2 <- sum(norm_prob_dezeroed * (logit_values ^ 2), na.rm = TRUE)
  data.frame(mean = mom_match1,
             sd =  sqrt(mom_match2 - (mom_match1 ^ 2)))
}

ltla_curr <- "Adur"
M_curr <- unlist(ltla_df[match(ltla_curr, ltla_df$ltla), "M"])
test <- logit_mom(x = ltla_prevalence[[ltla_curr]][[output_to_transform]], values = control$I_seq / M_curr)


d <- read.csv("/mnt/c/Temp/logit_moments_v2.csv")
d$mid_week <- format(strptime(d$mid_week, format = "%d/%m/%Y"), format = "%Y-%m-%d")

pdf("/mnt/c/Temp/logit_cross_check_means.pdf", 9, 9)
par(mfrow = c(4, 4))
for(ltla_curr in ltla_unique[1:16]) {
  plot(d[d$ltla == ltla_curr, "mean"], out_all[out_all$ltla == ltla_curr, "mean"], main = ltla_curr, xlab = "Brieuc", ylab = "George")
  abline(0, 1)
}
dev.off()

pdf("/mnt/c/Temp/logit_cross_check_sds.pdf", 9, 9)
par(mfrow = c(4, 4))
for(ltla_curr in ltla_unique[1:16]) {
  plot(d[d$ltla == ltla_curr, "sd"], out_all[out_all$ltla == ltla_curr, "sd"], main = ltla_curr, xlab = "Brieuc", ylab = "George")
  abline(0, 1)
}
dev.off()


head(d)
head(out_all)
str(d)
out_all <- out_all[order(out_all$ltla, out_all$mid_week), ]
d <- d[order(d$ltla, d$mid_week), ]
plot(d$mean, out_all$mean)
plot(d$sd, out_all$sd)

str(d)
str(out_all)

# Sanity check :)



out_all[1:30, ]

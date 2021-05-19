library(prevdebiasr)
library(dplyr)

pcr_infectious_df <- readr::read_csv("data/moment_match_infectious.csv")
region_df <- readr::read_csv("data/region.csv") %>%
  left_join(pcr_infectious_df, by = c("phe_region", "mid_week"))

imperfect <- TRUE
control <- get_control_parameters()
control$delta_AR_rho <- 0.99
control$delta_AR_sd <- 1
control$R_AR_sd <- 0.2

# Calculate quantiles for prevalence based on randomised testing data
I_post_quant <- prevdebiasr:::randomised_testing_prevalence(region_df,
                                                            control,
                                                            imperfect)

n_time <- nrow(region_df)
n_quant <- control$n_quant_approx_bias
n_quant_p1 <- 100

quant_approx <- seq(n_quant) / (n_quant + 1)
delquants <- NULL
nu_quant <- matrix(NA, n_time, n_quant)

nu <- boot::logit((region_df$Nt - region_df$nt) / region_df$M)
region_df$p1_low <- NA
region_df$p1_high <- NA
region_df$p2 <- boot::inv.logit(nu)

for (qnum in 1:n_quant) {
  I <- I_post_quant[ ,qnum]
  
  if (imperfect) {
    # Get range of feasible p1 = ilogit(nu + delta) (uniform prior over this range)
    for (i in 1:n_time) {
      # No false negatives
      true_pos_fn <- min(I[i], qbinom(1e-10,
                                      region_df$nt[i],
                                      1 - control$alpha_testing))
      region_df$p1_low[i] <- qbeta(1e-10,
                                 true_pos_fn + 1,
                                 I[i] - true_pos_fn + 1)
      
      # No false positives
      true_pos_fp <- min(I[i], region_df$nt[i] + qbinom(1 - 1e-10,
                                                      region_df$Nt[i] - region_df$nt[i],
                                                      control$beta_testing))
      region_df$p1_high[i] <- qbeta(1 - 1e-10,
                                  true_pos_fp + 1,
                                  I[i] - true_pos_fp + 1)
    }
    
    ll_targeted <- matrix(NA, n_time, n_quant_p1)
    p1_quants <- mapply(function(x, y) seq(x, y, length.out = n_quant_p1),
                        region_df$p1_low, region_df$p1_high)

    for (p1_qnum in 1:n_quant_p1) {

      region_df$p1 <- p1_quants[p1_qnum, ]

      ll_targeted[ ,p1_qnum] <- targeted_testing_loglik(region_df, I,
                                                        control, imperfect)
      
    }
    ll_prev <- apply(ll_targeted, 1, normalise_logprob)

    # Moment matching for beta shape parameters
    p1_mean <- colSums(exp(ll_prev) * p1_quants)
    p1_var <- colSums(exp(ll_prev) * (p1_quants ^ 2)) - (p1_mean ^ 2)
    p1_var <- pmax(1e-10, p1_var) # for numerical stability
    beta_shape1 <- p1_mean * (((p1_mean * (1 - p1_mean)) / p1_var) - 1)
    beta_shape2 <- (1 - p1_mean) * (((p1_mean * (1 - p1_mean)) / p1_var) - 1)
    
    # Get quantiles for gamma = nu + delta
    gamma_quant <- matrix(NA, n_time, n_quant)
    for (i in which(!is.na(region_df$nt) & (region_df$nt < ((1 - control$beta_testing) * I)))) { # Ignore I quantiles incompatible with nt here
      gamma_quant[i, ] <- boot::logit(stats::qbeta(p = quant_approx,
                                                   shape1 = beta_shape1[i],
                                                   shape2 = beta_shape2[i]))
    }
    delta_quant <- gamma_quant - nu
    delquants <- cbind(delquants, delta_quant)
    
  } else {
    # Below posterior is based on a uniform, i.e. beta(1,1), prior on
    # logit.inv(gamma) = logit.inv(nu + delta)
    beta_shape1 <- region_df$nt + 1
    beta_shape2 <- I - region_df$nt + 1
    
    gamma_quant <- matrix(NA, n_time, n_quant)
    for (i in which(!is.na(region_df$nt) & beta_shape2 > 0)) { # Ignore I quantiles incompatible with nt here
      gamma_quant[i, ] <- boot::logit(stats::qbeta(p = quant_approx,
                                                   shape1 = beta_shape1[i],
                                                   shape2 = beta_shape2[i]))
    }
    delta_quant <- gamma_quant - nu
    delquants <- cbind(delquants, delta_quant)
  }
}

delta_post_moment1 <- rowMeans(delquants, na.rm = TRUE)
delta_post_moment2 <- rowMeans(delquants^2, na.rm = TRUE)
delta_post_mean <- delta_post_moment1
delta_post_sd <- sqrt(delta_post_moment2 - (delta_post_moment1)^2)


# Examine approximation
id <- paste0("AR", control$delta_AR_rho,
             "sd", control$delta_AR_sd,
             "Rsd", control$R_AR_sd)
plot_dir <- file.path("plots", id, "PCR_positive")

for (phe_reg in c("London", "North West")) {
  indvc <- which(region_df$phe_region == phe_reg & 
                   region_df$mid_week >= as.Date("2020-11-23"))
  graphics.off()
  
  pdf(file.path(plot_dir, paste0("SI_moment_match_", phe_reg, ".pdf")), 8, 7)
  par(mfrow = c(3, 3))
  for (j in indvc) {
    hist(delquants[j, ], probability = TRUE, main = "", xlab = region_df$mid_week[j],
         breaks=seq(min(delquants[j, ], na.rm=T),max(delquants[j, ], na.rm=T), l = 9))
    xseq <- quantile(delquants[j, ], seq(0.0, 1, by = 0.01), na.rm = TRUE)
    lines(x = xseq, 
          y = dnorm(xseq, delta_post_mean[j], delta_post_sd[j]),
          col = "red")
  }
  mtext(phe_reg, side = 3, line = -3, outer = TRUE)
  dev.off()
}



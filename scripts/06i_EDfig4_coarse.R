library(prevdebiasr)
control <- get_control_parameters()
############################
# Import weekly Pillar 2 and REACT
############################

ltla_pop <- readr::read_csv("data/ltla_pop.csv")
react_date_df <- readr::read_csv("data/react_dates.csv")
ltla_df <- readr::read_csv("data/ltla.csv")
raw_pillar2_df <- readr::read_csv("data/raw_pillar2.csv")
react_ltla_df <- readr::read_csv("data/react_ltla.csv")


# Quantiles to plot
quant_plot <- c(0.025, 0.5, 0.975)

id <- "AR0.99sd1Rsd0.2"

out_dir <- file.path("output", id)
plot_dir <- file.path("plots", id)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

ltla_df <- readr::read_csv("data/ltla.csv")

debiased_ltla_strategy <- c("pcr_perfect", "national", "no_round8")
debiased_estimates <- list()
for (strategy in debiased_ltla_strategy) {
  
  debiased_ltla <- readRDS(file.path(out_dir, 
                                     paste0("ltla_prevalence_", strategy, ".RDS")))
  
  # Extract REACT round output
  react_ind <- which(unique(ltla_df$mid_week) %in% react_date_df$mid_week)
  debiased_ll <- lapply(debiased_ltla,
                        function(x) t(apply(x$norm_post[react_ind, ], 1, function(v) { 
                          if(any(is.na(v))) { 
                            return( rep(NA, length(quant_plot))) 
                          } else { 
                            return(control$I_seq[findInterval(quant_plot, 
                                                              cumsum(v)) + 1])
                          }
                        })))
  
  debiased_ltla_estimates <- abind::abind(
    debiased_ll, along = 3, new.names = list(as.character(7:11), 
                                             c("low", "mid", "up"),
                                             unique(raw_pillar2_df$ltla))
  ) %>%
    reshape2::melt(varnames = c("round", "est", "ltla")) %>%
    pivot_wider(names_from = "est") %>%
    left_join(ltla_pop, by = "ltla") %>%
    mutate(m = mid / M * 100,
           l = low / M * 100,
           u = up / M * 100)
  
  debiased_estimates[[strategy]] <- debiased_ltla_estimates
}
graphics.off()
cm2in <- 0.39
pdf(file.path(plot_dir, "SI_national_bias_correction_789.pdf"), width = 18 * cm2in, 18 * cm2in, pointsize = 7)

par(mfrow = c(4, 3), oma = c(2, 2, 2, 10), mar = c(4, 4, 2.5, 2.5))
max_prev <- 5
max_prev_uncorr <- 30
rounds_todo <- 7:9
panel_count <- 0
point_col <- rgb(red = .5, green = .4, blue = .8, alpha = .75)

for (meth in c("raw", debiased_ltla_strategy)) {
  for(this_round in rounds_todo) { 
    panel_count <- panel_count + 1
    
    ###########################################
    # Plot uncorrected Pillar 2 
    ########################################
    comp_1 <- react_ltla_df %>%
      ungroup() %>%
      filter(round == this_round) %>%
      select(l, m, u)
    
    if (meth == "raw") {
      comp_2 <- raw_pillar2_df %>%
        ungroup() %>%
        filter(round == this_round) %>%
        select(l, m, u)
      
      ylab <- "Uncorrected Pillar 1+2 prevalence (%)"
    } else {
      comp_2 <- debiased_estimates[[meth]] %>%
        ungroup() %>%
        filter(round == this_round) %>%
        select(l, m, u)
      
      ylab <- "Debiased prevalence (%)"
    }
    
    cor_est <- round(cor(comp_1$m, comp_2$m, use = "p", me = "sp"), 2)
    bias_mn <- mean(comp_2$m - comp_1$m, na.rm = T)
    bias_se <- sd(comp_2$m - comp_1$m, na.rm = T) / sqrt(nrow(comp_1))
    
    plot(comp_1$m, comp_2$m, main = "", xlim = c(0, max_prev), 
         ylim = c(0, ifelse(meth == "raw", max_prev_uncorr, max_prev)),
         xlab = "", ylab = "", ty = "n")
    title(ylab = ylab, xlab = "REACT prevalence (%)", line = 3)
    
    for (j in 1:nrow(comp_1)) {
      lines(x = comp_1[j, c("l", "u")], y = rep(comp_2$m[j], 2), col = grey(.8))
    }
    for (j in 1:nrow(comp_1)) {
      lines(x = rep(comp_1$m[j], 2), y = comp_2[j, c("l", "u")])
    }
    points(comp_1$m, comp_2$m, col = point_col, pch = 19, cex = .4)
    abline(0, 1)
    
    mtext(side = 3, 
          text = paste0("Bias = ", 
                        round(bias_mn, 2), 
                        "% (SE = ", 
                        round(bias_se, 2), 
                        ")"), 
          cex = .8, line = -2)
    
    mtext(side = 3, at = -.5, text = paste0("(", letters[panel_count], ")"), 
          cex = 1, line = .5)
  }
}
# mtext(side = 1, outer = T, text = "REACT prevalence in % (95% CIs)", line = 1)
# mtext(side = 2, outer = T, text = "Pillar 1+2 prevalence in % (95% CIs)", line = 1)
mtext(side = 4, outer = T, text = "Uncorrected", line = 0, at = .875, las = 1)
mtext(side = 4, outer = T, text = "PHE regional \ndelta", line = 0, at = .625, las = 1)
mtext(side = 4, outer = T, text = "National delta", line = 0, at = .375, las = 1)
mtext(side = 4, outer = T, text = "Without REACT \nround 8", line = 0, at = .125, las = 1)
mtext(side = 3, outer = T, text = "REACT round 7", line = -1, at = .15, las = 1)
mtext(side = 3, outer = T, text = "REACT round 8", line = -1, at = .5, las = 1)
mtext(side = 3, outer = T, text = "REACT round 9", line = -1, at = .85, las = 1)

dev.off()

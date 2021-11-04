### Extended Data Figure ###
library(tidyr)

control <- prevdebiasr::get_control_parameters()
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

debiased_ltla <- readRDS(file.path(out_dir, "ltla_prevalence_pcr_perfect.RDS"))

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

graphics.off()

pdf(file.path(plot_dir, "fig2_bias_correction_1011.pdf"), width = 16 * cm2in, 12 * cm2in, pointsize = 7)

par(mfrow = c(2, 2), oma = c(2, 2, 2, 2), mar = c(4, 4, 2.5, 2.5))
max_prev <- 1.5
max_prev_uncorr <- 8
rounds_todo <- 10:11
panel_count <- 0
point_col <- rgb(red = .5, green = .4, blue = .8, alpha = .75)

for (meth in c("raw", "debiased")) {
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
    }
    
    if (meth == "debiased") {
      comp_2 <- debiased_ltla_estimates %>%
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
    
    mtext(side = 3, at = -.2, text = paste0("(", letters[panel_count], ")"), 
          cex = 1, line = .5)
  }
}
# mtext(side = 1, outer = T, text = "REACT prevalence in % (95% CIs)", line = 1)
# mtext(side = 2, outer = T, text = "Pillar 1+2 prevalence in % (95% CIs)", line = 1)
# mtext(side = 4, outer = T, text = "Uncorrected", line = 0, at = .75, las = 1)
# mtext(side = 4, outer = T, text = "Corrected", line = 0, at = .25, las = 1)
mtext(side = 3, outer = T, text = "REACT round 10", line = -1, at = .25, las = 1)
mtext(side = 3, outer = T, text = "REACT round 11", line = -1, at = .75, las = 1)

dev.off()

raw_pillar2_df %>%
  group_by(round) %>%
  summarise(range(Nt)) %>%
  print(n=10)

react_ltla_df %>%
  group_by(round) %>%
  summarise(range(number_samples)) %>%
  print(n=10)

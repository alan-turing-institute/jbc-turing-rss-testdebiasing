###########################################
# Plot variation in delta across PHE regions 
###########################################

library(dplyr)
source("scripts/plot_utils.R")

plot_delta_curves <- function(dl, phe_regs, del_axis_side = 4,
                              ylim_use = c(1, 4.5), add_legend = F, legend_shift = 0,
                              add_t_cut = F, add_react_round_legend = F) {

  d <- filter(dl, phe_region == "London")
  regnames <- dl %>% pull(phe_region) %>% unique()
  plot_map <- data.frame(date = seq(min(as.Date(d$mid_week)) - 3, 
                                    max(as.Date(d$mid_week)) + 3, 
                                    by = 1))
  plot_map$day_index <- 1:nrow(plot_map)

  all_phe_regions_del_mean <- dl %>% 
    group_by(phe_region) %>% 
    group_map(~ .x$delta_prior_mean) %>% 
    do.call(what = cbind)
  colnames(all_phe_regions_del_mean) <- regnames

  all_phe_regions_del_sd <- dl %>% 
    group_by(phe_region) %>% 
    group_map(~ .x$delta_prior_sd) %>% 
    do.call(what = cbind)
  colnames(all_phe_regions_del_sd) <- regnames
  
  mid_pl <- all_phe_regions_del_mean[, phe_regs, drop = F]
  upp_pl <- mid_pl + 2 * all_phe_regions_del_sd[, phe_regs, drop = F]
  low_pl <- mid_pl - 2 * all_phe_regions_del_sd[, phe_regs, drop = F]
  xpl <- plot_map$day_index
  col_reg <- RColorBrewer::brewer.pal(9, "Set1")
  col_reg <- col_reg[c(9, 1:8)]
  names(col_reg) <- phe_regs
  matplot(x = 0, y = 0, ty = "n", ylim = ylim_use, xlim = range(xpl), 
          xaxt = "n", yaxt = "n", xlab = "", las = 2, ylab = "", xaxs = "i")
  xpl_week <- plot_map[match(d$mid_week, plot_map$date), "day_index"]
  for(phe_region in phe_regs) {
    colc <- col_reg[phe_region]
    lines(xpl_week, mid_pl[, phe_region], lty = 1, lwd = 3, col = colc)
    lines(xpl_week, upp_pl[, phe_region], lty = 1, lwd = 1, col = colc)
    lines(xpl_week, low_pl[, phe_region], lty = 1, lwd = 1, col = colc)
  }

  axis(side = del_axis_side, las = 2)
  mtext(side = del_axis_side, text = expression(delta), las = 1, line = 3, cex = 1.3)
  add_react_sampling_intervals_to_plot(d = plot_map, prev_max = ylim_use[2], 
                                       col_data = col_data, type = "daily",
                                       prev_min = ylim_use[1])
  annotate_months(plot_map, add_axis = T, shade_alpha = .2, 
                  for_presentation = T)
  leg_move_mult <- .8
  if (add_legend) {
    par(xpd = NA)
    legend(x = legend_shift, y = ylim_use[2],#x = min(xpl) - diff(range(xpl)) * leg_move_mult, y = ylim_use[2], 
           legend = phe_regs, 
           col = col_reg[phe_regs], lty = 1, lwd = 3, cex = 0.8, 
           bg = "white")
    par(xpd = F)
  }  
  if (add_react_round_legend) {
    par(xpd = NA)
    legend(x = legend_shift, y = ylim_use[1], # , x = min(xpl) - diff(range(xpl)) * leg_move_mult
           legend = "REACT rounds", 
           col = col_data["react_samples"], pch = 15, cex = 0.8, bg = "white", pt.cex = 1.6)
    par(xpd = F)
  }  
}

id <- "AR0.99sd1Rsd0.2"
type <- "Infectious"
out_dir <- file.path("output", id, type)
plot_dir <- file.path("plots", id, type)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

delta_df <- readr::read_csv(file.path(out_dir, "delta.csv"))
cm2in <- 0.39
pdf(paste0(plot_dir, "/delta_estimates.pdf"), width = 18 * cm2in, 7 * cm2in, pointsize = 7)
par(mfrow = c(1, 2), oma = c(1, 1, 1, 11), mar = c(3, 3, 2, 3))
phe_regs <- c("London", "Yorkshire and The Humber", "South West", 
              "East of England", "West Midlands", "North West", "South East", 
              "North East", "East Midlands")
plot_delta_curves(delta_df, phe_regs, del_axis_side = 2, add_legend = T, legend_shift = 950,
                  add_t_cut = F, add_react_round_legend = T, ylim_use = c(1, 4.5))
plot_delta_curves(delta_df, "London", del_axis_side = 2, add_legend = F, legend_shift = 0,
                  add_t_cut = F, add_react_round_legend = F, ylim_use = c(1, 4.5))
mtext(expression(bold("Log-odds ratio of being tested in the infected versus non-infected populations, " ~ delta)),
      line = 1, cex = 1.2, adj = 1.2)
dev.off()



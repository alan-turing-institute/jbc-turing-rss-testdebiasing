#############################################################
# Source all functions, get global parameters
#############################################################
source("scripts/plot_utils.R")
control <- prevdebiasr::get_control_parameters()

param_df <- expand.grid(delta_AR_rho = c(0.975, 0.99, 0.999),
                        delta_AR_sd = c(1, 2),
                        R_AR_sd = c(0.2),
                        type = c("PCR_positive"))

for (i in 1:nrow(param_df)) {
  #############################################################
  # Process debiased prevalence output
  #############################################################
  id <- paste0("AR", param_df$delta_AR_rho[i],
               "sd", param_df$delta_AR_sd[i],
               "Rsd", param_df$R_AR_sd[i])
  out_dir <- file.path("output", id)
  plot_dir <- file.path("plots", id)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  ltla_prevalence <- readRDS(file.path(out_dir, "ltla_prevalence_pcr_perfect.RDS"))
  
  norm_cdf_I_pl_all <- list()
  for(this_ltla in names(ltla_prevalence)) {
    this_M <- ltla_pop %>%
      filter(ltla == this_ltla) %>%
      pull(M)
    norm_cdf_bin_inds <- t(apply(ltla_prevalence[[this_ltla]]$norm_post, 1, 
                                 function(v) findInterval(qplot, cumsum(v), all.inside = T)))
    norm_cdf_I_pl <- cbind(control$bin.d$e.l[norm_cdf_bin_inds[, 1]], 
                           control$bin.d$mu[norm_cdf_bin_inds[, 2]],
                           control$bin.d$e.u[norm_cdf_bin_inds[, 3]]) / this_M * 100
    norm_cdf_I_pl_all[[this_ltla]] <- norm_cdf_I_pl
  }
  
  plot_map <- data.frame(date = seq(min(ltla_df$mid_week) - 3, 
                                    max(ltla_df$mid_week) + 3, 
                                    by = 1))
  plot_map$day_index <- 1:nrow(plot_map)
  xpl <- plot_map$day_index
  
  ###########################################################
  # Date range for R prediction in paper
  ###########################################################
  #date_range_out <- paste(range(ir$mid_week_cut), collapse = " - ")
  #write.table(date_range_out, file = paste(overleaf$numbers, "/R_pred_date_range.txt", sep = ""), 
  #            col.names = F, row.names = F, quote = F)
  
  
  ###########################################################
  #  R prediction Figure in paper
  ###########################################################
  graphics.off()
  pdf(file.path(plot_dir, "uncorrected_and_cross_sec_corr.pdf"), 12, 6)
  
  ltla_plot <- choose_focus_ltlas(ltla_df)
  par(mfrow = c(2, 5), mar = c(3, 3, 2, 1))
  pch_ests <- 19
  cex_ests <- 1
  col_curr <- list(corrected = rgb(0, 0, 1, alpha = .5), uncorrected = rgb(1, 0, 0, alpha = .5))
  plot(0, xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n", ty = "n")
  legend(x = "center", legend = c("Raw Pillar 1+2", "Corrected", "REACT estimate", "REACT rounds"),
         pch = c(rep(pch_ests, 3), 15), pt.cex = c(rep(cex_ests, 3), 1.6), 
         col = c(col_curr$uncorrected, col_curr$corrected, rep(col_data["react_samples"], 2)), cex = 1.2)
  
  for (this_ltla in ltla_plot) {
    colc <- 1
    d <- ltla_df %>%
      filter(ltla == this_ltla)
    d$raw_prop <- d$nt / d$Nt * 100
    raw_prop_bin_conf <- Hmisc::binconf(x = d$nt, n = d$Nt) * 100
    
    ymin <- .0
    ylim_use <- c(0, max(d$raw_prop))
    
    matplot(x = 0, y = 0, ty = "n", ylim = ylim_use, xlim = range(xpl), 
            xaxt = "n", yaxt = "s", xlab = "", las = 2, ylab = "", xaxs = "i", 
            main = this_ltla)
    mtext(side = 2, text = "%", line = 2, las = 2, cex = .8)
    abline(h = 5, col = "light grey", lty = 3)
    add_react_sampling_intervals_to_plot(d = plot_map, prev_max = ylim_use[2], 
                                         col_data = col_data, type = "daily")
    xpl_week <- plot_map[match(d$mid_week, plot_map$date), "day_index"]
    annotate_months(plot_map, add_axis = T, shade_alpha = .2, for_presentation = T)
    
    for(j in 1:length(xpl_week)) {
      points(x = xpl_week[j], y = norm_cdf_I_pl_all[[this_ltla]][j, 2], lty = 1, 
             lwd = 1, col = col_curr$corrected, pch = 19, cex = .8)
      lines(x = rep(xpl_week[j], 2), 
            y = pmax(ymin / 2, norm_cdf_I_pl_all[[this_ltla]][j, c(1, 3)]), 
            lty = 1, lwd = 1, col = col_curr$corrected)
      points(x = xpl_week[j], y = raw_prop_bin_conf[j, 1], lty = 1, lwd = 1, 
             col = col_curr$uncorrected, pch = 19, cex = .8)
      lines(x = rep(xpl_week[j], 2), 
            y = pmax(ymin / 2, raw_prop_bin_conf[j, 2:3]), 
            lty = 1, lwd = 1, col = col_curr$uncorrected)
    }
    
    # REACT plots
    for (this_round in c(7, 8)) {
      react_quants <- react_ltla_df %>%
        filter(ltla == this_ltla & round == this_round) %>%
        select(l,m,u) %>%
        as.matrix()
      
      mid_round_curr <- react_date_df %>%
        filter(round == this_round) %>%
        pull(mid_date)
      lines(x = rep(plot_map[match(mid_round_curr, plot_map$date), "day_index"], 2), 
            y = react_quants[c(1, 3)], col = col_data["react_samples"])
      points(x = plot_map[match(mid_round_curr, plot_map$date), "day_index"], 
             y = react_quants[2], col = col_data["react_samples"], pch = 19)
    }
    mtext(outer = T, side = 2, text = "Prevalence %")
  }
  
  dev.off()
  
  if (id == "AR0.99sd1Rsd0.2") {
    overleaf_dir <- "~/Dropbox/Apps/Overleaf/Estimating local prevalence from targeted testing data/figures"
    file.copy(paste0(plot_dir, "/uncorrected_and_cross_sec_corr.pdf"),
              paste0(overleaf_dir, "/uncorrected_and_cross_sec_corr.pdf"),
              overwrite = TRUE)
  }
}

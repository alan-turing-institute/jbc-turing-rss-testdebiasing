cut_date <- "2021-01-17" 
source("scripts/plot_utils.R")
source("scripts/SIR_utils.R")

param_df <- expand.grid(delta_AR_rho = c(0.975, 0.99, 0.999),
                        delta_AR_sd = c(1, 2),
                        R_AR_sd = c(0.2),
                        type = c("PCR_positive", "Infectious"))

for (i in 1:nrow(param_df)) {
  
  id <- paste0("AR", param_df$delta_AR_rho[i],
               "sd", param_df$delta_AR_sd[i],
               "Rsd", param_df$R_AR_sd[i])
  type <- param_df$type[i]
  
  out_dir <- file.path("output", id, type)
  plot_dir <- file.path("plots", id, type)
  
  out_files <- list.files(file.path(out_dir, "SIR"), full.names = TRUE)
  ltla_unique <- sub(".RDS", "", basename(out_files))
  SIR_model_results <- lapply(out_files, readRDS)
  names(SIR_model_results) <- sub(".RDS", "", basename(out_files))
  
  rd <- read.csv("data/UK_hotspot_Rt_estimates_Imperial.csv")
  
  its_keep <- (control_SIR$burn_in + 1):control_SIR$n_iters
  
  ltla_plot <- sort(choose_focus_ltlas(ltla_df))
  
  
  graphics.off()
  pdf(file.path(plot_dir, "R_comparison_for_paper.pdf"), 11, 7)
  
  par(mfrow = c(2, 5), mar = c(3, 3, 2, 1), oma = c(2, 4, 2, 2))
  pch_ests <- 19
  cex_ests <- 1
  col_curr <- list(corrected = rgb(0, 0, 1, alpha = .5), uncorrected = rgb(1, 0, 0, alpha = .5))
  plot(0, xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
  col_use <- list(debias = rgb(0, 0, 1, alpha = .5), imperial = rgb(1, 0, 0, alpha = .5))
  legend(x = "center", legend = c("De-biasing model", "Imperial model"),
         lty = 1, lwd = 3, col = unlist(col_use), cex = 1.2, bty = "n")
  
  for (ltla_curr in ltla_plot) {
    d <- ltla_df %>%
      filter(ltla == ltla_curr)
    
    plot_map <- data.frame(date = seq(min(as.Date(d$mid_week)) - 3, 
                                      max(as.Date(d$mid_week)) + 3, 
                                      by = 1))
    plot_map$day_index <- 1:nrow(plot_map)
    plot_map$date2 <- format(as.Date(plot_map$date), format = "%d/%m/%Y")
    
    xpl <- plot_map$day_index
    
    what_pl <- "R"
    ylim_use <- switch(what_pl, I = c(0, max_prev2), R = c(0, 2.5))
    matplot(x = 0, y = 0, ty = "n", ylim = ylim_use, xlim = range(xpl), 
            xaxt = "n", yaxt = "s", xlab = "", las = 2, ylab = "", xaxs = "i", 
            main = substr(ltla_curr, 1, 20))
    abline(h = 1)

    xpl_week <- plot_map[match(d$mid_week, plot_map$date), "day_index"]
    if(what_pl == "R") {
      SIR_model_results[[ltla_curr]]$R_quant <- t(apply(SIR_model_results[[ltla_curr]]$R[, its_keep], 1, 
                                                        function(v) quantile(v, qplot, na.rm = T)))
      
      matc <- SIR_model_results[[ltla_curr]]$R_quant
    }
    colc <- 1
    lines(xpl_week, matc[, 2], lty = 1, lwd = 4, col = col_use$debias)
    lines(xpl_week, matc[, 3], lty = 1, lwd = 1, col = col_use$debias)
    lines(xpl_week, matc[, 1], lty = 1, lwd = 1, col = col_use$debias)

    imperial_R <- rd %>%
      filter(area == ltla_curr & date %in% as.character(plot_map$date)) %>%
      select(CIlow, Rt, CIup)
    lines(plot_map$day_index, imperial_R[, 2], lty = 1, lwd = 4, col = col_use$imperial)
    lines(plot_map$day_index, imperial_R[, 3], lty = 1, lwd = 1, col = col_use$imperial)
    lines(plot_map$day_index, imperial_R[, 1], lty = 1, lwd = 1, col = col_use$imperial)
    annotate_months(plot_map, add_axis = T, shade_alpha = .2, for_presentation = T)
  }

  mtext(side = 2, outer = T, at = .25, text = "Effective R value")
  
  dev.off()

}

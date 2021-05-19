source("scripts/SIR_utils.R")

############################
# Import weekly Pillar 2 and REACT
############################

ltla_unique <- unique(ltla_df$ltla)
mid_week_unique <- unique(ltla_df$mid_week)
n_weeks <- length(mid_week_unique)
all_mid_week_cut <- mid_week_unique
date_recent <- max(mid_week_unique)

param_df <- expand.grid(delta_AR_rho = c(0.975, 0.99, 0.999),
                        delta_AR_sd = c(1, 2),
                        R_AR_sd = c(0.2),
                        type = c("PCR_positive", "Infectious"))

for (i in 1:nrow(param_df)) {
  
  id <- paste0("AR", param_df$delta_AR_rho[i],
               "sd", param_df$delta_AR_sd[i],
               "Rsd", param_df$R_AR_sd[i])
  
  type <- param_df$type[i]
  out_dir <- file.path("output", id, type, "SIR")
  plot_dir <- file.path("plots", id, type)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  out_files <- list.files(out_dir, full.names = TRUE)
  SIR_model_results <- lapply(out_files, readRDS)
  names(SIR_model_results) <- sub(".RDS", "", basename(out_files))
  
  R_all <- I_all <- data.frame()
  Rl <- Il <- list()
  
  quant_plot <- c(0.025, 0.5, 0.975)
  for (ltla_curr in ltla_unique) {
    
    this_M <- ltla_pop %>%
      filter(ltla == ltla_curr) %>%
      pull(M)
    
    saml_biased <- SIR_model_results[[ltla_curr]]
    
    its_keep <- (control_SIR$burn_in + 1):control_SIR$n_iters
    
    I_quant_curr <- t(apply(saml_biased$I[, its_keep], 1, function(v) quantile(v, quant_plot, na.rm = T))) / this_M * 100
    SIR_model_results[[ltla_curr]]$I_quant <- I_quant_curr
    SIR_model_results[[ltla_curr]]$R_quant <- t(apply(saml_biased$R[, its_keep], 1, 
                                                      function(v) quantile(v, quant_plot, na.rm = T)))
    SIR_model_results[[ltla_curr]]$R_quant[n_weeks, ]
    
    R_add <- as.data.frame(t(SIR_model_results[[ltla_curr]]$R_quant[n_weeks, ]))
    I_add <- as.data.frame(t(SIR_model_results[[ltla_curr]]$I_quant[n_weeks, ]))
    names(I_add) <- names(R_add) <- c("l", "m", "u")
    R_add$ltla <- I_add$ltla <- ltla_curr
    R_all <- rbind(R_all, R_add)
    I_all <- rbind(I_all, I_add)
  }
  rownames(R_all) <- R_all$ltla
  rownames(I_all) <- I_all$ltla
  
  
  ###########################################################
  # Prevalence R scatter Figure in paper
  ###########################################################
  pdf(paste0(plot_dir, "/prevalence_curves_new.pdf"), 10, 7)
  
  layout(mat = matrix(c(1, 1, 2, 3), 2, 2))
  par(mar = c(4, 4, 2, 2), oma = c(1, 1, 1, 1))
  lin_col <- gray(.8)

  max_R <- ceiling(max(R_all$u) * 10) / 10
  min_R <- 0
  point_col <- rgb(red = .5, green = .4, blue = .8, alpha = .75)
  R_all$norm_add <- runif(nrow(I_all), min = -.025, max = .025)
  R_all$m_jit <- R_all$m + R_all$norm_add
  
  plot(R_all$m_jit + R_all$norm_add, I_all[, "m"], xlab = "", ylab = "", 
       xlim = c(min_R, max_R), ylim = c(0, ceiling(max(I_all[, "m"]))),  ty = "n")
  mtext(side = 3, text = paste0("Prevalence vs R for week of ", date_recent), line = .5)
  abline(v = 1)
  
  for(i in 1:nrow(R_all)) {
    lines(x = rep(R_all$m_jit[i], 2), y = c(I_all$l[i], I_all$u[i]), col = lin_col) 
    lines(x = c(R_all$l[i], R_all$u[i]), y = rep(I_all$m[i], 2), col = lin_col) 
  }
  
  points(R_all$m_jit, I_all[, "m"], col = point_col, pch = 19, cex = .7)
  mtext(side = 3, at = min_R, text = "(a)", cex = 1, line = .5)
  
  R_label <- R_all[order(R_all$m_jit, decreasing = T)[1:5], "ltla"]
  I_label <- I_all[order(I_all$m, decreasing = T)[1:5], "ltla"]
  #dput(c(R_label, I_label))
  ltla_curve <- c("Bolsover", "West Devon", "Bassetlaw", "North Hertfordshire", "Plymouth", 
                  "Newham", "Ealing", "Hounslow", "Brent", "Slough", "Leicester")[c(6, 7, 1, 2, 11)]
  col_reg <- rainbow(length(ltla_curve), alpha = .75)
  names(col_reg) <- ltla_curve
  for(ltla_curr in ltla_curve) {
    text(x = R_all[ltla_curr, "m_jit"] - .01, y = I_all[ltla_curr, "m"] + .04, 
         labels = ltla_curr, pos = 4, cex = 1, col = col_reg[ltla_curr])
  }
  mtext(side = 2, text = "LTLA prevalence (%)", line = 3)
  mtext(side = 1, text = "Effective R value", line = 3)
  
  plot_map <- data.frame(date = seq(min(as.Date(ltla_df$mid_week)) - 3, 
                                    max(as.Date(ltla_df$mid_week)) + 3, 
                                    by = 1))
  plot_map$day_index <- 1:nrow(plot_map)
  xpl <- plot_map$day_index
  
  max_prev_I <- ceiling(max(sapply(SIR_model_results[ltla_curve], 
                                function(x) max(x$I_quant[ ,3]))))
  
  panel_count <- 1
  for (what_pl in c("I", "R")) {
    panel_count <- panel_count + 1
    ylim_use <- switch(what_pl, I = c(0, max_prev_I), R = c(0, 2.5))
    matplot(x = 0, y = 0, ty = "n", ylim = ylim_use, xlim = range(xpl), 
            xaxt = "n", yaxt = "s", xlab = "", las = 2, ylab = "", xaxs = "i")
    if(what_pl == "R")
      abline(h = 1)
    
    for (ltla_curr in ltla_curve) {
      d <- ltla_df %>% 
        filter(ltla == ltla_curr)
      xpl_week <- plot_map[match(d$mid_week, plot_map$date), "day_index"]
      if(what_pl == "I") {
        matc <- SIR_model_results[[ltla_curr]]$I_quant
      }
      if(what_pl == "R") {
        matc <- SIR_model_results[[ltla_curr]]$R_quant
      }
      colc <- col_reg[ltla_curr]
      lines(xpl_week, matc[, 2], lty = 1, lwd = 4, col = colc)
      lines(xpl_week, matc[, 3], lty = 1, lwd = 1, col = colc)
      lines(xpl_week, matc[, 1], lty = 1, lwd = 1, col = colc)
    }
    abline(v = match(date_recent, plot_map$date))
    axis(side = 3, at = match(date_recent, plot_map$date), labels = date_recent)
    mtext(side = 3, text = switch(what_pl, I = "Prevalence (%)", R = "Effective R value"), line = .5)
    
    annotate_months(plot_map, add_axis = T, shade_alpha = .2, for_presentation = T)
    if(what_pl == "I") {
      legend(x = "topleft", legend = ltla_curve, col = col_reg[ltla_curve], lty = 1, lwd = 3, cex = .8, bg = "white")
    }
    mtext(side = 3, at = -.5, text = paste0("(", letters[panel_count], ")"), cex = 1, line = .5)
  }
  
  dev.off()
  
}

############################
# Import weekly Pillar 2 and REACT
############################

source("scripts/plot_utils.R")
source("scripts/SIR_utils.R")

ltla_unique <- unique(ltla_df$ltla)
mid_week_unique <- unique(ltla_df$mid_week)
n_weeks <- length(mid_week_unique)
all_mid_week_cut <- mid_week_unique
date_recent <- max(mid_week_unique)

### DELTA HYPERPARAMETERS ###

param_df <- expand.grid(delta_AR_sd = c(1, 2),
                        delta_AR_rho = c(0.975, 0.99, 0.999))

sensitivity_plot_delta <- function(type, what_pl) {
  par(mfrow = c(3, 2))
  for (i in 1:nrow(param_df)) {
    
    id <- paste0("AR", param_df$delta_AR_rho[i],
                 "sd", param_df$delta_AR_sd[i],
                 "Rsd", 0.2)
    
    panel_label <- substitute(paste(psi, " = ", rho, ", ", sigma[epsilon], " = ", sig), 
                              list(rho = param_df$delta_AR_rho[i],
                                   sig = param_df$delta_AR_sd[i]))
    
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
      
      I_quant_curr <- t(apply(saml_biased$I[, its_keep], 1, 
                              function(v) quantile(v, quant_plot, na.rm = T))) / this_M * 100
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
    
    lin_col <- gray(.8)
    max_R <- ceiling(max(R_all$u) * 10) / 10
    min_R <- 0
    point_col <- rgb(red = .5, green = .4, blue = .8, alpha = .75)
    
    plot_map <- data.frame(date = seq(min(as.Date(ltla_df$mid_week)) - 3, 
                                      max(as.Date(ltla_df$mid_week)) + 3, 
                                      by = 1))
    plot_map$day_index <- 1:nrow(plot_map)
    xpl <- plot_map$day_index
    ltla_curve <- c("Bolsover", "West Devon", "Bassetlaw", "North Hertfordshire", "Plymouth", 
                    "Newham", "Ealing", "Hounslow", "Brent", "Slough", "Leicester")[c(6, 7, 1, 2, 11)]
    
    max_prev_I <- ceiling(max(sapply(SIR_model_results[ltla_curve], 
                                     function(x) max(x$I_quant[ ,3]))))
    
    col_reg <- rainbow(length(ltla_curve), alpha = .75)
    names(col_reg) <- ltla_curve
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

    mtext(side = 3, text = panel_label, 
          line = .5)
    
    annotate_months(plot_map, add_axis = T, shade_alpha = .2, for_presentation = T)
    legend(x = "topleft", legend = ltla_curve, col = col_reg[ltla_curve], lty = 1, lwd = 3, cex = .8, bg = "white")
    
  }
}

for (type in c("PCR_positive", "Infectious")) {
  for (whatpl in c("I", "R")) {
    pdf(paste0("plots/SI_delta_", type, "_", whatpl, ".pdf"), 10, 7)
    sensitivity_plot_delta(type, whatpl)
    dev.off()
    
    file.copy(paste0("plots/SI_delta_", type, "_", whatpl, ".pdf"),
              paste0(overleaf_dir, "/SI_delta_", type, "_", whatpl, ".pdf"),
              overwrite = TRUE)
  }
}

### SENS/SPEC HYPERPARAMETERS ###

param_df <- expand.grid(alpha = c(0.001, 0.003),
                        beta = c(0.05, 0.1, 0.3))

sensitivity_plot_sensspec <- function(type, what_pl) {
  par(mfrow = c(3, 2))
  for (i in 1:nrow(param_df)) {
    
    id <- paste0("alpha", param_df$alpha[i],
                 "beta", param_df$beta[i])
    
    panel_label <- substitute(paste(alpha, " = ", alp, ", ", beta, " = ", bet), 
                              list(alp = param_df$alpha[i],
                                   bet = param_df$beta[i]))
    
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
      
      I_quant_curr <- t(apply(saml_biased$I[, its_keep], 1, 
                              function(v) quantile(v, quant_plot, na.rm = T))) / this_M * 100
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
    
    lin_col <- gray(.8)
    max_R <- ceiling(max(R_all$u) * 10) / 10
    min_R <- 0
    point_col <- rgb(red = .5, green = .4, blue = .8, alpha = .75)
    
    plot_map <- data.frame(date = seq(min(as.Date(ltla_df$mid_week)) - 3, 
                                      max(as.Date(ltla_df$mid_week)) + 3, 
                                      by = 1))
    plot_map$day_index <- 1:nrow(plot_map)
    xpl <- plot_map$day_index
    ltla_curve <- c("Bolsover", "West Devon", "Bassetlaw", "North Hertfordshire", "Plymouth", 
                    "Newham", "Ealing", "Hounslow", "Brent", "Slough", "Leicester")[c(6, 7, 1, 2, 11)]
    
    max_prev_I <- ceiling(max(sapply(SIR_model_results[ltla_curve], 
                                     function(x) max(x$I_quant[ ,3]))))
    
    col_reg <- rainbow(length(ltla_curve), alpha = .75)
    names(col_reg) <- ltla_curve
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

    mtext(side = 3, text = panel_label, 
          line = .5)
    
    annotate_months(plot_map, add_axis = T, shade_alpha = .2, for_presentation = T)
    legend(x = "topleft", legend = ltla_curve, col = col_reg[ltla_curve], lty = 1, lwd = 3, cex = .8, bg = "white")
    
  }
}

for (type in c("PCR_positive", "Infectious")) {
  for (whatpl in c("I", "R")) {
    pdf(paste0("plots/SI_sensspec_", type, "_", whatpl, ".pdf"), 10, 7)
    sensitivity_plot_sensspec(type, whatpl)
    dev.off()
  }
}

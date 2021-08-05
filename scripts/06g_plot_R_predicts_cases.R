
source("scripts/plot_utils.R")
source("scripts/SIR_utils.R")

unique_ltla <- unique(ltla_df$ltla)
mid_week_unique <- sort(unique(ltla_df$mid_week))
cut_dates <- seq(as.Date("2020-10-18"), rev(mid_week_unique)[2], by = 7)

################################
# collect results from different cut dates
################################

id <- "AR0.99sd1Rsd0.2"
type <- "Infectious"
  
out_dir <- file.path("output", id, type)
plot_dir <- file.path("plots", id, type)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

ir <- data.frame()
for (mid_week_int in cut_dates) {
  mid_week_cut <- as.Date(mid_week_int, origin = "1970-01-01")
  results_file <- paste0(out_dir, "/cut_", mid_week_cut, ".RDS")
  SIR_model_results <- readRDS(results_file)
  
  for (ltla_curr in unique_ltla) {
    d_ltla <- ltla_df %>%
      mutate(nt_per = nt / (M / 1e5),
             Nt_per = Nt / (M / 1e5),
             nt_plus_1 = c(diff(log2(nt)), NA),
             nt_plus_2 = c(diff(log2(nt), lag = 2), NA, NA)) %>%
      filter(ltla == ltla_curr &
               mid_week <= mid_week_cut)
    
    saml_biased <- SIR_model_results[[ltla_curr]]
    its_keep <- (control_SIR$burn_in + 1):control_SIR$n_iters
    SIR_model_results[[ltla_curr]]$I_quant <- 
      t(apply(saml_biased$I[, its_keep], 1, 
              function(v) quantile(v, qplot, na.rm = T))) / d_ltla$M[1] * 100
    SIR_model_results[[ltla_curr]]$R_quant <- t(apply(saml_biased$R[, its_keep], 1, 
                                                      function(v) quantile(v, qplot, na.rm = T)))
    SIR_model_results[[ltla_curr]]$R_quant
    R_add <- as.data.frame(SIR_model_results[[ltla_curr]]$R_quant)
    I_add <- as.data.frame(SIR_model_results[[ltla_curr]]$I_quant)
    names(I_add) <- paste0("I_", c("l", "m", "u"))
    names(R_add) <- paste0("R_", c("l", "m", "u"))
    add <- d_ltla %>%
      mutate(mid_week_cut = mid_week_cut) %>%
      bind_cols(I_add, R_add)
    
    ir <- bind_rows(ir, add)
  }
}


###########################################################
#  R prediction Figure in paper
###########################################################  
graphics.off()
pdf(file.path(plot_dir, "validation_2_for_paper.pdf"), 11, 7)

ylims <- list(list(c(-3.2, 3.2), c(-2, 2), c(-1.5, 1.5)), list(c(-6, 6), c(-2.5, 2.5), c(-2, 2)))
case_ranges <- list(c(0, 201), c(201, 501), c(501, Inf))
n_ran <- length(case_ranges)
par(mfrow = c(2, n_ran), mar = c(3, 3, 2, 2), oma = c(4, 4, 7, 12))

case_labs <- sapply(case_ranges, function(x) paste(x[1], x[2] - 1, sep = " - "))
case_labs <- gsub("501 - Inf", "> 500", case_labs)
cor_res <- list()

for (n_week_ahead in 1:2) {
  
  cor_res[[n_week_ahead]] <- list()
  
  for(j in 1:length(case_ranges)) {
    
    ir_keep <- ir %>%
      filter(mid_week == mid_week_cut & 
               nt_per >= case_ranges[[j]][1] & 
               nt_per < case_ranges[[j]][2])
    xpl <- ir_keep$R_m
    ypl <- ir_keep[[paste0("nt_plus_", n_week_ahead)]]

    cor_out <- cor.test(xpl, ypl, use = "p", me = "sp")
    cor_res[[n_week_ahead]][[case_labs[j]]] <- cor_out
    plot(xpl, ypl, main = "", xlab = "", ylab = "", xlim = c(0, 3.1), 
         ylim = ylims[[n_week_ahead]][[j]], las = 2, ty = "n")
    if (j == n_ran) {
      axis(side = 4, at = c(-1, 1), labels = c("Halving", "Doubling"), las = 2, col = 1)
    }
    abline(h = 0, v = 1)
    abline(h = -1, col = "green", lty = 3)
    abline(h = 1, col = "red", lty = 3)
    points(xpl, ypl, cex = .8)
    if (n_week_ahead == 1) {
      mtext(side = 3, line = 3.5, text = case_labs[j], col = "blue")
    }
    mtext(side = 3, line = .5, text = paste("rho = ", round(cor_out$estimate, 2)), cex = .8)
  }
}
mtext(outer = T, side = 3, line = 3.5, text = "Weekly case numbers per 100,000 at baseline", col = "blue")
mtext(outer = T, side = 4, at = c(.75, .25), text = c("One week ahead", "Two weeks ahead"), las = 2)
mtext(outer = T, side = 4, at = c(.75, .25) - .03, text = rep("change in # cases", 2), las = 2)
mtext(outer = T, side = 1, line = 1, text = "Estimated R at baseline")
mtext(outer = T, side = 2, line = 1, text = "Future log2 change in case numbers")

dev.off()

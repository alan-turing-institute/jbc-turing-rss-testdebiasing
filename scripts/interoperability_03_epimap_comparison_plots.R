library(dplyr)
source("scripts/plot_utils.R")
source("scripts/SIR_utils.R")

# Run 00_download_data.R and 01_preprocess_data.R first
ltla_df <- readr::read_csv("data/ltla.csv")

# Link to file in "Interoperability of models" Overleaf, file "shared_data/IR_for_interop.csv"
IR <- read.csv("https://www.dropbox.com/s/43kjmxq9kqunp22/IR_for_interop.csv?dl=1", stringsAsFactors = F)

# Link to file in "Interoperability of models" Overleaf, file "shared_data/Rt_estimates_Epimap_combined.csv"
rd <- read.csv("https://www.dropbox.com/s/txo5s4wllnzht2o/Rt_estimates_Epimap_combined.csv?dl=1") %>%
  mutate(date = as.Date(Date)) %>%
  filter(date <= max(ltla_df$mid_week))

epimap_min_date <- min(rd$date)
epimap_max_date <- max(rd$date)

mid_week_unique <- unique(ltla_df$mid_week)
week_ind <- which(mid_week_unique >= epimap_min_date)
ltla_plot <- sort(choose_focus_ltlas(ltla_df))

# graphics.off()
# pdf(file.path(plot_dir, "R_comparison_epimap.pdf"), 11, 7)
par(mfrow = c(2, 5), mar = c(3, 3, 2, 1), oma = c(2, 4, 2, 2))
pch_ests <- 19
cex_ests <- 1
col_curr <- list(corrected = rgb(0, 0, 1, alpha = .5), uncorrected = rgb(1, 0, 0, alpha = .5))
plot(0, xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
col_use <- list(debias = rgb(0, 0, 1, alpha = .5), epimap = rgb(1, 0, 0, alpha = .5))
legend(x = "center", legend = c("De-biasing model", "Epimap model"),
       lty = 1, lwd = 3, col = unlist(col_use), cex = 1.2, bty = "n")

for (ltla_curr in ltla_plot) {
  d <- ltla_df %>%
    filter(ltla == ltla_curr & mid_week >= epimap_min_date & mid_week <= epimap_max_date)
  
  plot_map <- data.frame(date = seq(min(as.Date(d$mid_week)), 
                                    max(as.Date(d$mid_week)), 
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
  if (what_pl == "R") {
    which(IR$ltla == ltla_curr)
    which(IR$mid_week %in% d$mid_week)
    str(d)
    matc <- IR[IR$ltla == ltla_curr & IR$mid_week %in% as.character(d$mid_week), c("R_l", "R_m", "R_u")]
  }
  colc <- 1
  lines(xpl_week, matc[, 2], lty = 1, lwd = 4, col = col_use$debias)
  lines(xpl_week, matc[, 3], lty = 1, lwd = 1, col = col_use$debias)
  lines(xpl_week, matc[, 1], lty = 1, lwd = 1, col = col_use$debias)

  epimap_R <- rd %>%
    filter(area == ltla_curr & date %in% plot_map$date) %>%
    select(Rt_2_5, Rt_50, Rt_97_5)
  lines(plot_map$day_index, epimap_R[, 2], lty = 1, lwd = 4, col = col_use$epimap)
  lines(plot_map$day_index, epimap_R[, 3], lty = 1, lwd = 1, col = col_use$epimap)
  lines(plot_map$day_index, epimap_R[, 1], lty = 1, lwd = 1, col = col_use$epimap)
  annotate_months(plot_map, add_axis = T, shade_alpha = .2, for_presentation = T)
}
mtext(side = 2, outer = T, at = .25, text = "Effective R value")













# dev.off()

# id <- "AR0.99sd1Rsd0.2"
# imperfect <- F
# type <- c("Infectious", "PCR_positive")[2]
# type_in_file_path <- paste0(type, "_", ifelse(imperfect, "Imperfect", "Perfect"))
# 
# out_dir <- file.path("output", id, type_in_file_path, "SIR")
# plot_dir <- file.path("plots", id, type_in_file_path)
# dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
# output_plot_dir <- c("~/Downloads", plot_dir)[2]

# out_dir <- file.path("output", id, type)
# plot_dir <- file.path("plots", id, type)

# out_files <- list.files(file.path(out_dir), full.names = TRUE)
# ltla_unique <- sub(".RDS", "", basename(out_files))
# SIR_model_results <- lapply(out_files, readRDS)
# names(SIR_model_results) <- sub(".RDS", "", basename(out_files))


library(rgdal)
library(dplyr)
library(ggplot2)
library(patchwork)
source("scripts/plot_utils.R")
source("scripts/SIR_utils.R")
options(bitmapType = "cairo-png")

############################
# Import weekly Pillar 2 and REACT
##########################
ltla_df <- readr::read_csv("data/ltla.csv")
ltla_unique <- unique(ltla_df$ltla)
mid_week_unique <- unique(ltla_df$mid_week)
n_weeks <- length(mid_week_unique)
all_mid_week_cut <- mid_week_unique
date_recent <- max(mid_week_unique)

param_df <- expand.grid(delta_AR_rho = c(0.975, 0.99, 0.999),
                        delta_AR_sd = c(1, 2),
                        R_AR_sd = c(0.2),
                        type = c("PCR_positive", "Infectious"))
LTLA_shp_Reg <- get_ltla_shape_file()

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
  
  IR <- data.frame()
  R_all <- I_all <- data.frame()
  Rl <- Il <- list()
  
  for (ltla_curr in ltla_unique) {
    
    this_M <- ltla_pop %>%
      filter(ltla == ltla_curr) %>%
      pull(M)
    
    saml_biased <- SIR_model_results[[ltla_curr]]
    its_keep <- (control_SIR$burn_in + 1):control_SIR$n_iters
    I_quant_curr <- t(apply(saml_biased$I[, its_keep], 1, 
                            function(v) quantile(v, qplot, na.rm = T))) / this_M * 100
    SIR_model_results[[ltla_curr]]$I_quant <- I_quant_curr
    SIR_model_results[[ltla_curr]]$R_quant <- t(apply(saml_biased$R[, its_keep], 1, 
                                                      function(v) quantile(v, qplot, na.rm = T)))
    SIR_model_results[[ltla_curr]]$R_quant[n_weeks, ]
    R_add <- as.data.frame(SIR_model_results[[ltla_curr]]$R_quant)
    I_add <- as.data.frame(SIR_model_results[[ltla_curr]]$I_quant)
    rownames(R_add)
    names(R_add) <- paste0("R_", c("l", "m", "u"))
    names(I_add) <- paste0("I_", c("l", "m", "u"))
    add_all <- cbind(I_add, R_add)
    add_all$namdat <- rownames(add_all)
    add_all$ltla <- ltla_curr
    #add_all$mid_week <- sapply(strsplit(add_all$namdat, "_"), function(x) x[2])
    add_all$mid_week <- mid_week_unique
    IR <- rbind(IR, add_all)
  }
  
  n_week_plot <- 20
  week_plot <- rev(rev(mid_week_unique)[1 + 1 * (0:(n_week_plot - 1))])
  
  xmin <- 82672
  xmax <- 655653.8
  ymin <- 5342.7
  ymax <- 657536
  xoffset <- (1 / 5) * (xmax - xmin)
  yoffset <- (1 / 5) * (ymax - ymin)
  
  ### Prevalence maps ###
  
  plot_max <- ceiling(max(IR$I_m) * 2) / 2
  week_plot_curr = week_plot[1]
  
  IR_week <- IR[IR$mid_week == week_plot_curr, ]
  curr_y <- IR_week[match(LTLA_shp_Reg$LAD20NM, IR_week$ltla),"I_m"]
  curr_y[is.na(curr_y)] <- mean(curr_y, na.rm = T)
  
  LTLA_shp_Reg$lf_mean <- curr_y
  
  y_lnd <- ggplot(data = filter(LTLA_shp_Reg, RGN20NM == "London")) + 
    geom_sf(aes(fill = lf_mean ), color = NA, show.legend = FALSE) + 
    theme_void() + 
    ggtitle("London") + 
    scale_fill_viridis_c(option = "A", limits = c(0, plot_max), direction = -1) +
    theme(plot.title = element_text(size = 7)) 
  
  
  y_ne <- ggplot(data = filter(LTLA_shp_Reg, RGN20NM == "North East")) + 
    geom_sf(aes(fill = lf_mean ), color = NA, show.legend = FALSE) + 
    scale_fill_viridis_c(option = "A", limits = c(0, plot_max), direction = -1) +
    theme_void() +
    ggtitle("North East") + 
    theme(plot.title = element_text(size = 7)) 
  lnd_grob <- ggplotGrob(y_lnd)
  ne_grob <- ggplotGrob(y_ne)
  
  ltla_plot <- ggplot(data = LTLA_shp_Reg) + 
    geom_sf(aes(fill = lf_mean), colour=NA) +
    theme_void() +
    theme(plot.title = element_text(size = 7)) + 
    ggtitle(as.Date(week_plot_curr, origin = "1970-01-01")) +
    annotation_custom(grob = lnd_grob, xmin = xmin, xmax = xmin+0.35*(xmax-xmin), 
                      ymin = ymin+1.2*yoffset, ymax =ymin + 1.2*yoffset + 0.35*(ymax-ymin))# +
  #  annotation_custom(grob = ne_grob, xmin = xmin, xmax = xmin+0.35*(xmax-xmin), 
  #                    ymin = ymin+3.5*yoffset, ymax =ymin + 3.5*yoffset + 0.35*(ymax-ymin))
  
  for (week_plot_curr in week_plot[-1]) {
    IR_week <- IR[IR$mid_week == week_plot_curr, ]
    curr_y <- IR_week[match(LTLA_shp_Reg$LAD20NM, IR_week$ltla),"I_m"]
    curr_y[is.na(curr_y)] <- mean(curr_y, na.rm = T)
    LTLA_shp_Reg$lf_mean <- curr_y
    
    y_lnd <- ggplot(data = filter(LTLA_shp_Reg, RGN20NM == "London")) + 
      geom_sf(aes(fill = lf_mean ), color = NA, show.legend = FALSE) + 
      theme_void() +
      ggtitle("London") + 
      scale_fill_viridis_c(option = "A", limits = c(0, plot_max), direction = -1) +
      theme(plot.title = element_text(size = 7)) 
    
    y_ne <- ggplot(data = filter(LTLA_shp_Reg, RGN20NM == "North East")) + 
      geom_sf(aes(fill = lf_mean ), color = NA, show.legend = FALSE)  +
      theme_void() +
      ggtitle("North East") + 
      scale_fill_viridis_c(option = "A", limits = c(0, plot_max), direction = -1) +
      theme(plot.title = element_text(size = 7))
    
    lnd_grob <- ggplotGrob(y_lnd)
    ne_grob <- ggplotGrob(y_ne)
    
    ltla_plot <- ltla_plot + 
      ggplot(data = LTLA_shp_Reg) + 
      geom_sf(aes(fill = lf_mean), colour=NA, show.legend = FALSE) +
      theme_void() +
      theme(plot.title = element_text(size = 7)) + 
      ggtitle(as.Date(week_plot_curr, origin = "1970-01-01")) +
      annotation_custom(grob = lnd_grob, xmin = xmin, 
                        xmax = xmin+0.35*(xmax-xmin), 
                        ymin = ymin+1.2*yoffset, 
                        ymax= ymin + 1.2*yoffset + 0.35*(ymax-ymin))# +  
    #   annotation_custom(grob = ne_grob, xmin = xmin, xmax = xmin+0.35*(xmax-xmin), 
    #                     ymin = ymin+3.5*yoffset, ymax =ymin + 3.5*yoffset + 0.35*(ymax-ymin)) 
  }
  
  p_out <- ltla_plot +
    plot_layout(guides = "collect") & 
    scale_fill_viridis_c(option = "A", limits = c(0, plot_max), direction = -1) & 
    labs(fill="Prevalence (%)")
  
  ggsave(file.path(plot_dir, "prev_maps.png"), p_out, width = 12, height = 8, units = "in")
  
  
  ### Effective R maps ###
  max_R <- 2
  week_plot_curr = week_plot[1]
  
  IR_week <- IR[IR$mid_week == week_plot_curr, ]
  curr_y <- IR_week[match(LTLA_shp_Reg$LAD20NM, IR_week$ltla), "R_m"]
  curr_y[curr_y > max_R] <- max_R
  curr_y[is.na(curr_y)] <- mean(curr_y, na.rm = T)
  
  LTLA_shp_Reg$lf_mean <- curr_y
  
  y_lnd <- ggplot(data = filter(LTLA_shp_Reg, RGN20NM == "London")) + 
    geom_sf(aes(fill = lf_mean ), color = NA, show.legend = FALSE) + 
    theme_void() + 
    ggtitle("London") + 
    scale_fill_distiller(palette = "RdYlGn", limits =c(0, 2), 
                         breaks = seq(0, 2, 0.5), 
                         labels = c(seq(0, 1.5, 0.5), "2.0+")) +
    theme(plot.title = element_text(size = 7)) 
  
  
  y_ne <- ggplot(data = filter(LTLA_shp_Reg, RGN20NM == "North East")) + 
    geom_sf(aes(fill = lf_mean ), color = NA, show.legend = FALSE)  + 
    scale_fill_distiller(palette = "RdYlGn", limits =c(0, 2), 
                         breaks = seq(0, 2, 0.5), 
                         labels = c(seq(0, 1.5, 0.5), "2.0+")) +
    theme_void() +
    ggtitle("North East")+ 
    theme(plot.title = element_text(size = 7)) 
  lnd_grob <- ggplotGrob(y_lnd)
  ne_grob <- ggplotGrob(y_ne)
  
  ltla_plot <- ggplot(data = LTLA_shp_Reg) + 
    geom_sf(aes(fill = lf_mean), colour=NA) +
    theme_void() +
    theme(plot.title = element_text(size = 7)) + 
    ggtitle(as.Date(week_plot_curr, origin = "1970-01-01")) +
    annotation_custom(grob = lnd_grob, xmin = xmin, xmax = xmin+0.35*(xmax-xmin), 
                      ymin = ymin+1.2*yoffset, ymax =ymin + 1.2*yoffset + 0.35*(ymax-ymin))# +
  #  annotation_custom(grob = ne_grob, xmin = xmin, xmax = xmin+0.35*(xmax-xmin), 
  #                    ymin = ymin+3.5*yoffset, ymax =ymin + 3.5*yoffset + 0.35*(ymax-ymin))
  
  
  
  for (week_plot_curr in week_plot[-1]) {
    IR_week <- IR[IR$mid_week == week_plot_curr, ]
    curr_y <- IR_week[match(LTLA_shp_Reg$LAD20NM, IR_week$ltla), "R_m"]
    curr_y[curr_y > max_R] <- max_R
    curr_y[is.na(curr_y)] <- mean(curr_y, na.rm = T)
    LTLA_shp_Reg$lf_mean <- curr_y
    
    y_lnd <- ggplot(data = filter(LTLA_shp_Reg, RGN20NM == "London")) + 
      geom_sf(aes(fill = lf_mean ), color = NA, show.legend = FALSE) + 
      theme_void() +
      ggtitle("London") + 
      scale_fill_distiller(palette = "RdYlGn", limits =c(0, 2), 
                           breaks = seq(0, 2, 0.5), 
                           labels = c(seq(0, 1.5, 0.5), "2.0+")) +
      theme(plot.title = element_text(size = 7)) 
    
    y_ne <- ggplot(data = filter(LTLA_shp_Reg, RGN20NM == "North East")) + 
      geom_sf(aes(fill = lf_mean ), color = NA, show.legend = FALSE)  +
      theme_void() +
      ggtitle("North East") + 
      scale_fill_distiller(palette = "RdYlGn", limits =c(0, 2), 
                           breaks = seq(0, 2, 0.5), 
                           labels = c(seq(0, 1.5, 0.5), "2.0+")) + 
      theme(plot.title = element_text(size = 7)) 
    
    lnd_grob <- ggplotGrob(y_lnd)
    ne_grob <- ggplotGrob(y_ne)
    
    ltla_plot <- ltla_plot + 
      ggplot(data = LTLA_shp_Reg) + 
      geom_sf(aes(fill = lf_mean), colour=NA, show.legend = FALSE) +
      theme_void() +
      theme(plot.title = element_text(size = 7)) + 
      ggtitle(as.Date(week_plot_curr, origin = "1970-01-01")) +
      annotation_custom(grob = lnd_grob, xmin = xmin, xmax = xmin+0.35*(xmax-xmin), 
                        ymin = ymin+1.2*yoffset, ymax =ymin + 1.2*yoffset + 0.35*(ymax-ymin)) #+  
    #    annotation_custom(grob = ne_grob, xmin = xmin, xmax = xmin+0.35*(xmax-xmin), 
    #                      ymin = ymin+3.5*yoffset, ymax =ymin + 3.5*yoffset + 0.35*(ymax-ymin)) 
  }
  
  p_out <- ltla_plot + 
    plot_layout(guides = "collect") & 
    scale_fill_distiller(palette = "RdYlGn", limits = c(0, 2), 
                         breaks = seq(0, 2, 0.5), 
                         labels = c(seq(0, 1.5, 0.5), "2.0+")) & 
    labs(fill="Rt")
  
  ggsave(file.path(plot_dir, "R_maps.png"), p_out, width = 12, height = 8, units = "in")
  
}

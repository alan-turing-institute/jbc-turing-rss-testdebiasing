
n_week_plot <- 20
week_plot <- rev(rev(mid_week_unique)[1 + 2 * (0:(n_week_plot - 1))])

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
  scale_fill_distiller(palette = "RdYlBu", limits =c(0, 2), 
                       breaks = seq(0, 2, 0.5), 
                       labels = c(seq(0, 1.5, 0.5), "2.0+")) +
  theme(plot.title = element_text(size = 7)) 


y_ne <- ggplot(data = filter(LTLA_shp_Reg, RGN20NM == "North East")) + 
  geom_sf(aes(fill = lf_mean ), color = NA, show.legend = FALSE)  + 
  scale_fill_distiller(palette = "RdYlBu", limits =c(0, 2), 
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
    scale_fill_distiller(palette = "RdYlBu", limits =c(0, 2), 
                         breaks = seq(0, 2, 0.5), 
                         labels = c(seq(0, 1.5, 0.5), "2.0+")) +
    theme(plot.title = element_text(size = 7)) 
  
  y_ne <- ggplot(data = filter(LTLA_shp_Reg, RGN20NM == "North East")) + 
    geom_sf(aes(fill = lf_mean ), color = NA, show.legend = FALSE)  +
    theme_void() +
    ggtitle("North East") + 
    scale_fill_distiller(palette = "RdYlBu", limits =c(0, 2), 
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
  scale_fill_distiller(palette = "RdYlBu", limits = c(0, 2), 
                       breaks = seq(0, 2, 0.5), 
                       labels = c(seq(0, 1.5, 0.5), "2.0+")) & 
  labs(fill="Rt")

ggsave(file.path(plot_dir, "R_maps.png"), p_out, width = 12, height = 8, units = "in")


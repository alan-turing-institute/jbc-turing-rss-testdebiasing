library(ggplot2)
library(patchwork)
library(rgdal)
library(sf)
library(xtable)
source("scripts/plot_utils.R")
source("scripts/SIR_utils.R")

### WRANGLE DATA ###
ltla_df <- readr::read_csv("data/ltla.csv")
ltla_unique <- unique(ltla_df$ltla)
source("scripts/prep_BAME-IMD.R")

mid_week_unique <- unique(ltla_df$mid_week)
this_week <- max(mid_week_unique)
n_weeks <- length(mid_week_unique)

id <- "AR0.99sd1Rsd0.2"

type <- "Infectious"
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
  
  R_add <- as.data.frame(SIR_model_results[[ltla_curr]]$R_quant)
  I_add <- as.data.frame(SIR_model_results[[ltla_curr]]$I_quant)
  names(I_add) <- names(R_add) <- c("l", "m", "u")
  R_add$ltla <- I_add$ltla <- ltla_curr
  R_add$mid_week <- I_add$mid_week <- mid_week_unique
  R_all <- rbind(R_all, R_add)
  I_all <- rbind(I_all, I_add)
}

joint_df <- ltla_df %>%
  left_join(I_all, by = c("ltla", "mid_week")) %>%
  left_join(cov_data, by = "ltla") %>%
  mutate(nt_per_100k = 1e5 * nt / M,
         m_per_100k = 1e3 * m,
         l_per_100k = 1e3 * l,
         u_per_100k = 1e3 * u,
         missing_per_100k = m_per_100k - nt_per_100k,
         u_missing_per_100k = u_per_100k - nt_per_100k,
         l_missing_per_100k = l_per_100k - nt_per_100k,
         missing_prop = 1 -  (nt_per_100k / m_per_100k),
         u_missing_prop = 1 - (nt_per_100k / l_per_100k),
         l_missing_prop = 1 - (nt_per_100k / u_per_100k))

last_week_df <- joint_df %>%
  filter(mid_week == this_week) 

### Scatter plot ###

p1 <- ggplot(last_week_df, aes(m_per_100k, missing_prop)) + 
  geom_point(aes(color = bame_quint)) + 
  ylab("Estimated proportion of cases missing") +
  xlab("Estimated prevalence (per 100k)") +
  theme_minimal() +
  viridis::scale_color_viridis(discrete = TRUE)

ggsave("~/Downloads/missing_scatter.png", p1)

### Tables ###

# Counts
nround_dig <- 0
miss_cases_ci <- paste0("(", round(last_week_df$l_missing_per_100k, nround_dig), 
                        " - ", round(last_week_df$u_missing_per_100k, nround_dig), ")")
d <- data.frame(LTLA = last_week_df$ltla, 'Cases per 100K' = round(last_week_df$nt_per_100k, nround_dig), 
                'Missing Cases per 100K' =  round(last_week_df$missing_per_100k, nround_dig), 
                '95% CI' = miss_cases_ci)
d <- d[order(d$Missing.Cases, decreasing = T), ]
tab_out <- xtable(d[1:10, ], digits = nround_dig)
names(tab_out) <- c('LTLA', 'Cases', 'Estimated Missing Cases', '(95% CI)')
print(tab_out,  file = "~/Downloads/missing_cases_table.txt", include.rownames = FALSE)

# Proportions

nround_dig <- 2
miss_cases_ci <- paste0("(", round(last_week_df$l_missing_prop, nround_dig), 
                        " - ", round(last_week_df$u_missing_prop, nround_dig), ")")
d <- data.frame(LTLA = last_week_df$ltla, 'Cases per 100K' = round(last_week_df$nt_per_100k, 0), 
                'Proportion of missing cases' =  round(last_week_df$missing_prop, nround_dig), 
                '95% CI' = miss_cases_ci)
d <- d[order(d$Proportion.of.missing.cases, decreasing = T), ]
tab_out <- xtable(d[1:10, ], digits = c(0, 0, 0, nround_dig, nround_dig))
names(tab_out) <- c('LTLA', 'Cases', 'Estimated proportion of cases missing', '(95% CI)')
print(tab_out,  file = "~/Downloads/proportion_missing_cases_table.txt", include.rownames = FALSE)

#write_csv(last_week_df, "~/Downloads/prev_estimates_2021-06-20")

### PLOT DATA ###
options(bitmapType = "cairo-png", device = "X11")
LTLA_shp_Reg <- get_ltla_shape_file()

# Merge Northamptonshires

st_union_by = function(geo, group) {
  
  y2 = list()
  #loop over by groups and merge units
  for (i in unique(group)) {
    #which units
    z = geo[group == i]
    #merge
    y = Reduce(sf::st_union, z)
    y2[[i]] = y
  }
  
  y3 <- sf::st_sfc(y2)
  #
  sf::st_sf(data.frame(lad20cd = names(y3), geom = sf::st_sfc(y3)))
}

LTLA2Reg <- readr::read_csv("data/Local_Authority_District_to_Region__December_2020__Lookup_in_England.csv")
west_northamptonshire <- c("Daventry", "Northampton", "South Northamptonshire")
north_northamptonshire <- c("Corby", "East Northamptonshire", 
                           "Kettering", "Wellingborough")
LTLA_shp_Reg$lad20cd[LTLA_shp_Reg$LAD20NM %in% north_northamptonshire] <- "E06000061"
LTLA_shp_Reg$lad20cd[LTLA_shp_Reg$LAD20NM %in% west_northamptonshire] <- "E06000062"
LTLA_shp_EN <- st_union_by(LTLA_shp_Reg$geometry, LTLA_shp_Reg$lad20cd) %>%
  left_join(LTLA2Reg, by = c("lad20cd" = "LAD20CD"))
LTLA_shp_EN$LAD20NM[LTLA_shp_EN$lad20cd == "E06000061"] <- "North Northamptonshire"
LTLA_shp_EN$LAD20NM[LTLA_shp_EN$lad20cd == "E06000062"] <- "West Northamptonshire"

bar_width <- 0.35 
xmin <- 82672
xmax <- 655653.8
ymin <- 5342.7
ymax <- 657536
xoffset <- (1 / 5) * (xmax - xmin)
yoffset <- (1 / 5) * (ymax - ymin)

# Positive tests
curr_y <- last_week_df$nt_per_100k[match(LTLA_shp_EN$LAD20NM, last_week_df$ltla)]
LTLA_shp_EN$nt_per_100k <- curr_y

var_plt <- ggplot(data = LTLA_shp_EN) + 
  geom_sf(aes(fill = nt_per_100k), color = NA) + 
  scale_fill_viridis_c(limits = c(0, max(curr_y)),
                       guide = guide_colourbar(title.position = "top"),
                       na.value = "gray90") +
  theme_void()

var_lnd <- ggplot(data = filter(LTLA_shp_EN, RGN20NM == "London")) + 
  geom_sf(aes(fill = nt_per_100k), color = NA, show.legend = FALSE) + 
  scale_fill_viridis_c(limits = c(0, max(curr_y)),
                       na.value = "gray90") +
  theme_void() + 
  ggtitle("London") + 
  theme(plot.title = element_text(size = 5))

nt_plt <- var_plt +
  annotation_custom(grob = ggplotGrob(var_lnd), xmin = xmin, xmax = xmin+0.35*(xmax-xmin), 
                    ymin = ymin+3.5*yoffset, ymax =ymin + 3.5*yoffset + 0.35*(ymax-ymin)) + 
  theme_void() +
  ggtitle("Positive cases\n(per 100k)") +
  theme(plot.title = element_text(size = 8), 
        legend.key.size = unit(0.2, 'cm'), #change legend key size
        legend.key.height = unit(0.2, 'cm'), #change legend key height
        legend.key.width = unit(bar_width, 'cm'), #change legend key width
        legend.title = element_text(size=7), #change legend title font size
        legend.text = element_text(size=5))

# Median prevalence
curr_y <- last_week_df$m_per_100k[match(LTLA_shp_EN$LAD20NM, last_week_df$ltla)]
LTLA_shp_EN$m_per_100k <- curr_y

var_plt <- ggplot(data = LTLA_shp_EN) + 
  geom_sf(aes(fill = m_per_100k), color = NA) + 
  scale_fill_viridis_c(limits = c(0, max(curr_y)),
                       guide = guide_colourbar(title.position = "top"),
                       na.value = "gray90") +
  theme_void()

var_lnd <- ggplot(data = filter(LTLA_shp_EN, RGN20NM == "London")) + 
  geom_sf(aes(fill = m_per_100k), color = NA, show.legend = FALSE) + 
  scale_fill_viridis_c(limits = c(0, max(curr_y)),
                       na.value = "gray90") +
  theme_void() + 
  ggtitle("London") + 
  theme(plot.title = element_text(size = 5))

prev_plt <- var_plt +
  annotation_custom(grob = ggplotGrob(var_lnd), xmin = xmin, xmax = xmin+0.35*(xmax-xmin), 
                    ymin = ymin+3.5*yoffset, ymax =ymin + 3.5*yoffset + 0.35*(ymax-ymin)) + 
  theme_void() +
  ggtitle("Debiased prevalence\n(per 100k)") +
  theme(plot.title = element_text(size = 8), 
        legend.key.size = unit(0.2, 'cm'), #change legend key size
        legend.key.height = unit(0.2, 'cm'), #change legend key height
        legend.key.width = unit(bar_width, 'cm'), #change legend key width
        legend.title = element_text(size=7), #change legend title font size
        legend.text = element_text(size=5))

# Prevalence - positive tests
curr_y <- last_week_df$missing_per_100k[match(LTLA_shp_EN$LAD20NM, last_week_df$ltla)]
LTLA_shp_EN$missing_per_100k <- curr_y

var_plt <- ggplot(data = LTLA_shp_EN) + 
  geom_sf(aes(fill = missing_per_100k), color = NA) + 
  scale_fill_viridis_c(limits = c(0, max(curr_y)),
                       guide = guide_colourbar(title.position = "top"),
                       na.value = "gray90") +
  theme_void()

var_lnd <- ggplot(data = filter(LTLA_shp_EN, RGN20NM == "London")) + 
  geom_sf(aes(fill = missing_per_100k), color = NA, show.legend = FALSE) + 
  scale_fill_viridis_c(limits = c(0, max(curr_y)),
                       na.value = "gray90") +
  theme_void() + 
  ggtitle("London") + 
  theme(plot.title = element_text(size = 5))

missing_plt <- var_plt +
  annotation_custom(grob = ggplotGrob(var_lnd), xmin = xmin, xmax = xmin+0.35*(xmax-xmin), 
                    ymin = ymin+3.5*yoffset, ymax =ymin + 3.5*yoffset + 0.35*(ymax-ymin)) + 
  theme_void() +
  ggtitle("'Missing' cases\n(per 100k)") +
  theme(plot.title = element_text(size = 8), 
        legend.key.size = unit(0.2, 'cm'), #change legend key size
        legend.key.height = unit(0.2, 'cm'), #change legend key height
        legend.key.width = unit(bar_width, 'cm'), #change legend key width
        legend.title = element_text(size=7), #change legend title font size
        legend.text = element_text(size=5))

# Positive / prevalence tests
curr_y <- last_week_df$missing_prop[match(LTLA_shp_EN$LAD20NM, last_week_df$ltla)]
LTLA_shp_EN$missing_prop <- curr_y

var_plt <- ggplot(data = LTLA_shp_EN) + 
  geom_sf(aes(fill = missing_prop), color = NA) + 
  scale_fill_viridis_c(limits = c(0, 1),
                       guide = guide_colourbar(title.position = "top"),
                       na.value = "gray90") +
  theme_void()

var_lnd <- ggplot(data = filter(LTLA_shp_EN, RGN20NM == "London")) + 
  geom_sf(aes(fill = missing_prop), color = NA, show.legend = FALSE) + 
  scale_fill_viridis_c(limits = c(0, 1),
                       na.value = "gray90") +
  theme_void() + 
  ggtitle("London") + 
  theme(plot.title = element_text(size = 5))

prop_plt <- var_plt +
  annotation_custom(grob = ggplotGrob(var_lnd), xmin = xmin, xmax = xmin+0.35*(xmax-xmin), 
                    ymin = ymin+3.5*yoffset, ymax =ymin + 3.5*yoffset + 0.35*(ymax-ymin)) + 
  theme_void() +
  ggtitle("Proportion of cases\nmissing") +
  theme(plot.title = element_text(size = 8), 
        legend.key.size = unit(0.2, 'cm'), #change legend key size
        legend.key.height = unit(0.2, 'cm'), #change legend key height
        legend.key.width = unit(bar_width, 'cm'), #change legend key width
        legend.title = element_text(size=7), #change legend title font size
        legend.text = element_text(size=5))

out_plt <- nt_plt + prev_plt + missing_plt + prop_plt + plot_layout(nrow = 1, ncol = 4) &
  theme(legend.position = "bottom")#+
#  plot_annotation(
#    title = paste0("Estimates of missing cases by LTLA for week of ", this_week)
#  )

ggsave("~/Downloads/combi_plt.png", out_plt, width = 5.5, height = 3, dpi = 500)

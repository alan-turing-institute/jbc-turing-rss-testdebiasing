
library(rgdal); library(sf); library(dplyr); library(ggplot2)
library(gridExtra); library(grid); library(lattice)
source("scripts/plot_utils.R")
shape_df <- get_ltla_shape_file() %>%
  rename(`PHE region` = RGN20NM)

plot_dir <- "plots"
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

p1 <- ggplot(shape_df) +
  geom_sf(aes(fill = `PHE region`), size = 0.1) +
  theme_void(base_size = 7, base_family="Helvetica")
ggsave("plots/phe_reg_ltla.pdf", p1, width = 9, height = 7, units = "cm")

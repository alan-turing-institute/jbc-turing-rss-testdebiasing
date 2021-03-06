# PHE regions
phe_region_unique <- c("London", "Yorkshire and The Humber", 
                       "South West", "East of England",
                       "West Midlands", "North West", "South East", 
                       "North East", "East Midlands")

# Plotting parameters
qplot <- c(0.025, .5, .975)


choose_focus_ltlas <- function(ltla_df) {
  
  ltla_phe_unique <- unique(ltla_df[, c("phe_region", "ltla")])
  ltla_phe_unique <- ltla_phe_unique[order(ltla_phe_unique$ltla), ]
  ltla_plot <- ltla_phe_unique[match(phe_region_unique, 
                                     ltla_phe_unique$phe_region) + 0,]
  
  pull(ltla_plot, ltla)[1:9]
}


transpar_const <- .85
qpl <- c(0.025, .5, .975)
col_data <- c()
col_data["smooth_rand"] <- 1
col_data["filtered_targ"] <- rgb(red = 1, green = 0, blue = .15, alpha = transpar_const)
col_data["filtered_rand"] <- rgb(red = 0, green = .7, blue = .7, alpha = transpar_const)
col_data["bayes_factor"] <- 1
col_data["react_samples"] <- rgb(red = 1, green = .5, blue = 0, alpha = transpar_const)


annotate_months <- function(d, add_axis = T, shade_alpha = .2, for_presentation = F){
  #d$month <- sapply(strsplit(d$date, "-"), function(v) paste(v[1:2], collapse = "-"))
  d$month <- months(d$date)
  d$month <- format(d$date, "%Y-%m")
  month_unique <- unique(d$month)  
  month_name <- format(strptime(paste0(month_unique, "-01"), format = "%Y-%m-%d"), "%b")
  names(month_name) <- month_unique
  greys_shade <- grey(c(0.75, 1), alpha = shade_alpha)
  month_shade_col <- greys_shade[1:length(month_unique) %% 2 + 1]
  names(month_shade_col) <- month_unique
  month_mean_index <- sapply(month_unique, function(m) mean(which(d$month == m)))
  month_lower_divider <- sapply(month_unique, function(m) min(which(d$month == m)) - .5)
  month_upper_divider <- sapply(month_unique, function(m) max(which(d$month == m)) + .5)
  for(month_shade in month_unique){
    x_coords <- c(rep(month_lower_divider[month_shade], 2), rep(month_upper_divider[month_shade], 2), month_lower_divider[month_shade])
    big_num <- 1e10
    y_coords <- big_num * c(-1, 1, 1, -1, -1)
    polygon(x = x_coords, y = y_coords, col = month_shade_col[month_shade], border = NA)  
  }
  if (add_axis) {
    labs <- month_name[month_unique]
    ats <- month_mean_index
 #   if(for_presentation) {
#      labs[labs == "May"] <- ""
#    }
    axis(side = 1, at = ats, labels = labs, las = 2, tcl = 0)
  }
}

add_react_sampling_intervals_to_plot <- function(d, col_data, prev_max, type = c("daily", "weekly")[1], prev_min = 0){
  
  # REACT DATES
  react.datel <- list()
  react.datel$rd1 <- c("2020-05-01", "2020-06-01")
  react.datel$rd2 <- c("2020-06-19", "2020-07-07")
  react.datel$rd3 <- c("2020-07-24", "2020-08-11")
  react.datel$rd4 <- c("2020-08-20", "2020-09-08")
  react.datel$rd5 <- c("2020-09-18", "2020-10-05")
  react.datel$rd6 <- c("2020-10-16", "2020-11-02")
  react.datel$rd7 <- c("2020-11-13", "2020-12-03")
  react.datel$rd7a <- c("2020-11-13", "2020-11-24")
  react.datel$rd7b <- c("2020-11-25", "2020-12-03")
  react.datel$rd8a <- c("2021-01-06", "2021-01-15")
  react.datel$rd8 <- c("2021-01-06", "2021-01-22")
  react.datel$rd9 <- c("2021-02-04", "2021-02-23")
  react.datel$rd10 <- c("2021-03-11", "2021-03-30")
  react.datel$rd11 <- c("2021-04-15", "2021-05-03")
  react.datel$rd12 <- c("2021-05-20", "2021-06-07")
  
  for(j in 1:length(react.datel)){
    if (type == "daily") {
      x_lims <- match(react.datel[[j]], as.character(d$date))
    }
    if (type == "weekly") {
      x_lims <- match(react.datel[[j]], d$mid_week)
    }
    if (is.na(x_lims[1])) x_lims[1] <- 0
    if (is.na(x_lims[2])) x_lims[2] <- nrow(d)
    polygon(x = x_lims[c(1, 2, 2, 1, 1)], 
            y = prev_min + c(0, 0, -1, -1, 0) * prev_max / 40, 
            col = col_data["react_samples"], 
            border = NA)
  }
}



get_ltla_shape_file <- function(merge_Hackney_CoL = TRUE) {
  
  ## Local Authority Districts (same as LTLA, checked) to Regions lookup:
  LTLA2Reg <- readr::read_csv("data/Local_Authority_District_to_Region__December_2020__Lookup_in_England.csv")
  
  # import LTLA shapefile
  LTLA_shp = rgdal::readOGR("data/Local_Authority_Districts_(May_2020)_Boundaries_UK_BFE-shp")
  LTLA_shp <- sf::st_as_sf(LTLA_shp)
  
  # remove LTLAs for which we don't have testing data (i.e. all that are not in England)
  LTLA_shp_EN = LTLA_shp %>% 
    dplyr::filter(lad20cd %in% LTLA2Reg$LAD20CD)
  
  # remove island
  islands = c("Isle of Wight",   "Isles of Scilly")
  LTLA_shp_EN = LTLA_shp_EN[!LTLA_shp_EN$lad20nm %in% islands,]

  if (merge_Hackney_CoL) {
    
    # assign to City of London the code of Hackney and use it as
    # grouping variable
    LTLA_shp_EN$lad20cd[LTLA_shp_EN$lad20nm == "City of London"] <- "E09000012"
    
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
    
    LTLA_shp_EN <- st_union_by(LTLA_shp_EN$geometry, LTLA_shp_EN$lad20cd)
  }
  
  # add region lookup to LTLA shapefile
  LTLA_shp_Reg  <- dplyr::left_join(LTLA_shp_EN, LTLA2Reg, by = c("lad20cd" = "LAD20CD"))
  
  # reorder shapefile so that all LTLA belonging to the same region have neighboring indexes (this simplify indexing when summing stuff over region in the nimble code)
  LTLA_shp_Reg <- LTLA_shp_Reg[order(LTLA_shp_Reg$RGN20CD),]
  
  return(LTLA_shp_Reg)
}


########################################################################
# Calculate binomial confidence intervals using Wilson's method
# (Function adapted from Hmisc R package under GPL-3 License)
# (See https://cran.r-project.org/web/packages/Hmisc/index.html )
########################################################################
binconf <- function (x, n, alpha = 0.05) {
  
  bc <- function(x, n, alpha) {
    zcrit <- -qnorm(alpha/2)
    z2 <- zcrit * zcrit
    p <- x/n
    cl <- (p + z2/2/n + c(-1, 1) * zcrit * sqrt((p * (1 - 
                                                        p) + z2/4/n)/n))/(1 + z2/n)
    if (x == 1) 
      cl[1] <- -log(1 - alpha)/n
    if (x == (n - 1)) 
      cl[2] <- 1 + log(1 - alpha)/n
    
    c(x/n, cl)
  }
  
  mat <- matrix(ncol = 3, nrow = length(x))
  for (i in 1:length(x)) mat[i, ] <- bc(x[i], n[i], alpha = alpha)
  dimnames(mat) <- list(rep("", dim(mat)[1]), c("m", "l", "u"))
  
  as.data.frame(mat, row.names = NULL)
  
}

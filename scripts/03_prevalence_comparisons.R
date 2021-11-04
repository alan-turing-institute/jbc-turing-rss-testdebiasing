dir.create("output", showWarnings = FALSE)

####################
# Get REACT dates
####################

react_round_df <- readr::read_csv("data/react_round.csv")
ltla_df <- readr::read_csv("data/ltla.csv")

round_start_dates <- as.Date(c("2020-11-13", "2021-01-06", "2021-02-04",
                               "2021-03-11", "2021-04-15"))
round_end_dates <- as.Date(c("2020-12-03", "2021-01-22", "2021-02-23",
                             "2021-03-30", "2021-05-03"))

react_date_df <- tibble(round = 7:11, 
                        start_date = round_start_dates,
                        end_date = round_end_dates) %>%
  mutate(mid_date = as.Date(round(0.5 * (as.numeric(start_date) + 
                                     as.numeric(end_date))),
                            origin = "1970-01-01")) %>%
  right_join(distinct(ltla_df, mid_week), by = character()) %>%
  group_by(round) %>%
  filter(abs(mid_week - mid_date) == min(abs(mid_week - mid_date)))
  
readr::write_csv(react_date_df, "data/react_dates.csv")

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

######################################
# Calculate raw pillar 2 estimates
######################################

raw_pillar2_df <- ltla_df %>%
  inner_join(react_date_df, by = "mid_week") %>%
  bind_cols(binconf(.$nt, .$Nt) * 100)

readr::write_csv(raw_pillar2_df, "data/raw_pillar2.csv")

############################################
# Calculate REACT LTLA exact Binomial CIs
############################################

old_bucks_names <- c("Aylesbury Vale", "Chiltern", "South Bucks", "Wycombe")
new_bucks_name <- "Buckinghamshire"

old_north_names <- c("Corby", "East Northamptonshire", "Kettering", "Wellingborough")
new_north_name <- "North Northamptonshire"

old_west_names <- c("Daventry", "Northampton", "South Northamptonshire")
new_west_name <- "West Northamptonshire"

react_ltla_df <- react_round_df %>%
  mutate(ltla = case_when(ltla %in% old_bucks_names ~ new_bucks_name,
                          ltla %in% old_north_names ~ new_north_name,
                          ltla %in% old_west_names ~ new_west_name,
                          TRUE ~ ltla)) %>%
  group_by(ltla, round) %>%
  summarise(positive = sum(positive),
            number_samples = sum(number_samples), .groups = "drop") %>%
  bind_cols(binconf(.$positive, .$number_samples) * 100)

readr::write_csv(react_ltla_df, "data/react_ltla.csv")

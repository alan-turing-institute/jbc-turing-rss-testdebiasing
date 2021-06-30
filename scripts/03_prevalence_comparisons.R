dir.create("output", showWarnings = FALSE)

####################
# Get REACT dates
####################

react_round_df <- readr::read_csv("data/react_round.csv")

react7_start_date <- as.Date("2020-11-13")
react7_end_date <- as.Date("2020-12-03")
react7_pillar2_mid_week <- as.Date("2020-11-22")

react8_start_date <- as.Date("2021-01-06")
react8_end_date <- as.Date("2021-01-22")
react8_pillar2_mid_week <- as.Date("2021-01-17")

react_date_df <- tibble(round = c(7, 8), 
                        start_date = c(react7_start_date, react8_start_date),
                        end_date = c(react7_end_date, react8_end_date)) %>%
  mutate(mid_date = as.Date(0.5 * (as.numeric(start_date) + 
                                     as.numeric(end_date)),
                            origin = "1970-01-01"))

######################################
# Calculate raw pillar 2 estimates
######################################

raw_pillar2_df <- ltla_df %>%
  left_join(react_date_df, by = character()) %>%
  group_by(round) %>%
  filter(abs(mid_week - mid_date) == min(abs(mid_week - mid_date))) %>%
  bind_cols(as_tibble(Hmisc::binconf(.$nt, .$Nt) * 100)) %>%
  rename(m = PointEst, l = Lower, u = Upper)

readr::write_csv(raw_pillar2_df, "output/raw_pillar2.csv")

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
  group_by(ltla, round, mid_round_date) %>%
  summarise(positive = sum(positive),
            number_samples = sum(number_samples), .groups = "drop") %>%
  bind_cols(as_tibble(Hmisc::binconf(.$positive, .$number_samples) * 100)) %>%
  rename(m = PointEst, l = Lower, u = Upper)

readr::write_csv(react_ltla_df, "output/react_ltla.csv")

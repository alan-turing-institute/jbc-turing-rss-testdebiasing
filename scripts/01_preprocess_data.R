library(readr)
library(readODS)
library(readxl)
library(dplyr)
library(tidyr)
library(zoo)

### Import weekly LTLA-level Pillar 1+2 testing data ###
path_to_pillar12 <- "data/Demographic_LA_tables.ods"

pillar12_Nt <- read_ods(path_to_pillar12, sheet = "Table_6", skip = 2) %>%
  pivot_longer(-(1:6), names_to = "week", values_to = "Nt")
pillar12_nt <- read_ods(path_to_pillar12, sheet = "Table_7", skip = 2) %>%
  mutate_at(vars(-(1:6)), as.numeric) %>%
  pivot_longer(-(1:6), names_to = "week", values_to = "nt") %>%
  mutate(week = if_else(week == "19/11/202-25/11/20", "19/11/20-25/11/20", week),
         nt = replace_na(nt, 0))
join_var <- names(pillar12_Nt)[1:7]
pillar12_df <- pillar12_Nt %>%
  left_join(pillar12_nt, by = join_var) %>%
  filter(LTLA != "Unknown") %>%
  mutate(week = gsub(" ", "", week)) %>%
  separate(week, c("week_start", "week_end"), sep = "-") %>%
  mutate(week_start = as.Date(week_start, "%d/%m/%y"),
         week_end = as.Date(week_end, "%d/%m/%y"),
         mid_week = as.Date((as.numeric(week_start) + as.numeric(week_end)) / 2,
                            origin = "1970-01-01")
  )

week_df <- distinct(pillar12_df, week_start, week_end, mid_week)
region_code <- pillar12_df %>%
  distinct(Code = Region, phe_region = `Region Name`)
ltla_code <- pillar12_df %>%
  distinct(Code = LTLA, ltla = `LTLA Name`)


### Import ONS population estimates ###
path_to_ons_pop <- "data/ukmidyearestimates20192020ladcodes.xls"
pop_size <- readxl::read_excel(path_to_ons_pop, sheet = "MYE2 - Persons", skip = 3) %>%
  select(1:4) %>%
  filter(Code %in% c(pillar12_df$LTLA, pillar12_df$Region))

northamptonshire_pop <- readr::read_csv("data/northamptonshire_pop_nomis.csv", skip = 5) %>%
  rename(ltla = 1, M = 2)

ltla_pop <- pop_size %>%
  inner_join(ltla_code, by = "Code") %>%
  select(ltla, M = `All ages`) %>%
  filter(ltla != "Northamptonshire") %>%
  bind_rows(northamptonshire_pop)

phe_region_pop <- pop_size %>%
  inner_join(region_code, by = "Code") %>%
  select(phe_region, M = `All ages`)


### Import REACT PHE region totals and positives ###
# Weekly data
path_to_react_totals <- "data/react_total.csv"
path_to_react_positives <- "data/react_positive.csv"

react_Nr <- read_csv(path_to_react_totals) %>%
  rename(date = X1) %>%
  pivot_longer(-date, names_to = "phe_region", values_to = "Nr")
react_nr <- read_csv(path_to_react_positives) %>%
  rename(date = X1) %>%
  pivot_longer(-date, names_to = "phe_region", values_to = "nr")
react_df <- react_Nr %>%
  left_join(react_nr, by = c("date", "phe_region")) %>%
  mutate(nr = replace_na(nr, 0)) %>%
  full_join(week_df, by = character()) %>%
  filter(date <= week_end & date >= week_start) %>%
  group_by(phe_region, mid_week) %>%
  summarise(Nr = sum(Nr),
            nr = sum(nr), .groups = "drop")
readr::write_csv(react_df, "data/react.csv")

# Round data
path_to_react_7a <- "data/unwt_ordered_ltla_prev7a.csv"
path_to_react_7b <- "data/unwt_ordered_ltla_prev7b.csv"
react_7a <- read_csv(path_to_react_7a)
react_7b <- read_csv(path_to_react_7b)
react_7 <- bind_rows(react_7a, react_7b) %>%
  group_by(ltla) %>%
  summarise(positive = sum(positive), 
            number_samples = sum(number_samples), .groups = "drop") %>%
  mutate(round = 7)

react_round_df <- react_7
for (this_round in 8:11) {
  path_to_react <- paste0("data/unwt_ordered_ltla_prev", this_round, ".csv")
  this_react <- read_csv(path_to_react) %>%
    select(ltla, positive, number_samples) %>%
    mutate(round = this_round)
  react_round_df <- bind_rows(react_round_df, this_react)
}

readr::write_csv(react_round_df, "data/react_round.csv")

### Get regional counts ###

region_df <- pillar12_df %>%
  rename(phe_region = `Region Name`) %>%
  group_by(phe_region, mid_week) %>%
  summarise(Nt = sum(Nt),
            nt = sum(nt), .groups = "drop") %>%
  full_join(react_df, by = c("phe_region", "mid_week")) %>%
  mutate(Nr = replace_na(Nr, 0),
         nr = replace_na(nr, 0)) %>%
  left_join(phe_region_pop, by = "phe_region")

write_csv(region_df, "data/region.csv")

### Get LTLA counts ###

ltla_df <- pillar12_df %>%
  select(ltla = `LTLA Name`, phe_region = `Region Name`, mid_week, Nt, nt) %>%
  left_join(ltla_pop, by = "ltla")

write_csv(ltla_df, "data/ltla.csv")


### Get vaccination data ###
phe_region_to_region_df <- tibble(phe_region = c("East of England",
                                                 "London",
                                                 "West Midlands",
                                                 "East Midlands",
                                                 "North East",
                                                 "North West",
                                                 "South East",
                                                 "South West",
                                                 "Yorkshire and The Humber"),
                                  region = c("East Of England",
                                             "London",
                                             "Midlands",
                                             "Midlands",
                                             "North East And Yorkshire",
                                             "North West",
                                             "South East",
                                             "South West",
                                             "North East And Yorkshire"))
region_pop <- phe_region_pop %>%
  left_join(phe_region_to_region_df, by = "phe_region") %>%
  group_by(region) %>%
  summarise(M = sum(M), .groups = "drop")
path_to_vax <- "data/COVID-19-monthly-announced-vaccinations.xlsx"
path_to_weekly_vax <- "data/COVID-19-weekly-announced-vaccinations.xlsx"

vax_start_mid_week <- as.Date("2020-12-13") + (2 * 7)
vax_data_mid_week <- as.Date("2021-01-27")
latest_vax_date <- as.Date("2021-06-20")

weekly_vax_df <- readxl::read_excel(path_to_weekly_vax,
  sheet = "NHS Region", skip = 15, n_max = 7, col_names = FALSE
) %>%
  select(region = 1, Count = 17) %>%
  mutate(date = latest_vax_date)

monthly_vax_df <- readxl::read_excel(path_to_vax,
                     sheet = "Vaccination Date",
                     skip = 11, n_max = 141) %>%
  dplyr::select(date = 1, 3:9) %>%
    tidyr::pivot_longer(-date, names_to = "region", values_to = "Count") %>%
    mutate(date = as.Date(date), region = gsub("...\\d", "", region)) %>%
    bind_rows(weekly_vax_df)
    
vax_df <- region_pop %>%
  rename(M_reg = M) %>%
  full_join(tibble(date = seq(vax_start_mid_week - 1, latest_vax_date, by = 1)), by = character()) %>%
  left_join(monthly_vax_df, by = c("date", "region")) %>%
  group_by(region) %>%
  arrange(region, date) %>%
  mutate(Count = if_else(date == vax_start_mid_week - 1, 0, Count),
         count_fill = round(zoo::na.approx(Count)),
         prop_vax = count_fill / M_reg) %>%
  right_join(phe_region_to_region_df, by = "region") %>%
    right_join(ltla_df, by = c("phe_region", "date" = "mid_week")) %>%
    mutate(
      prop_vax = if_else(date < vax_start_mid_week, 0, prop_vax),
      V = round(prop_vax * M)
    ) %>%
    ungroup() %>%
  select(ltla, mid_week = date, V) %>%
  arrange(ltla, mid_week)

write_csv(vax_df, "data/vaccination.csv")

# max(vax_df$mid_week)
# dim(vax_df[vax_df$ltla == "Adur", "V"])
# plot(1:sum(vax_df$ltla == "Adur"), unlist(vax_df[vax_df$ltla == "Adur", "V"]))

### Get variant data ###
path_to_sanger <- "data/UK_variant_data_Sanger.tsv"
variants_in <- as.data.frame(read_tsv(path_to_sanger))
variants_in$mid_week <- as.character(as.Date(variants_in$WeekEndDate) - 3)
variants_in$ltla <- unlist(ltla_code[match(variants_in$LTLA, ltla_code$Code), "ltla"])
variants_in$name_date <- paste0(variants_in$ltla, "_", variants_in$mid_week)
name_date_unique <- unique(variants_in$name_date)
variants_out <- variants_in[match(name_date_unique, variants_in$name_date), c("mid_week", "ltla", "name_date")]
variants_out$n_delta <- sapply(name_date_unique, function(x) 
  sum(variants_in[variants_in$name_date == x & variants_in$Lineage == "B.1.617.2", "Count"]))
variants_out$n_tot <- sapply(name_date_unique, function(x) 
  sum(variants_in[variants_in$name_date == x, "Count"]))
variants_out$prop_delta <- variants_out$n_delta / variants_out$n_tot
variants_out <- variants_out[order(variants_out$ltla, variants_out$mid_week), ]
variants_out$week <- variants_out$mid_week
readr::write_csv(variants_out, "data/variants.csv")

### Get epimap Rt data ###

path_to_combined_epimap <- "data/Rt_estimates_Epimap_combined.csv"
r_in <- read.csv(file = path_to_combined_epimap, stringsAsFactors = F)  
r_in$name_date <- paste0(r_in$area, "_", r_in$Date)
mid_week_unique <- sort(unique(as.character(week_df$mid_week)))
path_to_imperial <- "data/UK_hotspot_Rt_estimates_Imperial.csv"
r_to_dup <- read.csv(file = path_to_imperial, stringsAsFactors = F)  
str(r_to_dup)
r_out <- r_in[r_in$Date %in% mid_week_unique, c("area", "Date")]
colnames(r_out) <- c("area", "date")
r_out[, c("CIlow", "Rt", "CIup")] <- r_in[match(paste0(r_out$area, "_", r_out$date), r_in$name_date), c("Rt_2_5", "Rt_50", "Rt_97_5")]
r_out$coverage <- .95
path_to_preproc_Epimap_Rt <- "data/Rt_estimates_Epimap_combined_preprocessed.csv"
write.csv(r_out, file = path_to_preproc_Epimap_Rt, row.names = F)





















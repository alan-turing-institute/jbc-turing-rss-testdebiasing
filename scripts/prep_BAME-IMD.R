library(dplyr)
library(readr)
library(readxl)

old_bucks_names <- c("Aylesbury Vale", "Chiltern", "South Bucks", "Wycombe")
new_bucks_name <- "Buckinghamshire"

old_north_names <- c("Corby", "East Northamptonshire", "Kettering", "Wellingborough")
new_north_name <- "North Northamptonshire"

old_west_names <- c("Daventry", "Northampton", "South Northamptonshire")
new_west_name <- "West Northamptonshire"

path_to_ons_pop <- "data/ukmidyearestimates20192020ladcodes.xls"
pop_size <- readxl::read_excel(path_to_ons_pop, sheet = "MYE2 - Persons", skip = 3) %>%
  select(ltla = Name, M = 4) %>%
  filter(ltla %in% c(ltla_df$ltla, old_bucks_names, old_north_names, old_west_names))

northamptonshire_pop <- readr::read_csv("data/northamptonshire_pop_nomis.csv", 
                                        skip = 5) %>%
  rename(ltla = 1, M = 2)

ltla_pop <- pop_size %>%
  filter(ltla != "Northamptonshire") %>%
  bind_rows(northamptonshire_pop)

bame_raw <- read_csv("data/census2011_ethnicity_LSOA_black_southasian.csv", 
                     skip = 7)
LSOA2LTLA <- read_csv("data/Output_Area_to_Lower_Layer_Super_Output_Area_to_Middle_Layer_Super_Output_Area_to_Local_Authority_District__December_2020__Lookup_in_England_and_Wales.csv", 
                      col_types = cols(FID = col_skip(), OA11CD = col_skip()))
LSOA2LTLA <- distinct(LSOA2LTLA)
bame_raw <- bame_raw %>% 
  rename(LSOA11CD = `2011 super output area - lower layer`) %>% 
  mutate(LSOA11CD = gsub(" :.*","", LSOA11CD))

bame_prop <-  LSOA2LTLA %>% 
  left_join(bame_raw, by = c("LSOA11CD" = "LSOA11CD" )) %>% 
  dplyr::select(-c("RGN20CD", "RGN20NM")) %>%
  mutate(LAD20NM = case_when(
    LAD20NM %in% old_bucks_names ~ new_bucks_name,
    LAD20NM %in% old_west_names ~ new_west_name,
    LAD20NM %in% old_north_names ~ new_north_name,
    TRUE ~ LAD20NM
  )) %>%
  group_by(LAD20NM) %>% 
  mutate(Black = `Black/African/Caribbean/Black British`, 
         South_Asian = `Asian/Asian British: Pakistani` + `Asian/Asian British: Indian` + `Asian/Asian British: Bangladeshi`,
         Other_BAME = `All usual residents` - White - `Black/African/Caribbean/Black British` - 
           (`Asian/Asian British: Pakistani` + `Asian/Asian British: Indian` + 
              `Asian/Asian British: Bangladeshi`)) %>%
  summarise_at(c("All usual residents", "Black", "South_Asian", "White", "Other_BAME"), sum) %>%
  mutate(Black_prop = Black/`All usual residents`, 
         South_Asian_prop = South_Asian/`All usual residents`,
         Other_BAME_prop = Other_BAME/`All usual residents`,
         BAME = 1 - White/`All usual residents`) %>% 
  dplyr::select(ltla  = LAD20NM, Black_prop,South_Asian_prop,Other_BAME_prop, BAME )

bame_prop <- bame_prop %>% 
  ungroup %>% 
  mutate(bame_quint = gtools::quantcut(BAME, q=5, na.rm=TRUE))

IMD_score <- readxl::read_excel("data/File_10_-_IoD2019_Local_Authority_District_Summaries__lower-tier__.xlsx", 
                               sheet = "IMD") %>%
  select(ltla = 2, imd = 5) %>%
  mutate(ltla = case_when(
    ltla %in% old_bucks_names ~ new_bucks_name,
    ltla %in% old_west_names ~ new_west_name,
    ltla %in% old_north_names ~ new_north_name,
    TRUE ~ ltla
  )) %>%
  group_by(ltla) %>%
  summarise(imd = mean(imd))

cov_data <- left_join(bame_prop, IMD_score, by = 'ltla') %>%
  mutate(imd_quint = gtools::quantcut(imd, q=5, na.rm=TRUE))


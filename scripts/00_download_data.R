dir.create("data", showWarnings = FALSE)

# Make sure to install prevdebiasr first
# (GN: I did this by creating local prevdebiasr.tar.gz and using install.packages())
library(prevdebiasr)

### Download weekly Pillar 1+2 data ###
# https://www.gov.uk/government/publications/nhs-test-and-trace-england-statistics-21-january-to-27-january-2021

list.files()
# url_to_pillar12 <- "https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/958348/Demographic_LA_tables_w35.ods"
url_to_pillar12 <- "https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/991016/Demographic_LA_tables_week52.ods"
path_to_pillar12 <- "data/Demographic_LA_tables_w35.ods"
download.file(url_to_pillar12, path_to_pillar12)


### Download REACT data ###
# url_to_react <- "https://raw.githubusercontent.com/mrc-ide/reactidd/bfc3a2577a692441a8de7a6f7aadbaf06741f463/inst/extdata"
url_to_react <- "https://github.com/mrc-ide/reactidd/raw/master/inst/extdata"
path_to_react_totals <- "data/react_total.csv"
download.file(file.path(url_to_react, "total.csv"), path_to_react_totals)
path_to_react_positives <- "data/react_positive.csv"
download.file(file.path(url_to_react, "positive.csv"), path_to_react_positives)

path_to_react_round7a <- "data/unwt_ordered_ltla_prev7a.csv"
download.file(file.path(url_to_react, "unwt_ordered_ltla_prev7a.csv"), 
              path_to_react_round7a)

path_to_react_round7b <- "data/unwt_ordered_ltla_prev7b.csv"
download.file(file.path(url_to_react, "unwt_ordered_ltla_prev7b.csv"), 
              path_to_react_round7b)

path_to_react_round8 <- "data/unwt_ordered_ltla_prev8.csv"
download.file(file.path(url_to_react, "unwt_ordered_ltla_prev8.csv"), 
              path_to_react_round8)

### Download ONS population estimates ###
# https://www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/populationestimates/datasets/populationestimatesforukenglandandwalesscotlandandnorthernireland

url_to_ons_pop <- "https://www.ons.gov.uk/file?uri=/peoplepopulationandcommunity/populationandmigration/populationestimates/datasets/populationestimatesforukenglandandwalesscotlandandnorthernireland/mid2019april2020localauthoritydistrictcodes/ukmidyearestimates20192020ladcodes.xls"
path_to_ons_pop <- "data/ukmidyearestimates20192020ladcodes.xls"
download.file(url_to_ons_pop, path_to_ons_pop)

### Download NHS vaccination data ###
# https://www.england.nhs.uk/statistics/statistical-work-areas/covid-19-vaccinations/covid-19-vaccinations-archive/ 
# url_to_vax <- "https://www.england.nhs.uk/statistics/wp-content/uploads/sites/2/2021/01/COVID-19-monthly-announced-vaccinations-14-January-2021.xlsx"
url_to_vax <- "https://www.england.nhs.uk/statistics/wp-content/uploads/sites/2/2021/05/COVID-19-monthly-announced-vaccinations-13-May-2021.xlsx"
path_to_vax <- "data/COVID-19-monthly-announced-vaccinations-14-January-2021.xlsx"
download.file(url_to_vax, path_to_vax)

### Download shape files for spatial plots ###
# https://geoportal.statistics.gov.uk/datasets/054349b09c094df2a97f8ddbd169c7a7_0
url_to_lookup <- "https://opendata.arcgis.com/datasets/054349b09c094df2a97f8ddbd169c7a7_0.csv"
path_to_lookup <- "data/Local_Authority_District_to_Region__December_2020__Lookup_in_England.csv"
download.file(url_to_lookup, path_to_lookup)

# https://geoportal.statistics.gov.uk/datasets/9d86e7bcf2864343ab72e2914756b86d_0
url_to_shapefile <- "https://opendata.arcgis.com/datasets/9d86e7bcf2864343ab72e2914756b86d_0.zip?outSR=%7B%22latestWkid%22%3A27700%2C%22wkid%22%3A27700%7D"
temp <- tempfile()
path_to_shapefile <- "data/Local_Authority_Districts_(May_2020)_Boundaries_UK_BFE-shp"
download.file(url_to_shapefile , temp)
unzip(temp, exdir = path_to_shapefile)
unlink(temp)


### Download Imperial College Rt estimates ###

url_to_imperial <- "https://imperialcollegelondon.github.io/covid19local/downloads/UK_hotspot_Rt_estimates.csv"
path_to_imperial <- "data/UK_hotspot_Rt_estimates_Imperial.csv"
download.file(url_to_imperial, path_to_imperial)

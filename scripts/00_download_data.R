dir.create("data", showWarnings = FALSE)

# Make sure to install prevdebiasr first
# (GN: I did this by creating local prevdebiasr.tar.gz and using install.packages())
library(prevdebiasr)

### Download weekly Pillar 1+2 data ###
# Procedure for updating Pillar 1+2 data links
# 1.	Go to https://www.gov.uk/government/collections/nhs-test-and-trace-statistics-england-weekly-reports
# 2.	Scroll down to “Latest Report” and click the link
# 3.	Scroll down to “Demographic and regional information for people tested and testing positive, STARTDATE to ENDDATE: data tables”
# 4.	Copy that link address to the url_to_pillar12 <- "PASTE ME HERE"
url_to_pillar12 <- "https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/998213/Demographic_output_w56.ods"
path_to_pillar12 <- "data/Demographic_LA_tables.ods"
download.file(url_to_pillar12, path_to_pillar12)

### Download REACT data ###
# Procedure for updating REACT data links
# 1.	Go to https://github.com/mrc-ide/reactidd/raw/master/inst/extdata
# 2.	Find additional files of the form unwt_ordered_ltla_prev*.csv and add on in below format
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

path_to_react_round9 <- "data/unwt_ordered_ltla_prev9.csv"
download.file(file.path(url_to_react, "unwt_ordered_ltla_prev9.csv"), 
              path_to_react_round9)

path_to_react_round10 <- "data/unwt_ordered_ltla_prev10.csv"
download.file(file.path(url_to_react, "unwt_ordered_ltla_prev10.csv"), 
              path_to_react_round10)

path_to_react_round11 <- "data/unwt_ordered_ltla_prev11.csv"
download.file(file.path(url_to_react, "unwt_ordered_ltla_prev11.csv"), 
              path_to_react_round11)

### Download NHS vaccination data ###
# Procedure for updating vaccination data links
# 1.	Go to https://www.england.nhs.uk/statistics/statistical-work-areas/covid-19-vaccinations/covid-19-vaccinations-archive/ 
# 2.	Scroll to bottom of page where the "Monthly Covid-19 vaccinations data archive" is
# 3.  Copy and paste most recent monthly data URL below to url_to_vax <- "HERE PLEASE"
url_to_vax <- "https://www.england.nhs.uk/statistics/wp-content/uploads/sites/2/2021/06/COVID-19-monthly-announced-vaccinations-10-June-2021.xlsx"
url_to_most_recent_weekly_vax <- "https://www.england.nhs.uk/statistics/wp-content/uploads/sites/2/2021/06/COVID-19-weekly-announced-vaccinations-24-June-2021.xlsx"
path_to_vax <- "data/COVID-19-monthly-announced-vaccinations.xlsx"
path_to_weekly_vax <- "data/COVID-19-weekly-announced-vaccinations.xlsx"
download.file(url_to_vax, path_to_vax)
download.file(url_to_most_recent_weekly_vax, path_to_weekly_vax)


### Download ONS population estimates ###
# https://www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/populationestimates/datasets/populationestimatesforukenglandandwalesscotlandandnorthernireland
url_to_ons_pop <- "https://www.ons.gov.uk/file?uri=/peoplepopulationandcommunity/populationandmigration/populationestimates/datasets/populationestimatesforukenglandandwalesscotlandandnorthernireland/mid2019april2020localauthoritydistrictcodes/ukmidyearestimates20192020ladcodes.xls"
path_to_ons_pop <- "data/ukmidyearestimates20192020ladcodes.xls"
download.file(url_to_ons_pop, path_to_ons_pop)

### NOTE: West and North Northamptonshire estimates available from
# https://www.nomisweb.co.uk/datasets/pestsyoala

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

url_to_imperial_old <- "https://zenodo.org/record/4668480/files/ImperialCollegeLondon/covid19local-v21.0.zip"
path_to_imperial_old <- "data/UK_hotspot_Rt_estimates_2021-04-07.csv"
temp <- tempfile()
download.file(url_to_imperial_old, temp)
data <- readr::read_csv(unz(temp, "ImperialCollegeLondon-covid19local-7d3b1fd/downloads/UK_hotspot_Rt_estimates.csv"))
readr::write_csv(data, path_to_imperial_old)
unlink(temp)
rm(data)

### Download Sanger variant data ###

url_to_sanger_variant_data <- "https://covid-surveillance-data.cog.sanger.ac.uk/download/lineages_by_ltla_and_week.tsv"
path_to_sanger <- "data/UK_variant_data_Sanger.tsv"
download.file(url_to_sanger_variant_data, path_to_sanger)

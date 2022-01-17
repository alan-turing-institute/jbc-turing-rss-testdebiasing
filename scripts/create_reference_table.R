# install.packages("fs")
# library(googledrive)
# library(help = googledrive)
# 
# drive_get(path = )
# 
# temp <- tempfile()
# download.file("https://docs.google.com/spreadsheets/d/1ekku1zeFb0Z2lPSu0-LdmNkMBVF5Aa15PZotIeZ08jM/edit?usp=sharing",
#               temp)
# dl <- drive_download(
#   as_id("1AiZda_1-2nwrxI8fLD0Y6e5rTg7aocv0"), path = temp, overwrite = TRUE)
# out <- unzip(temp, exdir = tempdir())
# bank <- read.csv(out[14], sep = ";")



library(gsheet)
library(xtable)
reftab <- gsheet::gsheet2tbl('docs.google.com/spreadsheets/d/1ekku1zeFb0Z2lPSu0-LdmNkMBVF5Aa15PZotIeZ08jM/edit?usp=sharing')
tabout <- print(xtable::xtable(reftab, 
                                label = "tab:reference_table", 
                                # align = rep("r", ncol(reftab) + 1),
                               align = c("p{1cm}", "p{3cm}", "|", "p{3cm}", "|", "p{6cm}", "|", "p{3cm}"),
                               caption = "Details of related work"),
                                caption.placement = "top", 
                                sanitize.text.function = function(x){x}, 
                                include.rownames = F,
                                hline.after = 0:(nrow(reftab) - 1))

drop_folder <- "C:/Users/nicho/Dropbox/Apps/Overleaf/Estimating local prevalence from targeted testing data"

cat(tabout, file = file.path(drop_folder, "reference_table.tex"))

?xtable

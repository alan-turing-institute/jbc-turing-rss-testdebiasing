d <- read.csv("data/pillar12_data_for_reporting_summary.csv")
d$f <- as.numeric(gsub(",", "", d$Female))
d$m <- as.numeric(gsub(",", "", d$Male))
d$both <- d$f + d$m
d$totperc <- round(100 * (d$both / sum(d$both)))
paste(paste0(d$Age, "yr: ", format(d$both, nsmall = 0, big.mark = ","), paste0(" (", d$totperc, "%)")), collapse = "; ")
m_tot <- sum(d$m)
f_tot <- sum(d$f)
round(100 * m_tot / (m_tot + f_tot))

paste0(paste0("Males: ", format(m_tot, nsmall = 0, big.mark = ","), " (", round(100 * m_tot / (m_tot + f_tot)), "%)"), 
       "; ",
        paste0("Females: ", format(f_tot, nsmall = 0, big.mark = ","), " (", round(100 * f_tot / (m_tot + f_tot)), "%)"))

sum(d$f)
d$m / 
str(d)


d1 <- read.csv("data/round_8_go.csv")
age_unique <- c("5-12", "13-17", "18-24", "25-34", "35-44", "45-54", "55-64", "65+")
agg_N <- sapply(age_unique, function(age) sum(d1[d1$age_group == age, "number_samples"]))
d2 <- data.frame(Age = age_unique, both = agg_N)
d2$totperc <- round(100 * (d2$both / sum(d2$both)))
paste(paste0(d2$Age, "yr: ", format(d2$both, nsmall = 0, big.mark = ","), paste0(" (", d2$totperc, "%)")), collapse = "; ")
sapply()
unique(d1$age_group)
dput(age_unique)
str(d1)
range(d1$react_week_start_date)

###########################################
# load in, summarize by sample/date and plot GHG data from NHC 2019-2020

# A carter
# 5/5/2020
setwd("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/Projects/NHC Metabolism 2019-2020")

library(readr)
library(dplyr)

dat <- read_csv("../field data/2019/NHC gas data/NHC_2019-2020_processed_GHGdata.csv")

dat.mean <- dat %>% group_by(Site, Date) %>%
  summarise(CH4.ugL = mean(CH4.ugL, na.rm=T),
            CO2.ugL = mean(CO2.ugL, na.rm=T), 
            N2O.ugL = mean(N2O.ugL, na.rm=T), 
            )
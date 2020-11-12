# Pull data from streampulse website
# AM Carter
# 2020.03.31
# Adapted from:

#Basic StreamPULSE data processing pipeline
#Updated 10/29/18
#Contact Mike Vlah (vlahm13@gmail.com) with questions or comments.

# Update the streamPULSE package from github
# The StreamPULSE package is in development and changes frequently!
# If something doesn't work as expected, first try reinstalling.

#library(devtools)
#install_github('streampulse/StreamPULSE', dependencies=TRUE)
setwd("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism/data/")

# Load packages.
library(StreamPULSE)
library(lubridate)
library(tidyverse)

sites <- read_csv("siteData/NHCsite_metadata.csv")
sites$startdate.EST <- as.POSIXct(sites$startdate.UTC, 
                                  format="%m/%d/%Y %H:%M", tz = "UTC") %>%
                        with_tz(tzone="EST")
sites$enddate.EST <- as.POSIXct(sites$enddate.UTC, 
                                  format="%m/%d/%Y %H:%M", tz = "UTC") %>%
                        with_tz(tzone="EST")

# subset for the sites we actually want:
nhc_sites <- sites[c(1:5,7),]
mud_sites <- sites[8:9,]
# # find date range of available data for each site
# for(i in 1:nrow(sites)){
#    tmp <- query_available_data("NC",sites[i,]$sitecode)
#    sites[i,9:10] <- tmp$datebounds
#  }
#  sites[,9:10]<- with_tz(sites[,9:10], tz="UTC")
#  write.csv(sites, file = "siteData/NHCsite_metadata.csv", row.names = F)

# Download data from streampulse for sampling sites.
# dateRange <- c(sites$startdate.EST[2], sites$enddate.EST[2])
  
# download all the data for each at site:
for(i in 1:nrow(nhc_sites)){
  vars <- as.character(unlist(
    query_available_data("NC",nhc_sites[i,]$sitecode)$variables))
  dat <- request_data(nhc_sites[i,]$siteID, 
                      # startdate=dateRange[1], 
                      # enddate = dateRange[2], 
                      variables=vars)
  dd <- dat$data %>% 
    mutate(value = ifelse(flagtype %in% c("Bad Data", "Questionable"),
                          NA, value)) %>%
    select(DateTime_UTC, site,value, variable) %>%
    pivot_wider(names_from = variable, values_from = value)
  write_csv(dd, paste0("metabolism/raw/", dd$site[1], "_", 
                   as.Date(dat$specs$enddate), ".csv"))  
  # write_rds(dat, path = paste0("metabolism/raw/",
  #                              nhc_sites[i,]$sitecode, 
  #                              "_2019-2020_raw.rds"))
}

# Mud creek site data download
for(i in 1:nrow(mud_sites)){
  vars <- as.character(unlist(
    query_available_data("NC",mud_sites[i,]$sitecode)$variables))
  dat <- request_data(sites[i,]$siteID, 
                      startdate=mud_sites$startdate.EST[1], 
                      enddate = mud_sites$enddate.EST[1], 
                      variables=vars)
  write_rds(dat, path = paste0("data/metabolism/raw/",
                               nhc_sites[i,]$sitecode, "_", 
                               as.character(date(mud_sites$enddate.EST[1])), 
                               "_raw.rds"))
}


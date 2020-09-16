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
setwd("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism")

# Load packages.
library(StreamPULSE)
library(lubridate)


sites <- read.csv("data/siteData/NHCsite_metadata.csv", header=T, stringsAsFactors = F)

# # find date range of available data for each site
# for(i in 1:nrow(sites)){
#    tmp <- query_available_data("NC",sites[i,]$sitecode)
#    sites[i,9:10] <- tmp$datebounds
#  }
#  sites[,9:10]<- with_tz(sites[,9:10], tz="UTC")
#  write.csv(sites, file = "siteData/NHCsite_metadata.csv", row.names = F)

# Download data from streampulse for sampling sites.
dateRange <- c(min(sites[sites$type=="synoptic"&sites$sitecode!="MC751",]$startdate.UTC), 
               max(sites[sites$type=="synoptic"&sites$sitecode!="MC751",]$enddate.UTC))
sites[sites$type=="core",]$startdate.UTC <- dateRange[1]
sites[sites$type=="core",]$enddate.UTC <- dateRange[2]
sites$startdate.UTC <- as.POSIXct(sites$startdate.UTC, format="%m/%d/%Y %H:%M")
sites$enddate.UTC<- as.POSIXct(sites$enddate.UTC, format="%m/%d/%Y %H:%M")
# List of all the variables at site:
for(i in 1:nrow(sites)){
  vars <- as.character(unlist(query_available_data("NC",sites[i,]$sitecode)$variables))
  dat <- request_data(sites[i,]$siteID, startdate=sites[i,]$startdate.UTC,enddate = sites[i,]$enddate.UTC, variables=vars)
  write.csv(dat$data, file = paste0("data/metabolism/raw/",sites[i,]$sitecode, "_",
                                    as.character(as.Date(sites[i,]$enddate.UTC)),"_raw.csv"), row.names=F)
}

# specific MC751 download:
sites <- sites[sites$sitename=="MC751",]
vars <- as.character(unlist(query_available_data("NC", "MC751")$variables))
sites$startdate.UTC <- as.POSIXct(sites$startdate.UTC, format="%m/%d/%Y %H:%M")
sites$enddate.UTC<- as.POSIXct(sites$enddate.UTC, format="%m/%d/%Y %H:%M")
dat <- request_data(sites$siteID, startdate=sites$startdate.UTC, enddate=(sites$enddate.UTC), variables=vars)

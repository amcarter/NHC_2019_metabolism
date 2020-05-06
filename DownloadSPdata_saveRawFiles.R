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
setwd("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/Projects/NHC Metabolism 2019-2020")

# Load packages.
library(StreamPULSE)
library(streamMetabolizer)
library(tidyr)
library(lubridate)



sites <- read.csv("data/siteData/NHCsite_coordinates.csv", header=T, stringsAsFactors = F)
sites <- sites[sites$type!="USGS", ]
sites$siteID <-  paste("NC", sites$sitecode, sep="_")
sites$startdate.UTC <- as.POSIXct(NA)
sites$enddate.UTC <- as.POSIXct(NA)

# # find date range of available data for each site
# for(i in 1:nrow(sites)){
#   tmp <- query_available_data("NC",sites[i,]$sitecode)
#   sites[i,9:10] <- tmp$datebounds
# }
# sites[,9:10]<- with_tz(sites[,9:10], tz="UTC")
# write.csv(sites, file = "siteData/NHCsite_metadata.csv", row.names = F)

# Download data from streampulse for synoptic sampling sites.
synoptic_sites <- sites[sites$type=="synoptic",]

# List of all the variables at site:
for(i in 1:nrow(synoptic_sites)){
  vars <- as.character(unlist(query_available_data("NC",synoptic_sites[i,]$sitecode)$variables))
  dat <- request_data(synoptic_sites[i,]$siteID, variables=vars)
  write.csv(dat$data, file = paste0("data/metabolism/raw/csv/",synoptic_sites[i,]$sitecode, "_",
                                    as.character(as.Date(synoptic_sites[i,]$enddate.UTC)),"_raw.csv"),
            row.names=F)
  saveRDS(dat, file=paste0("data/metabolism/raw/rds/",synoptic_sites[i,]$sitecode, "_",
                           as.character(as.Date(synoptic_sites[i,]$enddate.UTC)),"_raw.rds"))
}


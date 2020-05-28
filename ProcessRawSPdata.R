####################################
# Process raw NHC datafiles:
#   1. Remove bad data
#   3. calculate % sat using streampulse package
#   4. calculate discharge by interpolating data from NHC and UNHC based on watershed size

library(streamMetabolizer)
library(lubridate)
library(tidyverse)

sites <- read_csv("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism/data/siteData/NHCsite_metadata.csv")
Qdat <- read_csv("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism/data/siteData/interpolatedQ_allsites.csv")
filelist <- list.files("data/metabolism/raw")

# get rid of MC751 and Mud for now
filelist <- filelist[-c(2,3)]

#look at required inputs for a bayesian model in stream Metabolizer
metab_inputs(type="bayes", input="data")

for(i in 1:length(filelist)){
  dat <- read.csv(paste0("data/metabolism/raw/",filelist[i]), header=T, stringsAsFactors = F)
  dat$DateTime_UTC <- ymd_hms(dat$DateTime_UTC)
  sitename <- dat$site[1]
  lat <- sites$latitude[sites$sitecode==sitename]
  lon <- sites$longitude[sites$sitecode==sitename]
  dat$value[dat$flagtype=="Bad Data"|dat$flagtype=="Questionable"] <- NA
  dat <- select(dat,-site,-region )%>%
    pivot_wider(names_from = variable, values_from = value) 
  if("AirPres_kPa"%in% colnames(dat)){dat<- select(dat, -AirPres_kPa)}
  # remove leading and ending NAs
  w<- which(!is.na(dat$DO_mgL))
  dat <- dat[min(w):max(w),]
  
  #Load discharge data
  Qname <- paste(sitename, "Q", sep=".")
  Q <- Qdat[,c(1,2,which(colnames(Qdat)==Qname))]
  colnames(Q) <- c("DateTime_UTC","AirPres_kPa", "Discharge_m3s")
  
  dat<- left_join(dat, Q, by="DateTime_UTC")
  
  #Convert datetime to solar time 
  dat$DateTime_EST <- with_tz(dat$DateTime_UTC, tz="Etc/GMT+5") # convert to EST timezone
  dat$solar.time <- calc_solar_time(dat$DateTime_EST, longitude=lon)
  
  # Calculate DO saturation, depth, light
  dat$DO.sat <- calc_DO_sat(dat$WaterTemp_C, dat$AirPres_kPa*10)
  dat$depth <- calc_depth(dat$Discharge_m3s)
  dat$light <- calc_light(dat$solar.time, latitude=lat,longitude=lon)
  
  dat <- select(dat, -DateTime_EST)
  write_csv(dat, paste0("data/metabolism/processed/",sitename,".csv" ))
}



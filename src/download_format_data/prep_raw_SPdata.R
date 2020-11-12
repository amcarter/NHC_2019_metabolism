####################################
# Process raw NHC datafiles:
#   3. calculate % sat using streampulse package
#   4. calculate discharge by interpolating data from NHC and UNHC based on watershed size

library(streamMetabolizer)
library(lubridate)
library(tidyverse)

setwd("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism/data/")

sites <- read_csv("siteData/NHCsite_metadata.csv")
Qdat <- read_csv("siteData/interpolatedQ_allsites.csv")
Qdat[which(Qdat$AirPres_kPa>106), ]$AirPres_kPa <- NA
filelist <- list.files("metabolism/raw")

# get rid of MC751 and Mud for now
# filelist <- filelist[-c(2,3)]

#look at required inputs for a bayesian model in stream Metabolizer
metab_inputs(type="bayes", input="data")

for(i in 1:length(filelist)){
  dat <- read_csv(paste0("metabolism/raw/",filelist[i]))
  sitename <- dat$site[1]
  lat <- sites$latitude[sites$sitecode==sitename]
  lon <- sites$longitude[sites$sitecode==sitename]
  sensor_offset <- sites$sensor_offset_m[sites$sitecode==sitename]
  
  if("AirPres_kPa"%in% colnames(dat)){dat<- select(dat, -AirPres_kPa)}
  if("AirTemp_C"%in% colnames(dat)){dat<- select(dat, -AirTemp_C)}
  
  # remove leading and ending NAs
  w <- which(!is.na(dat$DO_mgL))
  dat <- dat[min(w):max(w), ]
  
  #Load discharge data
  Qname <- paste(sitename, "Q", sep=".")
  Q <- Qdat %>% select (DateTime_UTC, AirPres_kPa, AirTemp_C, 
                        discharge_m3s = all_of(Qname))
  
  dat <- left_join(dat, Q, by="DateTime_UTC")
  
  #Convert datetime to solar time 
  dat$DateTime_EST <- with_tz(dat$DateTime_UTC, tz="EST") # convert to EST timezone
  dat$solar.time <- calc_solar_time(dat$DateTime_EST, longitude=lon)
  
  # Calculate DO saturation, depth, light
  dat$AirPres_mbar <- dat$AirPres_kPa*10
  dat$DO.sat <- calc_DO_sat(dat$WaterTemp_C, dat$AirPres_mbar,
                            salinity.water = 0, model = "garcia-benson")
    # This number is probably not very representative and 
    # needs to be changed based on field data
  dat$depth <- calc_depth(dat$discharge_m3s)
  dat$light <- calc_light(dat$solar.time, latitude=lat,longitude=lon)
  
  # Calculate Level data from water pressure
  # Depth = pressure_Pa = kg/ms2/(density_kg/m3*gravity_m/s2)
  # density is temperature dependent, for now I am assuming it's just 998 kg/m3
  
  pressure_Pa <- (dat$WaterPres_kPa-dat$AirPres_kPa)*1000
  dat$level_m <- sensor_offset + pressure_Pa/(998*9.8)

  
  dat <- dat %>% 
    select(-DateTime_UTC, -AirPres_mbar, -AirPres_kPa, 
           -AirTemp_C, -WaterPres_kPa) %>%
    mutate(DO.obs = DO_mgL, 
           temp.water = WaterTemp_C, 
           discharge = discharge_m3s)
  
  write_csv(dat, paste0("metabolism/processed/",sitename,".csv"))
}



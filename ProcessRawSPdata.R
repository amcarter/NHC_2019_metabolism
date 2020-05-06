####################################
# Process raw NHC datafiles:
#   1. Remove bad data
#   2. Inspect data visually
#   3. calculate % sat using streampulse package
#   4. calculate discharge by interpolating data from NHC and UNHC based on watershed size

library(streamMetabolizer)
library(lubridate)
library(tidyverse)

sites <- read_csv("../data/siteData/NHCsite_metadata.csv")
Qdat <- read_csv("../data/siteData/interpolatedQ_allsites.csv")
filelist <- list.files("../data/metabolism/raw/rds")

# get rid of MC751 for now
filelist <- filelist[-2]

#look at required inputs for a bayesian model in stream Metabolizer
metab_inputs(type="bayes", input="data")

for(i in 1:length(filelist)){
  dat <- readRDS(paste0("../data/metabolism/raw/rds/",filelist[i]))
  sitename <- dat$data$site[1]
  lat <- sites$latitude[sites$sitecode==sitename]
  lon <- sites$longitude[sites$sitecode==sitename]
  dat <- dat$data
  dat$value[dat$flagtype=="Bad Data"] <- NA
  dat <- select(dat,DateTime_UTC, value, variable )%>%
    pivot_wider(names_from = variable, values_from = value)
  
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
  dat$DO.sat <- calc_DO_sat(dat$WaterTemp_C, dat$AirPres_kPa)
  dat$depth <- calc_depth(dat$Discharge_m3s)
  dat$light <- calc_light(dat$solar.time, latitude=lat,longitude=lon)
  
  # select and renames columns for streamMetabolizer
  metab_dat <- select(dat, solar.time, DO.obs=DO_mgL, DO.sat, depth, temp.water=WaterTemp_C, light, discharge=Discharge_m3s )
  metab_dat$DO.obs <- na.approx(metab_dat$DO.obs, maxgap=12)
  metab_dat$DO.sat <- na.approx(metab_dat$DO.sat, maxgap=12)
  metab_dat$depth <- na.approx(metab_dat$depth, maxgap=12, na.rm=FALSE)
  metab_dat$temp.water <- na.approx(metab_dat$temp.water, maxgap=12)
  metab_dat$discharge <- na.approx(metab_dat$discharge, maxgap=12, na.rm=FALSE)
  write_csv(metab_dat, paste0("../data/metabolism/processed/",sitename,".csv" ))
}

getSPdata <- function(site, start_date, end_date, variables, airP, sensorOffsets){
  site_code <- paste('NC',site,sep='_')
  dat <- request_data(site_code, start_date, end_date, variables)
  # remove bad data
  w <- which(dat$data$flagtype== "Bad Data")
  dat$data$value[w]<-NA
  dat<- dat$data[,c(1,4,5)]%>% spread(key = "variable",value = "value")
  # attach air pressure data to calculate percent sat and depth
  dat <- full_join(airP, dat, by = "DateTime_UTC")

  # Use function from Stream Metabolizer to calculate percent saturation
  if(!("satDO_mgL" %in% colnames(dat))){
    dat$satDO_mgL <- calc_DO_sat(dat$WaterTemp_C, dat$AirPres_kPa*10)
    dat$persatDO <- dat$DO_mgL/dat$satDO_mgL
  }

  # Determine if level data is already available, if not calculate it from
  # Water Pressure. Then add sensor offset values
  if("WaterPres_kPa" %in% colnames(dat)){
    dat$Level_m <- (dat$WaterPres_kPa-dat$AirPres_kPa)*0.10197 +# m of water per kPa
      sensorOffsets$offset.m[which(sensorOffsets$site==site)]
  }


  return(dat)
}

# create a dataframe for all fo the data together:
alldat <- data.frame()

for(site in NHCsites2018){
  dat <- getSPdata(site, start_date, end_date, variables, airP, sensorOffsets)
  write.csv(dat, file = paste0("data/raw/",site,".csv"),row.names = F)
  dat$site <- rep(site, nrow(dat))
  alldat <- bind_rows(alldat, dat)
}

write.csv(alldat, file = "data/raw/allSites_20180706.csv", row.names = F)




#Create Time series object of just the DO values to plot with DY graphs
# DOdat <- data.frame(alldat[,c(1,2,
#                               grep("*.DO_mgL",colnames(alldat)),
#                               grep("*.WaterTemp_C", colnames(alldat)),
#                               grep("*.Level_m",colnames(alldat)),
#                               grep("*.DO.sat", colnames(alldat)))])
#
#
#
# write.csv(DOdat, file = "data/allsites_NHC2018DOdata.csv", row.names=F)
#
# #DOdat.xts <- xts(alldat[,grep("*.DO_mgL", colnames(alldat))],order.by=alldat[,1])
#dygraph(DOdat.xts) %>% dyRangeSelector()


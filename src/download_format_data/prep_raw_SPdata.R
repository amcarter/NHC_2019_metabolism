# Process raw NHC datafiles:
# files downloaded from SP portal
# prep for running metabolism models

library(streamMetabolizer)
library(LakeMetabolizer)
library(lubridate)
library(tidyverse)
library(dygraphs)
library(xts)

setwd("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism/data/")
source("../src/helpers.R")

# 1. load raw data files and metadata ####
sites <- read_csv("siteData/NHCsite_metadata.csv")
Qdat <- read_csv("rating_curves/interpolatedQ_allsites_modified.csv", guess_max = 10000)

filelist <- list.files("metabolism/raw")

# get rid of MC751 and Mud for now
# filelist <- filelist[-c(2,3)]

#look at required inputs for a bayesian model in stream Metabolizer
# metab_inputs(type="bayes", input="data")

# 2. Pair with interpolated Q ####
get_Q <- function(dat, Qdat){
  Qname <- paste(dat$site[1], "Q", sep=".")
  Q <- Qdat %>% select (DateTime_UTC, AirPres_kPa,
                        discharge = all_of(Qname)) %>%
    right_join(dat, by="DateTime_UTC") %>%
    arrange(DateTime_UTC)
  
  return(Q)
}

# 3. Estimate depth based on Leopold and Maddock/Raymond.  ####
#     12/18/2020 still need to calibrate site specific parameters
DQ <- data.frame(sitename = c("NHC","PM","CBP","WB","WBP","PWC","UNHC"),
                 c_m = rep(0.409, 7),   # depth at unit discharge
                 f = rep(.294, 7))      # exponent in depth discharge relation
# coefficients are default values for Leopold and Maddock 1953 equation D=cQ^f
# defined in Raymond 2012.

# 4. Calculate Level data from water pressure ####
# Depth = pressure_Pa/density/acelleration due to gravity = 
#   P (Pa) = ro (kg/m3) * gravity (m/s2) * depth (m)
# density is temperature dependent, this equation in lake metabolizer
# accounts for that based on Martin and McCutcheon 1999
calc_water_level<- function(dat, sites){
  
  sensor_offset <- sites$sensor_offset_m[sites$sitecode==dat$site[1]]
  ro <- water.density(dat$WaterTemp_C)
  pressure_Pa <- (dat$WaterPres_kPa-dat$AirPres_kPa)*1000
  level_m <- sensor_offset + pressure_Pa/(ro * 9.8)
  
  return(level_m)
}


# 5. Function to prepare datafile for drift correction ####

prep_file <- function(filename, sites, Qdat, DQ){
  dat <- read_csv(paste0("metabolism/raw/",filename), guess_max = 10000)
  lat <- sites$latitude[sites$sitecode == dat$site[1]]
  lon <- sites$longitude[sites$sitecode == dat$site[1]]
  
  if("AirPres_kPa"%in% colnames(dat)){dat<- select(dat, -AirPres_kPa)}
  
  # remove leading and ending NAs
  w <- which(!is.na(dat$DO_mgL))
  dat <- dat[min(w):max(w), ]
  
  #Load discharge data
  dat <- get_Q(dat, Qdat) 
  
  #Convert datetime to solar time 
  dat$DateTime_EST <- with_tz(dat$DateTime_UTC, tz="EST") # convert to EST timezone
  dat$solar.time <- calc_solar_time(dat$DateTime_EST, longitude=lon)
  
  # Calculate DO saturation, depth, light
  dat$AirPres_mbar <- dat$AirPres_kPa*10
  dat$DO.sat <- calc_DO_sat(dat$WaterTemp_C, dat$AirPres_mbar,
                            salinity.water = 0, model = "garcia-benson")
  dat$light <- calc_light(dat$solar.time, latitude=lat,longitude=lon)
  
  # This number is probably not very representative and 
  # needs to be changed based on field data
  dat$depth <- calc_depth(dat$discharge, 
                          c = DQ$c_m[DQ$sitename == dat$site[1]],
                          f = DQ$f[DQ$sitename == dat$site[1]])
  
  # Calculate Level data from water pressure
  dat$level_m <- calc_water_level(dat, sites)
  
  # rename variables needed for metabolism model
  dat <- dat %>% 
    select(-DateTime_EST, -AirPres_mbar, -AirPres_kPa, 
           -WaterPres_kPa) %>%
    rename(DO.obs = DO_mgL, 
           temp.water = WaterTemp_C)
  
  return(dat)
}

nhc_level <- read_csv("rating_curves/NHC_UNHC_Q.csv", guess_max = 10000)
for(i in 1:length(filelist)){
  
  filename <- filelist[i]
  dat <- prep_file(filename, sites, Qdat, DQ)
  if(dat$site[1] == "NHC"){
  dat <- nhc_level %>%
    select(DateTime_UTC, NHC_level) %>%
    right_join(dat) %>%
    mutate(level_m = NHC_level) %>%
    select(-NHC_level ) %>%
    arrange(DateTime_UTC)
  }  
  if(dat$site[1] == "UNHC"){
  dat <- nhc_level %>%
    select(DateTime_UTC, UNHC_level) %>%
    right_join(dat) %>%
    mutate(level_m = UNHC_level) %>%
    select(-UNHC_level ) %>%
    arrange(DateTime_UTC)
  }  
  
  write_csv(dat, paste0("metabolism/processed/", dat$site[1], ".csv"))
}

# 6. Drift correction functions based on YSI data ####
# load ysi data
ysi <- read_csv("water_chemistry/nhc_ysi_data.csv") %>%
  mutate(DateTime_EST = round_date(ymd_hms(paste(date, time), tz = "EST"), 
                                   unit = '15 minutes')) 
ysi <- read_csv("water_chemistry/sp_nhc_ysi_data.csv") %>%
      bind_rows(ysi) %>%
      mutate(ysi_level = waterdepth_cm/100) %>% 
      select(site, DateTime_EST, time, ysi_temp = watertemp_C, 
             ysi_DO = DO_mgL, ysi_level)

correct_file <- function(dat, ysi){
  
  ysid <- ysi %>% 
    filter(site == dat$site[1]) %>%
    right_join(dat) %>%
    arrange(DateTime_EST)
  
  yd <- ysid %>%
    select(DateTime_EST, DO.obs, ysi_DO,
           temp.water, ysi_temp,
           level_m, ysi_level) %>%
    mutate(DO.obs = na.approx(ifelse(!is.na(ysi_DO), NA, DO.obs),
                              na.rm = F, maxgap = 12), 
           temp.water = na.approx(ifelse(!is.na(ysi_temp), NA, temp.water),
                                  na.rm = F, maxgap = 12),
           level_m = na.approx(ifelse(!is.na(ysi_temp), NA, level_m),
                               na.rm = F, maxgap = 12)) 
  # DO <- yd[,-1] %>%
  #   xts(order.by = yd$DateTime_EST)
  # 
  # dygraph(DO) %>% dyRangeSelector()
  dat$DO.corr <- drift_correct(yd, "DO.obs", "ysi_DO")
  dat$T.corr <- drift_correct(yd, "temp.water", "ysi_temp")
  
  return(dat)
}

# 7. Run for files and check output ####
# for(i in 1:length(filelist)){

i=2


  filename <- filelist[i]
  dat <- prep_file(filename, sites, Qdat, DQ)
  dat <- correct_file(dat, ysi)
# ysi %>% filter(site == dat$site[1]) %>% 
#   right_join(dat) %>%
#   arrange(DateTime_EST) %>%
#   select(level_m, ysi_level) %>%
#   xts(order.by = dat$DateTime_EST) %>%
#   dygraph() %>%
#   dyRangeSelector()
  # check corrections  
  cor <- ysi %>% 
    filter(site == dat$site[1]) %>%
    right_join(dat) %>%
    arrange(DateTime_EST)
  
  corDO <- cor %>% select(DO.obs, DO.corr, ysi_DO) %>%
    xts(order.by = dat$DateTime_EST)
  show(dygraph(corDO) %>% dyRangeSelector())
  
  corT <- cor %>%
    select(temp.water, T.corr, ysi_temp) %>%
    xts(order.by = dat$DateTime_EST)
  show(dygraph(corT) %>% dyRangeSelector())
  
  df <- dat %>%
    mutate(DO.obs = DO.corr,
           temp.water = T.corr) %>%
    select(-DO.corr, -T.corr)
  
  write_csv(df, paste0("metabolism/processed/drift_corrected/", dat$site[1],"_dfc.csv"))


# Correct level data from NHC sites
library(LakeMetabolizer)
# 1. Setup for working through sites manually ####
# this needs to be commented out unless it's in use so this file can be sourced
# setwd('C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/hall_50yl2/')
# source('NHC_2019_metabolism/src/helpers.R')
# library(lubridate)
# library(tidyverse)
# library(zoo)
# library(xts)
# library(dygraphs)
# 
# sites <- read_csv("NHC_2019_metabolism/data/siteData/NHCsite_metadata.csv") %>%
#   slice(c(1:7))

plot_pres <- function(NHC, waterpres = "level_m", airpres = "waterdepth_m",
                      extra = NA){
  if(is.na(extra)){
    NHC %>% select(all_of(waterpres), all_of(airpres)) %>%
      xts(order.by = NHC$DateTime_UTC) %>%
      dygraph() %>%
      dyRangeSelector()
  } else {
    NHC %>% select(all_of(waterpres), all_of(airpres), all_of(extra)) %>%
      xts(order.by = NHC$DateTime_UTC) %>%
      dygraph() %>%
      dyRangeSelector()
  }
}

# 2. Load Data ####
# get airpressure data:
# NOAA_airpres <- StreamPULSE:::FindandCollect_airpres(sites$latitude[1], 
#                                                      sites$longitude[1],
#                                                      ymd_hms("2016-07-14 00:00:00"), 
#                                                      ymd_hms("2021-01-01 00:00:00"))
# NHC <- read_csv("NHC_2019_metabolism/data/metabolism/raw/NHC_2021-01-12.csv") %>%
#   select(DateTime_UTC, AirPres_kPa) %>%
#   left_join(NOAA_airpres)
# NOAA_airpres <- NHC %>%
#   mutate(air_kPa = air_kPa - 
#            mean(NHC$air_kPa - NHC$AirPres_kPa, na.rm = T)) %>%
#   select(-AirPres_kPa)
# write_csv(NOAA_airpres, "NHC_2019_metabolism/data/siteData/NOAA_airpres.csv")
NOAA_airpres <- read_csv("NHC_2019_metabolism/data/siteData/NOAA_airpres.csv")

# load field notes
ysi <- read_csv("NHC_2019_metabolism/data/water_chemistry/all_nhc_ysi_data.csv") %>% 
  filter(!is.na(Date)) %>%
  mutate(waterdepth_m = waterdepth_cm/100) %>%
  select(site, DateTime_UTC, waterdepth_m, notes) 

# load data files
filelist <- list.files("NHC_2019_metabolism/data/metabolism/raw")
nhc <- read_csv("NHC_2019_metabolism/data/metabolism/corrected_level/NHC_lvl.csv") %>%
  select(DateTime_UTC, level_nhc = level_m)

for(f in 1:nrow(sites)){
ff <- filelist[grep(paste0('^', sites$sitecode[f], '_'), filelist)]
dat <- read_csv(paste0("NHC_2019_metabolism/data/metabolism/raw/",ff), 
                guess_max = 100000) %>%
  left_join(NOAA_airpres, by = "DateTime_UTC") 
if("AirTemp_C" %in% colnames(dat)){ 
  dat <- select(dat, -AirTemp_C) }
if("AirPres_kPa" %in% colnames(dat)){ 
  dat <- select(dat, -AirPres_kPa) }

# Calculate depth from water pressure and add sensor offset 
# Depth = pressure_Pa = kg/ms2/(density_kg/m3*gravity_m/s2)
dat <- dat %>% 
  mutate(pressure_Pa = (WaterPres_kPa - air_kPa) * 1000,
         level_m = pressure_Pa/(water.density(WaterTemp_C) * 9.8) + 
           sites$sensor_offset_m[f]) %>%
  select(-pressure_Pa) %>%
  left_join(ysi[ysi$site == dat$site[1],])
  
  if(dat$site[1] == "NHC"){
    dat$level_m[dat$level_m < 0.42] <- NA
  }
  if(dat$site[1] == "PM"){
    dat <- dat %>%
      left_join(nhc)
    plot_pres(dat, "level_m", "waterdepth_m", "level_nhc")
    dat$level_m[dat$level_m < 0.57] <- NA
    tmp <- data.frame(starts = c(which(dat$DateTime_UTC ==
                                         ymd_hms("2020-02-12 14:00:00"))
    ))
    tmp$stops <- tmp$starts +1
    
    gaps <- rle_custom(is.na(dat$level_m)) 
    if(gaps$values[1] ==1) { gaps <- gaps[-1,]}
    if(gaps$values[nrow(gaps)] == 1) { gaps <- gaps[-nrow(gaps),]}
    tmp <- gaps %>%
      filter(values == 1) %>%
      bind_rows(tmp) %>%
      mutate(datetime = as.Date(dat$DateTime_UTC[starts]),
             starts = starts - 1,
             stops = stops + 1,
             lvl_start = dat$level_m[starts],
             lvl_stop = dat$level_m[stops],
             jump = lvl_stop - lvl_start) %>%
      slice(c(3:8, 11, 12, 15)) %>%
      arrange(starts) 
    n <- nrow(dat)
    # dat$level_m1 -> dat$level_m
    dat$level_m1 <- dat$level_m
    for(i in 1:(nrow(tmp))){
      dat$level_m[tmp$stops[i]:n] <-  dat$level_m[tmp$stops[i]:n] - tmp$jump[i]
    }
    d = ymd_hms("2019-05-01 00:00:00")
    nhc_diff <- dat %>%
      filter(DateTime_UTC < d + 60*60*24*30*4) %>%
      mutate(diff = level_nhc - level_m,
             before = ifelse(DateTime_UTC > d, FALSE, TRUE)) %>%
      group_by(before) %>%
      summarize(diff = median(diff, na.rm = T))
      
    dat$level_m[dat$DateTime_UTC < d] <- dat$level_m[dat$DateTime_UTC < d] +
      nhc_diff$diff[2] - nhc_diff$diff[1]
    
    dat$level_m <- na.approx(dat$level_m,na.rm = F)
    dat$level_m <- drift_correct(dat, "level_m", "waterdepth_m")
    # plot_pres(dat, "level_nhc", "waterdepth_m", "level_m")
    dat <- select(dat, -level_nhc)
  }
  if(dat$site[1] == "CBP"){
    dat <- dat %>%
      left_join(nhc)
    # plot_pres(dat, "level_m", "waterdepth_m", "level_nhc")
    dat$level_m[dat$level_m < 0.28] <- NA
    tmp <- data.frame(starts = c(which(dat$DateTime_UTC ==
                                         ymd_hms("2019-03-19 19:45:00")),
                                 which(dat$DateTime_UTC ==
                                         ymd_hms("2019-07-02 15:15:00")),
                                 which(dat$DateTime_UTC ==
                                         ymd_hms("2019-08-02 16:00:00")),
                                 which(dat$DateTime_UTC ==
                                         ymd_hms("2020-02-18 15:45:00"))
    ))
    tmp$stops <- tmp$starts +1
    
    gaps <- rle_custom(is.na(dat$level_m)) 
    if(gaps$values[1] ==1) { gaps <- gaps[-1,]}
    if(gaps$values[nrow(gaps)] == 1) { gaps <- gaps[-nrow(gaps),]}
    tmp <- gaps %>%
      filter(values == 1) %>%
      slice(c(1,2,5,7,11:15)) %>%
      bind_rows(tmp) %>%
      mutate(datetime = as.Date(dat$DateTime_UTC[starts]),
             starts = starts - 1,
             stops = stops + 1,
             lvl_start = dat$level_m[starts],
             lvl_stop = dat$level_m[stops],
             jump = lvl_stop - lvl_start) %>%
      arrange(starts) 
    n <- nrow(dat)
    # dat$level_m1 -> dat$level_m
    dat$level_m1 <- dat$level_m
    for(i in 1:(nrow(tmp))){
      dat$level_m[tmp$stops[i]:n] <-  dat$level_m[tmp$stops[i]:n] - tmp$jump[i]
    }
    d = ymd_hms("2019-05-10 00:00:00")
    nhc_diff <- dat %>%
      filter(DateTime_UTC < d + 60*60*24*28*4) %>%
      mutate(diff = level_nhc - level_m,
             before = ifelse(DateTime_UTC > d, FALSE, TRUE)) %>%
      group_by(before) %>%
      summarize(diff = median(diff, na.rm = T))
      
    dat$level_m[dat$DateTime_UTC < d] <- dat$level_m[dat$DateTime_UTC < d] +
      nhc_diff$diff[2] - nhc_diff$diff[1]
    
    dat$level_m <- na.approx(dat$level_m,na.rm = F)
    dat$level_m <- drift_correct(dat, "level_m", "waterdepth_m")
    # plot_pres(dat, "level_m1", "waterdepth_m", "level_d")
    dat <- select(dat, -level_nhc)
  }
  if(dat$site[1] == "WB"){
    dat <- dat %>%
      left_join(nhc)
    # plot_pres(dat, "level_m", "waterdepth_m", "level_nhc")
    dat$level_m[dat$level_m < 0.28] <- NA
    tmp <- data.frame(starts = c(which(dat$DateTime_UTC ==
                                         ymd_hms("2019-07-19 15:45:00"))
    ))
    tmp$stops <- tmp$starts +1
    
    gaps <- rle_custom(is.na(dat$level_m)) 
    if(gaps$values[1] ==1) { gaps <- gaps[-1,]}
    if(gaps$values[nrow(gaps)] == 1) { gaps <- gaps[-nrow(gaps),]}
    tmp <- gaps %>%
      filter(values == 1) %>%
      slice(c(2,5:9,13, 14)) %>%
      bind_rows(tmp) %>%
      mutate(datetime = as.Date(dat$DateTime_UTC[starts]),
             starts = starts - 1,
             stops = stops + 1,
             lvl_start = dat$level_m[starts],
             lvl_stop = dat$level_m[stops],
             jump = lvl_stop - lvl_start) %>%
      arrange(starts) 
    n <- nrow(dat)
    # dat$level_m1 -> dat$level_m
    dat$level_m1 <- dat$level_m
    for(i in 1:(nrow(tmp))){
      dat$level_m[tmp$stops[i]:n] <-  dat$level_m[tmp$stops[i]:n] - tmp$jump[i]
    }
    
    
    dat$level_m <- na.approx(dat$level_m,na.rm = F)
    dat$level_m <- drift_correct(dat, "level_m", "waterdepth_m")
    # plot_pres(dat, "level_nhc", "waterdepth_m", "level_d")
    d = ymd_hms("2019-05-19 00:00:00")
    nhc_diff <- dat %>%
      filter(DateTime_UTC < d + 60*60*24*25*4) %>%
      mutate(diff = level_nhc - level_m,
             before = ifelse(DateTime_UTC > d, FALSE, TRUE)) %>%
      group_by(before) %>%
      summarize(diff = median(diff, na.rm = T))
      
    dat$level_m[dat$DateTime_UTC < d] <- dat$level_m[dat$DateTime_UTC < d] +
      nhc_diff$diff[2] - nhc_diff$diff[1]
    dat <- select(dat, -level_nhc)
  }
  if(dat$site[1] == "WBP"){
    dat <- dat %>%
      left_join(nhc)
    # plot_pres(dat, "level_m", "waterdepth_m", "level_nhc")
    dat$level_m[dat$level_m < 0.66] <- NA
    tmp <- data.frame(starts = c(which(dat$DateTime_UTC ==
                                         ymd_hms("2019-05-10 16:30:00")),
                                 which(dat$DateTime_UTC ==
                                         ymd_hms("2019-11-20 18:00:00")),
                                 which(dat$DateTime_UTC ==
                                         ymd_hms("2019-12-12 19:15:00"))
    ))
    tmp$stops <- tmp$starts +1
    
    gaps <- rle_custom(is.na(dat$level_m)) 
    if(gaps$values[1] ==1) { gaps <- gaps[-1,]}
    if(gaps$values[nrow(gaps)] == 1) { gaps <- gaps[-nrow(gaps),]}
    tmp <- gaps %>%
      filter(values == 1) %>%
      slice(15) %>%
      bind_rows(tmp) %>%
      mutate(datetime = as.Date(dat$DateTime_UTC[starts]),
             starts = starts - 1,
             stops = stops + 1,
             lvl_start = dat$level_m[starts],
             lvl_stop = dat$level_m[stops],
             jump = lvl_stop - lvl_start) %>%
      arrange(starts) 
    n <- nrow(dat)
    # dat$level_m1 -> dat$level_m
    dat$level_m1 <- dat$level_m
    for(i in 1:(nrow(tmp))){
      dat$level_m[tmp$stops[i]:n] <-  dat$level_m[tmp$stops[i]:n] - tmp$jump[i]
    }
    
    # plot_pres(dat, "level_nhc", "waterdepth_m", "level_d")
    d = ymd_hms("2019-05-08 00:00:00")
    d1 = ymd_hms("2020-06-18 00:00:00")
    nhc_diff <- dat %>%
      filter(DateTime_UTC < d + 60*60*24*30*2 | DateTime_UTC > d1) %>%
      mutate(diff = level_nhc - level_m,
             before = case_when(DateTime_UTC < d ~ 1,
                                DateTime_UTC > d1 ~ 3,
                                TRUE ~ 2)) %>%
      group_by(before) %>%
      summarize(diff = median(diff, na.rm = T))
      
    dat$level_m[dat$DateTime_UTC < d] <- dat$level_m[dat$DateTime_UTC < d] +
      nhc_diff$diff[1] - nhc_diff$diff[2]
    dat$level_m[dat$DateTime_UTC > d1] <- dat$level_m[dat$DateTime_UTC > d1] +
      nhc_diff$diff[3] - nhc_diff$diff[2]
    dat$level_m <- na.approx(dat$level_m,na.rm = F)
    dat$level_m <- drift_correct(dat, "level_m", "waterdepth_m")
    dat <- select(dat, -level_nhc)
  }
  if(dat$site[1] == "PWC"){
    dat <- dat %>%
      left_join(nhc)
    # # plot_pres(dat, "level_m", "waterdepth_m", "level_nhc")
    dat$level_m[dat$level_m < 0.69] <- NA
    tmp <- data.frame(starts = c(which(dat$DateTime_UTC ==
                                         ymd_hms("2019-05-30 21:00:00")),
                                 which(dat$DateTime_UTC ==
                                         ymd_hms("2019-07-19 17:00:00"))
    ))
    tmp$stops <- tmp$starts + c( 1,2)
    
    gaps <- rle_custom(is.na(dat$level_m)) 
    if(gaps$values[1] ==1) { gaps <- gaps[-1,]}
    if(gaps$values[nrow(gaps)] == 1) { gaps <- gaps[-nrow(gaps),]}
    tmp <- gaps %>%
      filter(values == 1) %>%
      slice(1,3,5) %>%
      bind_rows(tmp) %>%
      mutate(datetime = as.Date(dat$DateTime_UTC[starts]),
             starts = starts - 1,
             stops = stops + 1,
             lvl_start = dat$level_m[starts],
             lvl_stop = dat$level_m[stops],
             jump = lvl_stop - lvl_start) %>%
      arrange(starts) 
    n <- nrow(dat)
    # dat$level_m1 -> dat$level_m
    dat$level_m1 <- dat$level_m
    for(i in 1:(nrow(tmp))){
      dat$level_m[tmp$stops[i]:n] <-  dat$level_m[tmp$stops[i]:n] - tmp$jump[i]
    }
    
    # plot_pres(dat, "level_nhc", "waterdepth_m", "level_d")
    
    dat$level_m <- na.approx(dat$level_m,na.rm = F)
    dat$level_m <- drift_correct(dat, "level_m", "waterdepth_m")
    dat <- select(dat, -level_nhc)
  }

  # plot_pres(dat, "level_m", 'waterdepth_m', "level_m1")
  dat <- dat %>%
    mutate(#level_m = case_when(is.na(level_m1) ~ NA_real_,
           #                    TRUE ~ level_m),
           level_m = na.approx(level_m, maxgap = 96*3, na.rm = F),
           site = dat$site[1]) %>%
    # select(-level_m1, -WaterPres_kPa, -waterdepth_m, -notes) %>%
    group_by(DateTime_UTC, site) %>%
    summarize_all(mean, na.rm = T)

write_csv(dat, paste0("NHC_2019_metabolism/data/metabolism/corrected_level/", 
                      dat$site[1], "_lvl.csv"))

}

# 3. compile all levels ####
filelist <- list.files("NHC_2019_metabolism/data/metabolism/corrected_level/")
dd <- data.frame()
for(f in 1:length(filelist)){
  d <- read_csv(paste0("NHC_2019_metabolism/data/metabolism/corrected_level/", 
                       filelist[f])) %>% 
  select(DateTime_UTC, level_m, site)
  dd <- bind_rows(dd, d)
}

write_csv(dd, "NHC_2019_metabolism/data/rating_curves/all_sites_level_corrected.csv")

detach('package:LakeMetabolizer', unload = TRUE)

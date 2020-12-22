# code for handling field ysi data to 
# 12/21/2020
# A Carter
#
# 1. format ysi data for uploading to the portal
# 2. format SP ysi data so it is easier to use, pairing it with my data
# this references the datasheet I use to enter values after I go in the field
# it updates the correctly formated file with only the new data

library(tidyverse)
library(lubridate)
library(streamMetabolizer)

setwd(metab_projdir)
sites <- read_csv(file="data/siteData/NHCsite_metadata.csv")

# 1. format ysi data for upload ####

ysi <- read_csv("data/water_chemistry/NHC_YSI_data.csv")

# 2. format SP ysi data ####

sp <- read_csv("data/water_chemistry/StreamPULSE Probe Measurements.csv") %>%
  mutate(date = as.Date(Date, format = "%m.%d.%Y"),
         DateTime_EST = round_date(ymd_hms(paste(date, Time), tz = "EST"),
                                   unit = "15 minutes")) %>%
  select(site = Site, DateTime_EST, date, time = Time,
         watertemp_C = 'Temp (DegC)',
         DO_mgL = 'DO (mg/L)', 
         spc_uscm = 'Cond (uS/cm)',
         pH = 'pH (units)',
         NTU = 'Turb (NTU)',
         waterdepth_cm = 'Water Depth (cm)',
         stage_m = 'Stage at sensors',
         notes = Notes) %>%
  filter(site %in% c("NHC", "UNHC")) %>%
  mutate(waterdepth_cm = case_when(!is.na(waterdepth_cm) ~ waterdepth_cm,
                                   is.na(waterdepth_cm) ~ stage_m *100)) %>%
  select(-stage_m)

# add air pressure data
NOAA_airpres <- StreamPULSE:::FindandCollect_airpres(lat = sites$latitude[1],
                                                     long = sites$longitude[1],
                                                     min(sp$DateTime_EST, na.rm = T),
                                                     max(sp$DateTime_EST, na.rm = T))
sp <- sp %>% 
  mutate(DateTime_UTC = with_tz(DateTime_EST, tz = "UTC")) %>%
  left_join(NOAA_airpres) %>%
  mutate(airpres_mmHg = air_kPa * 7.50062,
         DO_sat = calc_DO_sat(temp.water = watertemp_C,
                              pressure.air = air_kPa * 10),
         DO_pctsat = DO_mgL / DO_sat) %>%
  select(-air_kPa, -DO_sat, -DateTime_UTC)
  
# save formatted SP file for NHC sites
write_csv(sp, "data/water_chemistry/sp_nhc_ysi_data.csv")
nhc <- sp %>%
  filter(site == "NHC")
unhc <- sp %>%
  filter(site == "UNHC")

# plot(nhc$date, nhc$stage_m)
# points(nhc$date, nhc$waterdepth_cm/100)
plot(unhc$date, unhc$waterdepth_cm)


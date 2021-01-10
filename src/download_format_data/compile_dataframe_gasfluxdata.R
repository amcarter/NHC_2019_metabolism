# Prep data frame with variables for NHC ghg analysis
# 2021 01 03

library(tidyverse)
library(lubridate)
library(zoo)

setwd("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism")

# load ghg data
ghg <- read_csv("data/gas_data/NHC_2019-2020_processed_GHGdata.csv") %>%
  rename(watertemp_C.ysi = watertemp_C, 
         airpres_mmHg.ysi = airpres_mmHg, 
         pH.ysi = pH)
# Combine duplicates examine leaky measures
# To Do: Better quantify measurement error 
# (GC error, how to get obs error from duplicates)
# Account for error where mass/temperature is estimated

ghg <- ghg %>%
  mutate(datetime = round_date(ymd_hms(paste(date, time)), 
         unit = "15 minute")) %>%
  select(-Notes, -date, -time) %>%
  group_by(site, datetime) %>%
  summarize(across(everything(), mean, na.rm = T)) %>%
  ungroup() %>%
  mutate(date = as.Date(datetime))


 
nuts <- read_csv("data/water_chemistry/water_chemistry_2019-2020_compiled.csv") %>%
  select(-sample_name, -time) %>%
  mutate(site = toupper(site),
         br_mgl = ifelse(br_mgl == "<0.03", 0.015, as.numeric(br_mgl)),
         nh4n_mgl = ifelse(nh4n_mgl == "<0.01", 0.005, as.numeric(nh4n_mgl)),
         po4p_mgl = ifelse(po4p_mgl == "<0.01", 0.005, as.numeric(po4p_mgl)))
ysi <- read_csv("data/water_chemistry/all_nhc_ysi_data.csv") %>%
  select(site, date, spc_uScm.ysi = spc_uscm) %>%
  mutate(date = as.Date(date, format = "%m/%d/%Y"))
site_dat <- read_csv("data/siteData/NHCsite_metadata.csv") %>%
  slice(c(1:5,7:9)) %>%
  select(site = sitecode, latitude, longitude, CRS, 
         distance_upstream_m = distance_m,
         width_march_m = width_mar_m,
         ws_area_km2 = ws_area.km2,
         slope_nhd = slope, habitat)

met <- readRDS("data/metabolism/compiled/raymond_met.rds")$preds %>%
  as_tibble() %>%
  select(date, site, GPP, ER, K600, 
         discharge_m3s = discharge.daily, watertemp_C = temp.water) %>%
  mutate(site = toupper(site),
         across(c(-site, -date), na.approx, na.rm = F)) 

air <- read_csv("data/siteData/interpolatedQ_allsites.csv") %>%
  select(DateTime_UTC, AirPres_kPa, AirTemp_C) %>%
  group_by(date = as.Date(with_tz(DateTime_UTC, tz = "EST"))) %>%
  summarize_all(mean, na.rm = T)

sites <- site_dat[1:6, 1, drop = T]
raw_dat <- data.frame()
for(site in sites){
  dat <- read_csv(paste0("data/metabolism/processed/", site, ".csv"), 
                  guess_max = 10000)
  dat$site <- site
  raw_dat = bind_rows(raw_dat, dat)
}

raw_daily <- raw_dat %>%
  as.tibble() %>%
  mutate(spc_uScm = ifelse(!is.na(SpecCond_uScm),
                           SpecCond_uScm, SpecCond_mScm/1000),
         date = as.Date(with_tz(DateTime_UTC, tz = "EST"))) %>%
  select(date, site, level_m, site, spc_uScm, pH, DO.obs) %>%
  group_by(date, site) %>%
  summarize_all(mean, na.rm = T)

raw_dat <- raw_dat %>%
  as.tibble() %>%
  mutate(spc_uScm = ifelse(!is.na(SpecCond_uScm),
                           SpecCond_uScm, SpecCond_mScm/1000),
         datetime = with_tz(DateTime_UTC, tz = "EST")) %>%
  select(datetime, site, level_m, site, spc_uScm, DO.obs) %>%
  mutate(across(c(-datetime, -site), na.approx, na.rm = F)) 

dat <- ghg %>% 
  left_join(ysi, by = c("date", "site")) %>%
  left_join(nuts, by = c("date", "site")) %>%
  left_join(site_dat, by = "site") %>%
  left_join(air, by = "date") %>%
  left_join(met, by = c("date", "site")) %>%
  left_join(raw_dat, by = c("datetime", "site")) %>%
  select(-error_pct, -DateTime_UTC)

# compile water chemistry ####
chem <- nuts %>%
  filter(site != "MC751") %>%
  # left_join(ysi, by = c("date", "site")) %>%
  left_join(site_dat, by = "site") %>%
  left_join(met, by = c("date", "site")) %>%
  left_join(raw_daily, by = c("date", "site")) %>%
  select(-width_march_m, -ws_area_km2, -habitat, -pH)
    
apply(chem, 2, function(x) sum(is.na(x)))/nrow(chem)
# compare ysi data:
ysicomp <- dat %>% 
  select(watertemp_C.ysi, watertemp_C, airpres_mmHg.ysi, AirPres_kPa, 
         pH.ysi, pH, spc_uScm.ysi, spc_uScm, airtemp_C.ysi, AirTemp_C)

dat <- dat %>%
  mutate(pH = pH.ysi,
         spc_uScm = ifelse(!is.na(spc_uScm), spc_uScm, spc_uScm.ysi),
         watertemp_C = watertemp_C.ysi) %>%
  select(-ends_with(".ysi"))

# calculate gas flux ####
library(marelac)
# K_gas = K600(Sc/600)^-0.5
tmp <- data.frame(watertemp_C = unique(dat$watertemp_C))
K <-  data.frame(gas_schmidt(tmp$watertemp_C, c("CH4","CO2", "N2O")))
K$watertemp_C <- tmp$watertemp_C
K <- dat %>%
  select(watertemp_C, K600) %>%
  left_join(K, by = "watertemp_C") %>%
  mutate(across(c(-watertemp_C, -K600), ~(./600)^(-0.5)*K600, 
                .names = "K_{col}")) %>%
  select(-CH4, -CO2, -N2O)
  
dat <- dat %>% 
  left_join(K, by = c("watertemp_C", "K600")) %>%
  mutate(CH4.flux_ugld = (CH4.ugL - CH4.sat) * K_CH4,
         CO2.flux_ugld = (CO2.ugL - CO2.sat) * K_CO2,
         N2O.flux_ugld = (N2O.ugL - N2O.sat) * K_N2O)

# still need to get average depth and width upstream
# write_csv(dat, "data/gas_data/compiled_ghg_dataframe_allvariables.csv")
dat <- read_csv("data/gas_data/compiled_ghg_dataframe_allvariables.csv")

# pare down to have fewer NA's
apply(dat, 2, function(x) sum(is.na(x)))/nrow(dat)
d2 <- dat %>%
  filter(site != "MC751", 
         !is.na(K600)) %>%
  select(-na_mgl, -mg_mgl, -ca_mgl, -br_mgl,
         -starts_with("K_"), -CRS, -latitude, -longitude, -pH, -spc_uScm,
         -nh4n_mgl, -po4p_mgl, -doc_mgl, -tdn_mgl, -cl_mgl,  -datetime, 
         -no3n_mgl, -so4_mgl, -AirTemp_C, -AirPres_kPa, -ends_with("sat"),
         -distance_upstream_m, -width_march_m, -ws_area_km2, -slope_nhd) 

apply(d2, 2, function(x) sum(is.na(x)))
# write_csv(d2, "data/gas_data/compiled_ghg_dataframe_short.csv")

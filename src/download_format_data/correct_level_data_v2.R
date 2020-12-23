# Correct level data from NHC sites and 
# calculate average depth as best as we can for now


# this process is idosincratic for each site and will be based
# on times that the sensor was moved or when there is bad data.

library(tidyverse)
library(xts)
library(dygraphs)
library(lubridate)
library(zoo)
setwd('C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism/data')
source('../src/helpers.R')
# 1. Load Data ####
# load field notes

ysi <- read_csv("water_chemistry/nhc_ysi_data.csv") %>% 
  filter(!is.na(Date)) %>%
  mutate(waterdepth_m = waterdepth_cm/100) %>%
  select(site, DateTime_UTC, waterdepth_m, notes) 

plot_pres <- function(NHC, waterpres = "level_m", airpres = "waterdepth_m"){
  NHC %>% select(all_of(waterpres), all_of(airpres)) %>%
    xts(order.by = NHC$DateTime_UTC) %>%
    dygraph() %>%
    dyRangeSelector()
}

# CBP ####
# cbp <- read_csv(paste0("metabolism/processed/CBP.csv")) %>%
#   select(DateTime_UTC, temp = temp.water, level_m, discharge) %>%
#   left_join(ysi[ysi$site == "CBP",]) 
# 
# plot_pres(cbp, "level_m", "waterdepth_m")
# # Fill in missing dates
# cbp <- data.frame(DateTime_UTC = seq(min(cbp$DateTime_UTC), 
#                                      max(cbp$DateTime_UTC),
#                                      by = '15 min')) %>%
#   left_join(cbp)
# dd <- ymd_hms("2020-02-12 16:45:00")
# w <- which(cbp$DateTime_UTC == dd)
# cbp$level_m[(w+1):(w+3)] <- NA
# # first, piece together places where chunks are clearly shifted
# # these jumps occur at: (dttm is last pt of data before jump)
# gaps <- rle_custom(is.na(cbp$level_m)) %>%
#   filter(values == 1) %>%
#   mutate(starts = starts - 1,
#          stops = stops + 1)
# 
# # these gaps look like snapping them together would be incorrect
# # gaps <- gaps %>% 
# #   mutate(lvl_start = cbp$level_m[starts],
# #          lvl_stop = cbp$level_m[stops],
# #          jump = lvl_stop - lvl_start) %>%
# #   filter(jump >= 0.01 | jump <= -.01) %>%# only adjust for gaps larger than 1 cm
# #   slice(c(-1,-2))
# n <- nrow(cbp)
# g <- nrow(gaps)
# cbp$level_m1 <- cbp$level_m
# cbp$level_m <- na.approx(cbp$level_m, maxgap = 50, na.rm = F)
# for(i in 1:(g+1)){
#   if(i == g+1){ end = n } else {end = gaps$starts[i] }
#   if(i == 1){ start = 1 } else {start = gaps$stops[i-1] }
#   delta <- cbp$waterdepth_m[start:end] - cbp$level_m[start:end]
#   if(sum(!is.na(delta)) > 0){
#     cbp$level_m[start:n] <- cbp$level_m[start:n] + mean(delta, na.rm = T)
#   }
# }  
# cbp$level_m[(w+1):(w+3)] <- NA
# dd <- ymd_hms(c("2019-10-29 14:00:00", "2019-08-21 19:45:00"))
# w <- which(cbp$DateTime_UTC %in% dd)
# cbp$level_m[w] <- NA
# plot_pres(cbp, "level_m", "waterdepth_m")
# # nhc$level_m <- na.approx(nhc$level_m, maxgap = 50, na.rm = F)
# # nhc$level_d <- drift_correct(nhc, "level_m", "waterdepth_m")
# 
# cbp <- cbp %>%
#   select(-level_m1) %>%
#   mutate(level_m = na.approx(level_m, maxgap = 3*96, na.rm = F))
# 
# plot(cbp$level_m, cbp$discharge, log = "xy")
# abline(v = 1.1)
# sum(cbp$level_m > 1.1, na.rm = T)/sum(!is.na(cbp$level_m))
# 
# dat <- read_csv("metabolism/processed/NHC.csv") %>% 
#   select(DateTime_UTC, nhc = level_m) %>%
#   right_join(cbp) %>%
#   arrange(DateTime_UTC)
# 
# plot(dat$nhc, dat$level_m)
# ggplot(dat, aes(nhc, level_m, col = DateTime_UTC)) +
#   geom_point()
# 
# write_csv(cbp, "metabolism/processed/level_v2/CBP_lvl.csv")


# PM  ####
dat <- read_csv(paste0("metabolism/processed/PM.csv")) %>%
  select(DateTime_UTC, temp = temp.water, level_m, discharge) %>%
  left_join(ysi[ysi$site == "PM",]) 

plot_pres(dat, "level_m", "waterdepth_m")
dat$level_m[which(dat$level_m<.58)] <- NA
# Fill in missing dates
dat <- data.frame(DateTime_UTC = seq(min(dat$DateTime_UTC), 
                                     max(dat$DateTime_UTC),
                                     by = '15 min')) %>%
  left_join(dat)

# first, piece together places where chunks are clearly shifted
# these jumps occur at: (dttm is last pt of data before jump)
gaps <- rle_custom(is.na(dat$level_m)) %>%
  filter(values == 1) %>%
  mutate(starts = starts - 1,
         stops = stops + 1)

n <- nrow(dat)
g <- nrow(gaps)
dat$level_m1 <- dat$level_m
dat$level_m <- na.approx(dat$level_m, maxgap = 50, na.rm = F)
for(i in 1:(g+1)){
  if(i == g+1){ end = n } else {end = gaps$starts[i] }
  if(i == 1){ start = 1 } else {start = gaps$stops[i-1] }
  delta <- dat$waterdepth_m[start:end] - dat$level_m[start:end]
  if(sum(!is.na(delta)) > 0){
    dat$level_m[start:n] <- dat$level_m[start:n] + mean(delta, na.rm = T)
  }
}  
plot_pres(dat, "level_m", "waterdepth_m")

dat <- dat %>%
  select(-level_m1) %>%
  mutate(level_m = na.approx(level_m, maxgap = 3*96, na.rm = F))

plot(dat$level_m, dat$discharge)
write_csv(dat, "metabolism/processed/level_v2/PM_lvl.csv")

# WB  ####
dat <- read_csv(paste0("metabolism/processed/WB.csv")) %>%
  select(DateTime_UTC, temp = temp.water, level_m, discharge) %>%
  left_join(ysi[ysi$site == "WB",]) 

plot_pres(dat, "level_m", "waterdepth_m")
dat$level_m[which(dat$level_m<.29)] <- NA
# Fill in missing dates
dat <- data.frame(DateTime_UTC = seq(min(dat$DateTime_UTC), 
                                     max(dat$DateTime_UTC),
                                     by = '15 min')) %>%
  left_join(dat)

dat$level_m[34413]<- NA
# first, piece together places where chunks are clearly shifted
# these jumps occur at: (dttm is last pt of data before jump)
gaps <- rle_custom(is.na(dat$level_m)) %>%
  filter(values == 1) %>%
  mutate(starts = starts - 1,
         stops = stops + 1)
dat$waterdepth_m[1000] <- dat$level_m[1000] + .53
dat$waterdepth_m[5000] <- dat$level_m[5000] + .31
#dat$waterdepth_m[34415] <- dat$level_m[34415] + .19
n <- nrow(dat)
g <- nrow(gaps)
dat$level_m1 <- dat$level_m
dat$level_m <- na.approx(dat$level_m, maxgap = 50, na.rm = F)
for(i in 1:(g+1)){
  if(i == g+1){ end = n } else {end = gaps$starts[i] }
  if(i == 1){ start = 1 } else {start = gaps$stops[i-1] }
  delta <- dat$waterdepth_m[start:end] - dat$level_m[start:end]
  if(sum(!is.na(delta)) > 0){
    dat$level_m[start:n] <- dat$level_m[start:n] + mean(delta, na.rm = T)
  }
}  
plot_pres(dat, "level_m", "waterdepth_m")

dat <- dat %>%
  select(-level_m1) %>%
  mutate(level_m = na.approx(level_m, maxgap = 3*96, na.rm = F))

plot(dat$level_m, dat$discharge, log = "y")

write_csv(dat, "metabolism/processed/level_v2/WB_lvl.csv")

# WBP  ####
dat <- read_csv(paste0("metabolism/processed/WBP.csv")) %>%
  select(DateTime_UTC, temp = temp.water, level_m, discharge) %>%
  left_join(ysi[ysi$site == "WBP",]) 

plot_pres(dat, "level_m", "waterdepth_m")
dat$level_m[which(dat$level_m<.67)] <- NA
# Fill in missing dates
dat <- data.frame(DateTime_UTC = seq(min(dat$DateTime_UTC), 
                                     max(dat$DateTime_UTC),
                                     by = '15 min')) %>%
  left_join(dat)

# first, piece together places where chunks are clearly shifted
# these jumps occur at: (dttm is last pt of data before jump)
gaps <- rle_custom(is.na(dat$level_m)) %>%
  filter(values == 1) %>%
  mutate(starts = starts - 1,
         stops = stops + 1)
dat$waterdepth_m[100] <- dat$level_m[100] + .29

n <- nrow(dat)
g <- nrow(gaps)
dat$level_m1 <- dat$level_m
dat$level_m <- na.approx(dat$level_m, maxgap = 50, na.rm = F)
for(i in 1:(g+1)){
  if(i == g+1){ end = n } else {end = gaps$starts[i] }
  if(i == 1){ start = 1 } else {start = gaps$stops[i-1] }
  delta <- dat$waterdepth_m[start:end] - dat$level_m[start:end]
  if(sum(!is.na(delta)) > 0){
    dat$level_m[start:n] <- dat$level_m[start:n] + mean(delta, na.rm = T)
  }
}  
plot_pres(dat, "level_m", "waterdepth_m")
dat <- dat %>%
  select(-level_m1) %>%
  mutate(level_m = na.approx(level_m, maxgap = 3*96, na.rm = F))

plot(dat$level_m, dat$discharge, log = "y")

write_csv(dat, "metabolism/processed/level_v2/WBP_lvl.csv")

# PWC  ####
dat <- read_csv(paste0("metabolism/processed/PWC.csv")) %>%
  select(DateTime_UTC, temp = temp.water, level_m, discharge) %>%
  left_join(ysi[ysi$site == "PWC",]) 

plot_pres(dat, "level_m", "waterdepth_m")
dat$level_m[which(dat$level_m<.7)] <- NA
# Fill in missing dates
dat <- data.frame(DateTime_UTC = seq(min(dat$DateTime_UTC), 
                                     max(dat$DateTime_UTC),
                                     by = '15 min')) %>%
  left_join(dat)

# first, piece together places where chunks are clearly shifted
# these jumps occur at: (dttm is last pt of data before jump)
gaps <- rle_custom(is.na(dat$level_m)) %>%
  filter(values == 1) %>%
  mutate(starts = starts - 1,
         stops = stops + 1)
dat$waterdepth_m[100] <- dat$level_m[100] + .18
dat$waterdepth_m[5471] <- .9
dat$waterdepth_m[5481] <- .9

n <- nrow(dat)
g <- nrow(gaps)
dat$level_m1 <- dat$level_m
dat$level_m <- na.approx(dat$level_m, maxgap = 50, na.rm = F)
for(i in 1:(g+1)){
  if(i == g+1){ end = n } else {end = gaps$starts[i] }
  if(i == 1){ start = 1 } else {start = gaps$stops[i-1] }
  delta <- dat$waterdepth_m[start:end] - dat$level_m[start:end]
  if(sum(!is.na(delta)) > 0){
    dat$level_m[start:n] <- dat$level_m[start:n] + mean(delta, na.rm = T)
  }
}  
plot_pres(dat, "level_m", "waterdepth_m")

dat <- dat %>%
  select(-level_m1) %>%
  mutate(level_m = na.approx(level_m, maxgap = 3*96, na.rm = F))

plot(dat$level_m, dat$discharge, log = "y")

write_csv(dat, "metabolism/processed/level_v2/PWC_lvl.csv")


# CBP2  ####
dat <- read_csv(paste0("metabolism/processed/CBP.csv")) %>%
  select(DateTime_UTC, temp = temp.water, level_m, discharge) %>%
  left_join(ysi[ysi$site == "CBP",]) 

plot_pres(dat, "level_m", "waterdepth_m")
dat$level_m[which(dat$level_m<.26)] <- NA
# Fill in missing dates
dat <- data.frame(DateTime_UTC = seq(min(dat$DateTime_UTC), 
                                     max(dat$DateTime_UTC),
                                     by = '15 min')) %>%
  left_join(dat)
dat$level_m[32934:32936] <- NA
# first, piece together places where chunks are clearly shifted
# these jumps occur at: (dttm is last pt of data before jump)
gaps <- rle_custom(is.na(dat$level_m)) %>%
  filter(values == 1) %>%
  mutate(starts = starts - 1,
         stops = stops + 1)
dat$waterdepth_m[100] <- dat$level_m[100] + .55
dat$waterdepth_m[32932] <- 1.24
dat$waterdepth_m[32938] <- 1.24

n <- nrow(dat)
g <- nrow(gaps)
dat$level_m1 <- dat$level_m
#dat$level_m <- na.approx(dat$level_m, maxgap = 50, na.rm = F)
for(i in 1:(g+1)){
  if(i == g+1){ end = n } else {end = gaps$starts[i] }
  if(i == 1){ start = 1 } else {start = gaps$stops[i-1] }
  delta <- dat$waterdepth_m[start:end] - dat$level_m[start:end]
  if(sum(!is.na(delta)) > 0){
    dat$level_m[start:n] <- dat$level_m[start:n] + mean(delta, na.rm = T)
  }
}  
plot_pres(dat, "level_m", "waterdepth_m")
# dat$level_m[(w+1):(w+3)] <- NA
# dd <- ymd_hms(c("2019-10-29 14:00:00", "2019-08-21 19:45:00"))
# w <- which(dat$DateTime_UTC %in% dd)
# dat$level_m[w] <- NA
# nhc$level_m <- na.approx(nhc$level_m, maxgap = 50, na.rm = F)
# nhc$level_d <- drift_correct(nhc, "level_m", "waterdepth_m")

dat <- dat %>%
  select(-level_m1) %>%
  mutate(level_m = na.approx(level_m, maxgap = 3*96, na.rm = F))

plot(dat$level_m, dat$discharge, log = "y")

write_csv(dat, "metabolism/processed/level_v2/CBP_lvl.csv")


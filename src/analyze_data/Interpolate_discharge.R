#
# Interpolate discharge to NHC sites without rating curves

# AMCarter
# 2020.03.31

# library(devtools)
# install_github('streampulse/StreamPULSE', dependencies=TRUE)
library(StreamPULSE)
library(LakeMetabolizer)
library(tidyverse)
library(lubridate)
library(zoo)
library(xts)
library(dygraphs)

setwd("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism/data")

# load metadata ####
wsAreas <- read_csv(file="siteData/NHCsite_watersheds.csv")
source("../src/helpers.R")

#Add WS area to sites list
# sites <- left_join(sites, wsAreas[,c(2,7)], by = "sitecode")


# find date range for which we need discharge for NHC sites

sites <- read_csv(file="siteData/NHCsite_metadata.csv")
# sites <- sites[1:7,] # get rid of MC751 and mud, they are not along same continuum
#dateRange <- c(min(sites$startdate.UTC), max(sites$enddate.UTC))
#DateTime_UTC <- seq(dateRange[1], dateRange[2], by = 15*60)

# load data from sp portal ####
# read in data from SP portal for NHC and UNHC to get discharge
NHC_dat <- request_data("NC_NHC",  
                        variables=c("AirPres_kPa", "AirTemp_C",
                                    "WaterPres_kPa", "WaterTemp_C"))

NHC <- NHC_dat$data %>%
  mutate(value = ifelse(flagtype %in% c("Bad Data", "Questionable"), 
                        NA, value)) %>% 
           select(DateTime_UTC, value, variable) %>%
  pivot_wider(names_from=variable, values_from=value) %>%
  select(DateTime_UTC, pressure_kPa=WaterPres_kPa, temp = WaterTemp_C, 
         AirTemp_C, AirPres_kPa)

plot_pres <- function(NHC, waterpres = "level_m", airpres = "waterdepth_m"){
  NHC %>% select(all_of(waterpres), all_of(airpres)) %>%
    xts(order.by = NHC$DateTime_UTC) %>%
    dygraph() %>%
    dyRangeSelector()
}

plot_pres(NHC, "pressure_kPa", "AirPres_kPa")
# remove NHC out of water pressure values (all below 103 from looking at graph)
NHC$pressure_kPa[which(NHC$pressure_kPa<102.5)]<-NA

# fill in missing airpressure data ####

NOAA_airpres <- StreamPULSE:::FindandCollect_airpres(sites$latitude[1], 
                                                     sites$longitude[1],
                                     min(NHC$DateTime_UTC), 
                                     max(NHC$DateTime_UTC))

NHC <- full_join(NHC, NOAA_airpres, by="DateTime_UTC")

plot_pres(NHC, "air_temp", "AirTemp_C")

# imported air pressure is systematically higher than SP air pressure
diff <- mean(NHC$air_kPa - NHC$AirPres_kPa, na.rm = T)
NHC <- NHC %>% 
  mutate(air_kPa = air_kPa - diff,
         AirPres_kPa = ifelse(is.na(AirPres_kPa), air_kPa, AirPres_kPa),
         AirTemp_C = air_temp) %>%
  select(-air_kPa, -air_temp)
# NHC <- filter(NHC, DateTime_UTC > ymd_hms("2018-01-01 00:00:00"))
# remove leading and ending NA's
w <- range(which(!is.na(NHC$pressure_kPa)), which(!is.na(NHC$temp)))
NHC <- NHC[(w[1]:w[2]),]

UNHC_dat <- request_data("NC_UNHC", variables=c("WaterPres_kPa", "WaterTemp_C"))
UNHC <- UNHC_dat$data %>%
  mutate(value = ifelse(flagtype %in% c("Bad Data", "Questionable"), 
                        NA, value)) %>%
  select(DateTime_UTC, value, variable) %>%
  pivot_wider(names_from=variable, values_from=value) %>%
  select(DateTime_UTC, UNHC.pressure_kPa=WaterPres_kPa, UNHC.temp = WaterTemp_C)
# remove leading and ending NA's
w <- range(which(!is.na(UNHC$UNHC.pressure_kPa)), which(!is.na(UNHC$UNHC.temp)))
UNHC <- UNHC[(w[1]:w[2]),]
# UNHC <- filter(UNHC, DateTime_UTC > ymd_hms("2018-01-01 00:00:00"))

# Calculate depth from water pressure and add sensor offset ####
# Depth = pressure_Pa = kg/ms2/(density_kg/m3*gravity_m/s2)
# density is temperature dependent, for now I am assuming it's just 998 kg/m3
sensor_offsets <- read_csv("siteData/sensor_offsets.csv")
qq <- NHC %>%
  mutate(nhc.pressure_Pa = (pressure_kPa - AirPres_kPa) * 1000,
         nhc.level_m = nhc.pressure_Pa / (water.density(temp) * 9.8) + 
           sensor_offsets[sensor_offsets$site == "NHC", ]$offset_cm / 100) %>%
  select(DateTime_UTC, nhc.temp = temp, nhc.level_m, 
         AirPres_kPa, AirTemp_C) 
qq <- UNHC %>%
  full_join(qq) %>%
  mutate(unhc.pressure_Pa = (UNHC.pressure_kPa - AirPres_kPa) * 1000,
         unhc.level_m = unhc.pressure_Pa / (water.density(UNHC.temp) * 9.8) +
           sensor_offsets[sensor_offsets$site == "UNHC", ]$offset_cm / 100) %>%
  select(DateTime_UTC, nhc.temp, nhc.level_m, 
         unhc.temp = UNHC.temp, unhc.level_m,
         AirPres_kPa, AirTemp_C) %>%
  arrange(DateTime_UTC)
# write_csv(qq, "rating_curves/NHC_UNHC_raw_level.csv")
# join with field depth measurements to correct level ####
qq <- read_csv("rating_curves/NHC_UNHC_raw_level.csv", guess_max = 10000)
  # nhc ####
ysi <- read_csv("water_chemistry/all_nhc_ysi_data.csv") %>%
  select(site, DateTime_UTC, watertemp_C, waterdepth_cm) 

nhc <- qq %>%
  select(DateTime_UTC, temp = nhc.temp, level_m = nhc.level_m) %>%
  left_join(ysi[ysi$site == "NHC",]) %>%
  mutate(waterdepth_m = waterdepth_cm / 100) %>%
  select(-site, -waterdepth_cm)

plot_pres(nhc, "level_m", "waterdepth_m")
# Fill in missing dates
nhc <- data.frame(DateTime_UTC = seq(min(nhc$DateTime_UTC), 
                                     max(nhc$DateTime_UTC),
                                     by = '15 min')) %>%
  left_join(nhc)

# first, piece together places where chunks are clearly shifted
# these jumps occur at: (dttm is last pt of data before jump)
gaps <- rle_custom(is.na(nhc$level_m)) %>%
  filter(values == 1, lengths < 50) %>%
  mutate(datetime = nhc$DateTime_EST[starts],
         starts = starts - 1,
         stops = stops + 1)
big_gaps <- rle_custom(is.na(nhc$level_m)) %>%
  filter(values == 1, lengths > 96*3) %>%
  mutate(starts = starts - 1)

# these gaps look like snapping them together would be incorrect
gaps <- gaps %>% 
  slice(-c(11, 16, 21, 66, 82, 84, 94, 99:nrow(gaps))) %>%
  bind_rows(big_gaps[-1,]) %>%
  mutate(lvl_start = nhc$level_m[starts],
         lvl_stop = nhc$level_m[stops],
         jump = lvl_stop - lvl_start) %>%
  filter(jump >= 0.01 | jump <= -.01) # only adjust for gaps larger than 1 cm

n <- nrow(nhc)
nhc$level_m1 <- nhc$level_m
nhc$level_m <- na.approx(nhc$level_m, maxgap = 50, na.rm = F)
for(i in 1:nrow(gaps)){
  end = gaps$starts[i]
  if(i == 1){ start = 1 } else {start = gaps$stops[i-1] }
  delta <- nhc$waterdepth_m[start:end] - nhc$level_m[start:end]
  if(sum(!is.na(delta)) > 0){
    nhc$level_m[start:n] <- nhc$level_m[start:n] + mean(delta, na.rm = T)
  }
}  
# for(i in 2:nrow(gaps)){
#   end_chunk = min(big_gaps$starts[which(big_gaps$starts >= gaps$stops[i])])
#   nhc$level_m[gaps$stops[i]:end_chunk] <- 
#     nhc$level_m[gaps$stops[i]:end_chunk] - gaps$jump[i]
# }
# nhc$level_m <- na.approx(nhc$level_m, maxgap = 50, na.rm = F)
#nhc$level_d <- drift_correct(nhc, "level_m", "waterdepth_m")
dd <- ymd_hms(c("2018-04-05 17:00:00", "2018-04-06 18:00:00"))
w <- which(nhc$DateTime_UTC %in% dd)
nhc$level_m[c((w[1]+1), (w[1]+2), (w[2]-1), (w[2]-2))] <- NA
nhc$level_m[(w[1]+1):(w[2]-1)] <- nhc$level_m[(w[1]+1):(w[2]-1)] + 
  mean(nhc$level_m[w]) - mean(nhc$level_m[c(w[1]+3, w[2]-3)])
# nhc$level_d[nhc$DateTime_UTC < dd] <- nhc$level_m[nhc$DateTime_UTC < dd] + 
#   mean(nhc$waterdepth_m[nhc$DateTime_UTC < dd] - 
#          nhc$level_m[nhc$DateTime_UTC < dd], na.rm = T)
# nhc$level_d[nhc$DateTime_UTC > dd] <- nhc$level_m[nhc$DateTime_UTC > dd] + 
#   mean(nhc$waterdepth_m[nhc$DateTime_UTC > dd] - 
#          nhc$level_m[nhc$DateTime_UTC > dd], na.rm = T)
plot_pres(nhc, "level_m1", "level_m")
plot_pres(nhc, "level_m", "waterdepth_m")

nhc <- nhc %>%
  select(-level_m1) %>%
  mutate(level_m = na.approx(level_m, maxgap = 96*3, na.rm = F)) %>%
  filter(DateTime_UTC < ymd_hms("2020-04-01 00:00:00"))#%>%

w <- range(which(!is.na(nhc$level_m)))
nhc <- nhc[w[1]:w[2],]
plot_pres(nhc)

  # unhc ####
# ysi <- read_csv("water_chemistry/all_nhc_ysi_data.csv") %>%
#   select(site, DateTime_UTC, watertemp_C, waterdepth_cm) 

unhc <- qq %>%
  select(DateTime_UTC, temp = unhc.temp, level_m = unhc.level_m) %>%
  left_join(ysi[ysi$site == "UNHC",]) %>%
  mutate(waterdepth_m = waterdepth_cm / 100) %>%
  select(-site, -waterdepth_cm)

#plotdd(unhc)

# Fill in missing dates
unhc <- data.frame(DateTime_UTC = seq(min(unhc$DateTime_UTC), 
                                     max(unhc$DateTime_UTC),
                                     by = '15 min')) %>%
  left_join(unhc)

# first, piece together places where chunks are clearly shifted
# these jumps occur at: (dttm is last pt of data before jump)
gaps <- rle_custom(is.na(unhc$level_m)) %>%
  filter(values == 1, lengths < 50) %>%
  mutate(datetime = unhc$DateTime_UTC[starts],
         starts = starts - 1,
         stops = stops + 1)
tmp <- data.frame(starts = c(which(unhc$DateTime_UTC ==
                                     ymd_hms("2017-05-01 17:30:00")),
                             which(unhc$DateTime_UTC == 
                                     ymd_hms("2017-04-11 14:30:00")),
                             which(unhc$DateTime_UTC == 
                                     ymd_hms("2017-05-16 13:45:00")),
                             which(unhc$DateTime_UTC == 
                                     ymd_hms("2017-07-12 15:15:00")),
                             which(unhc$DateTime_UTC == 
                                     ymd_hms("2017-08-09 16:15:00")),
                             which(unhc$DateTime_UTC == 
                                     ymd_hms("2017-08-30 15:15:00")),
                             which(unhc$DateTime_UTC == 
                                     ymd_hms("2018-06-19 15:15:00"))
                             ))
tmp$stops <- tmp$starts + 1
big_gaps <- rle_custom(is.na(unhc$level_m)) %>%
  filter(values == 1, lengths > 96*3) %>%
  mutate(starts = starts - 1)

# these gaps look like snapping them together would be incorrect
gaps <- gaps %>% 
  slice(-c(34,52)) %>%
  bind_rows(tmp) %>%
  bind_rows(big_gaps[-1,]) %>%
  mutate(lvl_start = unhc$level_m[starts],
         lvl_stop = unhc$level_m[stops],
         jump = lvl_stop - lvl_start) %>%
  filter(jump >= 0.01 | jump <= -.01) %>% # only adjust for gaps larger than 1 cm
  arrange(starts)

n <- nrow(unhc)
unhc$level_m1 <- unhc$level_m
unhc$level_m <- na.approx(unhc$level_m, maxgap = 50, na.rm = F)
for(i in 1:nrow(gaps)){
  end = gaps$starts[i]
  if(i == 1){ start = 1 } else {start = gaps$stops[i-1] }
  delta <- unhc$waterdepth_m[start:end] - unhc$level_m[start:end]
  if(sum(!is.na(delta)) > 0){
    unhc$level_m[start:n] <- unhc$level_m[start:n] + mean(delta, na.rm = T)
  }
}  
dd <- ymd_hms("2019-12-04 18:15:00")
w <- which(unhc$DateTime_UTC == dd)
unhc$level_m[(w+1):n] <- unhc$level_m[(w+1):n] + unhc$level_m[w]- unhc$level_m[w+3] 
dd <- ymd_hms("2019-10-02 16:15:00")
w <- which(unhc$DateTime_UTC == dd)
unhc$level_m[(w+1):n] <- unhc$level_m[(w+1):n] + unhc$level_m[w]- unhc$level_m[w+1] 
# for(i in 2:nrow(gaps)){
#   end_chunk = min(big_gaps$starts[which(big_gaps$starts >= gaps$stops[i])])
#   unhc$level_m[gaps$stops[i]:end_chunk] <- 
#     unhc$level_m[gaps$stops[i]:end_chunk] - gaps$jump[i]
# }
# unhc$level_m <- na.approx(unhc$level_m, maxgap = 50, na.rm = F)
# unhc$level_d <- drift_correct(unhc, "level_m", "waterdepth_m")

plot_pres(unhc, "level_m", "waterdepth_m")

unhc <- unhc %>%
  select(-level_m1) %>%
  mutate(level_m = na.approx(level_m, maxgap = 96*3, na.rm = F))
#         rename(level_m = level_d) #%>%
  plot_pres(unhc)

w <- range(which(!is.na(unhc$level_m)))
unhc <- unhc[w[1]:w[2],]

qq <- unhc %>%
  rename(temp_unhc = temp, level_unhc = level_m) %>%
  select(DateTime_UTC, temp_unhc, level_unhc) %>%
  full_join(nhc, by = "DateTime_UTC") %>%
  arrange(DateTime_UTC) %>%
  select(DateTime_UTC, temp_unhc, temp_nhc = temp,
         level_unhc, level_nhc = level_m) 

plot_pres(qq, "level_nhc", "level_unhc")
plot(qq$level_unhc, qq$level_nhc)

 # write_csv(qq, "rating_curves/NHC_UNHC_corrected_level.csv")
# Calculate discharge from rating curves ####
# Q = a * level ^ b
# ZQdat_sp <- read_csv(file="siteData/NC_streampulseZQ_data.csv")
ZQdat <- read_csv(file="rating_curves/modified_ZQ_curves.csv")
qq <- read_csv("rating_curves/NHC_UNHC_corrected_level.csv", guess_max = 10000)
# qq <- read_csv("rating_curves/NHC_UNHC_raw_level.csv", guess_max = 10000) %>%
#   rename(level_nhc = nhc.level_m, level_unhc = unhc.level_m)

qq <- qq %>%
  mutate(discharge_nhc = exp(ZQdat$a[1] + ZQdat$b[1] * log(level_nhc)),
         discharge_unhc = exp(ZQdat$a[2] + ZQdat$b[2] * log(level_unhc)))
         # discharge_unhc3 = exp(ZQdat$a[4] + ZQdat$b[4] * log(level_unhc)),
         # discharge_nhc2 = exp(ZQdat$a[3] + ZQdat$b[3] * log(level_nhc)))
par(mfrow = c(2,1), mar = c(2,5,1,2))
plot(qq$DateTime_UTC, qq$discharge_nhc,
     type = 'l', log = "y", ylab = "nhc Q")
# plot(qq$DateTime_UTC, qq$discharge_nhc, 
#      type = 'l', ylim = c(0,2), ylab = "nhc Q")
abline(h = c(ZQdat$max_Q[1], ZQdat$min_Q[1]), col = "brown3")#, lty = 2)
# lines(qq$DateTime_UTC, qq$discharge_unhc, col = 3)
plot(qq$DateTime_UTC, qq$discharge_unhc,
     type = 'l', log = "y", ylab = "unhc Q")
abline(h = c(ZQdat$max_Q[2], ZQdat$min_Q[2]), col = "brown3", lty = 2)
#lines(qq$DateTime_UTC, qq$discharge_unhc2, col = 5)
par(mfrow = c(1,1))
plot(qq$DateTime_UTC, qq$discharge_nhc, log = "y", type = "l")
lines(qq$DateTime_UTC, qq$discharge_unhc, col = 5)
# lines(qq$DateTime_UTC, qq$discharge_nhc2, col = 2)
# lines(qq$DateTime_UTC, qq$discharge_unhc3, col = 3)

qq %>%
  select(discharge_nhc,  discharge_unhc) %>%
  xts(order.by = qq$DateTime_UTC) %>%
  dygraph() %>% dyRangeSelector()
nnhc <- sum(!is.na(qq$discharge_nhc))
nunhc <- sum(!is.na(qq$discharge_unhc))

ZQdat <- data.frame(site = c("NHC","UNHC"),
           above_rc = 
             c(sum(qq$discharge_nhc > ZQdat$max_Q[1], na.rm = T)/nnhc,
               sum(qq$discharge_unhc > ZQdat$max_Q[2], na.rm = T)/nunhc),
           below_rc = 
             c(sum(qq$discharge_nhc < ZQdat$min_Q[1], na.rm = T)/nnhc,
               sum(qq$discharge_unhc < ZQdat$min_Q[2], na.rm = T)/nunhc)) %>% 
  left_join(ZQdat)

write_csv(ZQdat, "rating_curves/modified_ZQ_curves.csv")

par(mfrow = c(1,1))
plot(qq$discharge_unhc, qq$discharge_nhc,log = "xy",
     xlab = "unhc Q", ylab = "nhc Q", pch = 20)
# ggplot(qq, aes(log(level_unhc), log(level_nhc), col = DateTime_UTC)) +
#   geom_point() 
# the relationship looks like it changes from year to year, so I will use the 
# fit from one year at a time.
# fill missing NHC level based on UNHC #### 
qm <- qq %>%
  mutate(level_nhc_mod = NA, 
         level_unhc_mod = NA, 
         year = year(DateTime_UTC))
cordat <- data.frame()
for(i in 2017:2019){
  if(i == 2017) {
    years = c(2016, 2017)}
  if(i == 2018){
    years = 2018 }
  if(i == 2019){
    years = c(2019, 2020) }
  
  dat <- qm %>%
    filter(year %in% years)
  
  m = lm(log(level_nhc) ~ log(level_unhc), dat)$coefficients
  m2 = lm(log(level_unhc) ~ log(level_nhc), dat)$coefficients
  qm <- qm %>%
    mutate(level_nhc_mod = ifelse(year %in% years, 
                                  exp(m[1] + m[2] * log(level_unhc)),
                                  level_nhc_mod),
           level_unhc_mod = ifelse(year %in% years, 
                                   exp(m2[1] + m2[2] * log(level_nhc)),
                                   level_unhc_mod))
}
  
par(mfrow = c(1,2))
plot(qm$level_nhc, qm$level_nhc_mod, pch = 20, 
     xlab = "measured", ylab = "modeled", main = "NHC")
abline(0, 1, col = 2)
plot(qm$level_unhc, qm$level_unhc_mod, pch = 20, 
     xlab = "measured", ylab = "modeled", main = "UNHC")
abline(0, 1, col = 2)

qm <- qm %>% 
  mutate(level_nhc_mod = ifelse(is.na(level_nhc), level_nhc_mod, level_nhc),
         level_unhc_mod = ifelse(is.na(level_unhc), level_unhc_mod, level_unhc),
         notes = case_when(is.na(level_nhc) ~ "nhc modeled",
                           is.na(level_unhc) ~ "unhc modeled"))

# snap interpolated levels  to their neighbors ####

nhc_gaps <- rle_custom(is.na(qm$level_nhc))
unhc_gaps <- rle_custom(is.na(qm$level_unhc))


# fill in and plot modeled data
par(mfrow = c(1,1))
plot(qm$DateTime_UTC, qm$level_nhc_mod, log="y", main = "NHC", pch = 20)
points(qm$DateTime_UTC[qm$notes == "nhc modeled"], 
       qm$level_nhc_mod[qm$notes=="nhc modeled"], pch = 20, col = "red")
plot(qm$DateTime_UTC, qm$level_unhc_mod, log="y", main = "UNHC", pch = 20)
points(qm$DateTime_UTC[qm$notes == "unhc modeled"], 
       qm$level_unhc_mod[qm$notes=="unhc modeled"], pch = 20, col = "red")

# find endpoints of measured and modeled data in the gaps
nhc_gaps <- nhc_gaps[nhc_gaps$values==1,]
unhc_gaps <- unhc_gaps[unhc_gaps$values==1,]

# Don't allow interpolated discharge to be lower than NHC min flow
NHCmin <- min(qm$level_nhc, na.rm=T)
m <- min(qm$level_nhc_mod, na.rm=T)
 # it isnt.

for(i in 1:nrow(nhc_gaps)){
  a <- nhc_gaps[i,]$starts
  b <- nhc_gaps[i,]$stops
  
  if(a==1) next
  startdiff <- qm$level_nhc[a-1] - qm$level_nhc_mod[a]
  if(b==nrow(qm)){
    enddiff <- startdiff
  } else{
      enddiff <- qm$level_nhc[b+1] - qm$level_nhc_mod[b]
  }
  if(is.na(startdiff)||is.na(enddiff)) next
  
  diffQ <- seq(startdiff, enddiff, length.out=nhc_gaps[i,]$lengths)
  
  qm$level_nhc_mod[a:b]<- qm$level_nhc_mod[a:b] + diffQ
}

# snap UNHC gaps
# Don't allow interpolated discharge to be lower than NHC min flow
UNHCmin <- min(qm$level_unhc, na.rm=T)
m <- min(qm$level_unhc_mod, na.rm=T)
# it isnt

for(i in 1:nrow(unhc_gaps)){
  a <- unhc_gaps[i,]$starts
  b <- unhc_gaps[i,]$stops
  
  if(a==1) next
  startdiff <- qm$level_unhc[a-1] - qm$level_unhc_mod[a]
  if(b==nrow(qm)){
    enddiff <- startdiff
  } else{
      enddiff <- qm$level_unhc[b+1] - qm$level_unhc_mod[b]
  }
  if(is.na(startdiff)||is.na(enddiff)) next
  
  diffQ <- seq(startdiff, enddiff, length.out=unhc_gaps[i,]$lengths)
  
  qm$level_unhc_mod[a:b]<- qm$level_unhc_mod[a:b] + diffQ
}

# double check that everything looks okay
par(mfrow = c(1,1))
plot(qm$DateTime_UTC, qm$level_nhc_mod, log="y", main = "NHC", pch = 20)
points(qm$DateTime_UTC[qm$notes == "nhc modeled"], 
       qm$level_nhc_mod[qm$notes=="nhc modeled"], pch = 20, col = "red")
plot(qm$DateTime_UTC, qm$level_unhc_mod, log="y", main = "UNHC", pch = 20)
points(qm$DateTime_UTC[qm$notes == "unhc modeled"], 
       qm$level_unhc_mod[qm$notes=="unhc modeled"], pch = 20, col = "red")

qmm <- qm %>%
  mutate(discharge_nhc_mod = exp(ZQdat$a[1] + ZQdat$b[1] * log(level_nhc_mod)),
         discharge_unhc_mod = exp(ZQdat$a[2] + ZQdat$b[2] * log(level_unhc_mod)))

tmp <- read_csv("rating_curves/NHC_UNHC_raw_level.csv", guess_max = 10000) %>%
  select(DateTime_UTC, AirPres_kPa)

NHC_UNHC_Q_interp <- select(qmm, DateTime_UTC,
                            NHC_level = level_nhc_mod,
                            NHC_Q = discharge_nhc_mod,
                            UNHC_level = level_unhc_mod,
                            UNHC_Q = discharge_unhc_mod,
                            notes) %>%
  left_join(tmp)

write_csv(NHC_UNHC_Q_interp, "rating_curves/NHC_UNHC_Q.csv")

# Interpolate discharge ####
sites <- read_csv("siteData/NHCsite_metadata.csv") %>%
  slice(1:7)
newQdat <- data.frame(matrix(NA, nrow=nrow(qmm), ncol = (1+nrow(sites))))
colnames(newQdat)<- c("DateTime_UTC",paste(sites$sitecode, "Q", sep="."))

newQdat$DateTime_UTC <- qmm$DateTime_UTC
newQdat$NHC.Q <- qmm$discharge_nhc_mod
newQdat$UNHC.Q <- qmm$discharge_unhc_mod

for(i in which(!is.na(newQdat$NHC.Q))){
  df <- data.frame(Q = c(newQdat$NHC.Q[i], newQdat$UNHC.Q[i]), 
                   area = c(wsAreas$ws_area.km2[c(1,7)]))
  m <- lm(Q ~ area, data=df)
  a <- summary(m)$coefficients[,1]
  if(length(a) != 2){ next }
  Qnew <- a[1] + a[2] * sites$ws_area.km2[2:6]
  newQdat[i,3:7] <- Qnew
  if(i %% 5000 == 0){print(newQdat$DateTime_UTC[i])}
}

newQdat <-  NHC_UNHC_Q_interp %>%
  select(DateTime_UTC, AirPres_kPa, notes) %>%
  full_join(newQdat, by="DateTime_UTC")
write_csv(newQdat, path = "rating_curves/interpolatedQ_allsites_modified.csv")
plot(newQdat$DateTime_UTC, newQdat$NHC.Q, col = "grey80", type = "l", log = "y")
lines(newQdat$DateTime_UTC, newQdat$PM.Q, col = "grey60")
lines(newQdat$DateTime_UTC, newQdat$CBP.Q, col = "grey50")
lines(newQdat$DateTime_UTC, newQdat$WB.Q, col = "grey40")
lines(newQdat$DateTime_UTC, newQdat$WBP.Q, col = "grey35")
lines(newQdat$DateTime_UTC, newQdat$UNHC.Q, col = "grey20")


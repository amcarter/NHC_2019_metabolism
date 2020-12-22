#####################################
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
# ZQdat_sp <- read_csv(file="siteData/NC_streampulseZQ_data.csv")
ZQdat <- read_csv(file="rating_curves/modified_ZQ_curves.csv")
wsAreas <- read_csv(file="siteData/NHCsite_watersheds.csv")
sites <- read_csv(file="siteData/NHCsite_metadata.csv")
source("../src/helpers.R")

#Add WS area to sites list
# sites <- left_join(sites, wsAreas[,c(2,7)], by = "sitecode")


# find date range for which we need discharge for NHC sites

sites <- sites[1:7,] # get rid of MC751 and mud, they are not along same continuum
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

plot_pres <- function(NHC, waterpres = "pressure_kPa", airpres = "AirPres_kPa"){
  NHC %>% select(all_of(waterpres), all_of(airpres)) %>%
    xts(order.by = NHC$DateTime_EST) %>%
    dygraph() %>%
    dyRangeSelector()
}

plot_pres(NHC)
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

# join with field depth measurements to correct level ####
ysi <- read_csv("water_chemistry/all_nhc_ysi_data.csv") %>%
  mutate(DateTime_EST = with_tz(DateTime_EST, tz = "EST")) %>%
  select(site, DateTime_EST, watertemp_C, waterdepth_cm) 

nhc <- qq %>%
  mutate(DateTime_EST = with_tz(DateTime_UTC, tz = "EST")) %>%
  select(DateTime_EST, temp = nhc.temp, level_m = nhc.level_m) %>%
  left_join(ysi[ysi$site == "NHC",]) %>%
  mutate(waterdepth_m = waterdepth_cm / 100) %>%
  select(-site, -waterdepth_cm)

plotdd <- function(nhc){
  nhc %>%
    select(waterdepth_m, level_m) %>%
    xts(order.by = nhc$DateTime_EST) %>%
    dygraph() %>%
    dyRangeSelector()
}

plotdd(nhc)
# Fill in missing dates
nhc <- data.frame(DateTime_EST = seq(min(nhc$DateTime_EST), 
                                     max(nhc$DateTime_EST),
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
  mutate(lvl_start = nhc$level_m[starts],
         lvl_stop = nhc$level_m[stops],
         jump = lvl_stop - lvl_start) %>%
  filter(jump >= 0.01 | jump <= -.01) # only adjust for gaps larger than 1 cm

n <- nrow(nhc)
for(i in 2:nrow(gaps)){
  end_chunk = min(big_gaps$starts[which(big_gaps$starts >= gaps$stops[i])])
  nhc$level_m[gaps$stops[i]:end_chunk] <- 
    nhc$level_m[gaps$stops[i]:end_chunk] - gaps$jump[i]
}
nhc$level_m <- na.approx(nhc$level_m, maxgap = 50, na.rm = F)
nhc$level_d <- drift_correct(nhc, "level_m", "waterdepth_m")
plot_pres(nhc, "level_m", "level_d")

nhc <- nhc %>%
  select(-level_m) %>%
  rename(level_m = level_d) #%>%
  plotdd()

w <- range(which(!is.na(nhc$level_m)))
nhc <- nhc[w[1]:w[2],]
# uunhc####
ysi <- read_csv("water_chemistry/all_nhc_ysi_data.csv") %>%
  mutate(DateTime_EST = with_tz(DateTime_EST, tz = "EST")) %>%
  select(site, DateTime_EST, watertemp_C, waterdepth_cm) 

unhc <- qq %>%
  mutate(DateTime_EST = with_tz(DateTime_UTC, tz = "EST")) %>%
  select(DateTime_EST, temp = unhc.temp, level_m = unhc.level_m) %>%
  left_join(ysi[ysi$site == "UNHC",]) %>%
  mutate(waterdepth_m = waterdepth_cm / 100) %>%
  select(-site, -waterdepth_cm)

#plotdd(unhc)

# Fill in missing dates
unhc <- data.frame(DateTime_EST = seq(min(unhc$DateTime_EST), 
                                     max(unhc$DateTime_EST),
                                     by = '15 min')) %>%
  left_join(unhc)

# first, piece together places where chunks are clearly shifted
# these jumps occur at: (dttm is last pt of data before jump)
gaps <- rle_custom(is.na(unhc$level_m)) %>%
  filter(values == 1, lengths < 50) %>%
  mutate(datetime = unhc$DateTime_EST[starts],
         starts = starts - 1,
         stops = stops + 1)
tmp <- data.frame(starts = c(which(unhc$DateTime_EST ==
                                     ymd_hms("2017-05-12 12:30:00", 
                                             tz = "EST")),
                             which(unhc$DateTime_EST == 
                                     ymd_hms("2017-07-12 10:15:00"))))
tmp$stops <- tmp$starts - 1
big_gaps <- rle_custom(is.na(unhc$level_m)) %>%
  filter(values == 1, lengths > 96*3) %>%
  mutate(starts = starts - 1)

# these gaps look like snapping them together would be incorrect
gaps <- gaps %>% 
  slice(-c(34,52)) %>%
  bind_rows(tmp) %>%
  mutate(lvl_start = unhc$level_m[starts],
         lvl_stop = unhc$level_m[stops],
         jump = lvl_stop - lvl_start) %>%
  filter(jump >= 0.01 | jump <= -.01) # only adjust for gaps larger than 1 cm

n <- nrow(unhc)
for(i in 2:nrow(gaps)){
  end_chunk = min(big_gaps$starts[which(big_gaps$starts >= gaps$stops[i])])
  unhc$level_m[gaps$stops[i]:end_chunk] <- 
    unhc$level_m[gaps$stops[i]:end_chunk] - gaps$jump[i]
}
unhc$level_m <- na.approx(unhc$level_m, maxgap = 50, na.rm = F)
unhc$level_d <- drift_correct(unhc, "level_m", "waterdepth_m")
plot_pres(unhc, "level_m", "level_d")

unhc <- unhc %>%
  select(-level_m) %>%
  rename(level_m = level_d) #%>%
  plotdd()

w <- range(which(!is.na(unhc$level_m)))
unhc <- unhc[w[1]:w[2],]

qq <- unhc %>%
  rename(temp_unhc = temp, level_unhc = level_m) %>%
  select(DateTime_EST, temp_unhc, level_unhc) %>%
  full_join(nhc, by = "DateTime_EST") %>%
  arrange(DateTime_EST) %>%
  select(DateTime_EST, temp_unhc, temp_nhc = temp,
         level_unhc, level_nhc = level_m) 

plot(qq$level_unhc, qq$level_nhc)

# Calculate discharge from rating curves ####
# Q = a * level ^ b
# figure out where this ZQdat_sp is coming from

m <- lm(log(discharge_cms) ~ log(level_m),
      data = ZQdat_sp[ZQdat_sp$site == "NHC",])
ab.nhc_sp <- summary(m)$coefficients[,1]
ab.nhc <- as.data.frame(ZQdat[1,2:3])
plot(ZQdat_sp$level_m[ZQdat_sp$site == "NHC"],
     ZQdat_sp$discharge_cms[ZQdat_sp$site == "NHC"],
     xlim = c(.2,2), ylim = c(0,20), main = "NHC")
lines(seq(.5, 10, by = .01), 
      exp(ab.nhc_sp[1] + log(seq(.5, 10, by = .01)) * ab.nhc_sp[2]),
      lty = 2)
lines(seq(.5, 10, by = .01),
      exp(as.numeric(ab.nhc[1]) + 
            as.numeric(ab.nhc[2]) * log(seq(.5, 10, by = .01))))

m <- lm(log(discharge_cms) ~ log(level_m),
        data = ZQdat_sp[ZQdat_sp$site == "UNHC",])
ab.unhc_sp <- summary(m)$coefficients[,1]
ab.unhc <- as.data.frame(ZQdat[2,2:3])
plot(ZQdat_sp$level_m[ZQdat_sp$site =="UNHC"],
     ZQdat_sp$discharge_cms[ZQdat_sp$site == "UNHC"],
     xlim = c(.2, 2), ylim = c(0,20), main = "UNHC")
lines(seq(.1, 2, by = .01),
      exp(ab.unhc_sp[1] + log(seq(.1, 2, by = .01)) * ab.unhc_sp[2]),
          lty = 2)
lines(seq(.1, 2, by = .01),
      exp(as.numeric(ab.unhc[1]) + 
            as.numeric(ab.unhc[2]) * log(seq(.1, 2, by = .01))))
par(mfrow = c(1,2))
et_cm/100) %>%
         unhc.temp = UNHC.temp, unhc.level_m,
         AirPres_kPa)

sets$site == "NHC",]$offset_cm/100,

qq <- left_join(NHC, UNHC, by="DateTime_UTC")

Qdat <- qq %>%
  mutate(nhc.discharge = exp(as.numeric(ab.nhc[1])
                             + as.numeric(ab.nhc[2]) * log(nhc.level_m)),
         nhc.discharge = ifelse(nhc.discharge > 500 | nhc.discharge < .002, 
                                NA, nhc.discharge),
         unhc.discharge = exp(as.numeric(ab.unhc[1]) + 
                                as.numeric(ab.unhc[2]) * log(unhc.level_m)),
         unhc.discharge = ifelse(unhc.discharge >= 120, NA, unhc.discharge),
         nhc.discharge = na.approx(nhc.discharge, na.rm = F, maxgap = 12),
         unhc.discharge = na.approx(unhc.discharge, na.rm = F, maxgap = 12))

par(mfrow = c(1,1))
plot(Qdat$unhc.discharge, (Qdat$nhc.discharge),log = "xy",
     xlab = "unhc Q", ylab = "nhc Q")
#     col = alpha(1,.01), pch = 20)#, xlim = c(0,5), ylim = c(0,1.2))
m <- lm(unhc.discharge ~ nhc.discharge,
        Qdat)
mm <- summary(m)$coefficients[,1]
abline(mm, col = 3)

m <- lm(log(nhc.discharge) ~ log(unhc.discharge),
        Qdat)
mm <- summary(m)$coefficients[,1]

lines(seq(.01, 50.01, by = .1), 
      exp(mm[1] + mm[2] * log(seq(.01,50.01, by = .1))), col = 2)
abline(mm, col = 2)

write_csv(Qdat, "rating_curves/NHC_UNHC_Q.csv")

######################################################
# fill missing NHC Q based on UNHC. 
# This should be revisited, there is clear historesis in the Q-Q relationship
Qdat <- read_csv("rating_curves/NHC_UNHC_Q.csv")
Qdat <- 
  Qdat %>%
  mutate(predNHC_Q = exp(mm[1] + mm[2] * log(unhc.discharge)),
         predUNHC_Q = exp((log(nhc.discharge) - mm[1])/mm[2]),
         notes = case_when(is.na(nhc.discharge) & !is.na(predNHC_Q) ~ "modeled NHC Q",
                           is.na(unhc.discharge) & !is.na(predUNHC_Q) ~ "modeled UNHC Q"),
         modNHC_Q = ifelse(is.na(nhc.discharge), predNHC_Q, nhc.discharge),
         modUNHC_Q = ifelse(is.na(unhc.discharge), predUNHC_Q, unhc.discharge))
         

# snap NHC interpolated points to their neighbors

NHC_gaps <- rle_custom(is.na(Qdat$nhc.discharge))
UNHC_gaps <- rle_custom(is.na(Qdat$unhc.discharge))


# fill in and plot modeled data

plot(Qdat$DateTime_UTC, Qdat$modNHC_Q, log="y", main = "NHC")
points(Qdat$DateTime_UTC[Qdat$notes=="modeled NHC Q"], 
       Qdat$modNHC_Q[Qdat$notes=="modeled NHC Q"], col = "red")

plot(Qdat$DateTime_UTC, Qdat$modUNHC_Q, log="y", main = "UNHC")
points(Qdat$DateTime_UTC[Qdat$notes=="modeled UNHC Q"], 
       Qdat$modUNHC_Q[Qdat$notes=="modeled UNHC Q"], col = "red")

# find endpoints of measured and modeled data in the gaps
NHC_gaps <- NHC_gaps[NHC_gaps$values==1,]
UNHC_gaps <- UNHC_gaps[UNHC_gaps$values==1,]

# Don't allow interpolated discharge to be lower than NHC min flow
NHCmin <- min(Qdat$nhc.discharge, na.rm=T)
m<- min(Qdat$predNHC_Q, na.rm=T)
t <- which(Qdat$predNHC_Q==m)

for(i in 1:nrow(NHC_gaps)){
  a<- NHC_gaps[i,]$starts
  b <- NHC_gaps[i,]$stops
  
  if(a==1) next
  startdiff <- Qdat$nhc.discharge[a-1] - Qdat$predNHC_Q[a]
  if(b==nrow(Qdat)){
    enddiff <- startdiff
  } else{
      enddiff <- Qdat$nhc.discharge[b+1] - Qdat$predNHC_Q[b]
  }
  if(is.na(startdiff)||is.na(enddiff)) next
   diffQ <- seq(startdiff, enddiff, length.out=NHC_gaps[i,]$lengths)
  if(t %in% seq(a, b)){
   tmp1 <- seq(startdiff,(NHCmin-Qdat$predNHC_Q[t]), length.out=(t-a))
   tmp2 <- seq((NHCmin-Qdat$predNHC_Q[t]), enddiff, length.out=(b-t+1))
   diffQ <- c(tmp1,tmp2)
  }
  Qdat$modNHC_Q[a:b]<- Qdat$predNHC_Q[a:b]+diffQ
  
}

# snap UNHC gaps
# Don't allow interpolated discharge to be lower than NHC min flow
UNHCmin <- min(Qdat$unhc.discharge, na.rm=T)
m<- min(Qdat$predNHC_Q, na.rm=T)
t <- which(Qdat$predNHC_Q==m)

for(i in 1:nrow(UNHC_gaps)){
  a<- UNHC_gaps[i,]$starts
  b <- UNHC_gaps[i,]$stops
  startdiff <- Qdat$unhc.discharge[a-1] - Qdat$predUNHC_Q[a]
  if(b==nrow(Qdat)){
    enddiff <- startdiff
  } else{
    enddiff <- Qdat$unhc.discharge[b+1] - Qdat$predUNHC_Q[b]
  }
  if(is.na(startdiff)||is.na(enddiff)) next
  diffQ <- seq(startdiff, enddiff, length.out=NHC_gaps[i,]$lengths)
  
  Qdat$modUNHC_Q[a:b]<- Qdat$predUNHC_Q[a:b]+diffQ
  
}

# double check that everything looks okay
plot(Qdat$DateTime_UTC, Qdat$modNHC_Q, log="y", main = "NHC")
points(Qdat$DateTime_UTC[Qdat$notes=="modeled NHC Q"], 
       Qdat$modNHC_Q[Qdat$notes=="modeled NHC Q"], col = "red")

Qdat$modNHC_Q[Qdat$modNHC_Q<.01] <- NA

plot(Qdat$DateTime_UTC, Qdat$modUNHC_Q, log="y", main = "UNHC")
points(Qdat$DateTime_UTC[Qdat$notes=="modeled UNHC Q"], 
       Qdat$modUNHC_Q[Qdat$notes=="modeled UNHC Q"], col = "red")
diff <- Qdat$modUNHC_Q[Qdat$DateTime_UTC == ymd_hms("2020-03-18 19:15:00")]-
  Qdat$modUNHC_Q[Qdat$DateTime_UTC == ymd_hms("2020-03-18 17:30:00")]

Qdat$modUNHC_Q[Qdat$DateTime_UTC > ymd_hms("2020-03-18 17:30:00") &
                 Qdat$DateTime_UTC < ymd_hms("2020-07-15 21:00:00")] <-
  Qdat$modUNHC_Q[Qdat$DateTime_UTC > ymd_hms("2020-03-18 17:30:00") &
                   Qdat$DateTime_UTC < ymd_hms("2020-07-15 21:00:00")] - diff/2
Qdat$modUNHC_Q[Qdat$DateTime_UTC > ymd_hms("2020-03-18 17:30:00") &
                 Qdat$DateTime_UTC < ymd_hms("2020-06-15 21:00:00")] <-
  Qdat$modUNHC_Q[Qdat$DateTime_UTC > ymd_hms("2020-03-18 17:30:00") &
                   Qdat$DateTime_UTC < ymd_hms("2020-06-15 21:00:00")] - diff/4


Qdat$modUNHC_Q[Qdat$modUNHC_Q<.02] <- NA

NHC_UNHC_Q_interp <- select(Qdat, DateTime_UTC, AirPres_kPa,
                            NHC_Q = modNHC_Q,
                            UNHC_Q = modUNHC_Q, notes)

write_csv(NHC_UNHC_Q_interp, "rating_curves/NHC_UNHC_Q.csv")

#########################################################
# Add modeled Q for NHC sites to a dataframe with columns for each site to interpolate
newQdat <- data.frame(matrix(NA, nrow=nrow(Qdat), ncol = (1+nrow(sites))))
colnames(newQdat)<- c("DateTime_UTC",paste(sites$sitecode, "Q", sep="."))

newQdat$DateTime_UTC <- Qdat$DateTime_UTC
newQdat$NHC.Q <- Qdat$modNHC_Q
newQdat$UNHC.Q <- Qdat$modUNHC_Q

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

newQdat <- full_join(newQdat, Qdat[,c(1,6,11)], by="DateTime_UTC")
newQdat <- read_csv("siteData/interpolatedQ_allsites.csv")
plot(newQdat$DateTime_UTC, newQdat$NHC.Q, col = "grey80", type = "l", log = "y")
lines(newQdat$DateTime_UTC, newQdat$PM.Q, col = "grey60")
lines(newQdat$DateTime_UTC, newQdat$CBP.Q, col = "grey50")
lines(newQdat$DateTime_UTC, newQdat$WB.Q, col = "grey40")
lines(newQdat$DateTime_UTC, newQdat$WBP.Q, col = "grey35")
lines(newQdat$DateTime_UTC, newQdat$UNHC.Q, col = "grey20")

write_csv(newQdat, path = "rating_curves/interpolatedQ_allsites_modified.csv")

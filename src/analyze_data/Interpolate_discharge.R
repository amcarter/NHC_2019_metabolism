#####################################
# Interpolate discharge to NHC sites without rating curves

# AMCarter
# 2020.03.31

# library(devtools)
# install_github('streampulse/StreamPULSE', dependencies=TRUE)
library(StreamPULSE)
library(tidyverse)
library(lubridate)
library(zoo)

setwd("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism/data")

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


# read in data from SP portal for NHC and UNHC to get discharge
NHC_dat <- request_data("NC_NHC",  
                        variables=c("AirPres_kPa", "AirTemp_C",
                                    "WaterPres_kPa", "WaterTemp_C"))

NHC <- NHC_dat$data
NHC$value[NHC$flagtype=="Bad Data"|NHC$flagtype=="Questionable"]<- NA
NHC <- NHC %>% select(DateTime_UTC, value, variable)%>%
  pivot_wider(names_from=variable, values_from=value) %>%
  select(DateTime_UTC, pressure_kPa=WaterPres_kPa, temp = WaterTemp_C, 
         AirTemp_C, AirPres_kPa)

# remove NHC out of water pressure values (all below 103 from looking at graph)
NHC$pressure_kPa[which(NHC$pressure_kPa<102)]<-NA

# fill in missing airpressure data

NOAA_airpres <- StreamPULSE:::FindandCollect_airpres(sites$latitude[1], 
                                                     sites$longitude[1],
                                     ymd_hms("2016-07-12 02:00:00"), 
                                     ymd_hms("2020-11-01 05:00:00"))

NHC <- full_join(NHC, NOAA_airpres, by="DateTime_UTC")

NHC[which(is.na(NHC$AirPres_kPa)),"AirPres_kPa"] <- NHC[which(is.na(NHC$AirPres_kPa)),"air_kPa"]
NHC[which(is.na(NHC$AirTemp_C)),"AirTemp_C"] <- NHC[which(is.na(NHC$AirTemp_C)),"air_temp"]

NHC <- select(NHC, -air_kPa, -air_temp)
# NHC <- filter(NHC, DateTime_UTC > ymd_hms("2018-01-01 00:00:00"))

UNHC_dat <- request_data("NC_UNHC", variables=c("WaterPres_kPa", "WaterTemp_C"))
UNHC <- UNHC_dat$data
UNHC$value[UNHC$flagtype=="Bad Data"|UNHC$flagtype=="Questionable"]<- NA
UNHC <- UNHC %>% select(DateTime_UTC, value, variable)%>%
  pivot_wider(names_from=variable, values_from=value) %>%
  select(DateTime_UTC, UNHC.pressure_kPa=WaterPres_kPa, UNHC.temp = WaterTemp_C)

w <- which(UNHC$UNHC.pressure_kPa<101)
w <- w[which(w<100000)]
UNHC$UNHC.pressure_kPa[w]<-NA
UNHC <- filter(UNHC, DateTime_UTC > ymd_hms("2018-01-01 00:00:00"))


qq <- left_join(NHC, UNHC, by="DateTime_UTC")
  
# Calculate depth from water pressure and add sensor offset
# Depth = pressure_Pa = kg/ms2/(density_kg/m3*gravity_m/s2)
# density is temperature dependent, for now I am assuming it's just 998 kg/m3
sensor_offsets <- read_csv("siteData/sensor_offsets.csv")
qq <- qq %>%
  mutate(nhc.pressure_Pa = (pressure_kPa - AirPres_kPa)*1000,
         nhc.level_m = nhc.pressure_Pa/(998 * 9.8) +
           sensor_offsets[sensor_offsets$site == "NHC",]$offset_cm/100,
         unhc.pressure_Pa = (UNHC.pressure_kPa - AirPres_kPa) * 1000,
         unhc.level_m = unhc.pressure_Pa/(998 * 9.8) + 
           sensor_offsets[sensor_offsets$site == "UNHC",]$offset_cm/100) %>%
  mutate(unhc.level_m = ifelse(unhc.level_m < 0.24, NA, unhc.level_m)) %>%
  select(DateTime_UTC, nhc.temp = temp, nhc.level_m,
         unhc.temp = UNHC.temp, unhc.level_m,
         AirPres_kPa)

# Calculate discharge from rating curves ####
# Q = a * level ^ b
par(mfrow = c(1,2))
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

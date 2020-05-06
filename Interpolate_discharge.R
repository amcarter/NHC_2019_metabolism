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

setwd("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism")

ZQdat <- read_csv(file="data/siteData/NC_streampulseZQ_data.csv")
wsAreas <- read_csv(file="data/siteData/NHCsite_watersheds.csv")
sites <- read_csv(file="data/siteData/NHCsite_metadata.csv")
source("code/helpers.R")

#Add WS area to sites list
sites <- left_join(sites, wsAreas[,c(2,7)], by = "sitecode")


# find date range for which we need discharge for NHC sites

sites <- sites[2:6,c(2,9:10)]
dateRange <- c(min(sites$startdate.UTC), max(sites$enddate.UTC))
DateTime_UTC <- seq(dateRange[1], dateRange[2], by = 15*60)


# read in data from NHC and UNHC to get discharge

NHC <- read.csv("data/metabolism/raw/NHC_2020-03-20_raw.csv", stringsAsFactors = FALSE, header=TRUE)
NHC$value[NHC$flagtype=="Bad Data"|NHC$flagtype=="Questionable"]<- NA
NHC <- NHC %>% select(DateTime_UTC, value, variable)%>%
  pivot_wider(names_from=variable, values_from=value) %>%
  select(DateTime_UTC, pressure_kPa=WaterPres_kPa, temp = WaterTemp_C, AirPres_kPa)

# remove NHC out of water pressure values (all below 103 from looking at graph)
NHC$pressure_kPa[which(NHC$pressure_kPa<103)]<-NA

UNHC <- read.csv("data/metabolism/raw/UNHC_2020-03-20_raw.csv", stringsAsFactors = FALSE, header=TRUE)
UNHC$value[UNHC$flagtype=="Bad Data"|UNHC$flagtype=="Questionable"]<- NA
UNHC <- UNHC %>% select(DateTime_UTC, value, variable)%>%
  pivot_wider(names_from=variable, values_from=value) %>%
  select(DateTime_UTC, UNHC.pressure_kPa=WaterPres_kPa, UNHC.temp = WaterTemp_C)

UNHC <- left_join(NHC, UNHC, by="DateTime_UTC")

UNHC$DateTime_UTC<- ymd_hms(UNHC$DateTime_UTC)

# Calculate depth from water pressure and add sensor offset
# Depth = pressure_Pa = kg/ms2/(density_kg/m3*gravity_m/s2)
# density is temperature dependent, for now I am assuming it's just 998 kg/m3
sensor_offsets <- read_csv("data/siteData/sensor_offsets.csv")

UNHC$pressure_Pa <- (UNHC$pressure_kPa-UNHC$AirPres_kPa)*1000
UNHC$U.pressure_Pa<- (UNHC$UNHC.pressure_kPa-UNHC$AirPres_kPa)*1000
UNHC$level_m <- sensor_offsets[sensor_offsets$site=="NHC",]$offset_cm/100 +
  UNHC$pressure_Pa/(998*9.8)

UNHC$U.level_m <- sensor_offsets[sensor_offsets$site=="UNHC",]$offset_cm/100 +
  UNHC$U.pressure_Pa/(998*9.8)
UNHC$U.level_m[UNHC$U.level_m<0.24]<-NA


# Calculate discharge from rating curves
# level = a * Q ^ b
m<-nls(level_m ~ a*discharge_cms^b,
       data=ZQdat[ZQdat$site=="NHC",2:3],start=list(a=1,b=1))
a.NHC <- summary(m)$coefficients[1]
b.NHC <- summary(m)$coefficients[2]#Summary of the regression statistics
plot(ZQdat$discharge_cms[ZQdat$site=="NHC"], ZQdat$level_m[ZQdat$site=="NHC"], xlab="Q", ylab="z", main="NHC")
lines(seq(.1,7, by=.1), a.NHC*seq(.1,7, by=.1)^b.NHC)

m<-nls(level_m ~ a*discharge_cms^b,
       data=ZQdat[ZQdat$site=="UNHC",2:3],start=list(a=1,b=1))
a.UNHC <- summary(m)$coefficients[1]
b.UNHC <- summary(m)$coefficients[2]#Summary of the regression statistics
plot(ZQdat$discharge_cms[ZQdat$site=="UNHC"], ZQdat$level_m[ZQdat$site=="UNHC"], xlab="Q", ylab="z", main="UNHC")
lines(seq(.1,3, by=.1), a.UNHC*seq(.1,3, by=.1)^b.UNHC)

UNHC$NHC_Q_cms <- (UNHC$level_m/a.NHC)^(1/b.NHC)
UNHC$UNHC_Q_cms <- (UNHC$U.level_m/a.UNHC)^(1/b.UNHC)


Qdat <- UNHC %>% select(DateTime_UTC, NHC_Q_cms, UNHC_Q_cms)

Qdat$NHC_Q_cms <- na.approx(Qdat$NHC_Q_cms, na.rm=FALSE, maxgap=12)
Qdat$UNHC_Q_cms <- na.approx(Qdat$UNHC_Q_cms, na.rm=FALSE, maxgap=12)

plot(Qdat$UNHC_Q_cms, Qdat$NHC_Q_cms, xlab = "UNHC Q", ylab = "NHC Q", log = "xy")
mQ <- glm(log(NHC_Q_cms)~log(UNHC_Q_cms), data=Qdat)
a=summary(mQ)$coefficients[1,1]
b=summary(mQ)$coefficients[2,1]
abline(a,b , col = "red")


######################################################
# fill missing NHC Q based on UNHC. 
# This should be revisited, there is clear historesis in the Q-Q relationship

Qdat$predNHC_Q <- exp(a+b*log(Qdat$UNHC_Q_cms))
Qdat$predUNHC_Q <- exp((log(Qdat$NHC_Q_cms)-a)/b)

Qdat$notes <- as.character(NA)
Qdat$notes[is.na(Qdat$NHC_Q_cms)&!is.na(Qdat$predNHC_Q)]<- "modeled NHC Q"
Qdat$notes[is.na(Qdat$UNHC_Q_cms)&!is.na(Qdat$predUNHC_Q)]<- "modeled UNHC Q"

Qdat$NHC_Q_cms[is.na(Qdat$NHC_Q_cms)] <- Qdat$predNHC_Q[is.na(Qdat$NHC_Q_cms)]
Qdat$UNHC_Q_cms[is.na(Qdat$UNHC_Q_cms)] <- Qdat$predUNHC_Q[is.na(Qdat$UNHC_Q_cms)]

plot(Qdat$DateTime_UTC, Qdat$NHC_Q_cms, log="y", main = "NHC")
points(Qdat$DateTime_UTC[Qdat$notes=="modeled NHC Q"], Qdat$NHC_Q_cms[Qdat$notes=="modeled NHC Q"], col = "red")

plot(Qdat$DateTime_UTC, Qdat$UNHC_Q_cms, log="y", main = "UNHC")
points(Qdat$DateTime_UTC[Qdat$notes=="modeled UNHC Q"], Qdat$UNHC_Q_cms[Qdat$notes=="modeled UNHC Q"], col = "red")

# snap NHC interpolated points to their neighbors

NHC_gaps <- rle_custom(is.na(Qdat$NHC_Q_cms))
UNHC_gaps <- rle_custom(is.na(Qdat$UNHC_Q_cms))

# Add columns for Q for NHC sites
newQdat <- data.frame(matrix(NA, nrow=nrow(Qdat), ncol = (1+nrow(sites))))
colnames(newQdat)<- c("DateTime_UTC",paste(sites$sitecode, "Q", sep="."))

newQdat$DateTime_UTC <- Qdat$DateTime_UTC

for(i in which(!is.na(Qdat$NHC_Q_cms))){
  df <- data.frame(Q = c(Qdat[i,3], Qdat[i,4]), area = c(wsAreas$ws_area.km2[c(1,7)]))
  m <- glm(Q~area, data=df)
  Qnew <- summary(m)$coefficients[1,1]+summary(m)$coefficients[2,1]*sites$ws_area.km2
  newQdat[i,2:6] <-Qnew
}

Qdat <- full_join(Qdat, newQdat, by="DateTime_UTC")%>% select(-predNHC_Q, -predUNHC_Q)

write_csv(Qdat, path = "../data/siteData/interpolatedQ_allsites.csv")

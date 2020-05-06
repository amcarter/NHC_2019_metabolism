#####################################
# Interpolate discharge to NHC sites without rating curves

# AMCarter
# 2020.03.31

library(tidyverse)
library(StreamPULSE)

ZQdat <- read_csv(file="../data/siteData/NC_streampulseZQ_data.csv")
wsAreas <- read_csv(file="../data/siteData/NHCsite_watersheds.csv")
sites <- read_csv(file="../data/siteData/NHCsite_metadata.csv")

#Add WS area to sites list
sites <- left_join(sites, wsAreas[,c(2,7)], by = "sitecode")


# find date range for which we need discharge for NHC sites

sites <- sites[2:6,c(2,9:10)]
dateRange <- c(min(sites$startdate.UTC), max(sites$enddate.UTC))
DateTime_UTC <- seq(dateRange[1], dateRange[2], by = 15*60)


# read in data from NHC and UNHC from database based on this daterange

vars <- c("AirPres_kPa","WaterPres_kPa", "WaterTemp_C")
NHC <- request_data("NC_NHC", startdate=dateRange[1], enddate=dateRange[2],variables=vars)$data
NHC <- NHC %>% select(DateTime_UTC, value, variable)%>%
  spread(key=variable, value=value) %>%
  rename(pressure_kPa=WaterPres_kPa, temp = WaterTemp_C)

# remove NHC out of water pressure values (all below 103 from looking at graph)
NHC$pressure_kPa[which(NHC$pressure_kPa<103)]<-NA



UNHC <- request_data("NC_UNHC", startdate=dateRange[1], enddate=dateRange[2],variables=vars)$data
UNHC$value[which(UNHC$flagtype=="Bad Data")]<- NA
UNHC <- select(UNHC, DateTime_UTC, variable, value)%>%
  spread(variable, value) %>%
  rename(pressure_kPa=WaterPres_kPa, temp = WaterTemp_C)
UNHC<-left_join(UNHC, NHC[,1:2])

# Calculate depth from water pressure and add sensor offset
# Depth = pressure_Pa = kg/ms2/(density_kg/m3*gravity_m/s2)
# density is temperature dependent, for now I am assuming it's just 998 kg/m3
sensor_offsets <- read_csv("../data/siteData/sensor_offsets.csv")

NHC$pressure_Pa <- (NHC$pressure_kPa-NHC$AirPres_kPa)*1000
NHC$level_m <- sensor_offsets[sensor_offsets$site=="NHC",]$offset_cm/100 +
  NHC$pressure_Pa/(998*9.8)

UNHC$level_m <- sensor_offsets[sensor_offsets$site=="UNHC",]$offset_cm/100 +
  (UNHC$pressure_kPa-UNHC$AirPres_kPa)*1000/(998*9.8)
UNHC$level_m[UNHC$level_m<0.2]<-NA


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

NHC$NHC_Q_cms <- (NHC$level_m/a.NHC)^(1/b.NHC)
UNHC$UNHC_Q_cms <- (UNHC$level_m/a.UNHC)^(1/b.UNHC)


Qdat <- NHC %>% select(DateTime_UTC, AirPres_kPa, NHC_Q_cms)%>% 
  full_join(select(UNHC, DateTime_UTC, UNHC_Q_cms))

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
Qdat$notes[is.na(Qdat$NHC_Q_cms)]<- "modeled NHC Q"
Qdat$notes[is.na(Qdat$UNHC_Q_cms)]<- "modeled UNHC Q"

Qdat$NHC_Q_cms[is.na(Qdat$NHC_Q_cms)] <- Qdat$predNHC_Q[is.na(Qdat$NHC_Q_cms)]
Qdat$UNHC_Q_cms[is.na(Qdat$UNHC_Q_cms)] <- Qdat$predUNHC_Q[is.na(Qdat$UNHC_Q_cms)]

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

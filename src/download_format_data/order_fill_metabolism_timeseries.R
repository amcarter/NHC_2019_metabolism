#############################################################################
# format and fill annual metabolism time series for NHC 2019-2020 sites

library(dplyr)
library(tidyverse)
library(zoo)


setwd("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism")

NHCdat <- readRDS("data/metabolism/condensed/allNHCsites.rds")
sites <- c("UNHC","PWC", "WBP","WB","CBP","PM","NHC")

metab <- NHCdat$metab
metab$year <- as.numeric(format(metab$date, "%Y"))
metab$DOY <- as.numeric(format(metab$date, "%j"))

# Clean up data to remove negative GPP, positive ER and high flow

metab[!is.na(metab$GPP)& metab$GPP<0,c("GPP.upper","GPP.lower","GPP")]<-NA
metab[!is.na(metab$ER)& metab$ER>0,c("ER","ER.upper","ER.lower")]<-NA
Q_threshold=.95

for(site in sites){
  Q <- quantile(metab[metab$site==site,]$discharge.m3s, Q_threshold, na.rm=T)
  w<-which(metab$site==site & metab$discharge.m3s>Q)
  metab[w,2:7]<-NA
}

# create dataframe with filled DOY ordered metabolism estimates
metab_ordered <- data.frame()
metab_filled <- data.frame()

for(site in sites){
  tmp<- metab[metab$site==site,]
  #List of years 
  year_list <- unique(tmp$year)
  #Number of days per year
  num_days <- ifelse(unique(as.numeric(tmp$year))%%4 != 0, 365, 366)  
  
  #Generating a complete series of dates for this timeperiod
  #Making DOY information
  temp = 0
  for(j in 1:length(num_days)){
    temp[j] = list(1:num_days[j])
  }
  
  #Making a complete reference of Year and DOY
  complete_dates <- cbind(rep(year_list,num_days), unlist(temp))
  colnames(complete_dates) <- c("year", "DOY")
  
  #Placing the data in the complete record
  ts_full <- merge(complete_dates, tmp, by = c("year", "DOY"), all.x = TRUE)
  
  #Adding in newjd data 
  sum_days <- cumsum(num_days)
  
  #Cycling through and adding the number of days from the previous year to DOY 
  #to create newjd  
  ts_full$newjd <- ts_full$DOY
  if(length(unique(ts_full$year)) > 1){  
    for(k in 2:length(unique(as.numeric(ts_full$year)))){
      ts_full[ts_full$year == unique(as.numeric(ts_full$year))[k],]$newjd <- as.numeric(ts_full[ts_full$year == unique(as.numeric(ts_full$year))[k],]$newjd) + (sum_days[k-1])  
    }
  } #End if statement
  
  #Placing the data in the correct order
  ts_ordered <- ts_full[order(ts_full$year, ts_full$DOY), ]  %>% select(-msgs.fit, -warnings, -errors)
  ts_ordered$site <- site
  ts_ordered$date <- as.Date(paste(ts_ordered$year, ts_ordered$DOY), format='%Y %j')
  #Calculating average GPP timeseries
  avg_trajectory <- aggregate(ts_ordered, by = list(ts_ordered$DOY), FUN = mean, na.rm = TRUE)
  tmp<- select(avg_trajectory, DOY, GPP, GPP.upper, GPP.lower, ER, ER.upper, ER.lower, K600, K600.upper, K600.lower,discharge.m3s)
  tmp$site <- site
  tmp[,2:11]<- as.data.frame(apply(tmp[,2:11], 2, na.approx, na.rm=FALSE))
  
  metab_filled <- bind_rows(metab_filled,tmp)
  metab_ordered <- bind_rows(metab_ordered,ts_ordered)
}
str(metab_ordered)

write_csv(metab_ordered, "data/metabolism/ordered_NHC_sites_metabolism.csv")
write_csv(metab_filled, "data/metabolism/filled_NHC_sites_metabolism.csv")


#===============================================================================
#Diagnostic plots for streams
#Created 4/7/2017
#Modified 11/13/2017
#===============================================================================
library("ks") #For making kernel density plots
library("tidyr")
library("dplyr")
library("zoo")
library("scales")
library("readr")

setwd("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism")
ylims=c(-15,5)
kplot=6
cumlim = c(-800,400)
GPP.col="#74B16D"
ER.col="#B99D82"
PAR.col = "#FFE083"
hall_met <- read_csv("../Hall_50yearslater_NHC/hall_table_15.csv")
NHCdat <- readRDS("data/metabolism/condensed/allNHCsites.rds")
metab <- NHCdat$metab
metab[!is.na(metab$GPP)&metab$GPP<0,c("GPP","GPP.upper","GPP.lower")]<-NA
metab[!is.na(metab$ER)&metab$ER>0,c("ER","ER.upper","ER.lower")]<-NA

data <- NHCdat$data
sites <- c("UNHC","WBP","WB","CBP","PM","NHC")
sitematch <- data.frame(hall = c("Concrete", "Blackwood", "Wood Bridge"), 
                        sitecode= c("CBP", "UNHC", "WB"))

alldat <- ts_ordered[1,]

pdf("figures/NHC2019diagnostics3.pdf",width = 7, height = 5.83,onefile = TRUE)

  for(site in sites){
  
  #Setting up graphical parameters for a multipanel graph
  layout(matrix(c(1, 1,2, 3), 2, 2, byrow = TRUE))
    par(oma=c(0,0,2,0))
    #-------------------------------------------------
#Reading in and processing the data
#-------------------------------------------------
  
    #Reading in the merged timeseries data. Data should be formatted with 3 columns: date, GPP, ER
    ts <- metab[metab$site==site,] 
  # remove flow above 95%
    ts[!is.na(ts$discharge.m3s) & ts$discharge.m3s>quantile(ts$discharge.m3s, .95, na.rm=T),
       c("GPP", "GPP.upper","GPP.lower","ER","ER.upper","ER.lower")] <-NA
  #Making an uniterrupted timeseries and add in a newjd column
    #Divide date into year and day
      ts$Year <- as.numeric(format(ts$date, "%Y"))
      ts$DOY <- as.numeric(format(ts$date, "%j"))
    #List of years 
      year_list <- unique(ts$Year)
    #Number of days per year
      num_days <- ifelse(unique(as.numeric(ts$Year))%%4 != 0, 365, 366)  
      
    #Generating a complete series of dates for this timeperiod
      #Making DOY information
        temp = 0
        for(j in 1:length(num_days)){
          temp[j] = list(1:num_days[j])
        }
  
    #Making a complete reference of Year and DOY
      complete_dates <- cbind(rep(year_list,num_days), unlist(temp))
        colnames(complete_dates) <- c("Year", "DOY")
          
    #Placing the data in the complete record
      ts_full <- merge(complete_dates, ts, by = c("Year", "DOY"), all.x = TRUE)
     
    #Adding in newjd data 
      sum_days <- cumsum(num_days)
      
    #Cycling through and adding the number of days from the previous year to DOY 
    #to create newjd  
      ts_full$newjd <- ts_full$DOY
      if(length(unique(ts_full$Year)) > 1){  
        for(k in 2:length(unique(as.numeric(ts_full$Year)))){
          ts_full[ts_full$Year == unique(as.numeric(ts_full$Year))[k],]$newjd <- as.numeric(ts_full[ts_full$Year == unique(as.numeric(ts_full$Year))[k],]$newjd) + (sum_days[k-1])  
        }
      } #End if stantement
        
    #Placing the data in the correct order
      ts_ordered <- ts_full[order(ts_full$Year, ts_full$DOY), ]  
   
      
        #Calculating average GPP timeseries
          avg_trajectory <- aggregate(ts_ordered, by = list(ts_ordered$DOY), FUN = mean, na.rm = TRUE)
          sd_trajectory <- aggregate(ts_ordered, by = list(ts_ordered$DOY), FUN = sd, na.rm = TRUE)

          # #Plotting the average trajectory and GAM
          #   plot(avg_trajectory[, "Group.1"], avg_trajectory[, "GPP"], pch = 20, col = "gray60",xlab="DOY",ylab="GPP")
          #   points(avg_trajectory[, "DOY"], avg_trajectory[, "pred_gam"], pch = 20, col = "black")
     
# Plot annual metabolism graph
    avg_trajectory$NPP <- avg_trajectory$GPP + avg_trajectory$ER
    avg_trajectory$NPP.upper <- avg_trajectory$GPP.lower+avg_trajectory$ER.upper 
    avg_trajectory$NPP.lower <- avg_trajectory$GPP.upper+avg_trajectory$ER.lower 
    par(cex=.9, mar=c(4,4,1,1))
    plot(avg_trajectory$DOY, avg_trajectory$NPP, type = "l", col = "grey25", xlab = "DOY", 
         ylab = "NEP gO2/m2/d", ylim = ylims, lwd = 3)
    polygon(na.approx(c(avg_trajectory$DOY, rev(avg_trajectory$DOY)), na.rm=F), 
            na.approx(c(avg_trajectory$NPP.lower, rev(avg_trajectory$NPP.upper)), na.rm=FALSE),
            col=alpha("grey25", 0.3), border=NA)
    #get hall data
    if(site %in% sitematch$sitecode){
      hall <- sitematch$hall[sitematch$sitecode==site]
      hall_dat <- hall_met[hall_met$site==hall,]
      hall_dat$newdate <- base::as.Date(hall_dat$newdate, format="%m/%d/%Y")
      hall_dat$DOY <-as.numeric(format(hall_dat$newdate, "%j"))
      
      points(hall_dat$DOY, hall_dat$GPP_gO2m2d-hall_dat$ER_gO2m2d, col = "brown3", pch=20, cex=1.7)
      full_year <- data.frame(DOY = seq(1:366))
      hdat <- hall_dat %>% group_by(DOY)%>%
        summarize(ER=mean(ER_gO2m2d, na.rm=T),
                  GPP=mean(GPP_gO2m2d, na.rm=T))
      hdat<-  left_join(full_year, hdat, by="DOY")
      hdata<- hdatb<-hdat
      hdata$DOY <- seq(-365,0)
      hdatb$DOY <- seq(367,(367+365))
      hdat_wrap <- rbind(hdata, hdat, hdatb)
      hdat_wrap$ER <- na.approx(hdat_wrap$ER, na.rm=F)
      hdat_wrap$GPP <- na.approx(hdat_wrap$GPP, na.rm=F)
      hdat <- hdat_wrap[hdat_wrap$DOY %in% seq(1,366),]
    } else {hdat = NA}
    
# ######################################################################################
# # storm responses
#     #read in daily storm data
#     allsites_storms<- read_csv('data/stormmet.csv')
#     # how many days do you want to look at?
#     predays <- 4
#     postdays <- 6
#     
#     storms = which(allsites_storms$storm==1)
#     starts <- storms-predays
#     ends <- storms+postdays
#     allsites_storms$stormday <- NA
#     for(i in 1:length(starts)){
#       n <-starts[i]
#       allsites_storms$stormday[n:(n+predays+postdays)]<- seq(-2,5, by=1)
#     }
#     storms<- allsites_storms[!is.na(allsites_storms$stormday),]
#     
#     #plot storm occurance relative to metab series
#       subset <- storms[storms$site==site,]
#       stormstarts <- which(subset$storm==1)-2
#       plot( seq(-predays, postdays),subset$GPP[stormstarts[1]:(stormstarts[1]+predays+postdays)],ylim=c(-10,5), 
#             ylab="gO2/m2/d", xlab="days since storm", type="n")
#       mtext("10 largest storms", 3, -1, adj=.9, cex=.9)
#       abline(v=0, col="gray", lty=2, lwd=2)
#       abline(h=0)
#       for(j in 1:length(stormstarts)){
#         n <- stormstarts[j]
#         lines(seq(-predays,postdays), subset$GPP[n:(n+predays+postdays)], col = alpha(GPP.col,.7), lwd=2)
#         lines(seq(-predays,postdays), subset$ER[n:(n+predays+postdays)], col = alpha(ER.col,.7), lwd=2)
#       }
    
#-------------------------------------------------
#Cumulative flux plots
#-------------------------------------------------
  #Removing NA values
    datGPP<- select(avg_trajectory, DOY, GPP, GPP.lower, GPP.upper)
      na_rm <- na.omit(datGPP)

  #Performing cumulative sum by year
    na_rm$csum_gpp <- cumsum(na_rm[,"GPP"])
    na_rm$csum_gpp.lower <- cumsum(na_rm[,"GPP.lower"])
    na_rm$csum_gpp.upper <- cumsum(na_rm[,"GPP.upper"])
    plot(na_rm[, "DOY"],na_rm[, "csum_gpp"], ylim =cumlim, 
        ylab = "Cumulative flux (g O2)", xlab = "DOY", pch = 20, col = GPP.col)
    polygon(na.approx(c(na_rm$DOY, rev(na_rm$DOY)), na.rm=FALSE),
              na.approx(c(na_rm$csum_gpp.lower, rev(na_rm$csum_gpp.upper)), na.rm=FALSE),
                        col = alpha(GPP.col, 0.3), border=NA)
    
    datER<- select(avg_trajectory, DOY, ER, ER.lower, ER.upper)
    na_rm <- na.omit(datER)
    
    
    na_rm$csum_er <- cumsum(na_rm[,"ER"])
    na_rm$csum_er.lower <- cumsum(na_rm[,"ER.lower"])
    na_rm$csum_er.upper <- cumsum(na_rm[,"ER.upper"])
    
    points(na_rm[, "DOY"], na_rm[, "csum_er"], pch = 20, col = ER.col)
       polygon(na.approx(c(na_rm$DOY, rev(na_rm$DOY)), na.rm=FALSE),
              na.approx(c(na_rm$csum_er.lower, rev(na_rm$csum_er.upper)), na.rm=FALSE),
              col = alpha(ER.col, 0.3), border=NA)      
    abline(h=0)  
    
    if(!is.na(hdat)){
      hdat$csumGPP <- cumsum(hdat$GPP)
      points(hdat$DOY, hdat$csumGPP, pch=20, cex=.5,col = "grey25")
      hdat$csumER <- cumsum(hdat$ER)
      points(hdat$DOY, -hdat$csumER, cex=.5, pch=20, col = "grey25")
    }
      

#-------------------------------------------------
#Kernel density plots
#-------------------------------------------------
  #Kernel density calculation
    alldat <- bind_rows(alldat, ts_ordered)
    kernel <- kde(na.omit(ts_ordered[, c("GPP", "ER")]))

  #Plotting the kernel density plot
    plot(kernel, xlab = "GPP", ylab = "ER", ylim = c(-kplot, 0), xlim = c(0, kplot),
      display = "filled.contour2",cont=c(30,60,90), col=c(NA,"grey75","grey50","grey25"))
    
    #Add 1:1 line
      abline(0, -1)
    #Add halldat
      if(!is.na(hdat)){
        points(hall_dat$GPP_gO2m2d, -hall_dat$ER_gO2m2d, col="brown3", pch=20, cex=1.2)
      }      
mtext=mtext(paste0("Metabolism for ",site),3, 0, outer=TRUE)
  }

par(mfrow = c(1,1),
    mar = c(4,4,4,2), 
    oma = c(0,0,0,0))
kernel <- kde(na.omit(alldat[,c("GPP","ER")]))
plot(kernel, xlab = "GPP", ylab = "ER", ylim = c(-kplot, 0), xlim = c(0, kplot),
     display = "filled.contour2",cont=c(30,60,90), col=c(NA,"grey75","grey50","grey25"),
     main = "All NHC sites metabolism")
abline(0,-1)
points(hall_met$GPP_gO2m2d, -hall_met$ER_gO2m2d, col="brown3", pch=20, cex=1.7)


dev.off()        





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
sites <- c("UNHC","PWC","WBP", "WB", "CBP", "PM", "NHC")
ylims=c(-15,5)
kplot=6
cumlim = c(-800,400)
GPP.col="#74B16D"
ER.col="#B99D82"
PAR.col = "#FFE083"
pdf("figures/NHC2019diagnostics.pdf",width = 7, height = 5.83,onefile = TRUE)
site <- sites[1]

  for(site in sites){
  
  #Setting up graphical parameters for a multipanel graph
  layout(matrix(c(1, 1,1,2, 2,2,  3, 3,4,4,5,5), 2, 6, byrow = TRUE))
    par(oma=c(0,0,2,0))
    #-------------------------------------------------
#Reading in and processing the data
#-------------------------------------------------
  
    #Reading in the merged timeseries data. Data should be formatted with 3 columns: date, GPP, ER
    dat <- readRDS(paste0("data/metabolism/condensed/condensed_",site,".rds"))
    ts <- dat$metab 
  #Replacing all negative GPP values with NA for the moment
    ts[!is.na(ts$GPP) & ts$GPP < 0, "GPP"] <- NA
          
  #Replacing extremely high measurements (need more refined criteria)
    ts[!is.na(ts$GPP) & ts$GPP > 100, "GPP"] <- NA
      
  #Replacing all positive ER values with NA for the moment
    ts[!is.na(ts$ER) & ts$ER > 0, "ER"] <- NA
  
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
    par(cex=.9, mar=c(4,4,1,1))
    plot(avg_trajectory$DOY, avg_trajectory$GPP, type = "l", col = GPP.col, xlab = "DOY", 
         ylab = "gO2/m2/d", ylim = ylims, lwd = 3)
    lines(avg_trajectory$DOY, avg_trajectory$ER, col = ER.col, lwd = 3)
    polygon(na.approx(c(avg_trajectory$DOY, rev(avg_trajectory$DOY)), na.rm=F), 
            na.approx(c(avg_trajectory$GPP.lower, rev(avg_trajectory$GPP.upper)), na.rm=FALSE),
            col=alpha(GPP.col, 0.3), border=NA)
    polygon(na.approx(c(avg_trajectory$DOY, rev(avg_trajectory$DOY)), na.rm=F), 
            na.approx(c(avg_trajectory$ER.lower, rev(avg_trajectory$ER.upper)), na.rm=FALSE),
            col=alpha(ER.col, 0.3), border=NA)
    legend("bottomleft", c("GPP", "ER"), bty = "n", lty = c(1,1,1),lwd = 3, 
           col = c(GPP.col,  ER.col))
    
######################################################################################
# storm responses
    #read in daily storm data
    allsites_storms<- read_csv('data/stormmet.csv')
    # how many days do you want to look at?
    predays <- 4
    postdays <- 6
    
    storms = which(allsites_storms$storm==1)
    starts <- storms-predays
    ends <- storms+postdays
    allsites_storms$stormday <- NA
    for(i in 1:length(starts)){
      n <-starts[i]
      allsites_storms$stormday[n:(n+predays+postdays)]<- seq(-2,5, by=1)
    }
    storms<- allsites_storms[!is.na(allsites_storms$stormday),]
    
    #plot storm occurance relative to metab series
      subset <- storms[storms$site==site,]
      stormstarts <- which(subset$storm==1)-2
      plot( seq(-predays, postdays),subset$GPP[stormstarts[1]:(stormstarts[1]+predays+postdays)],ylim=c(-10,5), 
            ylab="gO2/m2/d", xlab="days since storm", type="n")
      mtext("10 largest storms", 3, -1, adj=.9, cex=.9)
      abline(v=0, col="gray", lty=2, lwd=2)
      abline(h=0)
      for(j in 1:length(stormstarts)){
        n <- stormstarts[j]
        lines(seq(-predays,postdays), subset$GPP[n:(n+predays+postdays)], col = alpha(GPP.col,.7), lwd=2)
        lines(seq(-predays,postdays), subset$ER[n:(n+predays+postdays)], col = alpha(ER.col,.7), lwd=2)
      }
    
    
    
    
    
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
      

#-------------------------------------------------
#Kernel density plots
#-------------------------------------------------
  #Kernel density calculation
    kernel <- kde(na.omit(ts_ordered[, c("GPP", "ER")]))

  #Plotting the kernel density plot
    plot(kernel, xlab = "GPP", ylab = "ER", ylim = c(-kplot, 0), xlim = c(0, kplot),
      display = "filled.contour2",cont=c(30,60,90), col=c(NA,"grey75","grey50","grey25"))
    
    #Add 1:1 line
      abline(0, -1)

      
###############################################
    # GPP vs PAR
      
      NCdat <- readRDS("data/light/NC_UNHC_predicted.rds")
      
      NCdat <- NCdat[NCdat$Year==2019,]
      NCdat<- NCdat[NCdat$DOY>64,]
      NCdat$PAR_stream <- na.approx(NCdat$PAR_stream, na.rm=FALSE)
      daily_light <- group_by(NCdat, DOY)%>%
        summarize(PAR =sum(PAR_stream))
      dat <- select(ts, DOY, GPP)
      dat <- left_join(dat, daily_light, by="DOY")
      plot(dat$PAR, dat$GPP, col = GPP.col, pch=20, ylim=c(0,ylims[2]),
           xlim = c(0,8000),xlab="PAR (daily cummulative)", ylab = "GPP (gO2/m2/d)")
mtext=mtext(paste0("Metabolism for ",site),3, 0, outer=TRUE)
  }
dev.off()        





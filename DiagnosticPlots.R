#===============================================================================
#Diagnostic plots for streams
#Created 4/7/2017
#Modified 11/13/2017
#===============================================================================
library("mgcv") #For GAMM model
library("ks") #For making kernel density plots
library("tidyr")
library("dplyr")

setwd("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism")
sites <- c("UNHC","PWC","WBP", "WB", "CBP", "PM", "NHC")
ylims=c(-15,5)
kplot=6
cumlim = c(-800,400)
pdf("figures/NHC2019diagnostics.pdf",width = 7, height = 5.83,onefile = TRUE)

  for(site in sites){
  
  #Setting up graphical parameters for a multipanel graph
  layout(matrix(c(1, 2, 2, 1, 3, 3), 2, 3, byrow = TRUE), widths = c(2,1))

#-------------------------------------------------
#Reading in and processing the data
#-------------------------------------------------
  
    #Reading in the merged timeseries data. Data should be formatted with 3 columns: date, GPP, ER
    dat <- readRDS(paste0("data/metabolism/condensed/condensed_",site,".rds"))
    ts <- dat$metab %>% select(date, GPP, ER)
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
   
      
#-------------------------------------------------
#Trying to extract seasonality with a generalized additive model
#http://stats.stackexchange.com/questions/9506/stl-trend-of-time-series-using-r
#-------------------------------------------------
  #TEMPORARY WORKAROUND TO ALLOW THE SCRIPT TO COMPLETE
    #if(length(unique(ts_ordered$Year))>1 & nrow(na.omit(ts_ordered))/(length(unique(ts_ordered$Year))*365) >0.5){  
    #Fit a model with trend and seasonal components --- warning this is slow:
      mod <- gamm(GPP ~ s(DOY, bs = "cc") + s(newjd, bs = "cr"),
        data = ts_ordered, method = "REML",
        correlation = corAR1(form = ~ 1 | Year),
        knots = list(DOY = c(0, 366)))
  
    #The fitted model
      summary(mod$gam)

    #Predicting values with the fitted model
      pred_gam <- predict(mod$gam, newdata = ts_ordered)
      ts_ordered$pred_gam <- pred_gam
        
    # #Plot the fitted model along with all of the data
    #   plot(ts_ordered[,"newjd"],ts_ordered[,"GPP"],pch=20,col="grey60",xlab="Time",ylab="GPP",main=paste(i))
    #   lines(pred_gam ~ newjd, data = ts_ordered, col = "red", lwd = 2)
  
      #Comparing to average trajectory
        #Calculating average GPP timeseries
          avg_trajectory <- aggregate(ts_ordered, by = list(ts_ordered$DOY), FUN = mean, na.rm = TRUE)
          sd_trajectory <- aggregate(ts_ordered, by = list(ts_ordered$DOY), FUN = sd, na.rm = TRUE)

          # #Plotting the average trajectory and GAM
          #   plot(avg_trajectory[, "Group.1"], avg_trajectory[, "GPP"], pch = 20, col = "gray60",xlab="DOY",ylab="GPP")
          #   points(avg_trajectory[, "DOY"], avg_trajectory[, "pred_gam"], pch = 20, col = "black")
     
# Plot annual metabolism graph
    avg_trajectory$NPP <- avg_trajectory$GPP + avg_trajectory$ER
    llim <- min(avg_trajectory$ER, na.rm=T)
    ulim <- max(avg_trajectory$GPP, na.rm=T)
    plot(avg_trajectory$DOY, avg_trajectory$GPP, type = "l", col = "forestgreen", xlab = "DOY", 
         ylab = "gO2/m2/d", ylim = ylims, lwd = 2)
    lines(avg_trajectory$DOY, avg_trajectory$ER, col = "brown3", lwd = 2)
    legend("bottomleft", c("GPP", "NEP","ER"), bty = "n", lty = c(1,1,1),lwd = 2, 
           col = c("forestgreen", "grey60", "brown3"))
    
#-------------------------------------------------
#Kernel density plots
#-------------------------------------------------
  #Kernel density calculation
    kernel <- kde(na.omit(ts_ordered[, c("GPP", "ER")]))

  #Plotting the kernel density plot
      #Determine axes limits
      ll <- -min(ts_ordered$ER, na.rm=T)
      
    plot(kernel, xlab = "GPP", ylab = "ER", ylim = c(-kplot, 0), xlim = c(0, kplot),
      display = "filled.contour2")
    
    #Add 1:1 line
      abline(0, -1)

#-------------------------------------------------
#Cumulative flux plots
#-------------------------------------------------
  #Removing NA values
    na_rm <- na.omit(avg_trajectory)

  #Performing cumulative sum by year
    na_rm$csum_gpp <- cumsum(na_rm[,"GPP"])
    na_rm$csum_er <- cumsum(na_rm[,"ER"])
    na_rm$csum_npp <- cumsum(na_rm[,"NPP"])
    
  #Making plots
    #Average seasonal
      #lim <- max(abs(na_rm[, c("csum_gpp", "csum_er")]))
      plot(na_rm[, "DOY"],na_rm[, "csum_gpp"], ylim =cumlim, 
        ylab = "Cumulative flux", xlab = "DOY", pch = 20, col = "forestgreen")
      
      points(na_rm[, "DOY"], na_rm[, "csum_er"], pch = 20, col = "brown3")
      points(na_rm[, "DOY"], na_rm[, "csum_npp"], pch = 20, col = "grey60")
      
      mtext=mtext(paste0("Metabolism for ",site),3, -2, outer=TRUE)
  }
dev.off()            

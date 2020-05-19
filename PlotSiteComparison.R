# Inter site comparison of NHC metabolism
# comparison to whole SP database

library(dplyr)
library(tidyverse)
library("ks")


setwd("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism")

alldat <- readRDS("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/data/streampulse/gapPhilled_data.rds")

alldat[!is.na(alldat$GPP_filled)& alldat$GPP_filled<0,"GPP"]<-NA
alldat[!is.na(alldat$ER_filled)& alldat$ER_filled>0,"ER"]<-NA

smry <- group_by(alldat, DOY) %>%
  summarize(GPP=median(GPP_filled, na.rm=T),
            GPP.upper=quantile(GPP_filled, .75, na.rm = T),
            GPP.lower=quantile(GPP_filled, .25, na.rm=T),
            ER=median(ER, na.rm=T),
            ER.upper=quantile(ER_filled, .75, na.rm = T),
            ER.lower=quantile(ER_filled, .25, na.rm=T),
            K600=median(K600, na.rm=T),
            K600.upper=quantile(K600, .75, na.rm = T),
            K600.lower=quantile(K600, .25, na.rm=T))

NHCdat <- readRDS("data/metabolism/condensed/allNHCsites.rds")
sites <- c("UNHC", "WBP","WB","CBP","PM","NHC")

metab <- NHCdat$metab
metab$DOY <- as.numeric(format(metab$date, "%j"))

# Clean up data to remove negative GPP, positive ER and high flow

metab[!is.na(metab$GPP)& metab$GPP<0,"GPP"]<-NA
metab[!is.na(metab$ER)& metab$ER>0,"ER"]<-NA
Q_threshold=.95

for(site in sites){
  Q <- quantile(metab[metab$site==site,]$discharge.m3s, Q_threshold, na.rm=T)
  w<-which(metab$site==site & metab$discharge.m3s>Q)
  metab[w,2:7]<-NA
}


#############################################################################
# metabolism comparison to all SP sites
ylims=c(-10, 5)

par(mfrow=c(1,1))
plot(smry$DOY, smry$GPP, ylab='', yaxs='i', type='n',
     bty='n', lwd=4, xlab='', ylim=ylims, xaxs='i', xaxt='n', yaxt='n')
polygon(x=c(smry$DOY, rev(smry$DOY)),
        y=c(smry$GPP.lower, rev(smry$GPP.upper)),
        border=NA, col=alpha('forestgreen', alpha=0.4))
polygon(x=c(smry$DOY, rev(smry$DOY)),
        y=c(smry$ER.lower, rev(smry$ER.upper)),
        border=NA, col=alpha('sienna', alpha=0.4))
axis(2, las=2, line=0, xpd=NA, tck=-.02, labels=FALSE,
     at=round(seq(ylims[1],ylims[2], by=2), 1))
axis(2, las=2, line=-0.5, xpd=NA, tcl=0, col='transparent',
     at=round(seq(ylims[1], ylims[2], by=2), 1))
abline(h=0, lty=1, lwd=1)
mtext(expression(paste("g O"[2]~"m"^"-2"~" d"^"-1")), side=2, line=2.5)
mtext('DOY', side=1, line=2, font=2)

cols <- colorRampPalette(c("grey60", "grey25"))(6)
for(i in 1:6){
  tmp <- metab[metab$site==sites[i],]
  tmp<- tmp[order(tmp$DOY),]
  lines(tmp$DOY, tmp$GPP, col = cols[i], lwd=2)
  lines(tmp$DOY, tmp$ER, col = cols[i], lwd=2)
}


#####################################
# kernel density plot
#Kernel density calculation
k_all <- kde(na.omit(alldat[, c("GPP_filled", "ER_filled")]))
k_nhc <- kde(na.omit(metab[,c("GPP","ER")]))
#Plotting the kernel density plot
plot(k_all, xlab = "GPP", ylab = "ER", ylim = c(-15, 0), xlim = c(0, 15),display = "slice")

par(new=T)
plot(k_nhc, xaxt='n',yaxt='n', ylab='',xlab='',ylim = c(-15, 0), xlim = c(0, 15),display = "filled.contour")
#Add 1:1 line
abline(0, -1)

#########################################################
# O2 exceedence plots
  par(mfrow=c(4,3))
  layout.matrix <- matrix(c(1,2,7,8,3,4,9,10,5,6,11,12), nrow = 4, ncol = 3)
  
  layout(mat = layout.matrix,
         heights = c(1,3,1,3), # Heights of the two rows
         widths = c(3,3,3)) # Widths of the two columns

for(site in sites){
  # Summarize by day
  datdaily <- NHCdat$data[NHCdat$data$site==site,] %>% 
    dplyr::group_by(date) %>% 
    dplyr::summarise(DOsat = mean(DO.obs/DO.sat, na.rm=T))
  datdaily$DOsat <- na.approx(datdaily$DOsat, na.rm=F)
  
  # Calculate exceedence curve for a site
  tmp <- datdaily$DOsat[which(!is.na(datdaily$DOsat))]
  tmp<- sort(tmp, decreasing=T)
  index <- seq(1:length(tmp))
  freq <- 100*(index/(1+length(tmp)))
  
  
  par(mar = c(0,4,3,1))
  specs <- NHCdat$specs[NHCdat$specs$site==site,]
  plot(exp(specs$K600_lnQ_nodes_centers), specs$K600, log='x',ylab = "K600",xlab='',
       xaxt='n',type='b',pch=19,ylim = range(NHCdat$specs$K600))
 # axis(3, labels=F)
  mtext("K 600 vs ln(Q)",3,-1.2,cex=.8)
  par(mar=c(4,4,0,1))
  plot(freq,tmp, lwd = 2, type="l", xlab='% exceedence',ylab="DO %sat")#, col = col, lty=lty)
  mtext(site, 3,-2,adj=.8)  

  par(new=T)
  tmp <- metab[metab$site==site,]
  tmp<- tmp[order(tmp$DOY),]
  plot(tmp$DOY, -tmp$ER, type="l",ylim=c(0,30),xlab='',ylab='',xaxt='n',yaxt='n',col = "sienna", lwd=2)
  
}

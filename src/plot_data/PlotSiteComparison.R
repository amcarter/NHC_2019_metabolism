# Inter site comparison of NHC metabolism
# comparison to whole SP database

library(dplyr)
library(tidyverse)
library(zoo)

setwd("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism")

alldat <- readRDS("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/data/streampulse/gapPhilled_data.rds")
# NHC_model<- readRDS("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/data/streampulse/metabolism/NHC_modeled.rds")$metab
# NHC_model$year<-(as.numeric(format(NHC_model$date, "%Y")))
# NHC_model<-NHC_model[-which(NHC_model$year==2016),]

# alldat[!is.na(alldat$GPP_filled)& alldat$GPP_filled<0,"GPP"]<-NA
# alldat[!is.na(alldat$ER_filled)& alldat$ER_filled>0,"ER"]<-NA

smry <- group_by(alldat, DOY) %>%
  summarize(GPP=median(GPP_filled, na.rm=T),
            GPP.upper=quantile(GPP_filled, .75, na.rm = T),
            GPP.lower=quantile(GPP_filled, .25, na.rm=T),
            ER=median(ER_filled, na.rm=T),
            ER.upper=quantile(ER_filled, .75, na.rm = T),
            ER.lower=quantile(ER_filled, .25, na.rm=T),
            K600=median(K600, na.rm=T),
            K600.upper=quantile(K600, .75, na.rm = T),
            K600.lower=quantile(K600, .25, na.rm=T))

# metab_ordered <- read_csv("data/metabolism/ordered_NHC_sites_metabolism.csv")
# metab_filled <- read_csv("data/metabolism/filled_NHC_sites_metabolism.csv")
metab_filled <- readRDS("data/metabolism/compiled/raymond_met.rds")$cumulative %>%
  mutate(DOY = format(date, "%j"))  
sites <- c("UNHC","WBP","WB","CBP","PM","NHC")

# compile NHC and UNHC data
# NHC_model <- select(NHC_model,  -site)
# NHC_model$DOY <- as.numeric(format(NHC_model$date, "%j"))
# NHC_new <- metab_ordered[metab_ordered$site=="NHC"&metab_ordered$year==2020,]
# NHC<- bind_rows(NHC_model, NHC_new)
NHC_smry <- metab_filled %>%
  filter(site =="nhc") %>%
  group_by(DOY) %>%
  summarize(GPP_m=mean(GPP, na.rm=T),
            GPP.75=quantile(GPP, .75, na.rm = T),
            GPP.25=quantile(GPP, .25, na.rm=T),
            ER_m=mean(ER, na.rm=T),
            ER.75=quantile(ER, .75, na.rm = T),
            ER.25=quantile(ER, .25, na.rm=T),
            GPP.upper=max(GPP, na.rm = T),
            GPP.lower=min(GPP, na.rm=T),
            ER.upper=max(ER, na.rm = T),
            ER.lower=min(ER, na.rm=T))
UNHC_smry <- metab_filled %>%
  filter(site =="unhc") %>%
  group_by(DOY) %>%
  summarize(GPP_m=mean(GPP, na.rm=T),
            GPP.75=quantile(GPP, .75, na.rm = T),
            GPP.25=quantile(GPP, .25, na.rm=T),
            ER_m=mean(ER, na.rm=T),
            ER.75=quantile(ER, .75, na.rm = T),
            ER.25=quantile(ER, .25, na.rm=T),
            GPP.upper=max(GPP, na.rm = T),
            GPP.lower=min(GPP, na.rm=T),
            ER.upper=max(ER, na.rm = T),
            ER.lower=min(ER, na.rm=T))
            # K600=median(K600, na.rm=T))

# NHC_smry$GPP.lower[which(is.infinite(NHC_smry$GPP.lower))]<- NA
# NHC_smry$GPP.upper[which(is.infinite(NHC_smry$GPP.upper))]<- NA
# NHC_smry$ER.lower[which(is.infinite(NHC_smry$ER.lower))]<- NA
# NHC_smry$ER.upper[which(is.infinite(NHC_smry$ER.upper))]<- NA

# 
# UNHC_sp <- SP_models[SP_models$site=="UNHC",]%>% 
#   select(-region, -site, -valid_day, -warnings, -errors)
# UNHC_sp$date <- as.Date(UNHC_sp$date)
# UNHC_sp$DOY <- as.numeric(format(UNHC_sp$date, "%j"))
# UNHC_sp<- select(UNHC_sp,-date)
# UNHC_new <- metab_ordered[metab_ordered$site=="UNHC"&metab_ordered$year==2020,]
# UNHC_new <- select(UNHC_new, -date)
# UNHC<- bind_rows(UNHC_sp, UNHC_new)
# UNHC_smry <- group_by(UNHC, DOY) %>%
#   summarize(GPP_m=mean(GPP, na.rm=T),
#           GPP.75=quantile(GPP, .75, na.rm = T),
#           GPP.25=quantile(GPP, .25, na.rm=T),
#           ER_m=mean(ER, na.rm=T),
#           ER.75=quantile(ER, .75, na.rm = T),
#           ER.25=quantile(ER, .25, na.rm=T),
#           GPP.upper=quantile(GPP, 1, na.rm = T),
#           GPP.lower=quantile(GPP, 0, na.rm=T),
#           ER.upper=quantile(ER, 1, na.rm = T),
#           ER.lower=quantile(ER, 0, na.rm=T),
#           K600=median(K600, na.rm=T))
# 
# create aggregated timeseries data from NHC sites
metab_smry <- group_by(metab_filled, DOY) %>%
  summarize(GPP_m=mean(GPP, na.rm=T),
            GPP.75=quantile(GPP, .75, na.rm = T),
            GPP.25=quantile(GPP, .25, na.rm=T),
            ER_m=mean(ER, na.rm=T),
            ER.75=quantile(ER, .75, na.rm = T),
            ER.25=quantile(ER, .25, na.rm=T),
            GPP.upper=max(GPP, na.rm = T),
            GPP.lower=min(GPP, na.rm=T),
            ER.upper=max(ER, na.rm = T),
            ER.lower=min(ER, na.rm=T))
            # K600=median(K600, na.rm=T))

# metab_smry$GPP.lower[which(is.infinite(metab_smry$GPP.lower))]<- NA
# metab_smry$GPP.upper[which(is.infinite(metab_smry$GPP.upper))]<- NA
# metab_smry$ER.lower[which(is.infinite(metab_smry$ER.lower))]<- NA
# metab_smry$ER.upper[which(is.infinite(metab_smry$ER.upper))]<- NA

#############################################################################
# metabolism comparison to all SP sites

png("figures/metabolism_comparison_allSPdata_NHCmodelrun.png", width=550, height=320)

  par(mfrow = c(1,2), mar=c(0,1,0,0), oma=c(5,5,2,2))
  
  ylims=c(-15, 5)
  
  # between site comparison
  plot(smry$DOY, smry$GPP, ylab='', yaxs='i', type='n',
       bty='n', lwd=4, xlab='', ylim=ylims, xaxs='i', xaxt='n', yaxt='n')
  polygon(x=c(smry$DOY, rev(smry$DOY)),
          y=c(smry$GPP.lower, rev(smry$GPP.upper)),
          border=NA, col=alpha('forestgreen', alpha=.5))
  polygon(x=c(smry$DOY, rev(smry$DOY)),
          y=c(smry$ER.lower, rev(smry$ER.upper)),
          border=NA, col=alpha('sienna', alpha=.5))
  axis(2, las=2, line=0, xpd=NA, tck=-.02, labels=FALSE,
       at=round(seq(-10,ylims[2], by=2), 1))
  axis(2, las=2, line=-0.5, xpd=NA, tcl=0, col='transparent',
       at=round(seq(-10, ylims[2], by=2), 1),)
  abline(h=0, lty=1, lwd=1)
  mtext(expression(paste("g O"[2]~"m"^"-2"~" d"^"-1")), side=2, line=2.5)
  mtext('day of year', side=1, line=-3, font=1)
  
  # polygon(c(metab_smry$DOY, rev(metab_smry$DOY)),
  #         na.approx(c(metab_smry$GPP.25, rev(metab_smry$GPP.75)),na.rm=F),
  #         col=alpha("black",.5), border=NA)
  polygon(c(metab_smry$DOY, rev(metab_smry$DOY)),
          y=c(metab_smry$GPP.lower[1], rollmean(na.approx(c(metab_smry$GPP.lower, rev(metab_smry$GPP.upper)),na.rm=F),3), metab_smry$GPP.upper[1]),
          col=alpha("black",.5), border=NA)
  
  # polygon(c(metab_smry$DOY, rev(metab_smry$DOY)),
  #         na.approx(c(metab_smry$ER.25, rev(metab_smry$ER.75)),na.rm=F),
  #         col=alpha("black",.5), border=NA)
  polygon(x=c(metab_smry$DOY, rev(metab_smry$DOY)),
          y=-c(metab_smry$ER.lower[1], 
              rollmean(na.approx(c(metab_smry$ER.lower, 
                                   rev(metab_smry$ER.upper)),na.rm=F),3), 
              metab_smry$ER.upper[1]),
          col=alpha("black",.5), border=NA)
  
  
  mtext("Intersite metabolism", 3,  line=.6, font=4, col="grey30")
  mtext("(6 sites)", 3, cex=.8, font=4, line=-.2, col="grey30")
  
  # 
  # cols <- colorRampPalette(c("grey60", "grey25"))(6)
  # for(i in 1:6){
  #   tmp <- metab[metab$site==sites[i],]
  #   tmp<- tmp[order(tmp$DOY),]
  #   lines(tmp$DOY, tmp$GPP, col = cols[i], lwd=2)
  #   lines(tmp$DOY, tmp$ER, col = cols[i], lwd=2)
  # }
  
  #Interannual variation
  plot(smry$DOY, smry$GPP, ylab='', yaxs='i', type='n',
       bty='n', lwd=4, xlab='', ylim=ylims, xaxs='i', xaxt='n', yaxt='n')
  polygon(x=c(smry$DOY, rev(smry$DOY)),
          y=c(smry$GPP.lower, rev(smry$GPP.upper)),
          border=NA, col=alpha('forestgreen', alpha=.5))
  polygon(x=c(smry$DOY, rev(smry$DOY)),
          y=c(smry$ER.lower, rev(smry$ER.upper)),
          border=NA, col=alpha('sienna', alpha=.5))
  axis(2, las=2, line=0, xpd=NA, tck=-.02, labels=FALSE,
       at=round(seq(-10,ylims[2], by=2), 1))
  #axis(2, las=2, line=-0.5, xpd=NA, tcl=0, col='transparent',
  #     at=round(seq(ylims[1], ylims[2], by=2), 1))
  abline(h=0, lty=1, lwd=1)
  #mtext(expression(paste("g O"[2]~"m"^"-2"~" d"^"-1")), side=2, line=2.5)
  mtext('day of year', side=1, line=-3)
  
  # polygon(c(NHC_smry$DOY, rev(NHC_smry$DOY)),
  #         na.approx(c(NHC_smry$GPP.25, rev(NHC_smry$GPP.75)),na.rm=F),
  #         col=alpha("black",.5), border=NA)
  polygon(c(NHC_smry$DOY, rev(NHC_smry$DOY)),
          c(NHC_smry$GPP.lower[1], rollmean(na.approx(c(NHC_smry$GPP.lower, rev(NHC_smry$GPP.upper)),na.rm=F),3), NHC_smry$GPP.upper[1]),
          col=alpha("black",.5), border=NA)
  
  # polygon(c(NHC_smry$DOY, rev(NHC_smry$DOY)),
  #         na.approx(c(NHC_smry$ER.25, rev(NHC_smry$ER.75)),na.rm=F),
  #         col=alpha("black",.5), border=NA)
  polygon(c(NHC_smry$DOY, rev(NHC_smry$DOY)),
          y=-c(NHC_smry$ER.lower[1], 
              rollmean(na.approx(c(NHC_smry$ER.lower, 
                                   rev(NHC_smry$ER.upper)),na.rm=F),3), 
              NHC_smry$ER.upper[1]),
          col=alpha("black",.5), border=NA)
  
  mtext("Interannual metabolism", 3,  line=.6, font=4, col="grey30")
  mtext("(3 years)", 3, cex=.8, font=4, line=-.2, col="grey30")
  
  
  par(new=T, mfrow=c(1,1), mar=c(0,0,0,0), oma=c(1,0,0,0))
  plot(1,1,xaxt="n", yaxt="n",xlab='',ylab="", type="n", bty="n")
   legend("bottom", legend=c("All Sites GPP (middle 50%)","All Sites ER (middle 50%)", "Local Sites (full range)"),
         fill =c(alpha("forestgreen",.5), alpha("sienna", .5), alpha("black",.5)), 
         border=NA, bty="n", cex=1, xpd=NA)
 dev.off()
  #Interannual variation
  plot(smry$DOY, smry$GPP, ylab='', yaxs='i', type='n',
       bty='n', lwd=4, xlab='', ylim=ylims, xaxs='i', xaxt='n', yaxt='n')
  polygon(x=c(smry$DOY, rev(smry$DOY)),
          y=c(smry$GPP.lower, rev(smry$GPP.upper)),
          border=NA, col=alpha('forestgreen', alpha=.5))
  polygon(x=c(smry$DOY, rev(smry$DOY)),
          y=c(smry$ER.lower, rev(smry$ER.upper)),
          border=NA, col=alpha('sienna', alpha=.5))
  axis(2, las=2, line=0, xpd=NA, tck=-.02, labels=FALSE,
       at=round(seq(ylims[1],ylims[2], by=2), 1))
  #axis(2, las=2, line=-0.5, xpd=NA, tcl=0, col='transparent',
  #     at=round(seq(ylims[1], ylims[2], by=2), 1))
  abline(h=0, lty=1, lwd=1)
  #mtext(expression(paste("g O"[2]~"m"^"-2"~" d"^"-1")), side=2, line=2.5)
  mtext('DOY', side=1, line=2, font=2)
  
  polygon(c(UNHC_smry$DOY, rev(UNHC_smry$DOY)),
          na.approx(c(UNHC_smry$GPP.25, rev(UNHC_smry$GPP.75)),na.rm=F),
          col=alpha("black",.5), border=NA)
  polygon(c(UNHC_smry$DOY, rev(UNHC_smry$DOY)),
          -na.approx(c(UNHC_smry$GPP.lower, rev(UNHC_smry$GPP.upper)),na.rm=F),
          col=alpha("black",.3), border=NA)
  
  polygon(c(UNHC_smry$DOY, rev(UNHC_smry$DOY)),
          na.approx(c(UNHC_smry$ER.25, rev(UNHC_smry$ER.75)),na.rm=F),
          col=alpha("black",.5), border=NA)
  polygon(c(UNHC_smry$DOY, rev(UNHC_smry$DOY)),
          -na.approx(c(UNHC_smry$ER.lower, rev(UNHC_smry$ER.upper)),na.rm=F),
          col=alpha("black",.3), border=NA)
  
  mtext("UNHC interannual (n=3)", 3, -2, font=4, col="grey30")

dev.off()
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

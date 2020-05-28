#############################################################################
# Plot comparison between pools and runs from NHC site ordered and filled metabolism

library(dplyr)
library(tidyverse)

setwd("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism")

metab_ordered <- read_csv("data/metabolism/ordered_NHC_sites_metabolism.csv")
metab_filled <- read_csv("data/metabolism/filled_NHC_sites_metabolism.csv")
sites <- c("UNHC","PWC", "WBP","WB","CBP","PM","NHC")


runs <- c("UNHC","WB","PM")
pools <- c("NHC","WBP","CBP")

# plot parameters:
pool.col = "steelblue"
run.col="grey20"
bor.col="grey30"

avg_pool <- metab_filled[metab_filled$site %in% pools, ] %>% group_by(DOY)%>%
  summarize(GPP=mean(GPP, na.rm=T),
            GPP_upper=mean(GPP.upper, na.rm=T),
            GPP_lower=mean(GPP.lower, na.rm=T),
            ER=mean(ER, na.rm=T),
            ER_upper=mean(ER.upper, na.rm=T),
            ER_lower=mean(ER.lower, na.rm=T))
cum_pool <- as.data.frame(apply(avg_pool, 2, cumsum))
cum_pool$DOY <- 1:366

avg_run <- metab_filled[metab_filled$site %in% runs, ] %>% group_by(DOY)%>%
  summarize(GPP=mean(GPP, na.rm=T),
            GPP_upper=mean(GPP.upper, na.rm=T),
            GPP_lower=mean(GPP.lower, na.rm=T),
            ER=mean(ER, na.rm=T),
            ER_upper=mean(ER.upper, na.rm=T),
            ER_lower=mean(ER.lower, na.rm=T))
cum_run <- as.data.frame(apply(avg_run, 2, cumsum))
cum_run$DOY <- 1:366

png("figures/pool_run_metab_comparison.png", width=600, height=240)
  par(mfrow = c(1,3), mar = c(0,0,0,0), oma=c(5,5,4,5))
  layout(t(c(1,2,3)),widths=c(1.5,1.5,1.2), heights=(1) )
    plot(metab_ordered$DOY, metab_ordered$GPP, ylim = c(-12,4), type="n", xlab="",ylab="",xaxt="n",yaxt="n",col.axis=bor.col)
    for(i in pools){
      tmp<- metab_ordered[metab_ordered$site==i,]
      lines(tmp$DOY, tmp$GPP, col=pool.col, lwd=1.5)
      lines(tmp$DOY, tmp$ER, col=pool.col, lwd=1.5)
      polygon(c(tmp$DOY, rev(tmp$DOY)), na.approx(c(tmp$GPP.lower, rev(tmp$GPP.upper)), na.rm=F), 
              col=alpha(pool.col, .3), border=F)
      polygon(c(tmp$DOY, rev(tmp$DOY)), na.approx(c(tmp$ER.lower, rev(tmp$ER.upper)), na.rm=F), 
              col=alpha(pool.col, .3), border=F)
    }  
      
      
    abline(h=0, col=bor.col )
    mtext("Pools", 3, -2, adj=.8, cex=1.2, col=bor.col)
    axis(2, las=2, line=0, xpd=NA, tck=-.02, labels=FALSE,
         at=round(seq(-12,4, by=2), 1))
    axis(2, las=2, line=-0.5, xpd=NA, tcl=0, col='transparent',
         at=round(seq(-12, 4, by=2), 1))
    mtext(expression(paste("g O"[2]~"m"^"-2"~" d"^"-1")), side=2, line=2.5)
    t<- seq.Date(as.Date("2019-02-01"), as.Date("2019-12-31"),by = "2 months")
    axis(1, at =as.numeric(format(t, "%j")) , labels=FALSE)
    text(x = as.numeric(format(t, "%j")), y = par("usr")[3]-.5 , labels = format(t, "%b"), 
         pos = 1, xpd=NA, cex=1.5)
    
    plot(metab_ordered$DOY, metab_ordered$GPP, ylim = c(-12,4), type="n", xlab="",ylab="",xaxt="n",yaxt="n",col.axis=bor.col)
    for(i in runs){
      tmp<- metab_ordered[metab_ordered$site==i,]
      lines(tmp$DOY, tmp$GPP, col=run.col, lwd=1.5)
      lines(tmp$DOY, tmp$ER, col=run.col, lwd=1.5)
      polygon(c(tmp$DOY, rev(tmp$DOY)), na.approx(c(tmp$GPP.lower, rev(tmp$GPP.upper)), na.rm=F), 
              col=alpha(run.col, .3), border=F)
      polygon(c(tmp$DOY, rev(tmp$DOY)), na.approx(c(tmp$ER.lower, rev(tmp$ER.upper)), na.rm=F), 
              col=alpha(run.col, .3), border=F)
    }  
      
      
    abline(h=0, col=bor.col )
    mtext("Runs", 3, -2, adj=.8, cex=1.2, col=bor.col)
    t<- seq.Date(as.Date("2019-02-01"), as.Date("2019-12-31"),by = "2 months")
    axis(1, at =as.numeric(format(t, "%j")) , labels=FALSE)
    text(x = as.numeric(format(t, "%j")), y = par("usr")[3]-.5 , labels = format(t, "%b"), 
         pos = 1, xpd=NA, cex=1.5)
    
    mtext("Annual metabolism by habitat type", 
          outer=T, font=4, adj=.27, col=bor.col, line=1)
    
  
    par(mar=c(0,1,0,0))
    plot(cum_run$DOY, cum_run$GPP, ylim = c(-1050,350),type="l", lwd=2, col=run.col,
         xlab="",ylab="",xaxt="n",yaxt="n", col.axis=bor.col)  
    lines(cum_run$DOY, cum_run$ER, col=run.col, lwd=2)
    polygon(c(cum_run$DOY, rev(cum_run$DOY)), c(cum_run$GPP_lower, rev(cum_run$GPP_upper)), 
            col=alpha(run.col, .3), border=F)
    polygon(c(cum_run$DOY, rev(cum_run$DOY)), c(cum_run$ER_lower, rev(cum_run$ER_upper)), 
            col=alpha(run.col, .3), border=F)
    lines(cum_pool$DOY, cum_pool$GPP, lwd=2, col=pool.col)  
    lines(cum_pool$DOY, cum_pool$ER, col=pool.col, lwd=2)
    polygon(c(cum_pool$DOY, rev(cum_pool$DOY)), c(cum_pool$GPP_lower, rev(cum_pool$GPP_upper)), 
            col=alpha(pool.col, .3), border=F)
    polygon(c(cum_pool$DOY, rev(cum_pool$DOY)), c(cum_pool$ER_lower, rev(cum_pool$ER_upper)), 
            col=alpha(pool.col, .3), border=F)
    
    abline(h=0, col=bor.col )
    legend("bottomleft", legend=c("Pool","Run"), col=c(pool.col,run.col),
           lwd=2, cex=1.5, bty="n", lty=c(1,1))
    
    t<- seq.Date(as.Date("2019-02-01"), as.Date("2019-12-31"),by = "3 months")
    axis(1, at =as.numeric(format(t, "%j")) , labels=FALSE)
    text(x = as.numeric(format(t, "%j")), y = par("usr")[3]-45 , labels = format(t, "%b"), 
         pos = 1, xpd=NA, cex=1.5)
    
    axis(4, las=2, line=0, xpd=NA, tck=-.02, labels=FALSE,
         at=round(seq(-900,300, by=200), 1))
    axis(4, las=2, line=-0.5, xpd=NA, tcl=0, col='transparent',
         at=round(seq(-900,300, by=200), 1))
    mtext(expression(paste("Cumulative g O"[2]~"m"^"-2")), side=4, line=3)
    mtext("Cumulative Metabolism", 
           font=4,  col=bor.col, line=1)
dev.off()    

  
pool_m <- metab_filled[metab_filled$site %in% pools,] %>% group_by(DOY)%>%
    summarize(GPP_m=mean(GPP, na.rm=T),
            GPP.975=quantile(GPP, .975, na.rm = T),
            GPP.75=quantile(GPP, .75, na.rm=T),
            GPP.25=quantile(GPP, .25, na.rm = T),
            GPP.025=quantile(GPP, .025, na.rm=T),
            ER_m=mean(ER, na.rm=T),
            ER.975=quantile(ER, .975, na.rm = T),
            ER.75=quantile(ER, .75, na.rm=T),
            ER.25=quantile(ER, .25, na.rm = T),
            ER.025=quantile(ER, .025, na.rm=T)
            )



run_m <- metab_filled[metab_filled$site %in% runs,] %>% group_by(DOY)%>%
  summarize(GPP_m=mean(GPP, na.rm=T),
            GPP.975=quantile(GPP, .975, na.rm = T),
            GPP.75=quantile(GPP, .75, na.rm=T),
            GPP.25=quantile(GPP, .25, na.rm = T),
            GPP.025=quantile(GPP, .025, na.rm=T),
            ER_m=mean(ER, na.rm=T),
            ER.975=quantile(ER, .975, na.rm = T),
            ER.75=quantile(ER, .75, na.rm=T),
            ER.25=quantile(ER, .25, na.rm = T),
            ER.025=quantile(ER, .025, na.rm=T)
            )


plot(pool_m$DOY, pool_m$GPP_m, ylim = c(-12,4), type="l", xaxt="n",yaxt="n", xlab="",ylab="")
lines(pool_m$DOY,pool_m$ER_m)
polygon(c(pool_m$DOY,rev(pool_m$DOY)), na.approx(c(pool_m$GPP.25, rev(pool_m$GPP.75))), col=alpha("black",.4), border=F)
polygon(c(pool_m$DOY,rev(pool_m$DOY)), na.approx(c(pool_m$GPP.025, rev(pool_m$GPP.975))), col=alpha("black",.3), border=F)

polygon(c(pool_m$DOY,rev(pool_m$DOY)), na.approx(c(pool_m$ER.25, rev(pool_m$ER.75))), col=alpha("black",.4), border=F)
polygon(c(pool_m$DOY,rev(pool_m$DOY)), na.approx(c(pool_m$ER.025, rev(pool_m$ER.975))), col=alpha("black",.3), border=F)


plot(run_m$DOY, run_m$GPP_m,ylim=c(-12,4), col="steelblue")
lines(run_m$DOY,run_m$ER_m, col="steelblue")
polygon(c(run_m$DOY,rev(run_m$DOY)), na.approx(c(run_m$GPP.25, rev(run_m$GPP.75))), col=alpha("steelblue",.4), border=F)
polygon(c(run_m$DOY,rev(run_m$DOY)), na.approx(c(run_m$GPP.025, rev(run_m$GPP.975))), col=alpha("steelblue",.3), border=F)

polygon(c(run_m$DOY,rev(run_m$DOY)), na.approx(c(run_m$ER.25, rev(run_m$ER.75))), col=alpha("steelblue",.4), border=F)
polygon(c(run_m$DOY,rev(run_m$DOY)), na.approx(c(run_m$ER.025, rev(run_m$ER.975))), col=alpha("steelblue",.3), border=F)


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

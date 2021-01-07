#############################################################################
# Plot comparison between pools and runs from NHC site ordered and filled metabolism

library(dplyr)
library(tidyverse)
library(zoo)

setwd("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism")

metab_ordered <- read_csv("data/metabolism/ordered_NHC_sites_metabolism.csv")
metab_filled <- read_csv("data/metabolism/filled_NHC_sites_metabolism.csv")
# sites <- c("UNHC","PWC", "WBP","WB","CBP","PM","NHC")
dat <- readRDS("data/metabolism/compiled/raymond_met.rds")
metab_ordered <- dat$preds %>%
  as_tibble() %>%
  mutate(DOY = format(date, "%j")) %>%
  rename(discharge.m3s = discharge.daily,
         K600.lower = K600_2.5,
         K600.upper = K600_97.5)
K <- select(metab_ordered, date, site, K600, K600.upper, K600.lower)
metab_filled <- dat$cumulative %>%
  as_tibble() %>%
  select(date, GPP, GPP.upper, GPP.lower,
         ER, ER.upper, ER.lower,
         discharge.m3s = discharge.daily, site) %>%
  left_join(K, by = c("site", "date")) %>%
  mutate(DOY = as.numeric(format(date, "%j"))) %>%
  group_by(site, DOY) %>%
  select(-date) %>%
  summarize_all(mean, na.rm = T)

metab_across <- dat$cumulative %>%
  as_tibble() %>%
  select(date, site, year, GPP) %>%
  mutate(DOY = as.numeric(format(date, "%j"))) %>%
  group_by(site,year, DOY) %>%
  select(-date) %>%
  summarize_all(mean, na.rm = T) %>%
  pivot_wider(id_cols = c("year", "DOY"), names_from = "site",
              values_from = "GPP" , names_prefix = "gpp_")
metab_across <- dat$cumulative %>%
  as_tibble() %>%
  select(date, site, year, ER) %>%
  mutate(DOY = as.numeric(format(date, "%j"))) %>%
  group_by(site,year, DOY) %>%
  select(-date) %>%
  summarize_all(mean, na.rm = T) %>%
  pivot_wider(id_cols = c("year", "DOY"), names_from = "site",
              values_from = "ER" , names_prefix = "er_") %>%
  left_join(metab_across) %>%
  filter(year == 2019) %>%arrange(DOY)

tmp <- dat$cumulative %>%
  filter(site %in% c("nhc","unhc")) %>% 
  as_tibble() %>% 
  select(date, GPP, ER, Q = discharge.daily, site) %>% 
  pivot_wider(names_from = "site", values_from = c("GPP","ER", "Q")) %>%
  arrange(date)

plot(tmp$date, -tmp$ER_nhc, type = 'l', ylim = c(-20, 4))
lines(tmp$date, tmp$GPP_nhc)
lines(tmp$date, tmp$GPP_unhc, col = "brown3")
lines(tmp$date, -tmp$ER_unhc, col = "brown3")
par(new = T)
plot(tmp$date, tmp$Q_nhc, ylim = c(.01, 500000), type = "l", log = "y")
plot(metab_across$DOY, metab_across$er_nhc, type = "l")#, ylim = c(-12,3.5))
lines(metab_across$DOY, metab_across$er_unhc, col = 2)
lot(metab_across, aes(x = er_nhc, col = DOY)) +
  geom_point(aes(y=er_wbp), size = 2) +
  xlim(0,10) +
  ylim(0,10)
  geom_point(aes(y = gpp_cbp), col = 2) +
  geom_point(aes(y = gpp_wb), col = 3) +
  geom_point(aes(y = gpp_wbp), col = 4) +
  geom_point(aes(y = gpp_unhc), col = 5) 


sites = c("nhc", "pm", "cbp","wb","wbp","unhc")

LAI <- read_csv("data/light/NC_UNHC_LAI_2017.csv")

runs <- c("unhc","wb","nhc")
pools <- c("pm","wbp","cbp")

# plot parameters:
pool.col = "#91C3DC"
run.col="#555555"
bor.col="grey30"

avg_pool <- metab_filled[metab_filled$site %in% pools, ] %>% group_by(DOY)%>%
  summarize(GPP=mean(GPP, na.rm=T),
            GPP_upper=mean(GPP.upper, na.rm=T),
            GPP_lower=mean(GPP.lower, na.rm=T),
            ER=mean(ER, na.rm=T),
            ER_upper=mean(ER.upper, na.rm=T),
            ER_lower=mean(ER.lower, na.rm=T)) 
cum_pool <- avg_pool %>%
  mutate(across(-DOY, cumsum))

avg_run <- metab_filled[metab_filled$site %in% runs, ] %>% group_by(DOY)%>%
  summarize(GPP=mean(GPP, na.rm=T),
            GPP_upper=mean(GPP.upper, na.rm=T),
            GPP_lower=mean(GPP.lower, na.rm=T),
            ER=mean(ER, na.rm=T),
            ER_upper=mean(ER.upper, na.rm=T),
            ER_lower=mean(ER.lower, na.rm=T))
cum_run <- avg_run %>%
  mutate(across(-DOY, cumsum))

data.frame(apply(avg_run, 2, cumsum))data.frame(apply(avg_pool, 2, cumsum)) %>% as_tibble()cum_pool$DOY <- 1:366cum_run$DOY <- 1:366
png("figures/pool_run_metab_comparison.png", width=600, height=240)
  metab_ordered <- metab_filled
  par(mfrow = c(1,3), mar = c(0,0,0,0), oma=c(5,5,4,5))
  layout(t(c(1,2,3)),widths=c(1.5,1.5,1.2), heights=(1) )
    plot(metab_ordered$DOY, metab_ordered$GPP, ylim = c(-12,4), 
         type="n", xlab="",ylab="",xaxt="n",yaxt="n",col.axis=bor.col)
    for(i in pools){
      tmp<- metab_ordered[metab_ordered$site==i,]
      lines(tmp$DOY, tmp$GPP, col=pool.col, lwd=1.5)
      lines(tmp$DOY, -tmp$ER, col=pool.col, lwd=1.5)
      polygon(c(tmp$DOY, rev(tmp$DOY)), 
              na.approx(c(tmp$GPP.lower, rev(tmp$GPP.upper)), na.rm=F), 
              col=alpha(pool.col, .3), border=F)
      polygon(c(tmp$DOY, rev(tmp$DOY)),
              na.approx(-c(tmp$ER.lower, rev(tmp$ER.upper)), na.rm=F), 
              col=alpha(pool.col, .3), border=F)
    }  
      
      
    abline(h=0, col=bor.col )
    mtext("Pools (n=3)", 3, -2, adj=.9, cex=1.2, col=bor.col)
    axis(2, las=2, line=0, xpd=NA, tck=-.02, labels=FALSE,
         at=round(seq(-12,4, by=2), 1))
    axis(2, las=2, line=-0.5, xpd=NA, tcl=0, col='transparent',
         at=round(seq(-12, 4, by=2), 1))
    mtext(expression(paste("g O"[2]~"m"^"-2"~" d"^"-1")), side=2, line=2.5)
    t<- seq.Date(as.Date("2019-02-01"), as.Date("2019-12-31"),by = "2 months")
    axis(1, at =as.numeric(format(t, "%j")) , labels=FALSE)
    text(x = as.numeric(format(t, "%j")), y = par("usr")[3]-.5 , labels = format(t, "%b"), 
         pos = 1, xpd=NA, cex=1.5)
    
    plot(metab_ordered$DOY, metab_ordered$GPP, ylim = c(-12,4), 
         type="n", xlab="",ylab="",xaxt="n",yaxt="n",col.axis=bor.col)
    for(i in runs){
      tmp<- metab_ordered[metab_ordered$site==i,]
      lines(tmp$DOY, tmp$GPP, col=run.col, lwd=1.5)
      lines(tmp$DOY, -tmp$ER, col=run.col, lwd=1.5)
      polygon(c(tmp$DOY, rev(tmp$DOY)), 
              na.approx(c(tmp$GPP.lower, rev(tmp$GPP.upper)), na.rm=F), 
              col=alpha(run.col, .3), border=F)
      polygon(c(tmp$DOY, rev(tmp$DOY)), 
              na.approx(-c(tmp$ER.lower, rev(tmp$ER.upper)), na.rm=F), 
              col=alpha(run.col, .3), border=F)
    }  
      
      
    abline(h=0, col=bor.col )
    mtext("Runs (n=3)", 3, -2, adj=.9, cex=1.2, col=bor.col)
    t<- seq.Date(as.Date("2019-02-01"), as.Date("2019-12-31"),by = "2 months")
    axis(1, at =as.numeric(format(t, "%j")) , labels=FALSE)
    text(x = as.numeric(format(t, "%j")), y = par("usr")[3]-.5 , 
         labels = format(t, "%b"), pos = 1, xpd=NA, cex=1.5)
    
    mtext("Annual metabolism by habitat type", 
          outer=T, font=4, adj=.27, col=bor.col, line=1)
    
  
    par(mar=c(0,1,0,0))
    plot(cum_run$DOY, cum_run$GPP, ylim = c(-1200,400),type="l", lwd=2, col=run.col,
         xlab="",ylab="",xaxt="n",yaxt="n", col.axis=bor.col)  
    lines(cum_run$DOY, -cum_run$ER, col=run.col, lwd=2)
    polygon(c(cum_run$DOY, rev(cum_run$DOY)), 
            na.approx(c(cum_run$GPP_lower, rev(cum_run$GPP_upper)), na.rm=F), 
            col=alpha(run.col, .3), border=F)
    polygon(c(cum_run$DOY, rev(cum_run$DOY)), 
            na.approx(-c(cum_run$ER_lower, rev(cum_run$ER_upper)), na.rm=F), 
            col=alpha(run.col, .3), border=F)
    lines(cum_pool$DOY, cum_pool$GPP, lwd=2, col=pool.col)  
    lines(cum_pool$DOY, -cum_pool$ER, col=pool.col, lwd=2)
    polygon(c(cum_pool$DOY, rev(cum_pool$DOY)), 
            na.approx(c(cum_pool$GPP_lower, rev(cum_pool$GPP_upper)),na.rm=F), 
            col=alpha(pool.col, .3), border=F)
    polygon(c(cum_pool$DOY, rev(cum_pool$DOY)), 
            -na.approx(c(cum_pool$ER_lower, rev(cum_pool$ER_upper)),na.rm=F), 
            col=alpha(pool.col, .3), border=F)
    
    abline(h=0, col=bor.col )
    legend("bottomleft", legend=c("Pool","Run"), col=c(pool.col,run.col),
           lwd=2, cex=1.5, bty="n", lty=c(1,1))
    
    t<- seq.Date(as.Date("2019-02-01"), as.Date("2019-12-31"),by = "3 months")
    axis(1, at =as.numeric(format(t, "%j")) , labels=FALSE)
    text(x = as.numeric(format(t, "%j")), y = par("usr")[3]-45 , 
         labels = format(t, "%b"), 
         pos = 1, xpd=NA, cex=1.5)
    
    axis(4, las=2, line=0, xpd=NA, tck=-.02, labels=FALSE,
         at=round(seq(-1200,400, by=200), 1))
    axis(4, las=2, line=-0.5, xpd=NA, tcl=0, col='transparent',
         at=round(seq(-1200,400, by=200), 1))
    mtext(expression(paste("Cumulative g O"[2]~"m"^"-2")), side=4, line=3)
          
    mtext("Cumulative Metabolism",  font=4,  col=bor.col , line=1)
dev.off()    

  
###################################################################################
# Cumulative panel only

cum_sites <- data.frame()

for(i in sites){
  sitedat<- metab_filled[metab_filled$site==i,]%>% select(-K600, -K600.upper, -K600.lower)
  tmp<- as.data.frame(apply(sitedat, 2, cumsum))
  tmp$DOY <- 1:nrow(tmp)
  tmp$site <- i
  cum_sites <- rbind(cum_sites, tmp)
}

png(filename="figures/cumulative_pool_run_plot.png", width=450, height=450)
  par(mar=c(3,4,2,4), mfrow=c(1,1))
  plot(1,1, ylim = c(0, 28), xlim = c(1,366),
       xlab="",ylab="",xaxt="n",yaxt="n", bty="n", type="n")
  polygon(c(LAI$DOY, 366,1), c(LAI$LAI_proc, .5,.5), border=NA, 
          col = alpha("forestgreen", .3))
  mtext("Leaf Area",side=1, line=-7,  adj=.53,col=alpha("forestgreen",.8 ))
  # axis(4, las=2, line=0, xpd=NA, tck=-.02, labels=FALSE,
  #      at=round(seq(0,100, by=20), 1))
  # axis(4, las=2, line=-0.5, xpd=NA, tcl=0, col='transparent',
  #      at=round(seq(0,100, by=20), 1))
  # mtext("LAI", side=4, line=2.5, adj=.87)
  
  
  par(new=T)
  plot(cum_sites$DOY, cum_sites$GPP, xlim=c(0,366),
       ylim = c(-1400,300),type="n",
       xlab="",ylab="",xaxt="n",yaxt="n", bty="n")  
  for(i in runs){
    tmp <- cum_sites[cum_sites$site==i,]
    lines(tmp$DOY, -tmp$ER, col=run.col, lwd=2)
    lines(tmp$DOY, tmp$GPP, col=run.col, lwd=2)
  }
  for(i in pools){
    tmp <- cum_sites[cum_sites$site==i,]
    lines(tmp$DOY, -tmp$ER, col=pool.col, lwd=2)
    lines(tmp$DOY, tmp$GPP, col=pool.col, lwd=2)
  }
    
  abline(h=0, col=bor.col)
  t<- seq.Date(as.Date("2019-01-01"), as.Date("2020-01-05"),by = "2 months")
  axis(1, at =c(as.numeric(format(t, "%j")),366) , labels=FALSE)
  text(x = c(as.numeric(format(t, "%j")),366), y = par("usr")[3]-45 , 
       labels = format(t, "%b"), 
        pos = 1, xpd=NA, cex=1)
  axis(2, las=2, line=0, xpd=NA, tck=-.02, labels=FALSE,
       at=round(seq(-800,300, by=200), 1))
  axis(2, las=2, line=-0.5, xpd=NA, tcl=0, col='transparent',
       at=round(seq(-800,300, by=200), 1))
  mtext(expression(paste("Cumulative g O"[2]~"m"^"-2")), side=2, line=2.5, adj=.7)
  
  legend(-15,-220, legend=c("Pool","Run"), col=c(pool.col,run.col),
         lwd=2, cex=1, bty="n", lty=c(1,1))
  
  par(new=T)
  plot(metab_filled[metab_filled$site=="nhc",]$DOY,
       metab_filled[metab_filled$site=="nhc",]$discharge.m3s,
       ylim = c(min(metab_filled[metab_filled$site=="nhc",]$discharge.m3s, na.rm=T), 1e17),
       log="y", xaxt="n", yaxt="n", xlab="", ylab="", bty="n", type="l")
  axis(2, las=2, line=0, xpd=NA, tck=-.02, labels=FALSE,
       at=c(0.01, 1, 100))
  text(x = -14, y = c(0.01, 1, 100) , labels = c("0.01", "1", "100"), 
       pos = 2, xpd=NA, cex=1)
  mtext(expression(paste("discharge m"^"3"~"/s")), side=2, line=2.5, adj=0)
  
  abline(v=c(260, 300), lty=2, col = "grey50")
  abline(v=c(75, 120), lty=2, col = "grey50")
  
  mtext("litter\nfall ", 3, adj=.765, line=-.8, cex=.8, col="grey30")
  mtext(" spring\nwindow", 3, adj=.26, line=-.8, cex=.8, col="grey30")

dev.off()





polygon(c(cum_run$DOY, rev(cum_run$DOY)), na.approx(c(cum_run$ER_lower, rev(cum_run$ER_upper)), na.rm=F), 
        col=alpha(run.col, .3), border=F)
lines(cum_pool$DOY, cum_pool$GPP, lwd=2, col=pool.col)  
lines(cum_pool$DOY, cum_pool$ER, col=pool.col, lwd=2)
polygon(c(cum_pool$DOY, rev(cum_pool$DOY)), na.approx(c(cum_pool$GPP_lower, rev(cum_pool$GPP_upper)),na.rm=F), 
        col=alpha(pool.col, .3), border=F)
polygon(c(cum_pool$DOY, rev(cum_pool$DOY)), na.approx(c(cum_pool$ER_lower, rev(cum_pool$ER_upper)),na.rm=F), 
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

mtext("Cumulative Metabolism",  font=4,  col=bor.col , line=1)




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


plot(pool_m$DOY, pool_m$GPP_m, ylim = c(-12,4), 
     type="l", xaxt="n",yaxt="n", xlab="",ylab="")
lines(pool_m$DOY,-pool_m$ER_m)
polygon(c(pool_m$DOY,rev(pool_m$DOY)), 
        na.approx(c(pool_m$GPP.25, rev(pool_m$GPP.75))), 
        col=alpha("black",.4), border=F)
polygon(c(pool_m$DOY,rev(pool_m$DOY)), 
        na.approx(c(pool_m$GPP.025, rev(pool_m$GPP.975))), 
        col=alpha("black",.3), border=F)

polygon(c(pool_m$DOY,rev(pool_m$DOY)), 
        -na.approx(c(pool_m$ER.25, rev(pool_m$ER.75))), 
        col=alpha("black",.4), border=F)
polygon(c(pool_m$DOY,rev(pool_m$DOY)), 
        -na.approx(c(pool_m$ER.025, rev(pool_m$ER.975))), 
        col=alpha("black",.3), border=F)


plot(run_m$DOY, run_m$GPP_m,ylim=c(-12,4), type = "l", col="steelblue")
lines(run_m$DOY,-run_m$ER_m, col="steelblue")
polygon(c(run_m$DOY,rev(run_m$DOY)), 
        na.approx(c(run_m$GPP.25, rev(run_m$GPP.75))), 
        col=alpha("steelblue",.4), border=F)
polygon(c(run_m$DOY,rev(run_m$DOY)), 
        na.approx(c(run_m$GPP.025, rev(run_m$GPP.975))), 
        col=alpha("steelblue",.3), border=F)

polygon(c(run_m$DOY,rev(run_m$DOY)), 
        -na.approx(c(run_m$ER.25, rev(run_m$ER.75))), 
        col=alpha("steelblue",.4), border=F)
polygon(c(run_m$DOY,rev(run_m$DOY)), 
        -na.approx(c(run_m$ER.025, rev(run_m$ER.975))), 
        col=alpha("steelblue",.3), border=F)


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
k_all <- kde(na.omit(dat$cumulative[, c("GPP", "ER")]))
k_nhc <- kde(na.omit(metab[,c("GPP","ER")]))
#Plotting the kernel density plot
plot(k_all, xlab = "GPP", ylab = "ER", ylim = c(0, 8), 
     xlim = c(0, 8),display = "slice")

par(new=T)
plot(k_nhc, xaxt='n',yaxt='n', ylab='',xlab='',ylim = c(-15, 0), xlim = c(0, 15),
     display = "filled.contour")
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

#########################################
# read in and plot metabolism for NHC sites
# AMCarter 4-17-20

library(streamMetabolizer)
library(dplyr)
library(scales)
library(zoo)

# load RDS from model output:

filelist <- list.files("../data/metabolism/modeled")

pdf("../data/figures/NHC2019metabolim2.pdf",width = 7, height = 5.83,onefile = TRUE)
for(j in c(1,3,5,7,9)){
  i <- filelist[j]
  dat <- readRDS(paste0("../data/metabolism/modeled/",i))
  dat_metab <- dat$predictions
  dat_daily <- dat$fit@data_daily
  dat_dat <- dat$fit@data
  
  # remove Metabolism data from days where flow is above 99 percentile:
  Qmax <- quantile(dat_daily$discharge.daily, .98, na.rm=T)
  w<- which(dat_daily$discharge.daily>=Qmax)
  dat_metab[w,2:7]<- NA
  dat_metab[w,9] <- "point removed for high Q"
  
  par(mfrow=c(2,1))
  layout.matrix <- matrix(c(1,2,3), nrow = 3, ncol = 1)
  
  layout(mat = layout.matrix,
         heights = c(2,4,1), # Heights of the two rows
         widths = 6) # Widths of the two columns
  
   par(mar = c(0,4,1,1))
  plot(dat_dat$date, dat_dat$DO.obs/dat_dat$DO.sat, ylim = c(0,1.2), ylab = "DO %sat", type = "l", xaxt="n")
  par(mar = c(0,4,0,1))
   plot(dat_metab$date, dat_metab$GPP, ylim = c(-8,4),#c(min(dat_metab$ER.lower, na.rm=T), max(dat_metab$GPP.upper, na.rm=T)),
       col = "forestgreen", lwd = 2, xlab = "", ylab = "(g O2/m2/d)", type="l", xaxt="n")
    polygon(c(dat_metab$date, rev(dat_metab$date)), na.approx(c(dat_metab$GPP.lower, rev(dat_metab$GPP.upper)), na.rm=FALSE), 
            col = alpha("forestgreen", .3), border=NA)
    lines(dat_metab$date, dat_metab$ER, col = "brown3", lwd = 2)          
    polygon(c(dat_metab$date, rev(dat_metab$date)), na.approx(c(dat_metab$ER.lower, rev(dat_metab$ER.upper)), na.rm=FALSE), 
            col = alpha("brown3", .3), border=NA)
    mtext(paste0("Metabolism for ",substr(i, 1, nchar(i)-4)),3, -2)
    abline(h=0)
  par(mar = c(3,4,0,1))
    plot(dat$fit@data$solar.time, dat$fit@data$discharge, log="y", ylab = "Q (m3s)", xlab = "date", type = "l", xaxt="n") 
    par(new=T)
    plot(dat_metab$date, dat_metab$GPP, type="n", axes=FALSE, ylab = "")
    t<- seq.Date(as.Date("2019-04-01"), dat_metab$date[nrow(dat_metab)],by = "month")
    axis(1, at =t , labels=FALSE)
    text(x = t, y = par("usr")[3]-.5 , labels = format(t, "%b"), 
         pos = 1, xpd=NA)
}
dev.off()

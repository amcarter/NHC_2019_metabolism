#########################################
# read in and plot metabolism for NHC sites
# AMCarter 4-17-20

library(streamMetabolizer)
library(dplyr)
library(scales)
library(zoo)
library(xts)
library(dygraphs)
library(readr)

source("code/helpers.R")

setwd("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism")

# load RDS from model output:
hall_met <- read_csv("../Hall_50yearslater_NHC/hall_table_15.csv")
sitematch <- data.frame(hall = c("Concrete", "Blackwood", "Wood Bridge"), 
                        sitecode= c("CBP.rds", "UNHC.rds", "WB.rds"))
filelist <- list.files("data/metabolism/modeled")
###############################
# DO, metabolism, and Q for each site
# plot all figures in a single PDF
pdf("figures/NHC2019metabolim2.pdf",width = 7, height = 5.83,onefile = TRUE)
for(j in c(1:5)){
  i <- filelist[j]
  dat <- readRDS(paste0("data/metabolism/modeled/",i))
  dat_metab <- dat@fit$daily%>% 
    select(date, GPP_mean, GPP_2.5pct, GPP_97.5pct, ER_mean, ER_2.5pct, ER_97.5pct)
  dat_daily <- dat@data_daily 
  dat_dat <- dat@data
  
 
  
  # remove Metabolism data from days where flow is above 99 percentile:
  Qmax <- quantile(dat_daily$discharge.daily, .98, na.rm=T)
  w<- which(dat_daily$discharge.daily>=Qmax)
  dat_metab[w,2:7]<- NA
  
  par(mfrow=c(2,1))
  layout.matrix <- matrix(c(1,2,3), nrow = 3, ncol = 1)
  
  layout(mat = layout.matrix,
         heights = c(2,4,2), # Heights of the two rows
         widths = 6) # Widths of the two columns
  
   par(mar = c(0,4,1,1))
  plot(dat_dat$date, dat_dat$DO.obs/dat_dat$DO.sat, ylim = c(0,1.2), ylab = "DO %sat", type = "l", xaxt="n")
  par(mar = c(0,4,0,1))
   plot(dat_metab$date, dat_metab$GPP_mean, ylim = c(-8,4),#c(min(dat_metab$ER.lower, na.rm=T), max(dat_metab$GPP.upper, na.rm=T)),
       col = "forestgreen", lwd = 2, xlab = "", ylab = "(g O2/m2/d)", type="l", xaxt="n")
    polygon(c(dat_metab$date, rev(dat_metab$date)), na.approx(c(dat_metab$GPP_2.5pct, rev(dat_metab$GPP_97.5pct)), na.rm=FALSE), 
            col = alpha("forestgreen", .3), border=NA)
    lines(dat_metab$date, dat_metab$ER_mean, col = "brown3", lwd = 2)          
    polygon(c(dat_metab$date, rev(dat_metab$date)), na.approx(c(dat_metab$ER_2.5pct, rev(dat_metab$ER_97.5pct)), na.rm=FALSE), 
            col = alpha("brown3", .3), border=NA)
    mtext(paste0("Metabolism for ",substr(i, 1, nchar(i)-4)),3, -2)
    abline(h=0)
  
    #get hall data
    if(i %in% sitematch$sitecode){
      hall <- sitematch$hall[sitematch$sitecode==i]
      hall_dat <- hall_met[hall_met$site==hall,]
      hall_dat$newdate <- base::as.Date(hall_dat$newdate, format="%m/%d/%Y")
      points(hall_dat$newdate, hall_dat$GPP_gO2m2d, col = "forestgreen", pch=20, cex=2)
      points(hall_dat$newdate, -hall_dat$ER_gO2m2d, col = "brown3", pch=20, cex=2)
    }
  par(mar = c(3,4,0,1))
    plot(dat_daily$date, dat_daily$discharge.daily, log="y", ylab = "Q (m3s)", xlab = "date", type = "l", xaxt="n") 
    par(new=T)
    plot(dat_metab$date, dat_metab$GPP_mean, type="n", axes=FALSE, ylab = "")
    t<- seq.Date(as.Date("2019-04-01"), dat_metab$date[nrow(dat_metab)],by = "month")
    axis(1, at =t , labels=FALSE)
    text(x = t, y = par("usr")[3]-.5 , labels = format(t, "%b"), 
         pos = 1, xpd=NA)
}
dev.off()
####################################################

#look at output from one model

mod <- readRDS("data/metabolism/modeled/NHC.rds")
plot_metab_preds(mod)
DOpreds <- predict_DO(mod)
Al_plot_DO_preds(DOpreds, y_var = "conc", style="dygraphs")


####################################################

# Four panel figure with:
# (1) Metabolism timeseries
# (2) Metabolism fingerprint
# (3) Characteristic storm response
# (4) GPP vs light

# specs for plotting
par(mfrow=c(2,2), mar = c(3,4,0,0), oma = c(1,1,1,1))
sites <- c("UNHC","PWC","WBP", "WB", "CBP", "PM", "NHC")
ylims=c(-15,5)
for(site in sites)
# Panel 1: Met timeseries 
metab_dat <- readRDS(paste0("data/metabolism/condensed/condensed_",site,".rds"))
dat <- metab_dat$metab

# Remove metab for top 5% of flow days
w<- which(dat$discharge.m3s>quantile(dat$discharge.m3s, .95, na.rm=TRUE))
dat[w,2:7]<- NA

plot(dat$date, dat$GPP, ylim = ylims, col = "forestgreen", lwd = 2, xlab = "", ylab = "(g O2/m2/d)", type="l")
polygon(c(dat$date, rev(dat$date)), na.approx(c(dat$GPP.lower, rev(dat$GPP.upper)), na.rm=FALSE), 
        col = alpha("forestgreen", .3), border=NA)
lines(dat$date, dat$ER, col = "brown3", lwd = 2)          
polygon(c(dat$date, rev(dat$date)), na.approx(c(dat$ER.lower, rev(dat$ER.upper)), na.rm=FALSE), 
        col = alpha("brown3", .3), border=NA)
abline(h=0)

#Panel 2: metab fingerprint

# Panel 3: storm response
stormdat <- read_csv('data/stormmet.csv')
subset <- stormdat[stormdat$site==site,]

# Remove metab for top 5% of flow days
# w<- which(subset$discharge_m3s>quantile(subset$discharge_m3s, .975, na.rm=TRUE))
# subset[w,7:12]<- NA
# how many days do you want to look at?
predays <- 2
postdays <- 5

storms = which(subset$storm==1)
starts <- storms-predays
subset$stormday <- NA
for(i in 1:length(starts)){
  n <-starts[i]
  subset$stormday[n:(n+predays+postdays)]<- seq(-2,5, by=1)
}
storms<- subset[!is.na(subset$stormday),]

plot( seq(-2, 5),storms$GPP[starts[1]:(starts[1]+7)],ylim=ylims, 
        ylab="gO2/m2/d", xlab="days since storm", type="n")
  abline(v=0, col="gray", lty=2, lwd=2)
  abline(h=0)
  for(j in 1:length(starts)){
    n <- starts[j]
    lines(seq(-2,5), subset$GPP[n:(n+7)], col = alpha("forestgreen",.5), lwd=2)
    lines(seq(-2,5), subset$ER[n:(n+7)], col = alpha("brown3",.5), lwd=2)
  }
}


# panel 4:GPP vs light
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


setwd("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism")

#source("code/helpers.R")
# load RDS from model output:
hall_met <- read_csv("../Hall_50yearslater_NHC/hall_table_15.csv")
sitematch <- data.frame(hall = c("Concrete", "Blackwood", "Wood Bridge"), 
                        sitecode= c("CBP", "UNHC", "WB"))
NHCdat <- readRDS("data/metabolism/condensed/allNHCsites.rds")
metab <- NHCdat$metab
metab[!is.na(metab$GPP)&metab$GPP<0,c("GPP","GPP.upper","GPP.lower")]<-NA
metab[!is.na(metab$ER)&metab$ER>0,c("ER","ER.upper","ER.lower")]<-NA

data <- NHCdat$data
sites <- c("UNHC","WBP","WB","CBP","PM","NHC")
###############################
# DO, metabolism, and Q for each site
# plot all figures in a single PDF
pdf("figures/NHC2019metabolim.pdf",width = 7, height = 5.83,onefile = TRUE)
for(site in sites){
  s_metab<- metab[metab$site==site,]
  s_metab[!is.na(s_metab$discharge.m3s)&s_metab$discharge.m3s>quantile(s_metab$discharge.m3s, .95, na.rm=T),
          c("GPP","GPP.upper","GPP.lower","ER","ER.upper","ER.lower")]<-NA
  sdat<- data[data$site==site,]
  
  par(mfrow=c(2,1))
  layout.matrix <- matrix(c(1,2), nrow = 2, ncol = 1)
  
  layout(mat = layout.matrix,
         heights = c(5,2), # Heights of the two rows
         widths = 6) # Widths of the two columns
  
   par(mar = c(0,4,1,4))
   plot(s_metab$date, s_metab$GPP, ylim = c(-15,15),#c(min(dat_metab$ER.lower, na.rm=T), max(dat_metab$GPP.upper, na.rm=T)),
       col = "grey30", lwd = 2, xlab = "", ylab = "(g O2/m2/d)", type="l", xaxt="n")
    polygon(c(s_metab$date, rev(s_metab$date)), na.approx(c(s_metab$GPP.lower, rev(s_metab$GPP.upper)), na.rm=FALSE), 
            col = alpha("forestgreen", .5), border=NA)
    lines(s_metab$date, s_metab$ER, col = "grey30", lwd = 2)          
    polygon(c(s_metab$date, rev(s_metab$date)), na.approx(c(s_metab$ER.lower, rev(s_metab$ER.upper)), na.rm=FALSE), 
            col = alpha("sienna", .5), border=NA)
    mtext(paste0("Metabolism for ",site),3,-1)
    abline(h=0)
  
    #get hall data
    if(site %in% sitematch$sitecode){
      hall <- sitematch$hall[sitematch$sitecode==site]
      hall_dat <- hall_met[hall_met$site==hall,]
      hall_dat$newdate <- base::as.Date(hall_dat$newdate, format="%m/%d/%Y")
      points(hall_dat$newdate, hall_dat$GPP_gO2m2d, bg = "forestgreen", col="black",pch=21, cex=1.5, lwd=2)
      points(hall_dat$newdate, -hall_dat$ER_gO2m2d, col= "black",bg= "sienna", pch=21, cex=1.5, lwd=2)
    } 
  
  par(new=T)
  plot(sdat$date, 100*sdat$DO.obs/sdat$DO.sat, ylim = c(-120,120), ylab = "",xlab="", type = "l", xaxt="n",yaxt="n", lwd=.8)
  axis(4, las=2, line=0, xpd=NA, tck=-.02, labels=FALSE,
       at=round(seq(0,120, by=20), 1))
  axis(4, las=2, line=-0.5, xpd=NA, tcl=0, col='transparent',
       at=round(seq(0, 120, by=20), 1))
  mtext("DO (%sat)", 4, 2, adj=.8)
  legend("bottomright", c("DO","GPP", "ER"), bty = "n", lty = c(1,1,1),lwd = c(.8,3,3), 
         col = c("black","forestgreen", "sienna"))
  
  par(mar = c(3,4,0,4))
    plot(sdat$date, sdat$discharge, log="y", ylab = "Q (m3s)", xlab = "date", type = "l", xaxt="n") 
    par(new=T)
    plot(s_metab$date, s_metab$GPP, type="n", axes=FALSE, ylab = "")
    t<- seq.Date(as.Date("2019-04-01"), s_metab$date[nrow(s_metab)],by = "month")
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
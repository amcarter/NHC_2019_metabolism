library(tidyverse)
library(lubridate)
library(zoo)
library(pracma)
setwd("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism/data")
source("../src/streamMetabolizer/inspect_model_fits.r")


dat <- readRDS("metabolism/compiled/raymond_met.rds")

dat$preds <- dat$preds %>%
  mutate(year = case_when(!(site %in% c("nhc","unhc")) ~ 2019,
                          date < ymd("2018-03-01") ~ 2017,
                          date < ymd("2019-03-01") ~ 2018,
                          TRUE ~ 2019))
sites = c("nhc", "pm", "cbp","wb","wbp","unhc")


pdf(file = "../figures/raymond_met_with_Q.pdf", onefile = T)
par(mar = c(3,4,0,4), oma = c(4,0,2,0),
    mfrow = c(3,1))
for(s in sites){
  if(s == "nhc"){
    for(y in c(2017, 2018, 2019, 2019)){
      preds <- dat$preds %>%
        filter(site == s,
               year == y)
      plot_hall_metab(preds, ylim = c(-30,8))
      par(new = T)
      plot(preds$date, preds$discharge.daily, log = "y", type = "l",
           ylim = c(0.001, 
                    max(dat$preds$discharge.daily, na.rm = T)*500),
           axes = F, xlab = "", ylab = "")
      axis(4)
      mtext("daily discharge (m3s)", 4, 2.5)   
      mtext(paste(s, y), cex = 1, line = -1.5)
    }
    next
  }
  if(s == "unhc"){
    for(y in c( 2019, 2017, 2018, 2019)){
      preds <- dat$preds %>%
        filter(site == s,
               year == y)
      plot_hall_metab(preds, ylim = c(-30,8))
      par(new = T)
      plot(preds$date, preds$discharge.daily, log = "y", type = "l",
           ylim = c(0.001, 
                    max(dat$preds$discharge.daily, na.rm = T)*500),
           axes = F, xlab = "", ylab = "")
      axis(4)
      mtext("daily discharge (m3s)", 4, 2.5)   
      mtext(paste(s, y), cex = 1, line = -1.5)
    }
    next
  }
  y = 2019
  preds <- dat$preds %>%
    filter(site == s,
           year == y)
  plot_hall_metab(preds, ylim = c(-30,8))
  par(new = T)
  plot(preds$date, preds$discharge.daily, log = "y", type = "l",
       ylim = c(min(dat$preds$discharge.daily, na.rm = T), 
                max(dat$preds$discharge.daily, na.rm = T)*100),
       axes = F, xlab = "", ylab = "")
  axis(4)
  mtext("daily discharge (m3s)", 4, 2.5)   
  mtext(paste(s, y), cex = 1, line = -1.5)
  # if(s %in% c("cbp", "wbp")){
  #   plot(1,type = "n",axes = F, xlab = "", ylab = "")
  # }
}
dev.off()

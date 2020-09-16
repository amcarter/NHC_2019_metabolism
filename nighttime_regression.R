# Calculate K and ER from nighttime regression on sites
# equation from Hall and Ulseth 2020, originally from Hornberger and Kelly 1975

# A Carter
# 2020-07-17

library(dplyr)
library(streamMetabolizer)
library(readr)
library(mgcv)
library(lubridate)

setwd(metab_projdir)

dat <- read_csv("data/metabolism/processed/NHC.csv", guess_max = 100000)
dat <- select(dat, -flagtype, -flagcomment)

dat$datetime <- with_tz(dat$DateTime_UTC, "EST")


##Convert KO2 to K600
K600fromO2<-function(temp, KO2) {
  ((600/(1800.6-(120.1*temp)+(3.7818*temp^2)-(0.047608*temp^3)))^-0.5)*KO2
  
}

nightreg<-function(o2file, date){
  #get Q
  Q <- mean(o2file[o2file$date==date, ]$discharge, na.rm=T)
  #find sunset:
  o2fil <- o2file[o2file$datetime>= ymd_hms(paste(date,"15:00:00")) & o2file$datetime<=ymd_hms(paste(date+1,"12:00:00")),]
  sunset<- o2fil$datetime[which(o2fil$light==0)[1]]
  o2fil <- o2fil[o2fil$datetime>=sunset,]  
  equil <- o2fil$datetime[which.min(o2fil$DO_mgL)[1]]
  o2fil <- o2fil[o2fil$datetime<=equil,]  
  
  
  temp<-o2fil$temp
  oxy<-o2fil$DO_mgL
  oxsat <- o2fil$DO_sat
  
  ##moving average on oxy data
  oxyf1<-stats::filter(oxy,rep(1/3,3), sides=2)
  satf1<-stats::filter(oxsat,rep(1/3,3), sides=2)
  
  #trim the ends of the oxy data
  oxyf2<- oxyf1[c(-1,-length(oxyf1))]
  satf2<- satf1[c(-1,-length(satf1))]
  
  ##calculate delO/delt
  deltaO2<-((oxyf2[-1]-oxyf2[-length(oxyf2)])/15)*1440
  
  #Trim the first two and last one from the temp data to match the filter oxy data
  temptrim<-temp[c(-2:-1,-length(temp))]
  #calc the dodef
  satdef<-satf2[-1]-oxyf2[-1]
  
  #calculate regression and plot
  nreg<-lm(deltaO2~satdef)
  # plot(satdef,deltaO2, main=date)
  # abline(nreg)
  
  coeff<-coef(nreg)
  ci <- confint(nreg, "satdef", 0.95)[2]-confint(nreg, "satdef", 0.95)[1]
  out<-list(coeff, 
            K600=K600fromO2(mean(temp, na.rm=T), coeff[2]), 
            Q_m3s=Q,
            n=length(satdef), 
            ci = ci,
            r2_adj=summary(nreg)$adj.r.squared,
            p=summary(nreg)$coefficients["satdef","Pr(>|t|)"])##gives both the regression coefficients and converts the slope to K600
  out
}

site <- "MC751"
dat <- read_csv(paste0("data/metabolism/processed/",site,".csv"), guess_max = 100000)
dat <- select(dat, -flagtype, -flagcomment)

dat$datetime <- with_tz(dat$DateTime_UTC, "EST")

o2file <- dat %>%
  select(datetime, DO_mgL, DO_sat=DO.sat, temp=WaterTemp_C, light, discharge=Discharge_m3s)

o2file$date <- as.Date(o2file$datetime, tz="EST")
dates <- unique(o2file$date)
nightreg_results<- data.frame()
for(i in 1:length(dates)){
  date <- dates[i]
  mm<- try(nightreg(o2file, date))
  if('try-error' %in% class(mm)) next
  tmp <- data.frame(date = date,
                    K600 = unname(mm$K600),
                    discharge_m3s = mm$Q_m3s,
                    r2_adj = mm$r2_adj,
                    p = mm$p, ci=mm$ci)
  nightreg_results <- bind_rows(nightreg_results, tmp)
}

# edit results to include only values where:
# K600>=0
# r2_adj>=0.3
# p<.1
# remove all points with >10% of confidence intervals
ci.max <- quantile(nightreg_results$ci, 0.9, na.rm=T)

nightreg_results <- nightreg_results %>%
  dplyr::filter(K600>=0,
                ci<ci.max,
                r2_adj>0)
#                p<0.1)

#plot(log(nightreg_results$discharge_m3s), nightreg_results$K600, pch=20, cex=2*nightreg_results$r2_adj)
# look at specific questionable points:
#dd<-as.Date("2020-02-08")
#nightreg(o2file, dd)
#plot(o2file[o2file$date%in%dd:(dd+1), ]$datetime, o2file[o2file$date%in%dd:(dd+1),]$DO_mgL)

 # calculate a spline with base stats
x<- range(log(nightreg_results$discharge_m3s))
lnQ_nodes<- seq(x[1],x[2], length.out=9) 
#lnQ_nodes <- lnQ_bins[-1]-diff(lnQ_bins)/2
Q <- data.frame(logQ=lnQ_nodes)
#calculate a spline with CI using mgcv
nightreg_results$logQ <- log(nightreg_results$discharge_m3s)
ss <- gam(K600 ~ s(logQ), data=nightreg_results)
p<- predict(ss, newdata=Q, se.fit=TRUE)

kk=predict(ss, newdata=data.frame(logQ=nightreg_results$logQ), se.fit=T)
nightreg_results$K600.pred <- kk$fit
nightreg_results$K600.se <- kk$se.fit


write_csv(nightreg_results, paste0("data/gas_data/night_regression/nightreg_", site,".csv"))

#png(width=7, height=6, units='in', filename=paste0("figures/KQ_nightreg/",site,".png"), type='cairo', res=300)
  plot(ss,axes=F, xlab = "log Q", ylab = "K600", main=site)
  with(nightreg_results, points(logQ, K600-mean(K600, na.rm=T), pch=20, col="grey70"))
  points(Q$logQ, p$fit-mean(nightreg_results$K600, na.rm=T), col="brown2", pch=19)
  axis(1)
  par(new=T)
  plot(c(Q$logQ,Q$logQ), c((p$fit+2*p$se.fit),(p$fit-2*p$se.fit)),
     type="n", axes=F, ylab="", xlab="")
  axis(2)
#dev.off()
# make a data frame of results
KQ <- data.frame(site=site,
                 n_nreg = nrow(nightreg_results),
                 lnQ_nodes=lnQ_nodes,
                 K600_lnQ_nodes=p$fit,
                 K600_lnQ_nodes_se=p$se.fit)

KQ_all <- bind_rows(KQ_all, KQ)
rownames(KQ_all) <- NULL

write_csv(KQ_all,"data/siteData/KQ_nightreg_priors.csv")


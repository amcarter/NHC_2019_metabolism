############################################################################
# Simulated data for testing SM assumptions

library(tidyverse)
library(lubridate)
library(zoo)
library(streamMetabolizer)
library(LakeMetabolizer)

setwd(metab_projdir)

dat <- read_csv("data/metabolism/processed/UNHC.csv", guess_max = 100000)
light <- readRDS("data/light/NC_UNHC_predicted.rds")
light$DateTime_UTC <- with_tz(light$local_time, "UTC")

light <- select(light, "DateTime_UTC","LAI", "PAR_bc")
dat <- left_join(dat, light)

startdate <- ymd_hms("2019-05-15 04:00:00")
enddate <- ymd_hms("2019-06-01 04:00:00")
dat$PAR_bc <- na.approx(dat$PAR_bc, na.rm=F)
dat$DO.obs[dat$solar.time==ymd_hms("2019-05-30 20:59:41")]<- NA
dat <- dat[dat$solar.time<=enddate&dat$solar.time>=startdate,]

dat <- select(dat, "solar.time", DO.obs="DO_mgL", "DO.sat", "depth", temp.water="WaterTemp_C", light="PAR_bc", discharge="Discharge_m3s")
dat$temp.water <- na.approx(dat$temp.water)
dat$DO.sat <- na.approx(dat$DO.sat)
dat$DO.obs <- na.approx(dat$DO.obs)

###############
# run stream metabolizer to get K values

bayes_name <- mm_name(type="bayes", pool_K600="binned", err_obs_iid = TRUE,
                      err_proc_iid = TRUE, err_proc_acor = FALSE)
bayes_specs <- specs(bayes_name)

Qrange <- c(min(log(dat$discharge), na.rm=T),quantile(log(dat$discharge), .98, na.rm=T))
bayes_specs$K600_lnQ_nodes_centers <- c(-3.4, -3.2, -3, -2.8, -2.6, -2.4)

mm <- metab(bayes_specs, data=dat)

data_daily <- mm@data_daily %>% select(date, discharge.m3s=discharge.daily)
K600 <- mm@fit$daily  %>% select(date, K600=K600_daily_50pct, K600.lower=K600_daily_2.5pct, K600.upper=K600_daily_97.5pct)
data_daily <- left_join(data_daily, K600, by="date")
metab <- predict_metab(mm)
metab <- left_join(metab, data_daily, by="date")

data <- predict_DO(mm)

mod_specs<- get_specs(mm)
mod_specs$K600 <- exp(get_fit(mm)$KQ_binned$lnK600_lnQ_nodes_50pct)
mod_specs$K600.lower <- exp(get_fit(mm)$KQ_binned$lnK600_lnQ_nodes_2.5pct)
mod_specs$K600.upper <- exp(get_fit(mm)$KQ_binned$lnK600_lnQ_nodes_97.5pct)
mod_specs$K600.sd <- exp(get_fit(mm)$KQ_binned$lnK600_lnQ_nodes_sd)

filter_vec = names(mod_specs) %in% c("K600_lnQ_nodes_centers", "K600", "K600.upper", "K600.lower", "K600.sd")
QvsK600 <- as.data.frame(mod_specs[filter_vec])

plot(QvsK600$K600_lnQ_nodes_centers, QvsK600$K600, type="b")
points(log(data_daily$discharge.m3s), data_daily$K600)

# Add K values
# 
# find_K <- function(x){
#   N<- which.min(abs(log(x)-QvsK600$K600_lnQ_nodes_centers))
#   k <- rnorm(1, QvsK600$K600[N], QvsK600$K600.sd[N])
#   return(k)
# }
# 
# for(i in 1:nrow(dat)){
#   dat$K600[i] <- find_K(dat$discharge[i])
# }  


# get values for process and observation error:

err <- get_fit(mm)$inst
sigma_proc <- mean(err$err_proc_iid_sd, na.rm=T)
sigma_obs <- mean(err$err_obs_iid_sd, na.rm=T)




#####################################################################
# data for simulating DO values:
dat_daily <- metab%>% select(date,  K600.daily=K600, GPP.daily=GPP, ER.daily=ER, discharge.daily=discharge.m3s)
dat_daily$err.obs.sigma= .1
dat_daily$err.obs.phi = 0
dat_daily$err.proc.sigma = 1
dat_daily$err.proc.phi = 0

sdat <- select(dat, -discharge)

# define simulation parameters
mm <- metab_sim(
  specs(mm_name('sim'), err_obs_sigma=sigma_obs, err_proc_sigma=sigma_proc,
        GPP_daily=NULL, ER_daily=NULL, K600_daily=NULL),
  data=dat, data_daily=dat_daily)

get_params(mm)

predict_metab(mm)
DO_sim <- predict_DO(mm)

# Low K 
dat_daily$K600.daily<- 0.5
mm_low_K<-metab_sim(
  specs(mm_name('sim'), err_obs_sigma=sigma_obs, err_proc_sigma=sigma_proc,
        GPP_daily=NULL, ER_daily=NULL, K600_daily=NULL),
  data=dat, data_daily=dat_daily)


get_params(mm_low_K)

predict_metab(mm_low_K)
DO_sim <- predict_DO(mm_low_K)

#################################
# metab_sim wrapper so that the last DO value 
# can be used as the first DO value

sim_DO_multiple_days<- function(dat, dat_daily, DO_initial){
  # mm<- metab_sim(
  #   specs(mm_name('sim'), err_obs_sigma=dat_daily$err.obs.sigma[1], err_proc_sigma=dat_daily$err.proc.sigma[1],
  #         GPP_daily=NULL, ER_daily=NULL, K600_daily=NULL),
  #   data=dat, data_daily=dat_daily)
  dat$date <- as.Date(dat$solar.time)
  dates <- unique(dat$date)
  params <- data.frame()
  met_preds <- data.frame()
  DO_preds <- data.frame()
  for(i in 1:(length(dates)-1)){
    startdate <- ymd_hms(paste(dates[i], "04:00:00"))
    enddate <- ymd_hms(paste(dates[i+1], "04:00:00"))
    tmp <- dat[dat$solar.time <=enddate &dat$solar.time>=startdate,]
    tmp <- tmp %>% select(-date)
    tmp$DO.obs[1] <- DO_initial
    mm <- metab_sim(
      specs(mm_name('sim'), GPP_daily=NULL, ER_daily=NULL, K600_daily=NULL),
      data=tmp, data_daily=dat_daily[i,])
    params <- bind_rows(params, get_params(mm))
    met_preds <- bind_rows(met_preds, predict_metab(mm))
    DO_preds <- bind_rows(DO_preds, predict_DO(mm))
    DO_initial <- DO_preds$DO.obs[nrow(DO_preds)]
  }
  DO_simulated <- list(params=params, met=met_preds, DO_preds=DO_preds)
  return(DO_simulated)
}

DO_simulated <- sim_DO_multiple_days(sdat, dat_daily, sdat$DO.obs[1])
dat_low_K <- dat
dat_low_K$DO.obs <- DO_simulated$DO_preds$DO.obs

mm <- metab(bayes_specs, data=dat_low_K)
plot(predict_metab(mm)$date, predict_metab(mm)$GPP, type="b")

######################
## Model Metabolism ##
######################
# install.packages("streamMetabolizer", dependencies=TRUE, 
#                  repos=c("http://owi.usgs.gov/R","https://cran.rstudio.com"))
# update.packages(oldPkgs=c("streamMetabolizer","unitted"), dependencies=TRUE, 
#                 repos=c("http://owi.usgs.gov/R", "https://cran.rstudio.com"))
# devtools::install_github("USGS-R/streamMetabolizer", ref="develop")
devtools::find_rtools()
Sys.getenv('PATH')

library(rstan)
library(tidyselect)
library(tidyr)
library(dplyr)
library(ggplot2)
library(streamMetabolizer)
library(lubridate)
library(dygraphs)
library(imputeTS)
library(parallel)

setwd("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism/data")

sites <- read_csv("siteData/NHCsite_metadata.csv")

## Read in Data ####
# select variables for metabolism
read_metdata <- function(site){
  MP <- read_csv(paste0("metabolism/processed/", site, ".csv"), guess_max = 100000) %>%
    select(solar.time, DO.obs, DO.sat, depth, temp.water, light, discharge)
  return(MP)
}

NHC <- read_metdata("NHC")
PM <- read_metdata("PM")
CBP <- read_metdata("CBP")
WB <- read_metdata("WB")
WBP <- read_metdata("WBP")
PWC <- read_metdata("PWC")
UNHC <- read_metdata("UNHC")


#Visualize the data #####
# dat <- UNHC
# 
# dat %>% unitted::v() %>%
#   mutate(DO.pctsat = 100 * (DO.obs / DO.sat)) %>%
#   select(solar.time, starts_with('DO')) %>%
#   gather(type, DO.value, starts_with('DO')) %>%
#   mutate(units=ifelse(type == 'DO.pctsat', 'DO\n(% sat)', 'DO\n(mg/L)')) %>%
#   ggplot(aes(x=solar.time, y=DO.value, color=type)) + geom_line() + 
#   facet_grid(units ~ ., scale='free_y') + theme_bw() +
#   scale_color_discrete('variable')
# 
# labels <- c(depth='depth\n(m)', temp.water='water temp\n(deg C)', 
#             light='PAR\n(umol m^-2 s^-1)', discharge='Q\n(cms)')
# dat %>% unitted::v() %>%
#   select(solar.time, depth, temp.water, light, discharge) %>%
#   gather(type, value, depth, temp.water, light, discharge) %>%
#   mutate(
#     type=ordered(type, levels=c('depth','temp.water','light','discharge')),
#     units=ordered(labels[type], unname(labels))) %>%
#   ggplot(aes(x=solar.time, y=value, color=type)) + geom_line() + 
#   facet_grid(units ~ ., scale='free_y') + theme_bw() +
#   scale_color_discrete('variable')


## Model the data #####
## Figure out range of log of daily discharge and reset parameters accordingly
## Recalc potential mean based on Raymond relationship
#K600_daily_meanlog_meanlog = 4.77 + 0.55*log(slope) + (-0.52*log(depth))
#K600_daily_meanlog_sdlog = 0.70

## Set bayes specs
bayes_name_new <- mm_name(type='bayes', pool_K600="binned", 
                          err_obs_iid=TRUE, err_proc_iid = TRUE, 
                          ode_method = "trapezoid", deficit_src='DO_mod', 
                          engine='stan')

bayes_name_linear <- mm_name(type='bayes', pool_K600="linear", 
                             err_obs_iid=TRUE, err_proc_iid = TRUE, 
                             ode_method = "trapezoid", deficit_src='DO_mod', 
                             engine='stan')


set_up_model <- function(dat, site, bayes_name){
    
  ## Set bayes specs
  bayes_specs <- specs(bayes_name)
  ## Based on range of log daily Q
  Qrange <- c(quantile(log(dat$discharge),.02, na.rm=T),
              quantile(log(dat$discharge), .98, na.rm=T))
  bayes_specs$K600_lnQ_nodes_centers <- seq(Qrange[1], Qrange[2], length=7)
  
  ## Based on Pete Raymond's data
  slope <- sites[sites$sitecode==site, ]$slope
  meanlog <- 4.77+0.55*log(slope)+(-0.52*(log(median(dat$depth, na.rm = T))))
  bayes_specs$K600_lnQ_nodes_meanlog <- c(rep(meanlog, 7))
  bayes_specs$K600_lnQ_nodes_sdlog <- c(rep(0.7, 7))
  ## Change sigma
  bayes_specs$K600_daily_sigma_sigma <- 0.05
  return(bayes_specs)
}


####################
## Run Models
###################
bayes_specs_NHC <- set_up_model(NHC, "NHC", bayes_name_new)
bayes_specs_PM <- set_up_model(PM, "PM", bayes_name_new)
bayes_specs_CBP <- set_up_model(CBP, "CBP", bayes_name_new)
bayes_specs_WB <- set_up_model(WB, "WB", bayes_name_new)
bayes_specs_WBP <- set_up_model(WBP, "WBP", bayes_name_new)
bayes_specs_PWC <- set_up_model(PWC, "PWC", bayes_name_new)
bayes_specs_UNHC <- set_up_model(UNHC, "UNHC", bayes_name_new)


fit_NHC_raymond <- metab(bayes_specs_NHC, data=NHC)
fit_PM_raymond <- metab(bayes_specs_PM, data=PM)
fit_CBP_raymond <- metab(bayes_specs_CBP, data=CBP)
fit_WB_raymond <- metab(bayes_specs_WB, data=WB)
fit_WBP_raymond <- metab(bayes_specs_WBP, data=WBP)
fit_PWC_raymond <- metab(bayes_specs_PWC, data=PWC)
fit_UNHC_raymond <- metab(bayes_specs_UNHC, data=UNHC)

####################
## Run in parallel -- (AJ you can ignore this was for comparing two levels of data cleaning)
####################
# list26 <- list(W26_con, W26_lim)
# ##apply
# cl <- makeCluster(mc <- getOption("cl.cores", 45))
# ## Export data to the cluster
# clusterExport(cl=cl, varlist=c("list26","bayes_specs_W26"))
# ## Run streamMetabolizer
# SM26_new <- parLapply(cl, list26, function(x) streamMetabolizer::metab(bayes_specs_W26, data=x[,-1]))
# ## stop
# stopCluster(cl)
# 
# fit_26_con <- SM26_new[[1]]
# fit_26_lim <- SM26_new[[2]]
# 



###############################################
## Save info & Check Binning
################################################
## Check binning
Binning <- function(Site){
  fit_Site <- get_fit(Site)
  
  SM_output <- fit_Site$daily
  SM_day <- get_data_daily(Site)
  SM_KQbin <- fit_Site$KQ_binned
  SM_specs <- get_specs(Site)
  
  day <- data.frame(SM_day$discharge.daily, SM_output$K600_daily_50pct, 
                    rep('daily', dim(SM_output)[1]))
  colnames(day)<-c('Q', 'K600', 'Group')
  
  nodes<-data.frame(exp(as.numeric(as.character(SM_specs$K600_lnQ_nodes_centers))), 
                    exp(SM_KQbin$lnK600_lnQ_nodes_50pct), 
                    rep('node', dim(SM_KQbin)[1]))
  colnames(nodes)<-c('Q', 'K600', 'Group')
  KQ_plot<-rbind(day,nodes)
  
  ggplot(data=KQ_plot, aes(x=log(Q), y=K600, group=Group, colour=Group)) + 
    geom_point(size=3) +
    #geom_line() + 
    scale_color_manual(name="K-Q",
                       breaks = c("daily", "node"),
                       values=c("grey", "purple"),
                       labels=c("Daily","Bin")) +
    ylab("K600") +
    xlab("logQ") +
    theme_bw() +
    theme(legend.position = "top")
}

Binning(fit_NHC_raymond)

## Visualize
plot_metab_preds(predict_metab(fit_13_lim))

## K600 vs ER
KvER <- get_fit(fit_79_lim)
plot(KvER$daily$K600_daily_mean, KvER$daily$ER_daily_mean)
cor.test(KvER$daily$K600_daily_mean, KvER$daily$ER_daily_mean)

## Write Files
writefiles <- function(mod){
  data <- get_fit(mod)
  for (i in seq_along(data)) {
    filename = paste(names(data)[i], ".csv")
    write.csv(data[[i]], filename)
  }
  write.csv(unlist(get_specs(mod)),"specs.csv")
  write.csv(get_data_daily(mod), "datadaily.csv")
  write.csv(get_data(mod),"mod_and_obs_DO.csv")
}

## Create new folder for site and write csv info
## Reset working directory!!
getwd()
writefiles(fit_207_lim)

## Save models
## Reset working directory!
#saveRDS(fit_22_lim, file = "W22_limclean_pp.rds")

############################################################################
## Pooling by hand for 13, 22, 130, and 207 (site with strong K600 v ER)
############################################################################

## First, pick out high GPP dates for each site
high_GPP_data <- function(fit, dat, thresh){
daily <- get_fit(fit)$daily
high_dates <- daily$date[which(daily$GPP_mean > thresh)]

dat$Date <- as.Date(dat$solar.time)
dat$Hour <- hour(dat$solar.time)
dat$Date_Mod <- ifelse(test= dat$Hour < 4, yes = ymd(dat$Date) - 1, no= ymd(dat$Date))
dat$Date_Mod <- as.Date(dat$Date_Mod, origin = "1970-01-01")

dat_high <- dat[dat$Date_Mod %in% high_dates,]
dat_high <- dat_high[,1:8]
return(dat_high)
}

W13_dat_high <- high_GPP_data(fit_13_lim, W13_lim, 0.5) ## AJ - the thresholds were chosen subjectively by examining fits
W22_dat_high <- high_GPP_data(fit_22_lim, W22_lim, 1)
W130_dat_high <- high_GPP_data(fit_130_lim, W130_lim, 0.25)
W207_dat_high <- high_GPP_data(fit_207_lim, W207_lim, 0.25)

## Run high days only model with partial pooling
fit_13_high <- metab(bayes_specs_W13, data=W13_dat_high[,-1])
fit_22_high <- metab(bayes_specs_W22, data=W22_dat_high[,-1])
fit_130_high <- metab(bayes_specs_W130, data=W130_dat_high[,-1])
fit_207_high <- metab(bayes_specs_W207, data=W207_dat_high[,-1])

## Save Models
saveRDS(fit_13_high, file = "W13_high_pp.rds")
saveRDS(fit_22_high, file = "W22_high_pp.rds")
saveRDS(fit_130_high, file = "W130_high_pp.rds")
saveRDS(fit_207_high, file = "W207_high_pp.rds")

## Check binning
Binning(fit_13_high)

## Visualize
plot_metab_preds(predict_metab(fit_207_high))

## K600 vs ER
KER <- function(x) {
  KvER <- get_fit(x)
  plot(KvER$daily$K600_daily_mean, KvER$daily$ER_daily_mean)
  cor.test(KvER$daily$K600_daily_mean, KvER$daily$ER_daily_mean)
}

KER(fit_13_high) #fine
KER(fit_22_high) #fine
KER(fit_130_high) #still bad
KER(fit_207_high) #good

## Write csv files for high models
writefiles(fit_207_high)



####################################################################################
## Extract node centers from high models and rerun partial pooling but with those
####################################################################################

## Rerun 13
bayes_specs_W13_v2 <- bayes_specs_W13
bayes_specs_W13_v2$K600_lnQ_nodes_meanlog <- get_fit(fit_13_high)$KQ_binned$lnK600_lnQ_nodes_mean
bayes_specs_W13_v2$K600_lnQ_nodes_sdlog <- get_fit(fit_13_high)$KQ_binned$lnK600_lnQ_nodes_sd
bayes_specs_W13_v2

## Rerun 22
bayes_specs_W22_v2 <- bayes_specs_W22
bayes_specs_W22_v2$K600_lnQ_nodes_meanlog <- get_fit(fit_22_high)$KQ_binned$lnK600_lnQ_nodes_mean
bayes_specs_W22_v2$K600_lnQ_nodes_sdlog <- get_fit(fit_22_high)$KQ_binned$lnK600_lnQ_nodes_sd
bayes_specs_W22_v2

## Rerun 130
bayes_specs_W130_v2 <- bayes_specs_W130
bayes_specs_W130_v2$K600_lnQ_nodes_meanlog <- get_fit(fit_130_high)$KQ_binned$lnK600_lnQ_nodes_mean
bayes_specs_W130_v2$K600_lnQ_nodes_sdlog <- get_fit(fit_130_high)$KQ_binned$lnK600_lnQ_nodes_sd
bayes_specs_W130_v2

## Rerun 207
bayes_specs_W207_v2 <- bayes_specs_W207
bayes_specs_W207_v2$K600_lnQ_nodes_meanlog <- get_fit(fit_207_high)$KQ_binned$lnK600_lnQ_nodes_mean
bayes_specs_W207_v2$K600_lnQ_nodes_sdlog <- get_fit(fit_207_high)$KQ_binned$lnK600_lnQ_nodes_sd
bayes_specs_W207_v2


## Run the Models
fit_13_lim_v2 <- metab(bayes_specs_W13_v2, data=W13_lim[,-1])
fit_22_lim_v2 <- metab(bayes_specs_W22_v2, data=W22_lim[,-1])
fit_130_lim_v2 <- metab(bayes_specs_W130_v2, data=W130_lim[,-1])
fit_207_lim_v2 <- metab(bayes_specs_W207_v2, data=W207_lim[,-1])

## Check binning
Binning(fit_207_lim_v2)

## K600 vs ER
KER(fit_13_lim_v2) #good
KER(fit_22_lim_v2) #good
KER(fit_130_lim_v2) #still bad
KER(fit_207_lim_v2) #fine

## Visualize
plot_metab_preds(predict_metab(fit_130_lim_v2))

## Write csv files for high models
writefiles(fit_13_lim_v2)


## Save Models
saveRDS(fit_13_lim_v2, file = "W13_lim_highpriors.rds")
saveRDS(fit_22_lim_v2, file = "W22_lim_highpriors.rds")
saveRDS(fit_130_lim_v2, file = "W130_lim_highpriors.rds")
saveRDS(fit_207_lim_v2, file = "W207_lim_highpriors.rds")















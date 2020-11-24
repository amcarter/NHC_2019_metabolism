## Model Metabolism #
# Adapted from JRB's workflow
# install.packages("streamMetabolizer", dependencies=TRUE, 
#                  repos=c("http://owi.usgs.gov/R","https://cran.rstudio.com"))
# update.packages(oldPkgs=c("streamMetabolizer","unitted"), dependencies=TRUE, 
#                 repos=c("http://owi.usgs.gov/R", "https://cran.rstudio.com"))
# devtools::install_github("USGS-R/streamMetabolizer", ref="develop")

library(rstan)
library(tidyverse)
library(ggplot2)
library(streamMetabolizer)
library(lubridate)
library(dygraphs)
# library(parallel)

setwd("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism/data")

## Read in Data ####
sites <- read_csv("siteData/NHCsite_metadata.csv")

# select variables for metabolism
read_metdata <- function(site){
  MP <- read_csv(paste0("metabolism/processed/", site, ".csv"), guess_max = 100000) %>%
    select(solar.time, DO.obs, DO.sat, depth, temp.water, light, discharge)
  return(MP)
}

NHC <- read_metdata("NHC")
# nhc19 <- NHC %>%
#   mutate(year = year(solar.time)) %>%
#   filter(year == 2019) %>%
#   select(-year)

PM <- read_metdata("PM")
CBP <- read_metdata("CBP_lvl")
WB <- read_metdata("WB")
WBP <- read_metdata("WBP")
PWC <- read_metdata("PWC")
UNHC <- read_metdata("UNHC")

# Visualize the data #####
# dat <- nhc19
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


## Set up Model Specs #####
## Figure out range of log of daily discharge and reset parameters accordingly
## Recalc potential mean based on Raymond relationship
#K600_daily_meanlog_meanlog = 4.77 + 0.55*log(slope) + (-0.52*log(depth))
#K600_daily_meanlog_sdlog = 0.70

## Set bayes specs
bayes_name <- mm_name(type='bayes', pool_K600="binned", 
                          err_obs_iid=TRUE, err_proc_iid = TRUE, 
                          ode_method = "trapezoid", deficit_src='DO_mod', 
                          engine='stan')

bayes_name_complete <- mm_name(type='bayes', pool_K600="normal_sdzero",
                             err_obs_iid=TRUE, err_proc_iid = TRUE,
                             ode_method = "trapezoid", deficit_src='DO_mod',
                             engine='stan')


kq_all <- read_csv("siteData/KQ_nightreg_priors.csv")
kq_hall <- read_csv("siteData/KQ_hall_prior.csv")

set_up_model <- function(dat, site, bayes_name,
                         version = "raymond", kq_all = NULL
                         ){
  ## Set bayes specs
  bayes_specs <- specs(bayes_name)
  bayes_specs$keep_mcmcs <- FALSE
  ## Based on range of log daily Q
  daily <- dat %>%
    mutate(date = as.Date(solar.time)) %>% group_by(date) %>%
    summarize(discharge = mean(discharge, na.rm = T))
  Qrange <- c(quantile(log(daily$discharge),.02, na.rm=T),
              quantile(log(daily$discharge), .98, na.rm=T))
  
  ## Based on Pete Raymond's data
  slope <- sites[sites$sitecode==site, ]$slope
  meanlog <- 4.77+0.55*log(slope)+(-0.52*(log(median(dat$depth, na.rm = T))))
  
  if(version == "raymond"){
    bayes_specs$K600_lnQ_nodes_centers <- seq(Qrange[1], Qrange[2], length=7)
    bayes_specs$K600_lnQ_nodes_meanlog <- c(rep(meanlog, 7))
    bayes_specs$K600_lnQ_nodes_sdlog <- c(rep(0.7, 7))
  }
  
  if(version == "nreg"){
    kq <- kq_all[kq_all$site == site,]
    nodes <- kq$lnQ_nodes   
    step <- diff(nodes)[1]
    low <- floor((nodes[1]-Qrange[1])/step)
    high <- floor((Qrange[2]-nodes[7])/step)
    if(high > 0){
      nodes <- c(nodes,
                 seq(from = (nodes[7] + step), 
                     by = step , 
                     length.out = high))
      for(i in 1:high){
        kq <- bind_rows(kq, kq[7,])}
    }
    if(low > 0){
      nodes <- c(seq(from = (nodes[1] - low * step), 
                     by = step , 
                     length.out = low),
                 nodes)
      for(i in 1:low){
        kq <- bind_rows(kq[1], kq)}
    }
  
    bayes_specs$K600_lnQ_nodes_centers <- nodes
    bayes_specs$K600_lnQ_nodes_meanlog <- (kq$lnK600_lnQ_nodes)
    #bayes_specs$K600_lnQ_nodes_sdlog <- kq$lnK600_lnQ_nodes_sd
    bayes_specs$K600_lnQ_nodes_sdlog <- rep(0.01, length(nodes))
    
  }
  if(version == "hall"){
    nodes <- kq_all %>%
      filter(!is.na(nodes) & nodes >= Qrange[1] & nodes <= Qrange[2])
    
    bayes_specs$K600_lnQ_nodes_centers <- nodes$nodes
    bayes_specs$K600_lnQ_nodes_meanlog <- log(nodes$K600)
    bayes_specs$K600_lnQ_nodes_sdlog <- c(rep(0.7, nrow(nodes)))
  }
  ## Change sigma
  bayes_specs$K600_daily_sigma_sigma <- 0.05
  
  return(bayes_specs)
}


# Run models with Kpriors ####
# First, try several different priors for K

# pooled raymond k600, default K600_sigma_sigma = 0.24, K600_sd = 0.7
specs_raymond <- set_up_model(nhc19, "NHC",bayes_name, "raymond")
fit_raymond_ss24 <- metab(specs_raymond, data = nhc19)
saveRDS(fit_raymond_ss24,"metabolism/modeled/fit_nhc19_raymond_ss24.rds")

# pooled raymond k600, small K600_sigma_sigma = 0.05, K600_sd = 0.7
specs_raymond$K600_daily_sigma_sigma <- 0.05
fit_raymond_ss05 <- metab(specs_raymond, data = nhc19)
saveRDS(fit_raymond_ss05,"metabolism/modeled/fit_nhc19_raymond_ss05.rds")

# pooled nightregression, default K600_sigma_sigma = 0.24
bayes_specs_nreg <- set_up_model(nhc19, "NHC", bayes_name, "nreg", kq_all)
fit_nreg_ss24 <- metab(bayes_specs_nreg, data = nhc19)
saveRDS(fit_nreg_ss24,"metabolism/modeled/fit_nhc19_nreg_ss24.rds")

# pooled nightregression, default K600_sigma_sigma = 0.05
bayes_specs_nreg$K600_daily_sigma_sigma <- 0.05
bayes_specs_nreg <- set_up_model(nhc19, "NHC", bayes_name, "nreg", kq_all)
fit_nreg_ss05 <- metab(bayes_specs_nreg, data = nhc19)
saveRDS(fit_nreg_ss05,"metabolism/modeled/fit_nhc19_nreg_ss05.rds")

# pooled nightregression, default K600_sigma_sigma = 0.05, K600_sd = 0.7
bayes_specs_nreg$K600_lnQ_nodes_sdlog <- rep(0.7, 8)
fit_nreg_sd7_ss05 <- metab(bayes_specs_nreg, data = nhc19)
saveRDS(fit_nreg_sd7_ss05,"metabolism/modeled/fit_nhc19_nreg_sd7_ss05.rds")


# Fixed K model ####
bayes_specs_fixedK <- set_up_model(nhc19, "NHC", bayes_name, "nreg", kq_all)
prep_fake_Q <- function(dat, bayes_specs_fixedK){
  dat <- dat %>% 
    mutate(date = as.Date(solar.time)) 
  daily <- dat %>%
    group_by(date) %>%
    summarize(discharge = mean(discharge, na.rm = T)) %>%
    mutate(logQ = log(discharge))
  
  nbins <-  length(bayes_specs_fixedK$K600_lnQ_nodes_centers)
  
  Qbreaks <- c(bayes_specs_fixedK$K600_lnQ_nodes_centers[1] -
               diff(bayes_specs_fixedK$K600_lnQ_nodes_centers)[1]/2,
               bayes_specs_fixedK$K600_lnQ_nodes_centers +
               diff(bayes_specs_fixedK$K600_lnQ_nodes_centers)[1]/2)
  rQ <- range(daily$logQ, na.rm = T)          
  if(Qbreaks[1] > rQ[1]){Qbreaks[1] <- rQ[1]}
  if(Qbreaks[nbins] < rQ[2]){Qbreaks[nbins] <- rQ[2]}
  
  daily <- daily %>%
    mutate(Q = cut(daily$logQ, Qbreaks, labels = 1:nbins)) %>%
    select(date, Q)
  daily$Q <- as.numeric(daily$Q)
  dat <- left_join(dat, daily, by = "date") %>%
    select(-discharge) %>%
    mutate(discharge = exp(Q)) %>%
    select(-Q,-date)
  return(dat)
}
dat <- prep_fake_Q(nhc19, bayes_specs_fixedK)
bayes_specs_fixedK$K600_lnQ_nodes_centers <- seq(1:nbins)
bayes_specs_fixedK$K600_lnQ_nodes_sdlog <- rep(.01,nbins)
fit_nreg_fixedK <- metab(bayes_specs_fixedK, data= dat) 
saveRDS(fit_nreg_fixedK,"metabolism/modeled/fit_nhc19_fixedK.rds")
rm(fit_nreg_fixedK)

bayes_specs_fixedK$K600_lnQ_nodediffs_sdlog <- 2
fit_nreg_fixedK <- metab(bayes_specs_fixedK, data= dat) 
saveRDS(fit_nreg_fixedK,"metabolism/modeled/fit_nhc19_fixedK2.rds")

# Fixed K model with Hall K600 ####
# CBP 
bayes_specs_Hall <- set_up_model(CBP, "CBP", bayes_name, "hall", kq_hall)
dat <- prep_fake_Q(CBP, bayes_specs_Hall)
bayes_specs_Hall$K600_lnQ_nodes_centers <- 
  seq(1:length(bayes_specs_Hall$K600_lnQ_nodes_centers))
fit <- metab(bayes_specs_Hall, dat)
saveRDS(fit, "metabolism/modeled/fit_cbp_fixed_hallK.rds")

# PM 
bayes_specs_Hall <- set_up_model(PM, "PM", bayes_name, "hall", kq_hall)
dat <- prep_fake_Q(PM, bayes_specs_Hall)
bayes_specs_Hall$K600_lnQ_nodes_centers <- 
  seq(1:length(bayes_specs_Hall$K600_lnQ_nodes_centers))
fit <- metab(bayes_specs_Hall, dat)
saveRDS(fit, "metabolism/modeled/fit_pm_fixed_hallK.rds")

# WB 
bayes_specs_Hall <- set_up_model(WB, "WB", bayes_name, "hall", kq_hall)
dat <- prep_fake_Q(WB, bayes_specs_Hall)
bayes_specs_Hall$K600_lnQ_nodes_centers <- 
  seq(1:length(bayes_specs_Hall$K600_lnQ_nodes_centers))
fit <- metab(bayes_specs_Hall, dat)
saveRDS(fit, "metabolism/modeled/fit_wb_fixed_hallK.rds")

# WBP 
bayes_specs_Hall <- set_up_model(WBP, "WBP", bayes_name, "hall", kq_hall)
dat <- prep_fake_Q(WBP, bayes_specs_Hall)
bayes_specs_Hall$K600_lnQ_nodes_centers <- 
  seq(1:length(bayes_specs_Hall$K600_lnQ_nodes_centers))
fit <- metab(bayes_specs_Hall, dat)
saveRDS(fit, "metabolism/modeled/fit_wbp_fixed_hallK.rds")

# NHC 
for(i in seq(2017:2019)){
  dat <- NHC %>% filter(year(solar.time) == i)
  bayes_specs_Hall <- set_up_model(dat, "NHC", bayes_name, "hall", kq_hall)
  dat <- prep_fake_Q(dat, bayes_specs_Hall)
  bayes_specs_Hall$K600_lnQ_nodes_centers <- 
    seq(1:length(bayes_specs_Hall$K600_lnQ_nodes_centers))
  fit <- metab(bayes_specs_Hall, dat)
  saveRDS(fit, paste0("metaboli  sm/modeled/fit_nhc",i,"_fixed_hallK.rds"))
}

# NHC 
for(i in seq(2017:2019)){
  dat <- UNHC %>% filter(year(solar.time) == i)
  bayes_specs_Hall <- set_up_model(dat, "UNHC", bayes_name, "hall", kq_hall)
  dat <- prep_fake_Q(dat, bayes_specs_Hall)
  bayes_specs_Hall$K600_lnQ_nodes_centers <- 
    seq(1:length(bayes_specs_Hall$K600_lnQ_nodes_centers))
  fit <- metab(bayes_specs_Hall, dat)
  saveRDS(fit, paste0("metabolism/modeled/fit_unhc",i,"_fixed_hallK.rds"))
}



#inspect fits ####
load_inspect_fit <- function(filename){
  fit <- readRDS(paste0("metabolism/modeled/", filename))
  mm<-plot_metab_preds(predict_metab(fit))
  print(mm)
  rr<-get_fit(fit)$overall %>%
    select(ends_with('Rhat'))
  print(rr)
  rh<-get_fit(fit)$daily %>%
    select(date, ends_with('Rhat')) %>%
    ggplot(aes(x = date, y = K600_daily_Rhat)) +
    geom_line() +
    geom_hline(yintercept = 1.05, lty = 2, col = "red")
  print(rh)
  fit@mcmc <- NULL
  return(fit)
}

nhc_fixed <- load_inspect_fit("fit_nhc19_fixedK.rds")
nhc_fixed2 <- load_inspect_fit("fit_nhc19_fixedK2.rds")
fit_raymond_ss05 <- load_inspect_fit("fit_nhc19_raymond_ss05.rds")
fit_nreg_ss24 <- load_inspect_fit("fit_nhc19_nreg_ss24.rds")
fit_nreg_fixedK1 <- load_inspect_fit("fit_nhc19_fixedK.rds")
fit_nreg_sd7_ss05 <- load_inspect_fit("fit_nhc19_nreg_sd7_ss05.rds")

# Based on results from testing this on NHC2019, I will move forward using
# night regression, not normalized to the raymond K600 estimates with 
# a K600_logsd of 0.7 and a K600_sigma_sigma of 0.05

# Run models for all sites #### 

# NHC (2016 - 2020)
for(i in 2016:2020){
  dat <- NHC %>%
    mutate(year = year(solar.time)) %>%
    filter(year == i) %>%
    select(-year)
  bayes_specs <- set_up_model(dat, "NHC", bayes_name, "nreg", kq_all)
  fit <- metab(bayes_specs, data=dat)
  saveRDS(fit, paste0("metabolism/modeled/nhc_", i, "_nreg_v1.rds"))
}

# UNHC (2016 - 2020)
for(i in 2016:2020){
  dat <- UNHC %>%
    mutate(year = year(solar.time)) %>%
    filter(year == i) %>%
    select(-year)
  bayes_specs <- set_up_model(dat, "UNHC", bayes_name, "nreg", kq_all)
  fit <- metab(bayes_specs, data=dat)
  saveRDS(fit, paste0("metabolism/modeled/unhc_", i, "_nreg_v1.rds"))
}

# PM
bayes_specs_PM <- set_up_model(PM, "PM", bayes_name, "nreg", kq_all)
fit_PM <- metab(bayes_specs_PM, data=PM)
saveRDS(fit_PM, "metabolism/modeled/pm_nreg_v1.rds")
rm(fit_PM)

# CBP
bayes_specs_CBP <- set_up_model(CBP, "CBP", bayes_name, "nreg", kq_all)
fit_CBP <- metab(bayes_specs_CBP, data=CBP)
saveRDS(fit_CBP, "metabolism/modeled/cbp_nreg_v1.rds")
rm(fit_CBP)

# WB
bayes_specs_WB <- set_up_model(WB, "WB", bayes_name, "nreg", kq_all)
fit_WB <- metab(bayes_specs_WB, data=WB)
saveRDS(fit_WB, "metabolism/modeled/wb_nreg_v1.rds")
rm(fit_WB)

# WBP
bayes_specs_WBP <- set_up_model(WBP, "WBP", bayes_name, "nreg", kq_all)
fit_WBP <- metab(bayes_specs_WBP, data=WBP)
saveRDS(fit_WBP, "metabolism/modeled/wbp_nreg_v1.rds")
rm(fit_WBP)

# PWC
bayes_specs_PWC <- set_up_model(PWC, "PWC", bayes_name, "nreg", kq_all)
fit_PWC <- metab(bayes_specs_PWC, data=PWC)
saveRDS(fit_PWC, "metabolism/modeled/pwc_nreg_v1.rds")
rm(fit_PWC)

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
Binning <- function(Site, thresh = 0.5){
  fit_Site <- get_fit(Site)
  
  SM_output <- fit_Site$daily
  SM_day <- get_data_daily(Site)
  SM_KQbin <-  fit_Site$KQ_binned
  SM_specs <- get_specs(Site)
  
  day <- data.frame(SM_day$discharge.daily, 
                    SM_output$K600_daily_50pct, 
                    SM_output$GPP_50pct,
                    SM_output$K600_daily_Rhat,
                    rep('daily', dim(SM_output)[1]))
  colnames(day)<-c('Q', 'K600', 'GPP','Rhat', 'Group')
  
  gg<-ggplot(day, aes(x=log(Q), y = GPP, col=Rhat))+
    geom_point() +
    geom_hline(yintercept = thresh)
  print(gg)
  nodes<-data.frame(exp(SM_specs$K600_lnQ_nodes_centers), 
                     exp(SM_KQbin$lnK600_lnQ_nodes_50pct), 
                    rep('node', dim(SM_KQbin)[1]))
  colnames(nodes)<-c('Q', 'K600', 'Group')
  nodes$K600_prior <- exp(SM_specs$K600_lnQ_nodes_meanlog)
  KQ_plot<-bind_rows(day,nodes)
  
  ggplot(data=KQ_plot, aes(x=log(Q), y=K600, group=Group, colour=Group)) + 
    geom_point(size=3) +
    #geom_line() + 
    scale_color_manual(name="K-Q",
                       breaks = c("daily", "node"),
                       values=c("grey", "purple"),
                       labels=c("Daily","Bin")) +
    geom_point(aes(y = K600_prior),
               size = 3, col = "purple", pch = 21) +
    ylab("K600") +
    xlab("logQ") +
    theme_bw() +
    theme(legend.position = "top")
}

Binning(fit_raymond_ss05)
Binning(fit_nreg_ss24)
Binning(fit_nreg_fixedK1)
Binning(fit_nreg_sd7_ss05)

## Visualize
plot_metab_preds(predict_metab(fit_nreg_sd7_ss05))
plot_metab_preds(predict_metab(fit_nreg_fixedK))

## K600 vs ER
KvER <- get_fit(fit_raymond_ss24)
KvER <- get_fit(fit_raymond_ss05)
KvER <- get_fit(fit_nreg_ss24)
KvER <- get_fit(fit_nreg_sd7_ss05)
KvER <- get_fit(fit_nreg_fixedK)
plot(KvER$daily$K600_daily_mean, KvER$daily$ER_daily_mean)
cor(KvER$daily$K600_daily_mean, 
    KvER$daily$ER_daily_mean, 
    use = "na.or.complete")

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

## AJ - the thresholds were chosen subjectively by examining fits
nhc_dat_high <- high_GPP_data(fit_raymond_ss24, nhc19, 0.5) 
# W22_dat_high <- high_GPP_data(fit_22_lim, W22_lim, 1)
# W130_dat_high <- high_GPP_data(fit_130_lim, W130_lim, 0.25)
# W207_dat_high <- high_GPP_data(fit_207_lim, W207_lim, 0.25)

## Run high days only model with partial pooling
fit_ryd_ss24_high <- metab(bayes_specs_W13, data=W13_dat_high[,-1])
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
plot_metab_preds(predict_metab(fit_NHC_raymond))

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















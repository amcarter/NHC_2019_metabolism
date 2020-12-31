## Model Metabolism #
# adapted from JRB script
# This version runs metabolism on NHC sites using the K600 values from Hall 1972

# update.packages(oldPkgs=c("streamMetabolizer","unitted"), dependencies=TRUE, 
#                 repos=c("http://owi.usgs.gov/R", "https://cran.rstudio.com"))
# devtools::install_github("USGS-R/streamMetabolizer", ref="develop")

library(rstan)
library(tidyverse)
library(ggplot2)
library(streamMetabolizer)
library(lubridate)
library(dygraphs)

setwd("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism/data")

## Read in Data ####
sites <- read_csv("siteData/NHCsite_metadata.csv") %>%
  slice(1:7)

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


## Set bayes specs #####
bayes_name <- mm_name(type='bayes', pool_K600="binned", 
                          err_obs_iid=TRUE, err_proc_iid = TRUE, 
                          ode_method = "trapezoid", deficit_src='DO_mod', 
                          engine='stan')

kq_hall <- read_csv("siteData/KQ_hall_prior_from_equation.csv")

set_up_model <- function(dat, bayes_name, nodes ){

  ## Set bayes specs
  bayes_specs <- specs(bayes_name)
  bayes_specs$keep_mcmc_data <- FALSE
  ## Based on range of log daily Q
  daily <- dat %>%
    mutate(date = as.Date(solar.time)) %>% group_by(date) %>%
    summarize(discharge = mean(discharge, na.rm = T))
  
  bayes_specs$K600_lnQ_nodes_centers <- nodes$lnQ
  bayes_specs$K600_lnQ_nodes_meanlog <- (nodes$lnK600)
  bayes_specs$K600_lnQ_nodes_sdlog <- c(rep(0.7, nrow(nodes)))

  ## Change sigma
  bayes_specs$K600_daily_sigma_sigma <- 0.05
  
  return(bayes_specs)
}


# Model Runs ####
# CBP 
bayes_specs_Hall <- set_up_model(CBP, bayes_name, 
                                 kq_hall[kq_hall$site =="CBP",])
fit <- metab(bayes_specs_Hall, CBP)
saveRDS(fit, "metabolism/modeled/fit_cbp_prior_churchillK.rds")
rm(fit)
gc()

bayes_specs_Hall$K600_lnQ_nodes_sdlog <- c(rep(0.05, 
                                          length(bayes_specs_Hall$K600_lnQ_nodes_sdlog)))
fit <- metab(bayes_specs_Hall, CBP)
saveRDS(fit, "metabolism/modeled/fit_cbp_fixed_churchillK.rds")
rm(fit)
gc()

# PM 
bayes_specs_Hall <- set_up_model(PM, bayes_name, 
                                 kq_hall[kq_hall$site =="PM",])
fit <- metab(bayes_specs_Hall, PM)
saveRDS(fit, "metabolism/modeled/fit_pm_prior_churchillK.rds")
rm(fit)
gc()

bayes_specs_Hall$K600_lnQ_nodes_sdlog <- c(rep(0.05, 
                                          length(bayes_specs_Hall$K600_lnQ_nodes_sdlog)))
fit <- metab(bayes_specs_Hall, PM)
saveRDS(fit, "metabolism/modeled/fit_pm_fixed_churchillK.rds")
rm(fit)
gc()

# WB 
bayes_specs_Hall <- set_up_model(WB, bayes_name, 
                                 kq_hall[kq_hall$site =="WB",])
fit <- metab(bayes_specs_Hall, WB)
saveRDS(fit, "metabolism/modeled/fit_wb_prior_churchillK.rds")
rm(fit)
gc()

bayes_specs_Hall$K600_lnQ_nodes_sdlog <- c(rep(0.05, 
                                          length(bayes_specs_Hall$K600_lnQ_nodes_sdlog)))
fit <- metab(bayes_specs_Hall, WB)
saveRDS(fit, "metabolism/modeled/fit_wb_fixed_churchillK.rds")
rm(fit)
gc()

# WBP 
bayes_specs_Hall <- set_up_model(WBP, bayes_name, 
                                 kq_hall[kq_hall$site =="WBP",])
fit <- metab(bayes_specs_Hall, WBP)
saveRDS(fit, "metabolism/modeled/fit_wbp_prior_churchillK.rds")
rm(fit)
gc()

bayes_specs_Hall$K600_lnQ_nodes_sdlog <- c(rep(0.05, 
                                          length(bayes_specs_Hall$K600_lnQ_nodes_sdlog)))
fit <- metab(bayes_specs_Hall, WBP)
saveRDS(fit, "metabolism/modeled/fit_wbp_fixed_churchillK.rds")
rm(fit)
gc()

# PWC 
bayes_specs_Hall <- set_up_model(PWC, bayes_name, 
                                 kq_hall[kq_hall$site =="PWC",])
fit <- metab(bayes_specs_Hall, PWC)
saveRDS(fit, "metabolism/modeled/fit_pwc_prior_churchillK.rds")
rm(fit)
gc()

bayes_specs_Hall$K600_lnQ_nodes_sdlog <- c(rep(0.05, 
                                          length(bayes_specs_Hall$K600_lnQ_nodes_sdlog)))
fit <- metab(bayes_specs_Hall, PWC)
saveRDS(fit, "metabolism/modeled/fit_pwc_fixed_churchillK.rds")
rm(fit)
gc()

# NHC 
for(i in 2016:2019){
  dd = ymd_hms(paste0((i),"-03-01 04:00:00"))
  dat <- NHC %>% 
    filter(solar.time >= dd,
           solar.time <= dd + 365*24*60*60)
  bayes_specs_Hall <- set_up_model(dat, bayes_name, 
                                   kq_hall[kq_hall$site =="NHC",])
  fit <- metab(bayes_specs_Hall, dat)
  saveRDS(fit, paste0("metabolism/modeled/fit_nhc_",i,"_prior_churchillK.rds"))
  rm(fit)
  gc()
  
  bayes_specs_Hall$K600_lnQ_nodes_sdlog <- c(rep(0.05, 
                                            length(bayes_specs_Hall$K600_lnQ_nodes_sdlog)))
  fit <- metab(bayes_specs_Hall, dat)
  saveRDS(fit, paste0("metabolism/modeled/fit_nhc_",i,"_fixed_churchillK.rds"))
  rm(fit)
  gc()
}

# UNHC 
for(i in 2016:2019){
  dd = ymd_hms(paste0(i,"-03-01 04:00:00"))
  dat <- UNHC %>% 
    filter(solar.time >= dd,
           solar.time <= dd + 365*24*60*60)
  bayes_specs_Hall <- set_up_model(dat, bayes_name, 
                                   kq_hall[kq_hall$site =="UNHC",])
  fit <- metab(bayes_specs_Hall, dat)
  saveRDS(fit, paste0("metabolism/modeled/fit_unhc_",i,"_prior_churchillK.rds"))
  rm(fit)
  gc()
  
  bayes_specs_Hall$K600_lnQ_nodes_sdlog <- c(rep(0.05, 
                                            length(bayes_specs_Hall$K600_lnQ_nodes_sdlog)))
  fit <- metab(bayes_specs_Hall, dat)
  saveRDS(fit, paste0("metabolism/modeled/fit_unhc_",i,"_fixed_churchillK.rds"))
  rm(fit)
  gc()
}

# #inspect fits ####
sites <- c("nhc2017",# "nhc2018", "nhc2019",
           #"unhc2017", "unhc2018", "unhc2019",
           "pm", "cbp", "wb")#, "wbp")
source("../src/streamMetabolizer/inspect_model_fits.r")

pdf("../figures/nhc_met_models_hall_v2.pdf", width = 9, height = 6)
  for(site in sites){
    fit <- readRDS(paste0("metabolism/modeled/fit_", site, "_fixed_hallK.rds"))
    plot_diagnostics(fit, site, ylim = c(-15, 7), lim = 7)
  }
dev.off()

# compile met data
all_met <- data.frame()
for(site in sites){
  fit <- readRDS(paste0("metabolism/modeled/", site, "_nreg_v2.rds"))
  met <- predict_metab(fit) %>%
    mutate(site = !!site)
  dat <- fit@data %>%
    select(date, DO_mgl = DO.obs, DO.sat, temp.water, depth, discharge) %>%
    group_by(date) %>%
    summarize_all(mean, na.rm = T)
  met <- fit@fit$daily %>% 
    select(date, K600 = K600_daily_50pct,
           K600.lower = K600_daily_2.5pct,
           K600.upper = K600_daily_97.5pct) %>%
    left_join(met, by = "date") %>%
    left_join(dat, by = "date")
  all_met <- bind_rows(all_met, met)
}

all_met <- all_met %>%
  mutate(site = case_when(
    grepl("^nhc", site) ~ "nhc",
    grepl("^unhc", site) ~ "unhc",
    TRUE ~ site
  ))

write_csv(all_met, 
          "C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/hall_50yl/code/data/NHC_metab_allsites_fixedHallK_v2.csv")
ss_met <- all_met %>% 
  mutate(doy = format(date, "%j")) %>%
  group_by(doy) %>%
  summarize(GPP.upper = quantile(GPP, .975, na.rm = T),
            GPP.lower = quantile(GPP, .025, na.rm = T),
            ER.upper = quantile(ER, .975, na.rm = T), 
            ER.lower = quantile(ER, .025, na.rm = T),
            GPP = median(GPP, na.rm = T),
            ER = median(ER, na.rm = T))

png("../figures/hall_met_comparison_hallK_v2.png", height = 4, width = 7.5, 
    units = "in", res = 300)
  m <- matrix(c(1,1,1,2,2), nrow = 1)
  layout(m)
  plot_hall_metab(ss_met, ylim = c(-15,10), doy = T)
  plot_kde_hall_metab(ss_met)
  mtext("All sites metabolism", outer = T, line = -3, cex = 1.2)
dev.off()

hall_met <- read_csv("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/hall_50yl/code/data/NHC_metab_allsites_fixedHallK_v2.csv")
     

png("../figures/metabolism_contours_K_estimates_v2.png", height = 5, width = 5, 
    units = "in", res = 300)
  plot_kde_metab(hall_met, col = "steelblue", lim = 3.5)
  par(new = T)
  plot_kde_metab(ss_met, col = "darkred", lim = 3.5)
  legend("topright", cex = 1.2,
         c("Hall 1970 K", "night regression K"),
         fill = c(alpha("steelblue", .75), alpha("darkred", .75)), 
         border = NA, bty = "n")
  mtext("NHC Metabolism Estimates (n = 1473)", cex = 1.2)
dev.off()
# ## Check binning
# Binning <- function(Site, thresh = 0.5){
#   fit_Site <- get_fit(Site)
#   
#   SM_output <- fit_Site$daily
#   SM_day <- get_data_daily(Site)
#   SM_KQbin <-  fit_Site$KQ_binned
#   SM_specs <- get_specs(Site)
#   
#   day <- data.frame(SM_day$discharge.daily, 
#                     SM_output$K600_daily_50pct, 
#                     SM_output$GPP_50pct,
#                     SM_output$K600_daily_Rhat,
#                     rep('daily', dim(SM_output)[1]))
#   colnames(day)<-c('Q', 'K600', 'GPP','Rhat', 'Group')
#   
#   gg<-ggplot(day, aes(x=log(Q), y = GPP, col=Rhat))+
#     geom_point() +
#     geom_hline(yintercept = thresh)
#   print(gg)
#   nodes<-data.frame(exp(SM_specs$K600_lnQ_nodes_centers), 
#                      exp(SM_KQbin$lnK600_lnQ_nodes_50pct), 
#                     rep('node', dim(SM_KQbin)[1]))
#   colnames(nodes)<-c('Q', 'K600', 'Group')
#   nodes$K600_prior <- exp(SM_specs$K600_lnQ_nodes_meanlog)
#   KQ_plot<-bind_rows(day,nodes)
#   
#   ggplot(data=KQ_plot, aes(x=log(Q), y=K600, group=Group, colour=Group)) + 
#     geom_point(size=3) +
#     #geom_line() + 
#     scale_color_manual(name="K-Q",
#                        breaks = c("daily", "node"),
#                        values=c("grey", "purple"),
#                        labels=c("Daily","Bin")) +
#     geom_point(aes(y = K600_prior),
#                size = 3, col = "purple", pch = 21) +
#     ylab("K600") +
#     xlab("logQ") +
#     theme_bw() +
#     theme(legend.position = "top")
# }
# 
# Binning(fit_raymond_ss05)
# Binning(fit_nreg_ss24)
# Binning(fit_nreg_fixedK1)
# Binning(fit_nreg_sd7_ss05)
# 
# ## Visualize
# plot_metab_preds(predict_metab(fit_nreg_sd7_ss05))
# plot_metab_preds(predict_metab(fit_nreg_fixedK))
# 
# ## K600 vs ER
# KvER <- get_fit(fit_raymond_ss24)
# KvER <- get_fit(fit_raymond_ss05)
# KvER <- get_fit(fit_nreg_ss24)
# KvER <- get_fit(fit_nreg_sd7_ss05)
# KvER <- get_fit(fit_nreg_fixedK)
# plot(KvER$daily$K600_daily_mean, KvER$daily$ER_daily_mean)
# cor(KvER$daily$K600_daily_mean, 
#     KvER$daily$ER_daily_mean, 
#     use = "na.or.complete")
# 
# ## Write Files
# writefiles <- function(mod){
#   data <- get_fit(mod)
#   for (i in seq_along(data)) {
#     filename = paste(names(data)[i], ".csv")
#     write.csv(data[[i]], filename)
#   }
#   write.csv(unlist(get_specs(mod)),"specs.csv")
#   write.csv(get_data_daily(mod), "datadaily.csv")
#   write.csv(get_data(mod),"mod_and_obs_DO.csv")
# }
# 
# ## Create new folder for site and write csv info
# ## Reset working directory!!
# getwd()
# writefiles(fit_207_lim)
# 






# Fixed K model ####

# the only place that Q enters into the model is to predict the K-Q relationship
# for this reason, we can use a fake Q as nodes and with our data to simplify
# using a fixed K/Q relationship as a prior. In this case, we are binning the 
# data and then replacing everything within each bin with a node, assigned to 
# the corresponding K derived from the Hall 1972 paper

# prep_fake_Q <- function(dat, bayes_specs_fixedK){
#   dat <- dat %>% 
#     mutate(date = as.Date(solar.time)) 
#   daily <- dat %>%
#     group_by(date) %>%
#     summarize(discharge = mean(discharge, na.rm = T)) %>%
#     mutate(logQ = log(discharge))
#   
#   nbins <-  length(bayes_specs_fixedK$K600_lnQ_nodes_centers)
#   
#   Qbreaks <- c(bayes_specs_fixedK$K600_lnQ_nodes_centers[1] -
#                diff(bayes_specs_fixedK$K600_lnQ_nodes_centers)[1]/2,
#                bayes_specs_fixedK$K600_lnQ_nodes_centers +
#                diff(bayes_specs_fixedK$K600_lnQ_nodes_centers)[1]/2)
#   rQ <- range(daily$logQ, na.rm = T)          
#   if(Qbreaks[1] > rQ[1]){Qbreaks[1] <- rQ[1]}
#   if(Qbreaks[nbins] < rQ[2]){Qbreaks[nbins] <- rQ[2]}
#   
#   daily <- daily %>%
#     mutate(Q = cut(daily$logQ, Qbreaks, labels = 1:nbins)) %>%
#     select(date, Q)
#   daily$Q <- as.numeric(daily$Q)
#   dat <- left_join(dat, daily, by = "date") %>%
#     select(-discharge) %>%
#     mutate(discharge = exp(Q)) %>%
#     select(-Q,-date)
# 
#   return(dat)
# 
# }
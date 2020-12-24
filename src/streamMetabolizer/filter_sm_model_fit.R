## Model Metabolism #
# adapted from JRB script
# This version runs metabolism on NHC sites using the K600 values from Hall 1972

# update.packages(oldPkgs=c("streamMetabolizer","unitted"), dependencies=TRUE, 
#                 repos=c("http://owi.usgs.gov/R", "https://cran.rstudio.com"))
# devtools::install_github("USGS-R/streamMetabolizer", ref="develop")

library(tidyverse)
library(ggplot2)
library(streamMetabolizer)
library(lubridate)
library(xts)
library(dygraphs)
# library(imputeTS)
# library(parallel)

setwd("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism/data")

## Read in and format Data ####
sites <- read_csv("siteData/NHCsite_metadata.csv") %>%
  slice(1:7)

fit <- readRDS("metabolism/modeled/fit_cbp_fixed_hallK.rds")
preds <- predict_metab(fit)
hist(preds$ER.lower)
plot_rhats(fit)

# filter high flow days and impossible values #
deltaQ_max = 2

RC <- read_csv("rating_curves/modified_ZQ_curves.csv")
nhcQ <- read_csv("rating_curves/NHC_UNHC_Q.csv", guess_max = 10000) %>%
  mutate(date = as.Date(with_tz(DateTime_UTC, tz = "EST"))) %>%
  group_by(date) %>%
  summarize(nhc_q = mean(NHC_Q, na.rm = T),
            unhc_q = mean(UNHC_Q, na.rm = T),
            deltaQ = max(NHC_Q, na.rm = T)/min(NHC_Q, na.rm = T),
            maxq = max(NHC_Q, na.rm = T),
            maxqu = max(UNHC_Q, na.rm = T),
            good_flow = ifelse(maxq <= RC$max_Q[1] &
                                 maxqu <= RC$max_Q[3] &
                                 deltaQ < deltaQ_max, 
                               TRUE, FALSE))

# plot(nhcQ$deltaQ, ylim = c(1, 2))  
# abline(h = 1.1)
# plot(nhcQ$maxqu, nhcQ$deltaQ, pch = 20, log = "xy")
# points(nhcQ$maxqu[nhcQ$maxq > RC$max_Q[1]], 
#        nhcQ$deltaQ[nhcQ$maxq > RC$max_Q[1]], pch = 20, col = 2)
# abline(v = RC$max_Q[3], h = 2)
# sum(nhcQ$deltaQ > 2 | nhcQ$maxq > RC$max_Q[1], na.rm = T)/sum(!is.na(nhcQ$deltaQ))
# plot(nhcQ$date, nhcQ$nhc_q, pch = 20, col = 2, log = "y")
# points(nhcQ$date[nhcQ$good_flow], nhcQ$nhc_q[nhcQ$good_flow], 
#       pch = 20, log = "y")

flow_dates <- nhcQ %>%
  select(date, nhc_q, deltaQ, good_flow) %>%
  filter(!is.na(good_flow))


plot_zoom(fit@data)#, colnames(fit@data)[-c(1,2,7,8)])

GPP_min = 0
ER_max = 0

filelist <- list.files("metabolism/modeled/nreg")
met_summary <- data.frame()
all_preds <- data.frame()
all_filled_preds <- data.frame()

for(file in filelist) {
  
  fit <- readRDS(paste0("metabolism/modeled/nreg/", file))
  tmp <- str_match(string = file, pattern = '^([a-z]+)([0-9]+)?_([a-z]+)')
  site <- tmp[2]
  method <- tmp[4]
  
  if(site %in% c("nhc", "unhc")) {
    year = as.numeric(tmp[3])
  } else { 
    year = 2019 
  }
  
  out <- filter_model(fit, flow_dates)
  preds <- out[[1]] 
  coverage <- data.frame(site = site,
                         year = year,
                         method = method) 
  coverage <- bind_cols(coverage, out[[2]])
  plot_metab(preds, main = paste(site, year))
  
  out <- fill_summarize_met(preds)
  cum <- out[[1]] %>%
    mutate(site = site,
           method = method)
  preds <- preds %>%
    mutate(site = site,
           method = method)
  met_sum <- bind_cols(coverage, out[[2]])
  
  met_summary <- bind_rows(met_summary, met_sum)
  all_preds <- bind_rows(all_preds, preds)
  all_filled_preds <- bind_rows(all_filled_preds, cum)
  
}


# #inspect fits ####
sites <- c("nhc2017", "nhc2018", "nhc2019",
           "unhc2017", "unhc2018", "unhc2019",
           "pm", "cbp", "wb", "wbp")
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
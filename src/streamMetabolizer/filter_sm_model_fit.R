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
source("../src/streamMetabolizer/inspect_model_fits.r")
## Filter bad flow days ####
sites <- read_csv("siteData/NHCsite_metadata.csv") %>%
  slice(1:7)

# fit <- readRDS("metabolism/modeled/fit_cbp_fixed_hallK.rds")
# preds <- predict_metab(fit)
# hist(preds$ER.lower)
# plot_rhats(fit)

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

write_csv(flow_dates, "rating_curves/flow_dates_filter.csv")
plot_zoom(fit@data)#, colnames(fit@data)[-c(1,2,7,8)])



# fix interpolated model outputs ####
# the most recent runs (12/27/2020) interpolated huge gaps, so I am going to 
# cut that and resave the fits

# # cbp
# fit <- readRDS(paste0("metabolism/modeled/churchill_fixed/", filelist[1]))
# fit <- readRDS(paste0("metabolism/modeled/churchill/", filelist[1]))
# fit <- readRDS(paste0("metabolism/modeled/raymond/", filelist[1]))
# file
# dd <- ymd_hms(c("2019-05-14 16:14:46",
#                 "2019-06-03 15:59:46",
#                 "2019-08-30 18:44:46",
#                 "2019-09-12 18:44:46"))
# dt <- as.Date(dd)
# fit@data[fit@data$solar.time >= dd[1] &
#            fit@data$solar.time <= dd[2],3:9] <- NA
# fit@data[fit@data$solar.time >= dd[3] &
#            fit@data$solar.time <= dd[4],3:9] <- NA
# fit@fit$daily[fit@fit$daily$date >= dt[1] &
#                 fit@fit$daily$date <= dt[2],2:71] <- NA
# fit@fit$daily[fit@fit$daily$date >= dt[1] &
#                 fit@fit$daily$date <= dt[2],74] <- "missing data"
# fit@fit$daily[fit@fit$daily$date >= dt[3] &
#                 fit@fit$daily$date <= dt[4],2:71] <- NA
# fit@fit$daily[fit@fit$daily$date >= dt[3] &
#                 fit@fit$daily$date <= dt[4],74] <- "missing data"
# fit@data_daily[fit@data_daily$date >= dt[1] &
#                  fit@data_daily$date <= dt[2],2:7] <- NA
# fit@data_daily[fit@data_daily$date >= dt[3] &
#                  fit@data_daily$date <= dt[4],2:7] <- NA
# saveRDS(fit, paste0("metabolism/modeled/churchill_fixed/", filelist[1]))
# saveRDS(fit, paste0("metabolism/modeled/churchill/", filelist[1]))
# saveRDS(fit, paste0("metabolism/modeled/raymond/", filelist[1]))

# # nhc 2016
# fit <- readRDS(paste0("metabolism/modeled/churchill_fixed/", filelist[2]))
# fit <- readRDS(paste0("metabolism/modeled/churchill/", filelist[2]))
# fit <- readRDS(paste0("metabolism/modeled/raymond/", filelist[2]))
# file
# dd <- ymd_hms(c("2016-11-22 09:44:51",
#                 "2016-12-08 08:59:51",
#                 "2016-12-11 04:14:51",
#                 "2016-12-28 07:14:51"))
# dt <- as.Date(dd)
# fit@data[fit@data$solar.time >= dd[1] &
#            fit@data$solar.time <= dd[2],3:9] <- NA
# fit@data[fit@data$solar.time >= dd[3] &
#            fit@data$solar.time <= dd[4],3:9] <- NA
# fit@fit$daily[fit@fit$daily$date >= dt[1] &
#                 fit@fit$daily$date <= dt[2],2:71] <- NA
# fit@fit$daily[fit@fit$daily$date >= dt[1] &
#                 fit@fit$daily$date <= dt[2],74] <- "missing data"
# fit@fit$daily[fit@fit$daily$date >= dt[3] &
#                 fit@fit$daily$date <= dt[4],2:71] <- NA
# fit@fit$daily[fit@fit$daily$date >= dt[3] &
#                 fit@fit$daily$date <= dt[4],74] <- "missing data"
# fit@data_daily[fit@data_daily$date >= dt[1] &
#                  fit@data_daily$date <= dt[2],2:7] <- NA
# fit@data_daily[fit@data_daily$date >= dt[3] &
#                  fit@data_daily$date <= dt[4],2:7] <- NA
# saveRDS(fit, "metabolism/modeled/churchill_fixed/fit_nhc_2016_fixed_churchillK.rds")
# saveRDS(fit, "metabolism/modeled/churchill/fit_nhc_2016_prior_churchillK.rds")
# saveRDS(fit, "metabolism/modeled/raymond/fit_nhc_2016_uninformed_raymondK.rds")

# # nhc 2017
# fit <- readRDS(paste0("metabolism/modeled/churchill_fixed/", filelist[3]))
# fit <- readRDS(paste0("metabolism/modeled/churchill/", filelist[3]))
# fit <- readRDS(paste0("metabolism/modeled/raymond/", filelist[3]))
# file
# dd <- ymd_hms(c("2017-05-16 08:29:51",
#                 "2017-05-30 10:59:51",
#                 "2017-06-21 09:29:51",
#                 "2017-06-27 16:14:51",
#                 "2017-08-01 09:44:51",
#                 "2017-08-09 10:44:51",
#                 "2017-08-17 11:59:51",
#                 "2017-08-18 12:29:51"))
# dt <- as.Date(dd)
# fit@data[fit@data$solar.time >= dd[1] &
#            fit@data$solar.time <= dd[2],3:9] <- NA
# fit@data[fit@data$solar.time >= dd[3] &
#            fit@data$solar.time <= dd[4],3:9] <- NA
# fit@data[fit@data$solar.time >= dd[5] &
#            fit@data$solar.time <= dd[6],3:9] <- NA
# fit@data[fit@data$solar.time >= dd[7] &
#            fit@data$solar.time <= dd[8],3:9] <- NA
# fit@fit$daily[fit@fit$daily$date >= dt[1] &
#                 fit@fit$daily$date <= dt[2],2:71] <- NA
# fit@fit$daily[fit@fit$daily$date >= dt[1] &
#                 fit@fit$daily$date <= dt[2],74] <- "missing data"
# fit@fit$daily[fit@fit$daily$date >= dt[3] &
#                 fit@fit$daily$date <= dt[4],2:71] <- NA
# fit@fit$daily[fit@fit$daily$date >= dt[3] &
#                 fit@fit$daily$date <= dt[4],74] <- "missing data"
# fit@fit$daily[fit@fit$daily$date >= dt[5] &
#                 fit@fit$daily$date <= dt[6],2:71] <- NA
# fit@fit$daily[fit@fit$daily$date >= dt[5] &
#                 fit@fit$daily$date <= dt[6],74] <- "missing data"
# fit@fit$daily[fit@fit$daily$date >= dt[7] &
#                 fit@fit$daily$date <= dt[8],2:71] <- NA
# fit@fit$daily[fit@fit$daily$date >= dt[7] &
#                 fit@fit$daily$date <= dt[8],74] <- "missing data"
# fit@data_daily[fit@data_daily$date >= dt[1] &
#                  fit@data_daily$date <= dt[2],2:7] <- NA
# fit@data_daily[fit@data_daily$date >= dt[3] &
#                  fit@data_daily$date <= dt[4],2:7] <- NA
# fit@data_daily[fit@data_daily$date >= dt[5] &
#                  fit@data_daily$date <= dt[6],2:7] <- NA
# fit@data_daily[fit@data_daily$date >= dt[7] &
#                  fit@data_daily$date <= dt[8],2:7] <- NA
# saveRDS(fit, "metabolism/modeled/churchill_fixed/fit_nhc_2017_fixed_churchillK.rds")
# saveRDS(fit, "metabolism/modeled/churchill/fit_nhc_2017_prior_churchillK.rds")
# saveRDS(fit, "metabolism/modeled/raymond/fit_nhc_2017_uninformed_raymondK.rds")

# # nhc 2018
# fit <- readRDS(paste0("metabolism/modeled/churchill_fixed/", filelist[4]))
# fit <- readRDS(paste0("metabolism/modeled/churchill/", filelist[4]))
# fit <- readRDS(paste0("metabolism/modeled/raymond/", filelist[4]))
# file
# dd <- ymd_hms(c("2018-04-30 11:14:51",
#                 "2018-05-16 18:14:51",
#                 "2018-09-11 09:29:51",
#                 "2018-09-19 09:29:51",
#                 "2018-12-19 10:29:51",
#                 "2019-01-07 13:29:51"
#                 ))
# dt <- as.Date(dd)
# fit@data[fit@data$solar.time >= dd[1] &
#            fit@data$solar.time <= dd[2],3:9] <- NA
# fit@data[fit@data$solar.time >= dd[3] &
#            fit@data$solar.time <= dd[4],3:9] <- NA
# fit@data[fit@data$solar.time >= dd[5] &
#            fit@data$solar.time <= dd[6],3:9] <- NA
# fit@fit$daily[fit@fit$daily$date >= dt[1] &
#                 fit@fit$daily$date <= dt[2],2:71] <- NA
# fit@fit$daily[fit@fit$daily$date >= dt[1] &
#                 fit@fit$daily$date <= dt[2],74] <- "missing data"
# fit@fit$daily[fit@fit$daily$date >= dt[3] &
#                 fit@fit$daily$date <= dt[4],2:71] <- NA
# fit@fit$daily[fit@fit$daily$date >= dt[3] &
#                 fit@fit$daily$date <= dt[4],74] <- "missing data"
# fit@fit$daily[fit@fit$daily$date >= dt[5] &
#                 fit@fit$daily$date <= dt[6],2:71] <- NA
# fit@fit$daily[fit@fit$daily$date >= dt[5] &
#                 fit@fit$daily$date <= dt[6],74] <- "missing data"
# fit@data_daily[fit@data_daily$date >= dt[1] &
#                  fit@data_daily$date <= dt[2],2:7] <- NA
# fit@data_daily[fit@data_daily$date >= dt[3] &
#                  fit@data_daily$date <= dt[4],2:7] <- NA
# fit@data_daily[fit@data_daily$date >= dt[5] &
#                  fit@data_daily$date <= dt[6],2:7] <- NA
# saveRDS(fit, "metabolism/modeled/churchill_fixed/fit_nhc_2018_fixed_churchillK.rds")
# saveRDS(fit, "metabolism/modeled/churchill/fit_nhc_2018_fixed_churchillK.rds")
# saveRDS(fit, "metabolism/modeled/raymond/fit_nhc_2018_fixed_churchillK.rds")

# # unhc 2017
# fit <- readRDS(paste0("metabolism/modeled/churchill_fixed/", filelist[7]))
# fit <- readRDS(paste0("metabolism/modeled/churchill/", filelist[7]))
# fit <- readRDS(paste0("metabolism/modeled/raymond/", filelist[7]))
# 
# dd <- ymd_hms(c("2017-05-30 08:59:41",
#                 "2017-06-13 11:44:41",
#                 "2017-08-01 10:29:41",
#                 "2017-08-07 14:29:41",
#                 ))
# dt <- as.Date(dd)
# fit@data[fit@data$solar.time >= dd[1] &
#            fit@data$solar.time <= dd[2],3:9] <- NA
# fit@data[fit@data$solar.time >= dd[3] &
#            fit@data$solar.time <= dd[4],3:9] <- NA
# fit@fit$daily[fit@fit$daily$date >= dt[1] &
#                 fit@fit$daily$date <= dt[2],2:71] <- NA
# fit@fit$daily[fit@fit$daily$date >= dt[1] &
#                 fit@fit$daily$date <= dt[2],74] <- "missing data"
# fit@fit$daily[fit@fit$daily$date >= dt[3] &
#                 fit@fit$daily$date <= dt[4],2:71] <- NA
# fit@fit$daily[fit@fit$daily$date >= dt[3] &
#                 fit@fit$daily$date <= dt[4],74] <- "missing data"
# fit@data_daily[fit@data_daily$date >= dt[1] &
#                  fit@data_daily$date <= dt[2],2:7] <- NA
# fit@data_daily[fit@data_daily$date >= dt[3] &
#                  fit@data_daily$date <= dt[4],2:7] <- NA
# saveRDS(fit, "metabolism/modeled/churchill_fixed/fit_unhc_2017_fixed_churchillK.rds")
# saveRDS(fit, "metabolism/modeled/churchill/fit_unhc_2016_prior_churchillK.rds")
# saveRDS(fit, "metabolism/modeled/raymond/fit_unhc_2016_uninformed_raymondK.rds")

# # unhc 2018
# fit <- readRDS(paste0("metabolism/modeled/churchill_fixed/", filelist[4]))
# fit <- readRDS(paste0("metabolism/modeled/churchill/", filelist[4]))
# fit <- readRDS(paste0("metabolism/modeled/raymond/", filelist[4]))
# file
# dd <- ymd_hms(c("2018-04-30 11:14:51",
#                 "2018-05-16 18:14:51",
#                 "2018-09-11 09:29:51",
#                 "2018-09-19 09:29:51",
#                 "2018-12-19 10:29:51",
#                 "2019-01-07 13:29:51"
#                 ))
# dt <- as.Date(dd)
# fit@data[fit@data$solar.time >= dd[1] &
#            fit@data$solar.time <= dd[2],3:9] <- NA
# fit@data[fit@data$solar.time >= dd[3] &
#            fit@data$solar.time <= dd[4],3:9] <- NA
# fit@data[fit@data$solar.time >= dd[5] &
#            fit@data$solar.time <= dd[6],3:9] <- NA
# fit@fit$daily[fit@fit$daily$date >= dt[1] &
#                 fit@fit$daily$date <= dt[2],2:71] <- NA
# fit@fit$daily[fit@fit$daily$date >= dt[1] &
#                 fit@fit$daily$date <= dt[2],74] <- "missing data"
# fit@fit$daily[fit@fit$daily$date >= dt[3] &
#                 fit@fit$daily$date <= dt[4],2:71] <- NA
# fit@fit$daily[fit@fit$daily$date >= dt[3] &
#                 fit@fit$daily$date <= dt[4],74] <- "missing data"
# fit@fit$daily[fit@fit$daily$date >= dt[5] &
#                 fit@fit$daily$date <= dt[6],2:71] <- NA
# fit@fit$daily[fit@fit$daily$date >= dt[5] &
#                 fit@fit$daily$date <= dt[6],74] <- "missing data"
# fit@data_daily[fit@data_daily$date >= dt[1] &
#                  fit@data_daily$date <= dt[2],2:7] <- NA
# fit@data_daily[fit@data_daily$date >= dt[3] &
#                  fit@data_daily$date <= dt[4],2:7] <- NA
# fit@data_daily[fit@data_daily$date >= dt[5] &
#                  fit@data_daily$date <= dt[6],2:7] <- NA
# saveRDS(fit, "metabolism/modeled/churchill_fixed/fit_unhc_2018_fixed_churchillK.rds")
# saveRDS(fit, "metabolism/modeled/churchill/fit_unhc_2018_fixed_churchillK.rds")
# saveRDS(fit, "metabolism/modeled/raymond/fit_unhc_2018_fixed_churchillK.rds")

# # pm
# fit <- readRDS(paste0("metabolism/modeled/churchill_fixed/", filelist[5]))
# fit <- readRDS(paste0("metabolism/modeled/churchill/", filelist[5]))
# fit <- readRDS(paste0("metabolism/modeled/raymond/", filelist[5]))
# file
# dd <- ymd_hms(c("2019-03-27 09:14:50",
#                 "2019-04-16 11:29:50",
#                 "2019-04-22 09:44:50",
#                 "2019-05-08 14:59:50"))
# dt <- as.Date(dd)
# fit@data[fit@data$solar.time >= dd[1] &
#            fit@data$solar.time <= dd[2],3:9] <- NA
# fit@data[fit@data$solar.time >= dd[3] &
#            fit@data$solar.time <= dd[4],3:9] <- NA
# fit@fit$daily[fit@fit$daily$date >= dt[1] &
#                 fit@fit$daily$date <= dt[2],2:71] <- NA
# fit@fit$daily[fit@fit$daily$date >= dt[1] &
#                 fit@fit$daily$date <= dt[2],74] <- "missing data"
# fit@fit$daily[fit@fit$daily$date >= dt[3] &
#                 fit@fit$daily$date <= dt[4],2:71] <- NA
# fit@fit$daily[fit@fit$daily$date >= dt[3] &
#                 fit@fit$daily$date <= dt[4],74] <- "missing data"
# fit@data_daily[fit@data_daily$date >= dt[1] &
#                  fit@data_daily$date <= dt[2],2:7] <- NA
# fit@data_daily[fit@data_daily$date >= dt[3] &
#                  fit@data_daily$date <= dt[4],2:7] <- NA
# saveRDS(fit, paste0("metabolism/modeled/churchill_fixed/", filelist[5]))
# saveRDS(fit, paste0("metabolism/modeled/churchill/", filelist[5]))
# saveRDS(fit, paste0("metabolism/modeled/raymond/", filelist[5]))

# Compile model outputs ####
filelist <- list.files("metabolism/modeled/churchill_fixed/")
GPP_min = 0
ER_max = 0

met_summary <- data.frame()
all_preds <- data.frame()
all_filled_preds <- data.frame()

file <- filelist[1]

 pdf("../figures/model_diagnostics_nightreg.pdf", width = 9, height = 6)

for(file in filelist) {
  # uninformed churchill ests
  # fit <- readRDS(paste0("metabolism/modeled/nreg/", file))
  # tmp <- str_match(string = file, 
  #                  pattern = '^([a-z]+)([0-9]+)?_([a-z]+_[a-z0-9]+)')
  # fit <- readRDS(paste0("metabolism/modeled/churchill_uninformed/", file))
  # fixed churchill ests
  fit <- readRDS(paste0("metabolism/modeled/churchill_fixed/", file))
  # churchill ests
  # fit <- readRDS(paste0("metabolism/modeled/churchill/", file))
  tmp <- str_match(string = file,
                   pattern = '^[a-z]+_([a-z]+)_?([0-9]+)?_([a-z]+_[a-zA-Z]+)')
  # nightreg
  # fit <- readRDS(paste0("metabolism/modeled/nreg/", file))
  # tmp <- str_match(string = file, pattern = '^([a-z]+)([0-9]+)?_([a-z]+)')
  # plot_zoom(fit@data)
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
  
  out <- fill_summarize_met(preds)
  cum <- out[[1]] %>%
    mutate(site = site,
           method = method)
  preds <- preds %>%
    mutate(site = site,
           method = method)
  
  bad <- unique(preds$date[is.na(preds$GPP) & is.na(preds$ER)])
  dat <- fit@data %>%
    filter(!(date %in% bad))
  plot_zoom(dat)
  mcmc <- get_mcmc(fit)
  rstan::traceplot(mcmc, pars = c("K600_daily[111]",
                                  "K600_daily[112]",
                                  "K600_daily[113]",
                                  "K600_daily[114]",
                                  "K600_daily[115]",
                                  "K600_daily[116]",
                                  "K600_daily[117]",
                                  "K600_daily[118]",
                                  "K600_daily[119]",
                                  "K600_daily[120]"), nrow = 5)
  # plot_metab(preds, main = paste(site, year))
  plot_diagnostics(fit, preds, paste(site, year, method),
                   ylim = c(-15, 7), lim = 7)
  met_sum <- bind_cols(coverage, out[[2]])
  
  met_summary <- bind_rows(met_summary, met_sum)
  all_preds <- bind_rows(all_preds, preds)
  all_filled_preds <- bind_rows(all_filled_preds, cum)
  
}
dev.off()


ch_un <- all_preds %>%
  as_tibble() %>%
  filter(method == "uninformed_churchillK", 
         !(site %in% c("wbp", "pwc")))
ch_fx <- all_preds %>%
  as_tibble() %>%
  filter(method == "fixed_churchillK", 
         !(site %in% c("wbp", "pwc")))

write_csv(all_preds, "metabolism/compiled/churchill_fixed_k_met.csv")
write_csv(all_filled_preds, "metabolism/compiled/churchill_fixed_k_filled_met.csv")
write_csv(met_summary, "metabolism/compiled/churchill_fixed_summary.csv")
write_csv(all_preds, "metabolism/compiled/churchill_k_met.csv")
write_csv(all_filled_preds, "metabolism/compiled/churchill_k_filled_met.csv")
write_csv(met_summary, "metabolism/compiled/churchill_summary.csv")

# #inspect fits ####
filelist <- list.files("metabolism/modeled/churchill_fixed")

pdf("../figures/model_diagnostics_fixed_churchill.pdf", width = 9, height = 6)
  for(file in filelist) {
    fit <- readRDS(paste0("metabolism/modeled/churchill_fixed/", file))
    site <- str_match(string = file, pattern = '^[a-z]+_([a-z]+)')[2]
    
    plot_diagnostics(fit, site, ylim = c(-15, 7), lim = 7)
  }
dev.off()
ss_met <- all_preds %>% 
  mutate(doy = format(date, "%j")) %>%
  group_by(doy) %>%
  summarize(GPP.upper = quantile(GPP, .975, na.rm = T),
            GPP.lower = quantile(GPP, .025, na.rm = T),
            ER.upper = quantile(ER, .975, na.rm = T), 
            ER.lower = quantile(ER, .025, na.rm = T),
            GPP = median(GPP, na.rm = T),
              ER = median(ER, na.rm = T))

png("../figures/hall_met_comparison_nightreg.png", height = 6, width = 7, 
    units = "in", res = 300)
  m <- matrix(c(1,2,1,2,1, 
                2,1,3,1,3), nrow = 2)
  layout(m)
  plot_hall_metab(ss_met, ylim = c(-15,10), doy = T)
  plot_kde_hall_metab(ss_met, lim = 5)
  plot_k(all_preds)
  mtext("night regression", outer = T, line = -3, cex = 1.2)
dev.off()


# Convert to Carbon ####
range(ch_un$GPP *14/32, na.rm = T)
mean(ch_un$GPP *14/32, na.rm = T)
range(ch_un$ER *14/32, na.rm = T)
mean(ch_un$ER *14/32, na.rm = T)
range(ch_fx$GPP *14/32, na.rm = T)
mean(ch_fx$GPP *14/32, na.rm = T)
range(ch_fx$ER *14/32, na.rm = T)
mean(ch_fx$ER *14/32, na.rm = T)

met <- read_csv("metabolism/compiled/churchill_summary.csv") %>%
  mutate(nep = -er_cum + gpp_cum)

met_x <- met %>%
  filter(method == "fixed_churchillK") %>%
  select(site, year, n_days = total_days, pctcoverage,
         gpp_median, gpp_max, er_median, er_max, 
         gpp_max10d, er_max10d, gpp_cum, er_cum, nep) %>%
  mutate(across(all_of(c("gpp_median", "er_median",
                         "gpp_max", "er_max",
                         "gpp_cum", "er_cum", "nep")), ~.*14/32))
met_u <- met %>%
  filter(method == "uninformed_churchillK")%>%
  select(site, year, n_days = total_days, pctcoverage,
         gpp_median, gpp_max, er_median, er_max, 
         gpp_max10d, er_max10d, gpp_cum, er_cum, nep) %>%
  mutate(across(all_of(c("gpp_median", "er_median",
                         "gpp_max", "er_max",
                         "gpp_cum", "er_cum", "nep")), ~.*14/32))

range(met_x$nep, na.rm = T)
range(met_u$nep * 14/32, na.rm = T)
range(met_x$er_cum * 14/32/1000, na.rm = T)
write_csv(met_x, "metabolism/compiled/fixedKmet_table.csv")
write_csv(met_u, "metabolism/compiled/freeKmet_table.csv")



# Compile Hall Data ####
hall <- read_csv("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/hall_50yl/data/hall/hall_table_15.csv")  %>%
  mutate(date = as.Date(date, format = "%m/%d/%Y")) %>%
  group_by(site, date) %>%
  summarize(gpp_gcm2d = mean(GPP_gO2m2d, na.rm = T) * 14 / 32,
            er_gcm2d = mean(ER_gO2m2d, na.rm = T) * 14 / 32) %>%
  ungroup()
hall <- read_csv("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/hall_50yl/data/hall/hall_table_15.csv")  %>%
  mutate(doy = as.numeric(format(as.Date(date, format = "%m/%d/%Y"), "%j"))) %>%
  group_by(doy) %>%
  summarize(gpp_gcm2d = mean(GPP_gO2m2d, na.rm = T) * 14 / 32,
            er_gcm2d = mean(ER_gO2m2d, na.rm = T) * 14 / 32) %>%
  ungroup()
n <- nrow(hall)
hallf <- data.frame(doy = 366,
                   gpp_gcm2d = 0.258,
                   er_gcm2d = 0.394) %>%
  bind_rows(hall)

hallf <- data.frame(doy = seq(1:366)) %>%
  full_join(hallf) %>%
  arrange(doy) %>%
  mutate(across(-doy, na.approx, na.rm = F)) %>%
  filter(doy != 366)  
ncon <- sum(hall$site =="Concrete")
nblk <- sum(hall$site =="Blackwood")
nwb <- sum(hall$site =="Wood Bridge")

dd <- as.Date("1968-04-14")
h68 <- hall %>%
  filter(date >= dd,
         date < dd + 365,
         site == "Concrete")
h69 <- hall %>%
  filter(date >= dd+ 365,
         date < dd + 2*365,
         site == "Concrete")
h_all <- hall %>%
  mutate(doy = format(date, "%j")) %>%
  group_by(doy) %>%
  summarize(gpp_gcm2d = mean(gpp_gcm2d, na.rm = T),
            er_gcm2d = mean(er_gcm2d, na.rm = T))
rd <- range(hall$date[hall$site =="Wood Bridge"])
dates <- data.frame(date = seq(rd[1], rd[2], by = "day"))
nrow(dates)
h_con <- hall %>%
  filter(site == "Concrete") %>%
  select(-site) %>%
  full_join(dates) %>%
  arrange(date) %>%
  mutate(across(-date, na.approx, na.rm = F))

h_c <- hall %>%
  filter(site == "Concrete") 
 
h_wb <- hall %>%
  filter(site == "Wood Bridge") %>%
  mutate(doy = format(date, "%j")) %>%
  group_by(doy) %>%
  summarize(gpp_gcm2d = mean(gpp_gcm2d, na.rm = T),
            er_gcm2d = mean(er_gcm2d, na.rm = T))

h68 <- h_con %>%
  filter(date >= dd ,
         date < dd + 365)
h68 <- h_con %>%
  filter(date >= dd +365,
         date < dd + 365*2)

h68 <- hall                    
met <- h68 %>%
    summarize(gpp_mean = mean(gpp_gcm2d, na.rm = T),
              gpp_median = median(gpp_gcm2d, na.rm = T),
              gpp_max = max(gpp_gcm2d, na.rm = T),
              er_mean = -mean(er_gcm2d, na.rm = T),
              er_median = -median(er_gcm2d, na.rm = T),
              er_max = -max(er_gcm2d, na.rm = T))
met$site = "all"
met$year = NA_real_
  
cum <- h68 %>%
 # select(-site) %>%
 # mutate(across(-date, na.approx, na.rm = F)) %>%
 mutate(across(-doy, cumsum, .names = "{col}_cum")) 

n <- nrow(cum)
l = (n-9)
weekly <- tibble(date = cum$date[1:l],
                 GPP_week = rep(NA_real_, l),
                 ER_week = rep(NA_real_, l))
for(i in 1:l){
  weekly$GPP_week[i] <- sum(cum$gpp_gcm2d[i:(i+9)]) 
  weekly$ER_week[i] <- sum(cum$er_gcm2d[i:(i+9)]) 
}
met$gpp_max10d <- weekly$date[which.max(weekly$GPP_week)]
met$er_max10d <- weekly$date[which.max(weekly$ER_week)]
met$gpp_cum <- cum$gpp_gcm2d_cum[n]*365/(n)
met$er_cum <- cum$er_gcm2d_cum[n]*365/(n)
met$daterange <- as.character(paste(cum$date[1], "-", cum$date[nrow(cum)]))
met$pctcoverage <- nrow(h68)/n

 # met_hall <- data.frame()
met_hall <- bind_rows(met_hall, met)

####

sites <- c("nhc2017", "nhc2018", "nhc2019",
           "unhc2017", "unhc2018", "unhc2019",
           "pm", "cbp", "wb", "wbp")

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
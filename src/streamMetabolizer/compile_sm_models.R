## Model Metabolism #
# adapted from JRB script
# This version runs metabolism on NHC sites using the K600 values from Hall 1972
source("NHC_2019_metabolism/src/streamMetabolizer/inspect_model_fits.r")
# Compile model outputs ####

flow_dates <- read_csv("NHC_2019_metabolism/data/rating_curves/flow_dates_filter.csv")

filelist <- list.files("NHC_2019_metabolism/data/metabolism/modeled/finalQ")#[c(-5, -7)]
GPP_min = 0
ER_max = 0

met_summary <- data.frame()
all_preds <- data.frame()
all_filled_preds <- data.frame()

 pdf("figures/model_diagnostics_raymond.pdf", width = 9, height = 6)

for(file in filelist) {
  # uninformed raymond ests
  fit <- readRDS(paste0("NHC_2019_metabolism/data/metabolism/modeled/finalQ/", file))
  tmp <- str_match(string = file,
                   pattern = '^[a-z]+_([a-z]+)_?([0-9]+)?_([a-z]+_[a-z]+)')
  site <- tmp[2]
  method <- tmp[4]
  year = 2019
  
  if(site %in% c("nhc", "unhc")) {
    year = as.numeric(tmp[3])
  } else { 
    year = 2019 
  }
  if(site == "wbp"){
    fit@data <- fit@data %>% filter(date <= as.Date("2020-03-20"))
    fit@fit$daily <- fit@fit$daily %>% filter(date <= as.Date("2020-03-20"))
    fit@data_daily <- fit@data_daily %>% filter(date <= as.Date("2020-03-20"))
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
           year = year,
           method = method)
  preds <- preds %>%
    mutate(site = site,
           year = year,
           method = method)
  
  # bad <- unique(preds$date[is.na(preds$GPP) & is.na(preds$ER)])
  # dat <- fit@data %>%
  #   filter(!(date %in% bad))
  # plot_zoom(dat)
  # mcmc <- get_mcmc(fit)
  # rstan::traceplot(mcmc, pars = c("K600_daily[11]",
  #                                 "K600_daily[12]",
  #                                 "K600_daily[13]",
  #                                 "K600_daily[14]",
  #                                 "K600_daily[15]",
  #                                 "K600_daily[16]",
  #                                 "K600_daily[17]",
  #                                 "K600_daily[18]",
  #                                 "K600_daily[19]",
  #                                 "K600_daily[20]"), nrow = 5)
  plot_diagnostics(fit, preds, paste(site, year, method),
                   ylim = c(-15, 7), lim = 7)
  met_sum <- bind_cols(coverage, out[[2]])
  
  met_summary <- bind_rows(met_summary, met_sum)
  all_preds <- bind_rows(all_preds, preds)
  all_filled_preds <- bind_rows(all_filled_preds, cum)
  
}
dev.off()

saveRDS(list(preds = all_preds, 
             summary = met_summary, 
             cumulative = all_filled_preds),
        "NHC_2019_metabolism/data/metabolism/compiled/raymond_met.rds")


# #inspect fits ####
# ray_met <- all_preds %>% 
#   mutate(doy = format(date, "%j")) %>%
#   group_by(doy) %>%
#   summarize(GPP.upper = quantile(GPP, .975, na.rm = T),
#             GPP.lower = quantile(GPP, .025, na.rm = T),
#             ER.upper = quantile(ER, .975, na.rm = T), 
#             ER.lower = quantile(ER, .025, na.rm = T),
#             GPP = median(GPP, na.rm = T),
#               ER = median(ER, na.rm = T))

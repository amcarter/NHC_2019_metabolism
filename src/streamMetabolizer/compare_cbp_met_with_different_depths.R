## Compare model results for metabolism run at CBP site with my depths
#   and with the median depth from Hall 1970

source("NHC_2019_metabolism/src/streamMetabolizer/inspect_model_fits.r")
# Compile model outputs ####

flow_dates <- read_csv("NHC_2019_metabolism/data/rating_curves/flow_dates_filter.csv")

filelist <- list.files("NHC_2019_metabolism/data/metabolism/modeled/finalQ")#[c(-5, -7)]
GPP_min = 0
ER_max = 0

met_summary <- data.frame()
all_preds <- data.frame()
all_filled_preds <- data.frame()


pdf("figures/model_diagnostics_CBP_depth_comparison.pdf", width = 9, height = 6)

  year = 2019
  site = 'cbp'
  # estimates based on my depth data
  fit <- readRDS("NHC_2019_metabolism/data/metabolism/modeled/finalQ/fit_cbp_uninformed_raymond_K.rds")
  out <- filter_model(fit, flow_dates)
  preds <- out[[1]] 
  coverage <- data.frame(site = site,
                         year = year,
                         depth = 'today') 
  coverage <- bind_cols(coverage, out[[2]])
  
  out <- fill_summarize_met(preds)
  cum <- out[[1]] %>%
    mutate(site = site,
           year = year,
           depth = 'today')
  preds_1 <- preds %>%
    mutate(site = site,
           year = year,
           depth = 'today')
  
  plot_diagnostics(fit, preds_1, paste(site, year, method),
                   ylim = c(-15, 7), lim = 7)
  met_sum <- bind_cols(coverage, out[[2]])
  
  met_summary <- bind_rows(met_summary, met_sum)
  all_preds <- bind_rows(all_preds, preds_1)
  all_filled_preds <- bind_rows(all_filled_preds, cum)

  # estimates based on Hall 1970 depth data
  fit <- readRDS("NHC_2019_metabolism/data/metabolism/modeled/fit_cbp_uninformed_raymond_K_median_depth_Hall.rds")
  out <- filter_model(fit, flow_dates)
  preds <- out[[1]] 
  coverage <- data.frame(site = site,
                         year = year,
                         depth = 'hall') 
  coverage <- bind_cols(coverage, out[[2]])
  
  out <- fill_summarize_met(preds)
  cum <- out[[1]] %>%
    mutate(site = site,
           year = year,
           depth = 'hall')
  preds_2 <- preds %>%
    mutate(site = site,
           year = year,
           depth = 'hall')
  
  plot_diagnostics(fit, preds_2, paste(site, year, method),
                   ylim = c(-15, 7), lim = 7)
  met_sum <- bind_cols(coverage, out[[2]])
  
  met_summary <- bind_rows(met_summary, met_sum)
  all_preds <- bind_rows(all_preds, preds_2)
  all_filled_preds <- bind_rows(all_filled_preds, cum)
  

dev.off()


select(met_summary, site, year, depth, gpp_median, er_median, gpp_cum, er_cum)

today <- read_csv("C:/Users/Alice Carter/git/nhc_50yl/NHC_2019_metabolism/data/metabolism/processed/CBP.csv")
hall <- read_csv('C:/Users/Alice Carter/git/nhc_50yl/hall_50yl/code/data/hall_discharge_temp_daily_corrected_dates.csv')

png('figures/cbp_metabolism_comparison_depth_adjustment.png',
    width = 6, height = 6, res = 300, units = 'in')
par(mfrow = c(2, 2), 
    oma = c(1, 0, 2, 0))
plot_kde_hall_metab(preds_1, lim = 7, site = "CBP")
mtext('Estimates based on 2019 depths', cex = 1.2, line = .25)
plot_kde_hall_metab(preds_2, lim = 7, site = "CBP")
mtext('Estimates adjusted to 1969 depths', cex = 1.2, line = .25)
plot(density(today$depth, na.rm = T), col = alpha('grey25', .75), lwd = 2, 
     main = '', xlab = 'average depth (m)', xlim = c(0,1))
lines(density(hall$depth_m, na.rm = T), col = alpha('brown3', 0.75), lwd = 2)
plot(density(today$depth + 0.1365, na.rm = T), col = alpha('grey25', .75), lwd = 2, 
     main = '', xlab = 'average depth (m)', xlim = c(0,1))
lines(density(hall$depth_m, na.rm = T), col = alpha('brown3', 0.75), lwd = 2)
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

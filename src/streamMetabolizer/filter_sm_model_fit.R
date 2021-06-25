## convert model results to Carbon and compile all methods/years
source("NHC_2019_metabolism/src/streamMetabolizer/inspect_model_fits.r")
# Compile model outputs ####
O2toC <- 12.0107/(2*15.999)

sm_fit <- readRDS("NHC_2019_metabolism/data/metabolism/compiled/raymond_met.rds")
dir_fit <- readRDS("NHC_2019_metabolism/data/metabolism/hall/hall_met_60min_2021_01.rds")

sm_preds <- sm_fit[[1]] %>%
  as.tibble() %>%
  select(-ends_with("Rhat")) %>%
  mutate(#across(starts_with("GPP" ), ~ . * O2toC),
         #across(starts_with("ER", ignore.case = FALSE), ~ . * O2toC),
         site = toupper(site), 
         year = year(date)) %>%
  filter(date <= as.Date("2020-03-20"))
dir_preds <- dir_fit[[1]] %>%
  as.tibble() %>%
  mutate(#GPP = GPP * O2toC,
         #ER = ER * O2toC,
         method = "direct_calculation")
dir_preds <- sm_preds %>% 
  select(date, site, temp.min) %>%
  right_join(dir_preds, by = c("site", "date")) 

hall_preds <- read_csv('hall_50yl/code/data/hall/hall_table_15.csv') %>%
  mutate(date = as.Date(date, format = "%m/%d/%Y"),
         site = case_when(site == "Concrete" ~ "CBP",
                          site == "Blackwood" ~ "BLK",
                          site == "Wood Bridge" ~ "WB")) %>%
  select(date, site, depth = depth_m, GPP = GPP_gO2m2d, ER = ER_gO2m2d) %>%
  group_by(site, date) %>%
  summarize_all(mean, na.rm = T) %>%
  ungroup()
hall_qt <- read_csv("hall_50yl/code/data/hall/hall_discharge_temp_daily.csv")
hall_preds <- hall_preds %>%
  left_join(hall_qt, by = "date") %>%
  mutate(level_m = stage_cm/100,
         method = "direct_calculation",
         # GPP = GPP * O2toC,
         # ER = -ER * O2toC,
         era = "then",
         year = year(date)) %>%
  select(-stage_cm) %>%
  rename(discharge = discharge_m3s, temp.water = water_temp_C)

all_preds <- dir_preds %>%
  mutate(era = "now") %>%
  bind_rows( hall_preds) %>%
  mutate(doy = format(date, "%j"), 
         month = paste(format(date, "%m"), format(date, "%b")),
         pr = GPP/ER) 

write_csv(all_preds, "NHC_2019_metabolism/data/metabolism/compiled/daily_preds_direct_calculation.csv")


sm_preds_sum <- sm_preds %>%
  mutate(doy = format(date, "%j")) %>%
  group_by(doy) %>% 
  summarize(GPP.upper = quantile(GPP, .975, na.rm = T),
            GPP.lower = quantile(GPP, .025, na.rm = T),
            ER.upper = quantile(ER, .975, na.rm = T),
            ER.lower = quantile(ER, .025, na.rm = T),
            GPP = median(GPP, na.rm = T),
            ER = median(ER, na.rm = T)) %>%
  ungroup() 
sm_preds_sum <- sm_preds %>%
  mutate(doy = format(date, "%j")) %>%
  group_by(year, doy) %>% 
  summarize(GPP.upper = quantile(GPP, .975, na.rm = T),
            GPP.lower = quantile(GPP, .025, na.rm = T),
            ER.upper = quantile(ER, .975, na.rm = T),
            ER.lower = quantile(ER, .025, na.rm = T),
            GPP = median(GPP, na.rm = T),
            ER = median(ER, na.rm = T)) %>%
  ungroup() %>%
  bind_rows(sm_preds_sum) %>%
  mutate(era = "now", 
         method = "uninformed_raymond")

dir_preds_sum <- dir_preds %>%
  mutate(doy = format(date, "%j")) %>%
  group_by(doy) %>% 
  summarize(GPP.upper = quantile(GPP, .975, na.rm = T),
            GPP.lower = quantile(GPP, .025, na.rm = T),
            ER.upper = quantile(ER, .975, na.rm = T),
            ER.lower = quantile(ER, .025, na.rm = T),
            GPP = median(GPP, na.rm = T),
            ER = median(ER, na.rm = T)) %>%
  ungroup() 
dir_preds_sum <- dir_preds %>%
  mutate(doy = format(date, "%j")) %>%
  group_by(year, doy) %>% 
  summarize(GPP.upper = quantile(GPP, .975, na.rm = T),
            GPP.lower = quantile(GPP, .025, na.rm = T),
            ER.upper = quantile(ER, .975, na.rm = T),
            ER.lower = quantile(ER, .025, na.rm = T),
            GPP = median(GPP, na.rm = T),
            ER = median(ER, na.rm = T)) %>%
  ungroup() %>%
  bind_rows(dir_preds_sum) %>%
  mutate(era = "now", 
         method = "direct_calculation")

hall_preds_sum <- hall_preds %>%
  mutate(doy = format(date, "%j")) %>%
  group_by(doy) %>% 
  summarize(GPP.upper = quantile(GPP, .975, na.rm = T),
            GPP.lower = quantile(GPP, .025, na.rm = T),
            ER.upper = quantile(ER, .975, na.rm = T),
            ER.lower = quantile(ER, .025, na.rm = T),
            GPP = median(GPP, na.rm = T),
            ER = median(ER, na.rm = T)) %>%
  ungroup() %>%
  mutate(era = "then", 
         method = "direct_calculation")
  

all_preds_sum <- bind_rows(sm_preds_sum, dir_preds_sum, hall_preds_sum)

all_preds_sum %>%
  filter(GPP <3) %>%
ggplot( aes(doy, GPP), color = "forestgreen") +
  geom_point() +
  geom_point(aes(y = ER), color = "sienna") +
  facet_wrap(era~method)



# calculate cumulative data for Hall ####
met <- hall_preds %>%
  summarize(gpp_mean = mean(GPP, na.rm = T),
            gpp_median = median(GPP, na.rm = T),
            gpp_max = max(GPP, na.rm = T),
            er_mean = -mean(ER, na.rm = T),
            er_median = -median(ER, na.rm = T),
            er_max = -min(ER, na.rm = T))
cum <- data.frame(date = seq(hall_preds$date[1], 
                             hall_preds$date[nrow(hall_preds)], 
                             by = "day")) %>%
  as_tibble() %>%
  left_join(hall_preds) %>%
  select(date, depth, GPP, ER) %>%
  mutate(across(-date, na.approx, na.rm = F)) %>%
  mutate(across(-date, cumsum, .names = "{col}_cum")) 

n <- nrow(cum)
l = (n-9)
weekly <- tibble(date = cum$date[1:l],
                 GPP_week = rep(NA_real_, l),
                 ER_week = rep(NA_real_, l))
for(i in 1:l){
  weekly$GPP_week[i] <- sum(cum$GPP[i:(i+9)]) 
  weekly$ER_week[i] <- sum(cum$ER[i:(i+9)]) 
}
met$gpp_max10d <- weekly$date[which.max(weekly$GPP_week)]
met$er_max10d <- weekly$date[which.max(weekly$ER_week)]
se <- sum(is.na(cum$GPP))
met$gpp_cum <- cum$GPP_cum[n-se]*365/(n-se)
se <- sum(is.na(cum$ER))
met$er_cum <- cum$ER_cum[n-se]*365/(n-se)
met$daterange <- as.character(paste(cum$date[1], "-", cum$date[nrow(cum)]))
met$total_days <- sum(!is.na(hall_preds$GPP))
met$pctcoverage <- sum(!is.na(hall_preds$GPP))/nrow(cum)
met$site = "all"

met_sum <- met

for(site in c("CBP", "WB", "BLK")){
  preds <- hall_preds %>%
    filter(site == !! site)
  met <- preds %>%
    summarize(gpp_mean = mean(GPP, na.rm = T),
              gpp_median = median(GPP, na.rm = T),
              gpp_max = max(GPP, na.rm = T),
              er_mean = -mean(ER, na.rm = T),
              er_median = -median(ER, na.rm = T),
              er_max = -min(ER, na.rm = T))
  cum <- data.frame(date = seq(preds$date[1], 
                             preds$date[nrow(preds)], 
                             by = "day")) %>%
    as_tibble() %>%
    left_join(preds) %>%
    select(date, depth, GPP, ER) %>%
    mutate(across(-date, na.approx, na.rm = F)) %>%
    mutate(across(-date, cumsum, .names = "{col}_cum"))  
  
  n <- nrow(cum)
  l = (n-9)
  weekly <- tibble(date = cum$date[1:l],
                   GPP_week = rep(NA_real_, l),
                   ER_week = rep(NA_real_, l))
  for(i in 1:l){
    weekly$GPP_week[i] <- sum(cum$GPP[i:(i+9)]) 
    weekly$ER_week[i] <- sum(cum$ER[i:(i+9)]) 
  }
  met$gpp_max10d <- weekly$date[which.max(weekly$GPP_week)]
  met$er_max10d <- weekly$date[which.max(weekly$ER_week)]
  se <- sum(is.na(cum$GPP))
  met$gpp_cum <- cum$GPP_cum[n-se]*365/(n-se)
  se <- sum(is.na(cum$ER))
  met$er_cum <- cum$ER_cum[n-se]*365/(n-se)
  met$daterange <- as.character(paste(cum$date[1], "-", cum$date[nrow(cum)]))
  met$total_days <- sum(!is.na(preds$GPP))
  met$pctcoverage <- sum(!is.na(preds$GPP))/nrow(cum)
  met$site = site
  
  met_sum <- bind_rows(met_sum, met)
}

# Group summary data####

sm_fit$summary <- sm_fit$summary %>% mutate(site = toupper(site))
summary <- bind_rows(sm_fit$summary, dir_fit$summary, met_sum)

write_csv(summary, "NHC_2019_metabolism/data/metabolism/compiled/metabolism_summary_table_2021_01.csv")

ggplot(summary, aes(gpp_median, color = method))+
  geom_bar()
png("figures/tmp/er_mintemp_by_month.png", width = 7, height = 5, 
    res = 300, units = "in")
all_preds %>%
  filter(GPP < 3) %>%
  ggplot(aes(temp.min, ER, color = era)) +
  geom_point() +
  theme_bw() +
  facet_wrap(.~month, scales = "free")
dev.off()

all_preds %>%
  filter(GPP < 3) %>%
  ggplot(aes(temp.water, GPP, color = era)) +
  geom_point() +
  theme_bw() +
  facet_wrap(.~month, scales = "free")
mmarize(GPP.upper = quantile(GPP, .975, na.rm = T),
            GPP.lower = quantile(GPP, .025, na.rm = T),
            ER.upper = quantile(ER, .975, na.rm = T),
            ER.lower = quantile(ER, .025, na.rm = T),
            GPP = median(GPP, na.rm = T),
            ER = median(ER, na.rm = T))




  mutate(doy = format(date, "%j")) %>%
  group_by(doy) %>%
  summarize(GPP.upper = quantile(GPP, .975, na.rm = T),
            GPP.lower = quantile(GPP, .025, na.rm = T),
            ER.upper = quantile(ER, .975, na.rm = T), 
            ER.lower = quantile(ER, .025, na.rm = T),
            GPP = median(GPP, na.rm = T),
              ER = median(ER, na.rm = T))
hall_met <- all_preds %>% 
  mutate(doy = format(date, "%j")) %>%
  group_by(doy) %>%
  summarize(GPP.upper = quantile(GPP, .975, na.rm = T),
            GPP.lower = quantile(GPP, .025, na.rm = T),
            ER.upper = quantile(ER, .975, na.rm = T), 
            ER.lower = quantile(ER, .025, na.rm = T),
            GPP = median(GPP, na.rm = T),
              ER = median(ER, na.rm = T))

# png("../figures/hall_met_comparison_nightreg.png", height = 6, width = 7, 
#     units = "in", res = 300)
#   m <- matrix(c(1,2,1,2,1, 
#                 2,1,3,1,3), nrow = 2)
#   layout(m)
#   plot_hall_metab(hall_met, ylim = c(-15,10), doy = T)
#   plot_kde_hall_metab(hall_met, lim = 4)
#   plot_k(all_preds)
#   mtext("stream metabolizer", outer = T, line = -3, cex = 1.2)
# dev.off()


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
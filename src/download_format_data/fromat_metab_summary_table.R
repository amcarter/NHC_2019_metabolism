# Format and compile summary metab data from NHC raymond runs.
# 2021 01 05

library(tidyverse)
library(lubridate)
library(zoo)

setwd("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism")

# load site and met data ####
site_dat <- read_csv("data/siteData/NHCsite_metadata.csv") %>%
  slice(c(1:5,7)) %>%
  select(site = sitecode, distance_m, width_mar_m, slope)
site_dat <- read_csv("data/rating_curves/calculated_channel_dimensions_maroct.csv")

O2toC = 12.0107/(2*15.999)
met <- readRDS("data/metabolism/compiled/raymond_met.rds")
sum <- met$summary %>%
  as_tibble() %>%
  mutate(days = round(total_days * pctcoverage, 0)) %>%
  select(site, year, days, total_days, pctcoverage, 
         gpp_mean, gpp_max, peaktime_gpp = gpp_max10d, 
         er_mean, er_max, peaktime_er = er_max10d,
         gpp_cum, er_cum) %>%
  mutate(across(starts_with("er"), ~ - . * O2toC),
         across(starts_with("gpp"), ~ . * O2toC),
         nep_cum = gpp_cum + er_cum)
m <- met$preds %>%
  as_tibble() %>%
  select(GPP, ER) %>%
  summarize(gpp_mean = mean(GPP, na.rm = T),
            gpp_max = max(GPP, na.rm = T), 
            er_mean = mean(ER, na.rm = T), 
            er_max = min(ER, na.rm = T)) %>%
  mutate(across(.fns = ~ .* O2toC))

m$days <- sum(sum$days)
m$total_days <- sum(sum$total_days)
m$site <- "all"
sum <- bind_rows(sum, m)
sum$method <- "inverse modeling"
  
# compile hall method 2019 ####

hmet <- readRDS("data/metabolism/hall/hall_met_60min.rds")
hsum <- hmet$summary %>%
  as_tibble() %>%
  mutate(days = round(total_days * pctcoverage, 0)) %>%
  select(site, year, days, total_days, pctcoverage, 
         gpp_mean, gpp_max, peaktime_gpp = gpp_max10d, 
         er_mean, er_max, peaktime_er = er_max10d,
         gpp_cum, er_cum) %>%
  mutate(across(starts_with("er"), ~ - . * O2toC),
         across(starts_with("gpp"), ~ . * O2toC),
         nep_cum = gpp_cum + er_cum)

m <- hmet$preds %>%
  as_tibble() %>%
  select(GPP, ER) %>%
  summarize(gpp_mean = mean(GPP, na.rm = T),
            gpp_max = max(GPP, na.rm = T), 
            er_mean = mean(ER, na.rm = T), 
            er_max = min(ER, na.rm = T)) %>%
  mutate(across(.fns = ~ .* O2toC))

m$days <- sum(hsum$days)
m$total_days <- sum(hsum$total_days)
m$site <- "all"

cum <- hmet$preds %>%
  as_tibble() %>%
  group_by(doy = as.numeric(format(date, "%j"))) %>%
  select(doy, GPP, ER) %>%
  summarize_all(mean, na.rm = T) %>%
  ungroup() %>%
  mutate(across(-doy, ~ . * O2toC)) %>%
  full_join(data.frame(doy = 1:365)) %>%
  arrange(doy) %>%
  mutate(across(-doy, cumsum, .names = "{col}_cum"))

mms <- cum %>%
  left_join(data.frame(date = seq(as.Date("2019-01-01"), 
                                  as.Date("2019-12-31"), by = "day"),
                       doy = 1:365), by = "doy") %>%
  group_by(month = substr(date, 6, 7)) %>%
  summarize(er = mean(ER, na.rm = T),
            gpp = mean(GPP, na.rm = T)) %>%
  ungroup() 
  
m$peak_gpp <- month(as.Date(
  paste0("2019-", mms$month[which.max(mms$gpp)], "-01")), label = T, abbr = T)
m$peak_er = month(as.Date(
  paste0("2019-", mms$month[which.min(mms$er)],"-01")), label = T, abbr = T)
m$gpp_cum <- cum$GPP_cum[365]
m$er_cum <- cum$ER_cum[365]
m$nep_cum <- m$gpp_cum + m$er_cum
m$year <- 2017

hsum$method <- "direct calculation"
# compile and format Hall 1972 data ####
hall <- read_csv("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/hall_50yl/data/hall/hall_table_15.csv") %>%
  mutate(date = as.Date(date, format = "%m/%d/%Y"),
         site = case_when(site == "Concrete" ~ "CBP",
                          site == "Blackwood" ~ "BLK", 
                          site == "Wood Bridge" ~ "WB")) %>%
  group_by(site, date) %>%
  summarize(gpp = mean(GPP_gO2m2d, na.rm = T) * O2toC,
            er = -mean(ER_gO2m2d, na.rm = T) * O2toC) %>%
  ungroup() %>%
  arrange(date)

  s = "all"
  y = "all"
sumhall <- function(hall, s, y){
  mrow <- hall %>%
    summarize(gpp_mean = mean(gpp, na.rm = T),
              gpp_max = max(gpp, na.rm = T),
              er_mean = mean(er, na.rm = T),
              er_max = max(er, na.rm = T)) %>%
    mutate(site = s,
           year = y,
           days = min(length(unique(hall$date)), sum(!is.na(hall$gpp))),
           total_days = as.numeric(hall$date[nrow(hall)] - hall$date[1]) +1,
           pctcoverage = days/total_days)
  
  
  cum <- data.frame(date = seq(hall$date[1], hall$date[nrow(hall)], by = 'days')) %>%
    left_join(hall, by = "date") %>% 
    select(-site) %>%
    group_by(date) %>%
    summarize_all(mean, na.rm = T) %>%
    mutate(across(-date, na.approx),
           across(-date, cumsum, .names = "{col}_cum")) %>%
    as_tibble()
  n <- nrow(cum)
  l = (n-9)
  weekly <- tibble(date = cum$date[1:l],
                   GPP_week = rep(NA_real_, l),
                   ER_week = rep(NA_real_, l))
  for(i in 1:l){
    weekly$GPP_week[i] <- sum(cum$gpp[i:(i+9)]) 
    weekly$ER_week[i] <- sum(cum$er[i:(i+9)]) 
  }
  mrow$peaktime_gpp <- weekly$date[which.max(weekly$GPP_week)]
  mrow$peaktime_er <- weekly$date[which.min(weekly$ER_week)]
  mrow$gpp_cum <- cum$gpp_cum[n]*365/(n)
  mrow$er_cum <- cum$er_cum[n]*365/(n)
  mrow$nep_cum <- mrow$gpp_cum + mrow$er_cum

  return(mrow)
}
  
hall_dat <- sumhall(hall, "all", NA_real_)

hall_c <- hall %>%
  mutate(doy = as.numeric(format(date, "%j"))) %>%
  group_by(doy) %>%
  select(-date, -site) %>%
  summarize_all(mean, na.rm = T) %>%
  ungroup() %>%
  arrange(doy)
hall_c <- data.frame(date = seq(as.Date("1969-01-01"), as.Date("1970-01-01"), by = "day"),
           doy = c(1:365, 1)) %>%
  left_join(hall_c, by = "doy") %>%
  select(-doy) %>%
  mutate(site = "all")

mrow <- sumhall(hall_c, "all_doy", NA_real_)
hall_dat <- bind_rows(hall_dat, mrow)

for(s in unique(hall$site)) {
  halls <- hall %>%
    filter(site == s) %>%
    arrange(date)
  
  mrow <- sumhall(halls, s, NA_real_)
  hall_dat <- bind_rows(hall_dat, mrow)
}

for(i in 1:2){
  w1 <- max(which(hall$date <= hall$date[1] + (i-1) *365))
  w2 <- min(which(hall$date >= hall$date[1] + i*365))
  halls <- hall %>%
    filter(site == "CBP",
           date >= hall$date[w1],
           date <= hall$date[w2]) %>%
    arrange(date)
  mrow <- sumhall(halls, "CBP", 1967+i)
  hall_dat <- bind_rows(hall_dat, mrow)
}

hall_dat$method <- "direct calculation"
# Summarize metabolism by month ####
preds <- met$preds %>%
  as_tibble() %>%
  select(-starts_with("K600"), -ends_with("Rhat"),
         -errors,-good_flow, -method) %>%
  mutate(month = substr(date, 6, 7),
         across(starts_with(c("GPP", "ER")), ~ . * O2toC)) 
# month_preds <- preds %>%
#   filter(year == 2019) %>%
#   select(-date, -year) %>%
#   mutate(site = case_when(site == "nhc" ~ "NHC",
#                    site == "cbp" ~ "CBP",
#                    site == "pm" ~ "PM", 
#                    site == "wb" ~ "WB",
#                    site == "wbp" ~ "WBP",
#                    site == "unhc" ~ "UNHC")) %>%
#   group_by(site, month) %>%
#   summarize(across(.fns = list(mean = ~mean(.,na.rm = T), 
#                                sd = ~sd(.,na.rm = T)),
#                    .names = "{col}_{fn}")) %>%
#   ungroup() 
month_preds <- preds %>%
  select(-date) %>%
  mutate(site = case_when(site == "nhc" ~ "NHC",
                   site == "cbp" ~ "CBP",
                   site == "pm" ~ "PM",
                   site == "wb" ~ "WB",
                   site == "wbp" ~ "WBP",
                   site == "unhc" ~ "UNHC")) %>%
  group_by(site, year, month) %>%
  summarize(across(.fns = list(mean = ~mean(.,na.rm = T),
                               sd = ~sd(.,na.rm = T)),
                   .names = "{col}_{fn}")) %>%
  ungroup()

sum_month <- month_preds %>%
  left_join(site_dat, by = "site")
write_csv(sum_month, 
          "data/metabolism/compiled/raymond_summarized_monthly.csv")
# month_preds <- preds %>%
#   select(-date, -site, -year) %>% 
#   group_by(month) %>%
#   summarize(across(.fns = list(mean = ~mean(.,na.rm = T), 
#                                sd = ~sd(.,na.rm = T)),
#                    .names = "{col}_{fn}")) %>%
#   ungroup() %>%
#   mutate(site = "all") %>%
#   bind_rows(month_preds)
#   


# m_preds <- month_preds %>%
#   group_by(site, year) %>%
#   summarize(peak_gpp = month(as.Date(
#               paste0("2019-",month[which.max(GPP)],"-01")), 
#               label = T, abbr = T),
#             peak_er = month(as.Date(
#               paste0("2019-", month[which.min(ER)], "-01")),
#               label = T, abbr = T))

hpreds <- hmet$preds %>%
  as_tibble() %>%
  select(-K600, -good_flow, -depth) %>%
  mutate(month = substr(date, 6, 7),
         across(starts_with(c("GPP","ER")), ~ . * O2toC))

hmonth_preds <- hpreds %>%
  select(-date) %>%
  group_by(site, year, month) %>%
  summarize(across(.fns = list(mean = ~mean(.,na.rm = T), 
                               sd = ~sd(.,na.rm = T)),
                   .names = "{col}_{fn}")) %>%
  ungroup()

hmonth_preds <- hpreds %>%
  select(-date, -year, -site) %>%
  group_by(month) %>%
  summarize(across(.fns = list(mean = ~mean(.,na.rm = T), 
                               sd = ~sd(.,na.rm = T)),
                   .names = "{col}_{fn}")) %>%
  ungroup() %>%
  mutate(site = "all") %>%
  bind_rows(hmonth_preds)
hmonth_preds <- hpreds %>%
  filter(site == "CBP") %>%
  select(-date, -year, -site) %>%
  group_by(month) %>%
  summarize(across(.fns = list(mean = ~mean(.,na.rm = T), 
                               sd = ~sd(.,na.rm = T)),
                   .names = "{col}_{fn}")) %>%
  ungroup() %>%
  mutate(site = "CBP") %>%
  bind_rows(hmonth_preds)


# hm_preds <- hmonth_preds %>%
#   group_by(site, year) %>%
#   summarize(peak_gpp = month(as.Date(
#               paste0("2019-", month[which.max(GPP)], "-01")),
#               label = T, abbr = T),
#             peak_er = month(as.Date(
#               paste0("2019-", month[which.min(ER)],"-01")),
#               label = T, abbr = T))

# hallQ <- read_csv("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/hall_50yl/data/hall/hall_figure26_digitized_dailystage.csv")
# hallQT <- read_csv("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/hall_50yl/data/hall/hall_figure27_digitized_mean_daily_temp.csv") %>%
#   full_join(hallQ, by = "date") %>%
#   group_by(date) %>%
#   summarize(across(-notes, mean, na.rm = T)) %>%
#   ungroup()
#   
# hallQT <- data.frame(date = seq(min(hallQT$date), 
#                                 max(hallQT$date), by = "day")) %>%
#   left_join(hallQT, by = "date") %>%
#   mutate(across( -date, na.approx, na.rm = F),
#          notes = case_when(stage_cm < min(hall_rc$stage_cm) ~ "below RC",
#                            stage_cm > max(hall_rc$stage_cm) ~ "above RC"))
# write_csv(hallQT, "C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/hall_50yl/data/hall/hall_discharge_temp_daily.csv")
hallQT <- read_csv("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/hall_50yl/data/hall/hall_discharge_temp_daily.csv")

hall_month <- hall %>% 
  mutate(month = substr(date, 6, 7),
         year = as.numeric(substr(date, 1, 4))) %>%
  left_join(hallQT[,c(1,2,4)], by = "date") %>%
  select(-date) %>%
  group_by(site, year, month) %>%
  summarize(across(.fns = list(mean = ~mean(.,na.rm = T), 
                               sd = ~sd(.,na.rm = T)),
                   .names = "{col}_{fn}")) %>%  
  ungroup() 
hall_month <- hall %>% 
  mutate(month = substr(date, 6, 7)) %>%
  left_join(hallQT[,c(1,2,4)], by = "date") %>%
  select(-date, -site) %>%
  group_by(month) %>%
  summarize(across(.fns = list(mean = ~mean(.,na.rm = T), 
                               sd = ~sd(.,na.rm = T)),
                   .names = "{col}_{fn}")) %>%  
  ungroup() %>%
  mutate(site = "all") %>%
  bind_rows(hall_month)
hall_month <- hall %>% 
  filter(site == "CBP") %>%
  mutate(month = substr(date, 6, 7)) %>%
  left_join(hallQT[,c(1,2,4)], by = "date") %>%
  select(-date, -site) %>%
  group_by(month) %>%
  summarize(across(.fns = list(mean = ~mean(.,na.rm = T), 
                               sd = ~sd(.,na.rm = T)),
                   .names = "{col}_{fn}")) %>%  
  ungroup() %>%
  mutate(site = "CBP") %>%
  bind_rows(hall_month)

  
# compile all metabolism summaries ####
sum <- left_join(sum, m_preds, by = c("site", "year")) %>%
  select(-peaktime_er, -peaktime_gpp)


hsum <- left_join(hsum, hm_preds, by = c("site", "year")) %>%
  select(-peaktime_er, -peaktime_gpp)
hsum <- bind_rows(hsum, m)

hall_dat <- hall_dat %>%
  mutate(peak_gpp = month(peaktime_gpp, label = T, abbr = T),
         peak_er = month(peaktime_er, label = T, abbr = T)) %>%
  select(-peaktime_er, -peaktime_gpp)
  
compiled <- sum %>%
  mutate(site = as.character(site),
         site = case_when(site == "nhc" ~ "NHC",
                          site == "cbp" ~ "CBP",
                          site == "pm" ~ "PM", 
                          site == "wb" ~ "WB",
                          site == "wbp" ~ "WBP",
                          site == "unhc" ~ "UNHC")) %>%
  bind_rows(hsum, hall_dat) %>%
  left_join(site_dat)

write_csv(compiled, "data/metabolism/compiled/metabolism_summary_table_gC.csv")
comp <- compiled %>%
  mutate(group = case_when(year >2000 & method == "inverse modeling" ~ 
                             "raymond",
                           year > 2000 & method == "direct calculation" ~
                             "hall_method",
                           TRUE ~ "hall_data"),
         distance_m = ifelse(site == "BLK", 10000, distance_m)) %>%
  filter(!is.na(distance_m)) %>%
  filter(!(year %in% 2017:2018))

png("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/hall_50yl/figures/longitudinal_met_all_methods.png",
    width = 6, height = 4, res = 300, units = "in")
  ggplot(comp, aes(distance_m, gpp_cum))+
    geom_line() +
    geom_line(aes(y = er_cum)) +
    geom_line(aes(y = gpp_cum + er_cum), lty = 2) +
    facet_grid(group~.) +
    labs(title = "cumulative metabolism along 10 km",
           y = "GPP, ER, and NEP (gC/m2/y)")
  
dev.off()

# met vs geomorph plots ####
png("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/hall_50yl/figures/Q_by_month.png",
    width = 6, height = 4, res = 300, units = "in")
  ggplot(sum_month, aes(x = as.numeric(month), y = log(discharge.daily_mean), col = as.factor(year))) +
    geom_line() + 
    facet_wrap(.~site)
dev.off()

png("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/hall_50yl/figures/Q_by_month.png",
    width = 6, height = 4, res = 300, units = "in")
  sum_month %>% filter(month %in% )
  ggplot(sum_month, aes(x = as.numeric(month), y = GPP_mean, col = width_mar_m)) +
    geom_point(size = 2) 
    , col = as.factor(year))) +
    facet_wrap(.~site)
dev.off()


# plot metabolism by month ####

month_preds <- month_preds %>%
  mutate(month = as.numeric(month))
hall_month <- hall_month %>%
  mutate(month = as.numeric(month))

allm <- month_preds %>% 
  filter(site != "all")%>%
  group_by(site, month) %>%
  summarize(across(-year, mean, na.rm = T)) %>%
  ungroup()
all <- allm %>% 
  group_by(month) %>%
  summarize(across(-site, mean, na.rm = T)) %>%
  ungroup() %>%
  mutate(site = "all")

allm <- month_preds
png("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/hall_50yl/figures/2019monthly_avg_met_all_sites_raymond.png",
    width = 7.5, height = 5, res = 300, units = "in")

ggplot(allm, aes(x = as.numeric(month), y = GPP_mean)) +
  geom_line() +
  geom_ribbon(aes(ymin = GPP_mean - GPP_sd, 
                  ymax = GPP_mean + GPP_sd),
              fill = alpha("forestgreen", .3), col = NA) +
  geom_line(aes(y = ER_mean)) +
  geom_ribbon(aes(ymin = ER_mean - ER_sd, 
                  ymax = ER_mean + ER_sd),
              fill = alpha("sienna", .3), col = NA) +
  facet_wrap(.~ site) +
  ylim(-4,1.2) +
  geom_hline(yintercept = 0, col = "grey") +
  labs(title = "2019 monthly average metabolism from StreamMetabolizer",
       x = "month", y = "gC/m2/d")
dev.off()
png("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/hall_50yl/figures/2019monthly_PR_all_sites_raymond.png",
    width = 7.5, height = 5, res = 300, units = "in")

ggplot(allm, aes(x = as.numeric(month), y = -GPP_mean/ER_mean)) +
  geom_line() +
  facet_wrap(.~ site, scales = "free_y") +
  geom_hline(yintercept = 1, col = "grey30", lty = 2) +
  labs(title = "2019 Monthly Productivity:Respiration from StreamMetabolizer",
       x = "month", y = "GPP/ER")
dev.off()

ggplot(hall_month, aes(month, -GPP/ER, col = year)) +
  geom_line() +
  facet_wrap(.~site)+ 
  ylim(0,1.5) +
  theme_minimal()

png("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/hall_50yl/figures/ER_v_QT_raymond_met.png",
    width = 6, height = 4, res = 300, units = "in")
  month_preds %>%
    filter(site != "all") %>%
  ggplot(aes(log(discharge.daily_mean), ER_mean, 
                          col = (temp.water_mean))) +
    geom_point(size = 1.5) +
    facet_wrap(.~site, scales = "free") +
    geom_errorbar(aes(ymin = ER.lower_mean, ymax = ER.upper_mean))
dev.off()

png("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/hall_50yl/figures/annual_NEP_nowthen.png",
    width = 7.5, height = 6, res = 300, units = "in")
  par(mar = c(0,4,4,3), mfrow = c(2,1))
  plot(1, type = 'n', xlim = c(1,12), ylim = c(0,3), ylab = "gC/m2/d", 
       xlab = "month", xaxt = "n", 
       main = "Respiration in excess of GPP")    
  polygon(c(1:12,12:1), 
            c(-hmonth_preds$ER_mean[hmonth_preds$site == "CBP"],
              rev(hmonth_preds$GPP_mean[hmonth_preds$site == "CBP"])), 
          col = alpha("grey30", .3), border = "grey30", lwd = 2)
  polygon(c(1:8,10:12, 12:10,8:1),
            c(-hall_month$er_mean[hall_month$site == "CBP"],
              rev(hall_month$gpp_mean[hall_month$site == "CBP"])),
          col = alpha("steelblue", .3), border = "steelblue", lwd = 2)
  arrows(x0 = c(1:12) + .01,
         y0 = -hmonth_preds$ER_mean[hmonth_preds$site == "CBP"] +
           hmonth_preds$ER_sd[hmonth_preds$site == "CBP"],
         x1 = c(1:12)+.01,
         y1 = -hmonth_preds$ER_mean[hmonth_preds$site == "CBP"] -
           hmonth_preds$ER_sd[hmonth_preds$site == "CBP"],
         length = 0, col = "grey30", lwd = 2) 
  arrows(x0 = c(1:12) + .03,
         y0 = hmonth_preds$GPP_mean[hmonth_preds$site == "CBP"] +
           hmonth_preds$GPP_sd[hmonth_preds$site == "CBP"],
         x1 = c(1:12)+.03,
         y1 = hmonth_preds$GPP_mean[hmonth_preds$site == "CBP"] -
           hmonth_preds$GPP_sd[hmonth_preds$site == "CBP"],
         length = 0, col = "grey30", lwd = 2, lty = 2) 
  arrows(x0 = c(1:8, 10:12) - .03,
         y0 = hall_month$gpp_mean[hall_month$site == "CBP"] +
           hall_month$gpp_sd[hall_month$site == "CBP"],
         x1 = c(1:8, 10:12) -.03,
         y1 = hall_month$gpp_mean[hall_month$site == "CBP"] -
           hall_month$gpp_sd[hall_month$site == "CBP"],
         length = 0, col = "steelblue", lwd = 2, lty = 2) 
  arrows(x0 = c(1:8, 10:12) - .01,
         y0 = -hall_month$er_mean[hall_month$site == "CBP"] +
           hall_month$er_sd[hall_month$site == "CBP"],
         x1 = c(1:8, 10:12) -.01,
         y1 = -hall_month$er_mean[hall_month$site == "CBP"] -
           hall_month$er_sd[hall_month$site == "CBP"],
         length = 0, col = "steelblue", lwd = 2) 
  legend("topright",
         legend = c("Now", "Then"),
         fill = c(alpha("grey30", .3), alpha("steelblue", .3)),
         bty = "n", cex = 1.3, ncol = 2)
  mtext("CBP site",line =  -1.4, adj = 0.05, cex = 1.2)  
  
  par(mar = c(4,4,0,3))
  plot(1, type = 'n', xlim = c(1,12), ylim = c(0,3), ylab = "gC/m2/d", 
       xlab = "month", xaxt = "n")
  axis(1, at = seq(1, 12, by = 1), labels = month.abb[seq(1,12, by = 1)])
  polygon(c(1:12,12:1), 
            c(-hmonth_preds$ER_mean[hmonth_preds$site == "all"],
              rev(hmonth_preds$GPP_mean[hmonth_preds$site == "all"])), 
          col = alpha("grey30", .3), border = "grey30", lwd = 2)
  polygon(c(1:8,10:12, 12:10,8:1),
            c(-hall_month$er_mean[hall_month$site == "all"],
              rev(hall_month$gpp_mean[hall_month$site == "all"])),
          col = alpha("steelblue", .3), border = "steelblue", lwd = 2)
  arrows(x0 = c(1:12) + .01,
         y0 = -hmonth_preds$ER_mean[hmonth_preds$site == "all"] +
           hmonth_preds$ER_sd[hmonth_preds$site == "all"],
         x1 = c(1:12)+.01,
         y1 = -hmonth_preds$ER_mean[hmonth_preds$site == "all"] -
           hmonth_preds$ER_sd[hmonth_preds$site == "all"],
         length = 0, col = "grey30", lwd = 2) 
  arrows(x0 = c(1:12) + .03,
         y0 = hmonth_preds$GPP_mean[hmonth_preds$site == "all"] +
           hmonth_preds$GPP_sd[hmonth_preds$site == "all"],
         x1 = c(1:12)+.03,
         y1 = hmonth_preds$GPP_mean[hmonth_preds$site == "all"] -
           hmonth_preds$GPP_sd[hmonth_preds$site == "all"],
         length = 0, col = "grey30", lwd = 2, lty = 2) 
  arrows(x0 = c(1:8, 10:12) - .03,
         y0 = hall_month$gpp_mean[hall_month$site == "all"] +
           hall_month$gpp_sd[hall_month$site == "all"],
         x1 = c(1:8, 10:12) -.03,
         y1 = hall_month$gpp_mean[hall_month$site == "all"] -
           hall_month$gpp_sd[hall_month$site == "all"],
         length = 0, col = "steelblue", lwd = 2, lty = 2) 
  arrows(x0 = c(1:8, 10:12) - .01,
         y0 = -hall_month$er_mean[hall_month$site == "all"] +
           hall_month$er_sd[hall_month$site == "all"],
         x1 = c(1:8, 10:12) -.01,
         y1 = -hall_month$er_mean[hall_month$site == "all"] -
           hall_month$er_sd[hall_month$site == "all"],
         length = 0, col = "steelblue", lwd = 2) 
  mtext("All sites",line =  -1.4, adj = 0.05, cex = 1.2)  
dev.off()

png("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/hall_50yl/figures/PR_cbp_nowthen.png",
    width = 6, height = 4, res = 300, units = "in")

plot(c(1:8, 10:12), -hall_month$gpp_mean[hall_month$site =="CBP"]/
       hall_month$er_mean[hall_month$site =="CBP"], type = "l", 
     col = "steelblue", lwd = 2, ylim = c(0.1, 2.1), xaxt = "n",
     xlab = "month", ylab = "P:R ratio", 
     main = "Productivity to Respiration at CBP")
axis(1, at = seq(1, 12, by = 1), labels = month.abb[seq(1,12, by = 1)])

# arrows(x0 = c(1:8, 10:12, 12:10, 8:1),
#        y0 = -hall_month$gpp_mean[hall_month$site =="CBP"]/
#   hall_month$er_mean[hall_month$site =="CBP"] +
#     sqrt(hall_month$gpp_sd[hall_month$site =="CBP"]^2 +
#     hall_month$er_sd[hall_month$site =="CBP"]^2),
#   c(1:8, 10:12, 12:10, 8:1),
#   -hall_month$gpp_mean[hall_month$site =="CBP"]/
#     hall_month$er_mean[hall_month$site =="CBP"] -
#     sqrt(hall_month$gpp_sd[hall_month$site =="CBP"]^2 +
#            hall_month$er_sd[hall_month$site =="CBP"]^2),
#   col = "steelblue", length = 0)
# lines(c(1:8, 10:12), -hall_month$gpp_mean[hall_month$site =="CBP"]/
#   hall_month$er_mean[hall_month$site =="CBP"],
#   col = "steelblue", lwd = 2, lty = 2)
# lines(1:12, -hmonth_preds$GPP_mean[hmonth_preds$site == "all"]/
#         hmonth_preds$ER_mean[hmonth_preds$site == "all"], 
#       lwd = 2, col = "grey30")
lines(1:12, -hmonth_preds$GPP_mean[hmonth_preds$site == "CBP"]/
        hmonth_preds$ER_mean[hmonth_preds$site == "CBP"], 
      lwd = 2, col = "grey30")
# arrows(x0 = c(1:8, 10:12, 12:10, 8:1),
#        y0 = -hmonth_preds$GPP_mean[hmonth_preds$site =="CBP"]/
#          hmonth_preds$ER_mean[hmonth_preds$site =="CBP"] +
#          sqrt(hmonth_preds$GPP_sd[hmonth_preds$site =="CBP"]^2 +
#                 hmonth_preds$ER_sd[hmonth_preds$site =="CBP"]^2),
#        c(1:8, 10:12, 12:10, 8:1),
#        -hmonth_preds$GPP_mean[hmonth_preds$site =="CBP"]/
#          hmonth_preds$ER_mean[hmonth_preds$site =="CBP"] -
#          sqrt(hmonth_preds$GPP_sd[hmonth_preds$site =="CBP"]^2 +
#                 hmonth_preds$ER_sd[hmonth_preds$site =="CBP"]^2),
#        col = "grey30", length = 0)

legend("topright",
       legend=c("Now","Then"),
       col = c("grey30", "steelblue"),
       lwd = 2, cex = 1.2, bty = "n")
    
dev.off()   
      


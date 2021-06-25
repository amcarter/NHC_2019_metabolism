# Format and compile summary metab data from NHC raymond runs.
# 2021 01 05

library(tidyverse)
library(lubridate)
library(zoo)

# load site and met data ####
site_dat <- read_csv("NHC_2019_metabolism/data/siteData/NHCsite_metadata.csv") %>%
  slice(c(1:5,7)) %>%
  select(site = sitecode, distance_m, width_mar_m, slope)
site_dat <- read_csv("NHC_2019_metabolism/data/rating_curves/calculated_channel_dimensions_maroct.csv") %>%
  right_join(site_dat, by = "distance_m")

O2toC = 12.0107/(2*15.999)
sm_met <- read_csv("NHC_2019_metabolism/data/metabolism/compiled/daily_preds_stream_metabolizer.csv")
d_met <- read_csv("NHC_2019_metabolism/data/metabolism/compiled/daily_preds_direct_calculation.csv")
all_sum <- read_csv("NHC_2019_metabolism/data/metabolism/compiled/metabolism_summary_table_2021_01.csv")

# compile and format Hall 1972 data ####
hall <- read_csv("hall_50yl/code/data/hall/hall_table_15.csv") %>%
  mutate(date = as.Date(date, format = "%m/%d/%Y"),
         site = case_when(site == "Concrete" ~ "CBP",
                          site == "Blackwood" ~ "BLK", 
                          site == "Wood Bridge" ~ "WB")) %>%
  group_by(site, date) %>%
  summarize(gpp = mean(GPP_gO2m2d, na.rm = T),
            er = -mean(ER_gO2m2d, na.rm = T)) %>%
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
write_csv(hall_dat, "NHC_2019_metabolism/data/metabolism/compiled/metabolism_summary_table_then_2021_01_v2.csv")

# Summarize metabolism by month ####
d_met <- d_met %>%
  mutate(year = case_when(era == "then" ~ 1969,
                          TRUE ~ year)) %>%
  filter(GPP < 3)

monthly <- sm_met %>%
  rename(discharge = discharge.daily) %>%
  bind_rows(d_met) %>%
  select(-starts_with("K600"), -ends_with("Rhat"),
         -errors,-good_flow, -notes, -level_m, -doy) %>%
  mutate(month = as.numeric(format(date, "%m"))) %>%
  select(-date) %>%
  group_by(site, year, month, era, method) %>%
  summarize(across(.fns = list(mean = ~mean(.,na.rm = T),
                               sd = ~sd(.,na.rm = T)),
                   .names = "{col}_{fn}")) %>%
  ungroup()

cbp_monthly <- monthly %>%
  filter(site == "CBP")
geo_month <- left_join(monthly, site_dat, by = "site")

# png("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/hall_50yl/figures/longitudinal_met_all_methods.png",
#     width = 6, height = 4, res = 300, units = "in")
#   
#   all_sum %>% 
#     filter(is.na(method)) %>%
#   ggplot( aes(site, gpp_cum, color = year))+
#     geom_bar() 
#     geom_bar(aes(y = er_cum)) +
#     geom_bar(aes(y = er_cum + gpp_cum)) +
#     # facet_grid(era~.) +
#     labs(title = "cumulative metabolism along 10 km",
#            y = "GPP, ER, and NEP (gC/m2/y)")
#   
# dev.off()

# met vs geomorph plots ####
# png("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/hall_50yl/figures/Q_by_month.png",
#     width = 6, height = 4, res = 300, units = "in")
  geo_month %>%
    filter(method == "direct_calculation",
           era == "now") %>%
  ggplot(aes(x = month, y = GPP_mean, col = width_m)) +
    geom_point(size = 3)
    
png("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/hall_50yl/figures/Q_by_month.png",
    width = 6, height = 4, res = 300, units = "in")
  
ggplot(geo_month, aes(x = month, y = GPP_mean, col = width_m)) +
    geom_point(size = 2) 
    
dev.off()


# plot metabolism by month ####
png("figures/2019monthly_avg_pr_all_sites_directcalc.png",
    width = 7.5, height = 5, res = 300, units = "in")
monthly %>%
  filter(era == "now",
         year == 2019,
         method == "direct_calculation") %>%
ggplot(aes(x = month, y = -pr_mean)) +
  geom_line() +
  facet_wrap(.~ site) +
  geom_hline(yintercept = 1, col = "grey") +
  labs(title = "2019 monthly average P/R from direct calculation",
       x = "month", y = "log(P/R)")
dev.off()

png("figures/2019monthly_avg_met_all_sites_directcalc.png",
    width = 7.5, height = 5, res = 300, units = "in")
monthly %>%
  filter(era == "now",
         year == 2019,
         method == "direct_calculation") %>%
ggplot(aes(x = month, y = GPP_mean)) +
  geom_line() +
  geom_ribbon(aes(ymin = GPP_mean - GPP_sd, 
                  ymax = GPP_mean + GPP_sd),
              fill = alpha("forestgreen", .3), col = NA) +
  geom_line(aes(y = ER_mean)) +
  geom_ribbon(aes(ymin = ER_mean - ER_sd, 
                  ymax = ER_mean + ER_sd),
              fill = alpha("sienna", .3), col = NA) +
  facet_wrap(.~ site) +
  # ylim(-4,1.2) +
  geom_hline(yintercept = 0, col = "grey") +
  labs(title = "2019 monthly average metabolism from direct calculation",
       x = "month", y = "gC/m2/d")
dev.off()

png("figures/2019monthly_avg_met_all_years_directcalc.png",
    width = 7.5, height = 5, res = 300, units = "in")
monthly %>%
  filter(era == "now",
         site %in% c("NHC", "UNHC"),
         method == "direct_calculation") %>%
ggplot(aes(x = month, y = GPP_mean)) +
  geom_line() +
  geom_ribbon(aes(ymin = GPP_mean - GPP_sd, 
                  ymax = GPP_mean + GPP_sd),
              fill = alpha("forestgreen", .3), col = NA) +
  geom_line(aes(y = ER_mean)) +
  geom_ribbon(aes(ymin = ER_mean - ER_sd, 
                  ymax = ER_mean + ER_sd),
              fill = alpha("sienna", .3), col = NA) +
  facet_grid(site ~ year) +
  geom_hline(yintercept = 0, col = "grey") +
  labs(title = "Monthly average metabolism from direct calculation",
       x = "month", y = "gC/m2/d")
dev.off()

png("figures/2019monthly_avg_met_all_sites_SM.png",
    width = 7.5, height = 5, res = 300, units = "in")
monthly %>%
  filter(year == 2019,
         method == "uninformed_raymond") %>%
ggplot(aes(x = month, y = GPP_mean)) +
  geom_line() +
  geom_ribbon(aes(ymin = GPP_mean - GPP_sd, 
                  ymax = GPP_mean + GPP_sd),
              fill = alpha("forestgreen", .3), col = NA) +
  geom_line(aes(y = ER_mean)) +
  geom_ribbon(aes(ymin = ER_mean - ER_sd, 
                  ymax = ER_mean + ER_sd),
              fill = alpha("sienna", .3), col = NA) +
  facet_wrap(.~ site) +
  # ylim(-4,1.2) +
  geom_hline(yintercept = 0, col = "grey") +
  labs(title = "2019 monthly average metabolism from StreamMetabolizer",
       x = "month", y = "gC/m2/d")
dev.off()

png("figures/2019monthly_avg_met_all_years_SM.png",
    width = 7.5, height = 5, res = 300, units = "in")
monthly %>%
  filter(site %in% c("NHC", "UNHC"),
         method == "uninformed_raymond") %>%
ggplot(aes(x = month, y = GPP_mean)) +
  geom_line() +
  geom_ribbon(aes(ymin = GPP_mean - GPP_sd, 
                  ymax = GPP_mean + GPP_sd),
              fill = alpha("forestgreen", .3), col = NA) +
  geom_line(aes(y = ER_mean)) +
  geom_ribbon(aes(ymin = ER_mean - ER_sd, 
                  ymax = ER_mean + ER_sd),
              fill = alpha("sienna", .3), col = NA) +
  facet_grid(site ~ year) +
  geom_hline(yintercept = 0, col = "grey") +
  labs(title = "Monthly average metabolism from StreamMetabolizer",
       x = "month", y = "gC/m2/d")
dev.off()


now_19 <- monthly %>%
  filter(method == "direct_calculation",
         era == "now",
         year == 2019) 
now_y <- monthly %>%
  filter(method == "direct_calculation",
         era == "now",
         site %in% c("NHC", "UNHC")) 
then <- monthly %>% 
  filter(method == "direct_calculation",
         era == "then",
         site == "CBP")

png("figures/Met_within_across_thennow.png",
    width = 10, height = 4, res = 300, units = "in")
ylim = c(-4,2.5)
par(mfrow = c(1,3))
plot(now_19$month, now_19$GPP_mean, type = 'n', 
     ylab = "gO2/m2/d", xlab = "month", main = "across sites",
     ylim = ylim)
for(s in unique(now_19$site)){
  lines(now_19$month[now_19$site == s], now_19$GPP_mean[now_19$site == s],
        col = "forestgreen", lwd = 2 )
  lines(now_19$month[now_19$site == s], now_19$ER_mean[now_19$site == s],
        col = "sienna", lwd = 2)
}
abline(h = 0)
plot(now_y$month, now_y$GPP_mean, type = 'n', 
     ylab = "gO2/m2/d", xlab = "month", main = "across years",
     ylim = ylim)
for(s in unique(now_y$site)){
  now_yy <- now_y %>% 
    filter(site == s)
  for(y in unique(now_yy$year)){
    lines(now_yy$month[now_yy$year == y], now_yy$GPP_mean[now_yy$year == y],
          col = "forestgreen", lwd = 2)
    lines(now_yy$month[now_yy$year == y], now_yy$ER_mean[now_yy$year == y],
          col = "sienna", lwd = 2)
  }
}
abline(h = 0)
plot(now_19$month[now_19$site == "CBP"], now_19$GPP_mean[now_19$site == "CBP"], 
     type = 'l',ylab = "gO2/m2/d", xlab = "month", main = "across decades",
     ylim = ylim, col = "forestgreen", lwd = 2)
abline(h = 0)

lines(now_19$month[now_19$site == "CBP"], now_19$ER_mean[now_19$site == "CBP"],
      col = "sienna", lwd = 2)
lines(then$month, then$GPP_mean, col = "forestgreen", lwd = 2)
lines(then$month, then$ER_mean, col = "sienna", lwd = 2)

dev.off()

png("figures/PR_within_across_thennow.png",
    width = 10, height = 4, res = 300, units = "in")
par(mfrow = c(1,3))
ylim = c(0,1.5)
plot(now_19$month, -now_19$pr_mean, type = 'n', 
     ylab = "gO2/m2/d", xlab = "month", main = "across sites",
     ylim = ylim)
for(s in unique(now_19$site)){
  lines(now_19$month[now_19$site == s], -now_19$pr_mean[now_19$site == s],
        lwd = 2 )
}
abline(h = 1)
plot(now_y$month, -now_y$pr_mean, type = 'n', 
     ylab = "gO2/m2/d", xlab = "month", main = "across years",
     ylim = ylim)
for(s in unique(now_y$site)){
  now_yy <- now_y %>% 
    filter(site == s)
  for(y in unique(now_yy$year)){
    lines(now_yy$month[now_yy$year == y], -now_yy$pr_mean[now_yy$year == y],
          lwd = 2)
  }
}
abline(h = 1)
plot(now_19$month[now_19$site == "CBP"], -now_19$pr_mean[now_19$site == "CBP"], 
     type = 'l',ylab = "gO2/m2/d", xlab = "month", main = "across decades",
     ylim = ylim, lwd = 2)
abline(h = 1)
lines(then$month, -then$pr_mean, lwd = 2)

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
      


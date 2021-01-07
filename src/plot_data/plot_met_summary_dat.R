# Plot summarized metabolism data for Hall comparison

library(tidyverse)
library(lubridate)
library(ggplot2)
setwd("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism/data")

dat <- read_csv("metabolism/compiled/raymond_summary.csv")
dat <- read_csv("metabolism/compiled/metabolism_summary_tables_v2.csv") %>%
  mutate(method = factor(method, levels = c("hall", "hall_now", "raymond")))

# remove duplicate hall buisiness
ldat <- dat[c(-21, -22,-25),] %>%
  select(-data_days, -days, -pct_cov, -gpp_peak, -er_peak, -nep_cum_gC) %>%
  mutate(across(starts_with("er"),~ -.)) %>%
  pivot_longer(cols = c(-site, -year, -method),
               names_to = "variable",
               values_to = "gC_m2_time")
  
               
labs = c("median GPP daily", "max GPP daily", "median ER daily",
          "max ER daily",  "GPP cumulative annual", "ER cumulative annual")
                                             
names(labs) <- ldat$variable[1:6, drop = T]

png("../figures/metabolism_summary_comparison.png", 
    width = 6, height = 4, res = 300, units = "in")
ggplot(ldat, aes(x = method, y = gC_m2_time, fill = method)) +
  geom_boxplot() +
  facet_wrap(~ variable, scales = "free_y",
             labeller = labeller(variable = labs))+
  ggtitle("Metabolism summarized by site years")

dev.off()



lldat <- filter(ldat, method =="hall_now") %>%
  mutate(siteyear = paste0(site, year)) 
ggplot(lldat, aes(x = site, y = gC_m2_time, fill = factor(year))) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(variable ~., scales = "free_y")

rdat <- filter(dat, method == "raymond")

png("../figures/interannual_vs_intersite_met_distribution_pdfs.png",
    width = 6, height = 5, units = "in", res = 300)
par(mfrow = c(2,2), mar = c(4,4,2,2))
plot(density(rdat$gpp_median_gC), xlab = "median GPP", main = "")
abline(v = rdat$gpp_median_gC[rdat$site == "nhc", drop = T], lwd = 1.5, 
       col = "steelblue")
abline(v = rdat$gpp_median_gC[rdat$site == "unhc", drop = T], lwd = 1.5, 
       col = "steelblue", lty = 2)
abline(v = rdat$gpp_median_gC[rdat$year == 2019, drop = T]+0.004, 
       lwd = 1.5, col = "sienna3")
plot(density(rdat$er_median_gC), xlab = "median ER", main = "")
abline(v = rdat$er_median_gC[rdat$site == "nhc", drop = T], lwd = 1.5, 
       col = "steelblue")
abline(v = rdat$er_median_gC[rdat$site == "unhc", drop = T], 
       lwd = 1.5, col = "steelblue", lty = 2)
abline(v = rdat$er_median_gC[rdat$year == 2019, drop = T]+0.015, 
       lwd = 1.5, col = "sienna3")
plot(density(rdat$gpp_cum_gC), xlab = "cummulative GPP", main = "")
abline(v = rdat$gpp_cum_gC[rdat$site == "nhc", drop = T], lwd = 1.5, 
       col = "steelblue")
abline(v = rdat$gpp_cum_gC[rdat$site == "unhc", drop = T], lwd = 1.5, 
       col = "steelblue", lty = 2)
abline(v = rdat$gpp_cum_gC[rdat$year == 2019, drop = T]+1, 
       lwd = 1.5, col = "sienna3")
plot(density(rdat$er_cum_gC), xlab = "cummulative ER", main = "")
abline(v = rdat$er_cum_gC[rdat$site == "nhc", drop = T], lwd = 1.5, 
       col = "steelblue")
abline(v = rdat$er_cum_gC[rdat$site == "unhc", drop = T], 
       lwd = 1.5, col = "steelblue", lty = 2)
abline(v = rdat$er_cum_gC[rdat$year == 2019, drop = T]+ 4, 
       lwd = 1.5, col = "sienna3")
par(new = T, mfrow = c(1,1), mar = c(0,0,0,0))
plot(1,1, type = "n")
legend("top",
       c("2017-2019 NHC", "2017-2019 UNHC", "2019 all sites"),
       col = c("steelblue", "steelblue","sienna3"),
       lty = c(1,2,1), lwd = 1.5, ncol = 3, bty = "n")
dev.off()

var.test(rdat$gpp_cum_gC[rdat$site =="nhc"], 
         rdat$gpp_cum_gC[rdat$year == 2019], 
         alternative = "two.sided")


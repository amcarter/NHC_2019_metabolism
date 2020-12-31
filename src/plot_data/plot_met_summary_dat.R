# Plot summarized metabolism data for Hall comparison

library(tidyverse)
library(lubridate)
library(ggplot2)
setwd("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism/data")

dat <- read_csv("metabolism/compiled/metabolism_summary_tables_R.csv") %>%
  mutate(method = factor(method, levels = c("hall", "fixedk", "modeledk")))

ldat <- dat[c(-24,-27),] %>%
  select(-data_days, -days, -pct_cov, -gpp_peak, -er_peak, -nep_cum_gC) %>%
  mutate(across(starts_with("er"),~ -.)) %>%
  pivot_longer(cols = c(-site, -year, -method),
               names_to = "variable",
               values_to = "gC_m2_time")
  
               
labs = c("GPP median daily", "GPP max daily", "ER median daily",
          "ER max daily",  "GPP cumulative annual", "ER cumulative annual")
                                             
names(labs) <- ldat$variable[1:6, drop = T]

png("../figures/metabolism_summary_comparison.png", 
    width = 6, height = 4, res = 300, units = "in")
ggplot(ldat, aes(x = method, y = gC_m2_time, fill = method)) +
  geom_boxplot() +
  facet_wrap(~ variable, scales = "free_y",
             labeller = labeller(variable = labs))+
  ggtitle("Metabolism summarized by site years")

dev.off()


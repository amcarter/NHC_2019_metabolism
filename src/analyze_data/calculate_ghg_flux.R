#######################################
# Calculate GHG concentrations and fluxes 
# A Carter
# 2020 12 10

setwd(metab_projdir)

library(tidyverse)
library(ggplot2)

# read in data file
# This has been processed from raw GC data to dissolved concentrations using A Helton's headspace calcs worksheet

# NEED TO add calculation of variability here when summarizing duplicates
# might also be worth figuring out how to incorporate the error percent 
#   (based on times when I had to estimate vial mass)
gas <- read_csv("data/gas_data/NHC_2019-2020_processed_GHGdata.csv") %>%
  mutate(date = as.Date(Date, format="%m/%d/%Y")) %>%
  group_by(Site, date) %>%
  summarize_if(is.numeric, mean, na.rm = T)

gas <- gas[gas$Site!="MC751",]
# gas<- gas %>% select(-Date) %>%
#   left_join(k[, 1:3], by="date")

# plot(k$date, k$K600)
# par(new=T)
# plot(k$date, k$discharge_m3s, log="y", type="l")

ggplot(gas) +
 aes(x = CO2.ugL, y = CH4.ugL, size = `water_temp_C`) +
 geom_point() +
 theme_minimal()


gas <- pivot_longer(data = gas, 
                    cols = c("CH4.ugL","CO2.ugL","N2O.ugL"),
                    names_to = "gas",
                    values_to = "gas_ugL")
gas$site <- factor(gas$Site, levels=c("UNHC","PWC","WBP","WB","CBP","PM","NHC"))
gas %>% select(gas, -Site)

ggplot(gas) +
  aes(x=date, y=gas_ugL, color=Site)+
  geom_point()+
  facet_grid(gas~., scales="free")


site <- "NHC"
k <- read_csv(paste0("data/gas_data/night_regression/nightreg_", site,".csv"))
k$site <- site

sites <- c( "UNHC", "PM","WB","WBP","CBP")

for(site in sites){
  kk <- read_csv(paste0("data/gas_data/night_regression/nightreg_", site,".csv"))
  dates <- data.frame(date = seq(kk$date[1], kk$date[nrow(kk)], by = "1 day"))
  kk <- kk %>% right_join(dates, by = "date") %>%
    arrange(date) %>% 
    select(date, discharge_m3s, lnK600 = lnK600.pred, lnK600.sd) %>%
    mutate(lnK600 = na.approx(lnK600, na.rm = F),
           lnK600.sd = na.approx(lnK600.sd, na.rm = F))
  kk$site <- site
  k<- bind_rows(k, kk)
}
  #k[k$K600<0,"K600"] <- 0
# K600 is the value from each day's night regression. 
# lnK600.pred is the one selected based on the KQ relationship
k <- k %>% 
  select(date, site, lnK600, lnK600.sd)


gas_flux <- left_join(gas, k, by=c("date", "site")) %>%
  mutate(K600 = exp(lnK600),
         flux_ugLd = gas_ugL*K600,
         flux_ugLd.sd = gas_ugL*exp(lnK600.sd),
         site = factor(site, levels=c("UNHC","PWC","WBP","WB","CBP","PM","NHC"))) %>%
  select(site, date, gas, water_temp_C, air_pres_mmHg, 
         K600, flux_ugLd, flux_ugLd.sd)


ggplot(gas_flux) +
  aes(x=date, y=K600, color=Site)+
  geom_point()


ggplot(gas_flux) +
  aes(x=date, y=flux_ugLd, color=Site, ylab = "ug/L/d")+
  geom_point()+
  facet_grid(gas~., scales="free")



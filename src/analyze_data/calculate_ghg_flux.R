#######################################
# Calculate GHG concentrations and fluxes 
# A Carter
# 2020 08 10

setwd(metab_projdir)

library(tidyverse)
library(ggplot2)

# read in data file
# This has been processed from raw GC data to dissolved concentrations using A Helton's headspace calcs worksheet

gas <- read_csv("data/gas_data/NHC_2019-2020_processed_GHGdata.csv")

gas$date <- as.Date(gas$Date, format="%m/%d/%Y")
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
gas %>% select(gas, -Site, -Date)

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
  kk$site <- site
  k<- bind_rows(k, kk)
}
  #k[k$K600<0,"K600"] <- 0
k <-k %>% 
  select(date, site, K600.pred, K600.se)

gas <- left_join(gas, k, by=c("date", "site"))
ggplot(gas) +
  aes(x=date, y=K600.pred, color=Site)+
  geom_point()

gas$flux_ugLd <- gas$gas_ugL*gas$K600.pred
gas$site <- factor(gas$Site, levels=c("UNHC","PWC","WBP","WB","CBP","PM","NHC"))


ggplot(gas) +
  aes(x=date, y=flux_ugLd, color=Site, ylab = "ug/L/d")+
  geom_point()+
  facet_grid(gas~., scales="free")



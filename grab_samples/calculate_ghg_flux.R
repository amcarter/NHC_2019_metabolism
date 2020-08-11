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



ggplot(gas) +
 aes(x = CO2.ugL, y = CH4.ugL, colour = Site, size = `water_temp_C`) +
 geom_point() +
 scale_color_hue() +
 theme_minimal()

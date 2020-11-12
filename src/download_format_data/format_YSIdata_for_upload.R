##############################
# code to format ysi data for uploading to the portal
# this references the datasheet I use to enter values after I go in the field
# it creates a new correctly formated file with only the data after a specified date until now


setwd(metab_projdir)

ysi <- read_csv("data/water_chemistry/NHC_YSI_data.csv")
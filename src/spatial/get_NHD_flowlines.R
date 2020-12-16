################################
# Get NHD flowlines and data for NHC
# A Carter
# 12/8/2020

# Consult these docs from this Download link
# (ftp://ftp.horizon-systems.com/NHDplus/NHDPlusV21/Documentation/NHDPlusV2_User_Guide.pdf) 
# to understand NHD column values, their units etc.
# This is long, so just control-f on the information you need.
# 
# This map has all USGS gages
# (https://maps.waterdata.usgs.gov/mapper/index.html) 
# 
# This tool is nice for visualizing HUCs
# (http://waterqualityexplorer.rc.duke.edu:3838/explorer/) 


library(tidyverse)
library(nhdplusTools)
library(sf)
library(tmap)
library(ggplot2)
tmap_mode("view") #makes map plots interactive

setwd("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism")

# Get NHD flowlines ####
# Get sample station locations
sites <- read_csv("data/siteData/NHCsite_coordinates.csv")
metab_sites <- sites[c(1:5,7),]
NHC_sites_sf <- sites[c(1:5,7),2:4] %>% 
  st_as_sf(coords=c("longitude","latitude"),remove=F, crs=4326)

# get all flowlines upstream of Blands USGS station
NHC_gage_id <- "USGS-02097314"
NHC_site <- list(featureSource = "nwissite",  #prepping below web quirie below, this just indicates the kind of data we have
                 featureID = NHC_gage_id) 
NHC_line <- navigate_nldi(NHC_site, "UT", "", 100) #returns all upstream flowlines

NHC_gages_sf <- navigate_nldi(NHC_site, "UT", "nwissite", 100) #returns all upstream USGS gages

# now map flowlines and gages
tm_shape(NHC_line) + tm_lines() +
  tm_shape(NHC_sites_sf) + tm_dots(col = "purple") +  
  tm_shape(NHC_gages_sf) + tm_dots(col = "red")


# Download everything from NHD! ####
#Flowlines above are just the ID's . Download everything!
save_file_name <- "data/siteData/NHC_NHD_subset.gpkg" #gkpg is an open source geospatial format
subset_nhdplus(NHC_line$nhdplus_comid, 
               save_file_name, 
               "download")
NHC_nhd <- st_read(save_file_name)
glimpse(NHC_nhd)
 

# attach comids to each sample station
comid_points<- rep(NA, nrow(NHC_sites_sf))
for(i in 1:nrow(NHC_sites_sf)){
  comid_points[i]<- discover_nhdplus_id(NHC_sites_sf[i,])  
}

sample_site_NHD_reaches <- tibble(comid=comid_points) %>%
  left_join(NHC_nhd)%>% st_as_sf()

siteDat <- sample_site_NHD_reaches %>% 
  select(comid, streamorder = streamorde, 
         pathlength, slope, 
         watershed_area_km2 = totdasqkm )
siteDat$site <- NHC_sites_sf$site

siteDat %>% 
  select(pathlength, streamorder, slope, watershed_area_km2) %>%
  as_tibble() %>%
  pivot_longer(cols = c(2,3,4), names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = pathlength, y = value)) +
    geom_line() +
    facet_wrap(.~variable, scales = "free_y", ncol = 1)


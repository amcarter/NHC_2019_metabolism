################################
# Make Map of NHC watershed with different subsets of sites for presentation 
# A Carter
# 3.27.2021

library(nhdplusTools)
library(tidyverse)
library(sf)
library(tmap)
library(rgdal)
library(maps)
library(rgee)

setwd('C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/hall_50yl2/NHC_2019_metabolism')

# Get sample station locations
syn_sites <- read.csv("data/map_files/NC_synoptic_site_coordinates.csv", 
                  header=T, stringsAsFactors = F)
upper_sites <- read.csv('data/map_files/NHCsite_50yl_coordinates.csv', 
                        header = T, stringsAsFactors = F)
syn_sites_sf <- st_as_sf(syn_sites[c(1:13),], coords=c("Long","Lat"),
                         remove=F, crs=4326)
wwtp_sf <- st_as_sf(syn_sites[15,], coords=c("Long","Lat"),remove=F, crs=4326)
up_sites_sf <- st_as_sf(upper_sites, coords = c('longitude', 'latitude'),
                        remove = F, crs = 4326)
carter_sites <- filter(up_sites_sf, type == 'now')
hall_sites <- filter(up_sites_sf, type == 'hall')

# get streamlines, etc

# nhc_lines = nhdplusTools::navigate_nldi(
#   nldi_feature = list(featureSource = "comid",
#                       featureID = as.character(comid)),
#   mode = 'UT', data_source = '', distance_km = 100)
# 
# st_write(st_combine(nhc_lines),
#          dsn = 'data/map_files',
#          layer = 'stream_lines',
#          driver = 'ESRI shapefile', 
#          delete_layer = TRUE)
# 
# hall_study_reach = filter(nhc_lines,
#                           nhdplus_comid %in% c(8895490, 8895362, 8895420, 8895440)) %>%
#   st_combine()
# 
# st_write(hall_study_reach,
#          dsn = 'data/map_files',
#          layer = 'hall_study_reach',
#          driver = 'ESRI shapefile',
#          delete_layer = TRUE)
# 
# nhc_ripar = hall_study_reach %>%
#   st_transform(crs = 4326) %>%
#   st_buffer(dist = 250) %>%
#   st_transform(crs = 4326)
# 
# st_write(nhc_ripar,
#          dsn = 'data/map_files',
#          layer = 'riparian',
#          driver = 'ESRI shapefile',
#          # driver = 'GeoJSON')
#          delete_layer = TRUE)

#synoptic jobbies
cur_nhd <- st_read("data/map_files/NHC_NHD_subset.gpkg") 
NHC_mainstem <- cur_nhd %>% 
  filter(gnis_name=="Mud Creek" | gnis_name=="New Hope Creek"|comid==8888400)
UNHC_NHD <- NHC_mainstem %>%
  filter(hydroseq > 250019295 &gnis_name=="New Hope Creek")
UMC_NHD <- NHC_mainstem %>%
  filter(hydroseq > 250087220 & gnis_name=="Mud Creek" )
longitudinal_transect <- NHC_mainstem %>%
  filter(!comid %in% UNHC_NHD$comid)%>%
  filter(!comid %in% UMC_NHD$comid)

# load shape files
duke_forest_boundary <- st_read("data/map_files/2019_boundary.shp")
korstian_div <- filter(duke_forest_boundary, DIVISION == 'Korstian') %>%
  sf::st_transform(crs=4326)
stream_line <- st_read("data/map_files/stream_lines.shp")
riparian_boundary <- st_read("data/map_files/riparian.shp")
study_reaches_line <- st_read("data/map_files/hall_study_reach.shp")
watershed_boundary <- st_read('data/map_files/nhcwsboundary.shp') 
mud_ws_boundary <- st_read('data/map_files/mudwsboundary.shp') 
upper_ws_boundary <- st_read("data/map_files/nhc_wb_streamstats.shp")


# make maps
# tmap_mode('view')
tmap_mode('plot')
basic_map = tm_shape(watershed_boundary) + 
  tm_polygons(alpha=0, border.col="black", lwd=1.5) +
  tm_shape(stream_line) + tm_lines(col='grey45', lwd = 1.2) +
  tm_compass(type="arrow", position=c("right", "bottom", show.labels=3),
             size=3, text.size=1) +
  tm_scale_bar(text.size = 1, position="right", width = .2) +
  tm_style(style='white') +
  tm_layout(frame=FALSE, bg.color="white") 

tmap_save(basic_map, filename="figures/map/basic_watershed.png",
          bg="white", dpi = 300)

synoptic = tm_shape(watershed_boundary) + 
  tm_polygons(alpha=0, border.col="black", lwd=1.5) +
  tm_shape(stream_line) + tm_lines(col='grey45',lwd=1.2) +
  # tm_shape(longitudinal_transect) + tm_lines(col='black', lwd=2) +
  tm_shape(syn_sites_sf) + tm_dots( col="brown3", size=0.4) +
  tm_compass(type="arrow", position=c("right", "bottom", show.labels=3),
             size=3, text.size=1) +
  tm_scale_bar(text.size = 1, position="right", width = .2) +
  tm_style(style='white') 


tmap_save(synoptic, filename="figures/map/synoptic_watershed.png",
          bg="white", dpi = 300)

longitudinal = tm_shape(watershed_boundary) + 
  tm_polygons(alpha=0, border.col="black", lwd=1.5) +
  tm_shape(stream_line) + tm_lines(col='grey45',lwd=1.2) +
  tm_shape(longitudinal_transect) + tm_lines(col='steelblue', lwd=3) +
  tm_shape(syn_sites_sf) + tm_dots( col="brown3", size=0.4) +
  tm_compass(type="arrow", position=c("right", "bottom", show.labels=3),
             size=3, text.size=1) +
  tm_scale_bar(text.size = 1, position="right", width = .2) +
  tm_style(style='white') 


tmap_save(longitudinal, filename="figures/map/longitudinal_watershed.png",
          bg="white", dpi = 300)


metab = tm_shape(watershed_boundary) + 
  tm_polygons(alpha=0, border.col="black", lwd=1.5) +
  tm_shape(korstian_div) + tm_polygons(alpha=0.5, col = 'springgreen3',
                                       border.col="transparent", lwd=.5) +
  # tm_shape(riparian_boundary) + tm_polygons(alpha=0, col="black", lwd=1.5,
  #                                           border.col='steelblue3', border.alpha=0.8) +
  # tm_shape(study_reaches_line) + tm_lines(col='steelblue3', lwd=2.5) +
  tm_shape(stream_line) + tm_lines(col='grey45',lwd=1.2) +
  # tm_shape(longitudinal_transect) + tm_lines(col='black', lwd=2) +
  tm_shape(carter_sites) + tm_dots(col="brown3", size=0.4) +
  # tm_shape(hall_sites) + tm_symbols(shape=3, col="black", size=0.6, border.lwd=2) +
  tm_compass(type="arrow", position=c("right", "bottom", show.labels=3),
             size=3, text.size=1) +
  tm_scale_bar(text.size = 1, position="right", width = .2) +
  tm_style(style='white') +
  tm_layout(frame=FALSE, bg.color="white") +
  tm_add_legend(type='fill', labels = 'Duke Forest', col = 'springgreen3', alpha=0.1,
                border.col='transparent') +
  tm_legend(show=TRUE, position=c('right', 'top'), outside=FALSE, bg.color='gray97',
            frame=FALSE, text.size=1.1)

tmap_save(metab, filename="figures/map/metab_watershed.png",
          bg="white", dpi = 300)

hall = tm_shape(watershed_boundary) + 
  tm_polygons(alpha=0, border.col="black", lwd=1.5) +
  tm_shape(korstian_div) + tm_polygons(alpha=0.5, col = 'springgreen3',
                                       border.col="transparent", lwd=.5) +
  # tm_shape(riparian_boundary) + tm_polygons(alpha=0, col="black", lwd=1.5,
  #                                           border.col='steelblue3', border.alpha=0.8) +
  # tm_shape(study_reaches_line) + tm_lines(col='steelblue3', lwd=2.5) +
  tm_shape(stream_line) + tm_lines(col='grey45',lwd=1.2) +
  # tm_shape(longitudinal_transect) + tm_lines(col='black', lwd=2) +
  tm_shape(hall_sites) + tm_dots( col="black", size=0.01) +
  tm_shape(carter_sites) + tm_dots(col="brown3", size=0.4) +
  tm_compass(type="arrow", position=c("right", "bottom", show.labels=3),
             size=3, text.size=1) +
  tm_scale_bar(text.size = 1, position="right", width = .2) +
  tm_style(style='white') +
  tm_layout(frame=FALSE, bg.color="white") +
  tm_add_legend(type='fill', labels = 'Duke Forest', col = 'springgreen3', alpha=0.1,
                border.col='transparent') +
  tm_add_legend(type='symbol', labels = '    1969 sites', col = 'black',
                size=0.01, shape=3, border.lwd=2) +
  tm_legend(show=TRUE, position=c('right', 'top'), outside=FALSE, bg.color='gray97',
            frame=FALSE, text.size=1.1)

tmap_save(hall, filename="figures/map/hall_watershed.png",
          bg="white", dpi = 300)
df = tm_shape(watershed_boundary) + 
  tm_polygons(alpha=0, border.col="black", lwd=1.5) +
  tm_shape(korstian_div) + tm_polygons(alpha=0.5, col = 'springgreen3',
                                       border.col="transparent", lwd=.5) +
  tm_shape(stream_line) + tm_lines(col='grey45',lwd=1.2) +
  # tm_shape(hall_sites) + tm_dots( col="black", size=0.01) +
  tm_compass(type="arrow", position=c("right", "bottom", show.labels=3),
             size=3, text.size=1) +
  tm_scale_bar(text.size = 1, position="right", width = .2) +
  tm_style(style='white') +
  tm_layout(frame=FALSE, bg.color="white") +
  tm_add_legend(type='fill', labels = 'Duke Forest', col = 'springgreen3', alpha=0.1,
                border.col='transparent') +
  tm_add_legend(type='symbol', labels = '    1969 sites', col = 'black',
                size=0.01, shape=3, border.lwd=2) +
  tm_legend(show=TRUE, position=c('right', 'top'), outside=FALSE, bg.color='gray97',
            frame=FALSE, text.size=1.1)

tmap_save(df, filename="figures/map/df_watershed.png",
          bg="white", dpi = 300)
metab = tm_shape(watershed_boundary) + 
  tm_polygons(alpha=0, border.col="black", lwd=1.5) +
  tm_shape(korstian_div) + tm_polygons(alpha=0.3, col = 'springgreen3',
                                       border.col="transparent", lwd=.5) +
  # tm_shape(riparian_boundary) + tm_polygons(alpha=0, col="black", lwd=1.5,
  #                                           border.col='steelblue3', border.alpha=0.8) +
  # tm_shape(study_reaches_line) + tm_lines(col='steelblue3', lwd=2.5) +
  tm_shape(stream_line) + tm_lines(col='grey45',lwd=1.2) +
  # tm_shape(longitudinal_transect) + tm_lines(col='black', lwd=2) +
  tm_shape(carter_sites) + tm_dots(col="brown3", size=0.4) +
  # tm_shape(hall_sites) + tm_symbols(shape=3, col="black", size=0.6, border.lwd=2) +
  tm_compass(type="arrow", position=c("right", "bottom", show.labels=3),
             size=3, text.size=1) +
  tm_scale_bar(text.size = 1, position="right", width = .2) +
  tm_style(style='white') +
  tm_layout(frame=FALSE, bg.color="white") +
  tm_add_legend(type='symbol', labels = '  Study sites', col = 'red2', size = 0.7,
                shape=1) +
  tm_add_legend(type='symbol', labels = '  Hall 1972 sites', col = 'black',
                size=0.5, shape=3, border.lwd=2) +
  tm_add_legend(type='line', labels = 'Study reach', col = 'steelblue3', lwd = 2.5) +
  tm_add_legend(type='line', labels = 'Riparian zone', col = 'steelblue3', lwd = 1) +
  tm_add_legend(type='fill', labels = 'Duke Forest', col = 'springgreen3', alpha=0.3,
                border.col='transparent') +
  tm_legend(show=TRUE, position=c('right', 'top'), outside=FALSE, bg.color='gray97',
            frame=FALSE, text.size=1.1)

tmap_save(basic_map, filename="figures/map/basic_watershed.png",
          bg="white", dpi = 300)

map_without_colon = tm_shape(watershed_boundary) + tm_polygons(alpha=0, border.col="black", lwd=1) +
  tm_shape(korstian_div) + tm_polygons(alpha=0.3, col = 'springgreen3',
                                       border.col="transparent", lwd=.5) +
  tm_shape(study_reaches_line_shortened) + tm_lines(col='steelblue3', lwd=2.5) +
  tm_shape(stream_line_shortened) + tm_lines(col='black', alpha=0.5, lwd=0.5) +
  tm_shape(carter_sites) + tm_symbols(shape=1, col="red2", size=0.6, border.lwd=2) +
  tm_scale_bar(text.size = 1, position="left") +
  tm_compass(type="arrow", position=c("right", "bottom", show.labels=3),
             size=5, text.size=1) +
  tm_style(style='white') +
  tm_layout(frame=TRUE, bg.color="white") +
  tm_add_legend(type='symbol', labels = '  Study sites', col = 'red2', size = 0.7,
                shape=1) +
  tm_add_legend(type='line', labels = 'Study reach', col = 'steelblue3', lwd = 2.5) +
  tm_add_legend(type='fill', labels = 'Duke Forest', col = 'springgreen3', alpha=0.3,
                border.col='transparent') +
  tm_legend(show=TRUE, position=c('left', 'top'), outside=FALSE, bg.color='gray97',
            frame=TRUE, text.size=1.1)

tmap_save(map_without_colon, filename="figs/map_without_colon.png",
          bg="white", dpi = 300)

tmap_save(map_with_colon, filename="figs/map_with_colon.png",
          bg="white", dpi = 300)
# plot
tmap_mode("view")

tm_shape(nhc_ws)+tm_polygons(alpha=0, border.col="black", lwd=.5)+
    tm_shape(mud_ws)+tm_polygons(alpha=1, border.col="black",lwd=.5)+
  tm_shape(cur_nhd)+tm_lines(col = "grey60") +
  tm_shape(longitudinal_transect) + tm_lines(lwd=2)+
  tm_shape(sites_sf)+tm_dots(col="brown3", size=.05)+
    tm_shape(wwtp_sf)+tm_markers(shape=3, col="lightblue",size=.05)+
  tm_shape(pt_sf) + tm_dots()+
  tm_scale_bar(text.size = 1, position = "left") 

tmap_mode("plot")
par(bg=NA)
map<-tm_shape(nhc_ws)+tm_polygons(alpha=0, border.col="black", lwd=.5)+
#  tm_shape(mud_ws)+tm_polygons(alpha=0, border.col="black",lwd=.5)+
  tm_shape(cur_nhd)+tm_lines(col = "grey60") +
  tm_shape(longitudinal_transect) + tm_lines(lwd=2)+
  tm_shape(sites_sf)+tm_dots(col="brown3", size=.05)+
#  tm_shape(wwtp_sf)+tm_markers(shape=3, col="lightblue",size=.05)+
  tm_scale_bar(text.size = 1, position = "left") +
  tm_compass(type="arrow",position=c("right","bottom", show.labels=3))+
  tm_layout(frame=FALSE, bg.color="transparent")
tmap_save(map, filename="NHCmap_scalebar.eps", bg="transparent", dpi = 1200, 
          )


# Plot of longitudinal transect
long_sites_sf <- sites_sf[sites_sf$site!="MC751",]
tmap_mode("view")

tm_shape(cur_nhd)+tm_lines(col = "grey80") +
  tm_shape(longitudinal_transect) + tm_lines(lwd=2)+
  tm_shape(long_sites_sf)+tm_dots(col="brown3", size=.05)+
  #  tm_shape(wwtp_sf)+tm_markers(shape=3, col="lightblue",size=.05)+
  tm_scale_bar(text.size = 1, position = "left") 

# Plot North carolina with piedmont shape


pied <- readOGR(dsn="ncpiedmont_shape",
                layer="Piedmont_shape")

par(bg=NA)
png("NCmap.png",bg="transparent", type="windows")
map('state',region='North Carolina',fill=TRUE, col="white",bg="transparent",lwd=2)
  plot(pied,
     add=TRUE, 
     col="grey90")
  points(wwtp_sf$Long, wwtp_sf$Lat, col="brown3", pch=22, cex=3)
dev.off()
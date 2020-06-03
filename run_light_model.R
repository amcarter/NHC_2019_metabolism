#===============================================================================
#Workflow for getting light estimates for sites in NHC watershed
#modified from Audrey's code
#Created 5/14/2020
#===============================================================================
 # library(devtools)
 # devtools::install_github("psavoy/StreamLightUtils")
 # devtools::install_github("psavoy/StreamLight", auth_token = "7a3e9b10f997aeada28437a6ec760bff3a807a2e") 
library(StreamLight)
library(StreamLightUtils)
library(dplyr)
library(lubridate)

#-------------------------------------------------
#Downloading NLDAS data
#-------------------------------------------------
  #Getting site information
    sites <- read.csv("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism/data/siteData/NHCsite_metadata.csv", header = TRUE)
    startdate <- "2019-01-01"
    sites <- select(sites, Site_ID=sitecode, Lat=latitude, Lon=longitude)
    sites$Site_ID<- paste0("NC_", as.character(sites$Site_ID))
    write.csv(sites, "C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism/data/light/NHC_site_requests.csv", row.names = FALSE, col.names=FALSE)

  #Downloading NLDAS data for all sites except Walker Branch
    save_dir<-DL_NLDAS_dir <- "C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism/data/light/NLDAS"
    
  NLDAS_DL_AMC<-function (save_dir, Site_ID, Lat, Lon, startDate) {
    http_string <- paste("http://hydro1.gesdisc.eosdis.nasa.gov/daac-bin/access/timeseries.cgi?variable=NLDAS:NLDAS_FORA0125_H.002:DSWRFsfc&location=GEOM:POINT(")
    start_split <- strsplit(startDate, "-")[[1]]
    url <- paste(http_string, Lon, ", ", Lat, ")&startDate=", 
                 start_split[1], "-", start_split[2], "-", 
                 start_split[3], "T00", "&type=asc2", sep = "")
    destfile <- paste(save_dir, "/", Site_ID, "_NLDAS.asc", 
                      sep = "")
    try_result <- try(download.file(url, destfile, method = "wininet"), 
                      silent = FALSE)
    if (class(try_result) == "try-error") {
      file.remove(destfile)
    }
  }
    
    
  mapply(NLDAS_DL_AMC, save_dir = DL_NLDAS_dir, sites[, "Site_ID"], sites[, "Lat"], 
      sites[, "Lon"], startdate)
 
    #List of successfully downloaded sites
    NLDAS_list <- str_sub(list.files("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism/data/light/NLDAS"),
                          1, -11)   
    
    #Processing the downloaded NLDAS data
    NLDAS_proc_rd <- "C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism/data/light/NLDAS"
    NLDAS_processed <- NLDAS_proc(NLDAS_proc_rd, NLDAS_list)  
    
#-------------------------------------------------
#MODIS data
#-------------------------------------------------
  #Get list of Site_ID's used for the request
    setwd("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism/data/light/AppEEARS_MODIS_LAI")
    request_sites <- read.csv("../NHC_site_requests.csv", header = TRUE)[, "Site_ID"] 
    
  #Directory where the zipped data is located
    zip_dir <- "C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism/data/light/AppEEARS_MODIS_LAI"

  #Request zip file
    zip_file <- "nhc_points.zip"

  #Unpack the MODIS AppEEARS request
    MOD_unpack <- AppEEARS_unpack(zip_file, zip_dir, request_sites)
    
  #Applying the function
    MOD_processed <- lapply(MOD_unpack, FUN = AppEEARS_proc, proc_type = "spline")
    
#-------------------------------------------------
#Create model driver files
#-------------------------------------------------
  #Applying the function
    driver_sd <- "C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism/data/light/driver_files"
    make_driver(sites[, c("Site_ID", "Lat", "Lon")], NLDAS_processed, 
      MOD_processed, save_dir = driver_sd)  

#-------------------------------------------------
#Extracting tree height
#-------------------------------------------------
 # extract_height(sites[, "Site_ID"], sites[, "Lat"], sites[, "Lon"])
sites <- read.csv("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism/data/light/NHC_parameters.csv", header = T)
#-------------------------------------------------
#Running the model    
#-------------------------------------------------
  #Function for batching over multiple sites
    batch_model <- function(Site, read_dir, save_dir){
      #Get the model driver
        driver_file <- readRDS(paste(read_dir, "/", Site, "_driver.rds", sep = ""))

      #Get model parameters for the site
        site_p <- params[params[, "Site_ID"] == Site, ]

      #Run the model
        modeled <- stream_light(driver_file, Lat = site_p[, "Lat"], Lon = site_p[, "Lon"],
          stream_azimuth = site_p[, "Azimuth"], bottom_width = site_p[, "Width"], BH = 0.6,
          BS = 100, WL = site_p[, "WL"], TH = site_p[, "TH"], overhang = site_p[, "overhang"],
          overhang_height = site_p[, "overhang_height"], x_LAD = site_p[, "x"])

      #Save the output
        saveRDS(modeled, paste(save_dir, "/", Site, "_predicted.rds", sep = ""))

    } #End batch_model
    
  #Applying the model to all sites
    model_rd <- "C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism/data/light/driver_files"
    model_sd <- "C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism/data/light/predicted"

  #Get the parameters for all sites
    params <- read.csv("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism/data/light/NHC_parameters.csv",
      colClasses = c("character", rep("numeric", 5)), header = TRUE)

  #Running the model
    lapply(params[, "Site_ID"], FUN = batch_model, read_dir = model_rd,
      save_dir = model_sd)
    
#-------------------------------------------------
#Summarize the outputs as daily observations
#-------------------------------------------------
  #Set the working directory for the predicted values
    predicted_dir <- "C:/Users/Phil/Desktop/HBEF_Audrey/predicted"
    
  #Just visualize the results
    lapply(params[, "Site_ID"], FUN = daily_light, pred_dir = predicted_dir, 
      plot = TRUE, return = FALSE)
    
  #If you want to store the results you can wrap it in another function
    save_fun <- function(Site_ID, pred_dir, save_dir){
      daily <- daily_light(pred_dir, Site_ID, plot = FALSE, return = TRUE)
      
      #Save the outputs
        setwd(save_dir)
        saveRDS(daily, paste(Site_ID, "_daily.rds", sep = ""))
    } #End save_fun
    
  #Run the function to store your results
    save_dir <- "C:/Users/Phil/Desktop"
    lapply(params[, "Site_ID"], FUN = save_fun, pred_dir = predicted_dir,
      save_dir = save_dir)
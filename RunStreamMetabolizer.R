#########################################
# Run Stream Metabolizer on NHC data
# AMCarter 4-17-20

 # library(remotes)
 # remotes::install_github('appling/unitted')
 # remotes::install_github("USGS-R/streamMetabolizer")

library(streamMetabolizer)
library(dplyr)
library(readr)

setwd("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism")

# load RDS from model output:

filelist <- list.files("data/metabolism/processed")
bayes_name <- mm_name(type="bayes", pool_K600="binned", err_obs_iid = TRUE,
                      err_proc_iid = TRUE, err_proc_acor = FALSE)
bayes_specs <- specs(bayes_name)

for(i in 2:length(filelist)){
  dat <- read_csv(paste0("data/metabolism/processed/",filelist[i]))
  # set nodes for stream metabolizer to bin discharge around
  Qrange <- c(min(log(dat$discharge), na.rm=T),quantile(log(dat$discharge), .98, na.rm=T))
  bayes_specs$K600_lnQ_nodes_centers <- seq(Qrange[1], Qrange[2], length=7)
  sitename <- substr(filelist[i], 1, nchar(filelist[i])-4)
  sitename
  #fit metab model
  mm <- metab(bayes_specs, data=dat)
  saveRDS(mm, paste0("data/metabolism/modeled/",sitename,".rds"))
}


#############################
# extract data from rds
for(i in 1:length(filelist)){
    sitename <- substr(filelist[i], 1, nchar(filelist[i])-4)
    mod <- readRDS(paste0("data/metabolism/modeled/",sitename,".rds"))
    data_daily <- mod@data_daily %>% select(date, discharge.m3s=discharge.daily)
    K600 <- mod@fit$daily %>% select(date, K600=K600_daily_50pct, K600.lower=K600_daily_2.5pct, K600.upper=K600_daily_97.5pct)
    data_daily <- left_join(data_daily, K600, by="date")
    metab <- predict_metab(mod)
    metab <- left_join(metab, data_daily, by="date")
    
    data <- predict_DO(mod)
    
    mod_specs<- get_specs(mod)
    mod_specs$K600 <- exp(get_fit(mod)$KQ_binned$lnK600_lnQ_nodes_50pct)
    mod_specs$K600.lower <- exp(get_fit(mod)$KQ_binned$lnK600_lnQ_nodes_2.5pct)
    mod_specs$K600.upper <- exp(get_fit(mod)$KQ_binned$lnK600_lnQ_nodes_97.5pct)
    
    
    filter_vec = names(mod_specs) %in% c("K600_lnQ_nodes_centers", "K600", "K600.upper", "K600.lower")
    QvsK600 <- mod_specs[filter_vec]
    
    
    metab_dat <- list(data=data, metab=metab, QvsK600=QvsK600)
    
    saveRDS(metab_dat, paste0("data/metabolism/condensed/condensed_",sitename,".rds"))
}

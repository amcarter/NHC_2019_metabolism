#########################################
# Run Stream Metabolizer on NHC data
# AMCarter 4-17-20

 # library(remotes)
 # remotes::install_github('appling/unitted')
 # remotes::install_github("USGS-R/streamMetabolizer")

library(streamMetabolizer)
library(dplyr)
library(readr)

setwd(metab_projdir)

# load RDS from model output: discharge
nreg<- read_csv("data/siteData/KQ_nightreg_priors.csv")
filelist <- list.files("data/metabolism/processed")
bayes_name <- mm_name(type="bayes", pool_K600="binned", err_obs_iid = TRUE,
                      err_proc_iid = TRUE, err_proc_acor = FALSE)
bayes_specs <- specs(bayes_name)


all_data <- data.frame()
all_metab <- data.frame()
all_specs <- data.frame()

for(i in 1:length(filelist)){
  dat <- read_csv(paste0("data/metabolism/processed/",filelist[i]))
  dat <- select(dat, solar.time, DO.obs=DO_mgL, DO.sat, depth, temp.water=WaterTemp_C, light, discharge=Discharge_m3s)
  site <- substr(filelist[i],1,(nchar(filelist[i])-4))
  # set nodes for stream metabolizer to bin discharge around
  #Qrange <- c(min(log(dat$discharge), na.rm=T),quantile(log(dat$discharge), .98, na.rm=T))
  #bayes_specs$K600_lnQ_nodes_centers <- seq(Qrange[1], Qrange[2], length=7)
  
  #get nodes from KQ relationships determined from nighttime regression
  KQ <- nreg[nreg$site==site,]
  bayes_specs$K600_lnQ_nodes_centers <- KQ$lnQ_nodes
  bayes_specs$K600_lnQ_nodes_meanlog <- log(KQ$K600_lnQ_nodes)
  bayes_specs$K600_lnQ_nodes_sdlog <- sqrt(abs(log(KQ$K600_lnQ_nodes_se)))
  bayes_specs$K600_daily_sigma_sigma <- .15 
  #fit metab model
  mm <- metab(bayes_specs, data=dat)

  data_daily <- mm@data_daily %>% select(date, discharge.m3s=discharge.daily)
  K600 <- mm@fit$daily  %>% select(date, K600=K600_daily_50pct, K600.lower=K600_daily_2.5pct, K600.upper=K600_daily_97.5pct)
    data_daily <- left_join(data_daily, K600, by="date")
    metab <- predict_metab(mm)
    metab <- left_join(metab, data_daily, by="date")
    metab$site <- sitename
    all_metab <- rbind(all_metab, metab)

    data <- predict_DO(mm)
    data$site <- sitename
    all_data <- rbind(all_data, data) 
    
    mod_specs<- get_specs(mm)
    mod_specs$K600 <- exp(get_fit(mm)$KQ_binned$lnK600_lnQ_nodes_50pct)
    mod_specs$K600.lower <- exp(get_fit(mm)$KQ_binned$lnK600_lnQ_nodes_2.5pct)
    mod_specs$K600.upper <- exp(get_fit(mm)$KQ_binned$lnK600_lnQ_nodes_97.5pct)
    
    filter_vec = names(mod_specs) %in% c("K600_lnQ_nodes_centers", "K600", "K600.upper", "K600.lower")
    QvsK600 <- as.data.frame(mod_specs[filter_vec])
    QvsK600$site <- sitename
    all_specs<- rbind(all_specs, QvsK600)
    
    metab_dat <- list(data=data, metab=metab, specs=QvsK600)
    
   # saveRDS(metab_dat, paste0("metabolism/", filelist[i],".rds"))
    
}

metab_dat <- list(data=all_data, metab=all_metab, specs=all_specs)
    
saveRDS(metab_dat,"metabolism/NHC_modeled.rds")
  




#############################
# extract data from rds
all_metab <- data.frame()
all_data <- data.frame()
all_specs <- data.frame()
for(i in 1:length(filelist)){
    sitename <- substr(filelist[i], 1, nchar(filelist[i])-4)
    mod <- readRDS(paste0("data/metabolism/modeled/",sitename,".rds"))
   
}

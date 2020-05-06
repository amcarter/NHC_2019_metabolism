#########################################
# Run Stream Metabolizer on NHC data
# AMCarter 4-17-20


library(streamMetabolizer)
library(dplyr)

# load RDS from model output:

filelist <- list.files("../data/metabolism/processed")
bayes_name <- mm_name(type="bayes", pool_K600="binned", err_obs_iid = TRUE,
                      err_proc_iid = TRUE, err_proc_acor = FALSE)
bayes_specs <- specs(bayes_name)

for(i in 1:length(filelist)){
  dat <- read_csv(paste0("../data/metabolism/processed/",filelist[1]))
  # set nodes for stream metabolizer to bin discharge around
  Qrange <- c(min(log(dat$discharge), na.rm=T),quantile(log(dat$discharge), .98, na.rm=T))
  bayes_specs$K600_lnQ_nodes_centers <- seq(Qrange[1], Qrange[2], length=7)
  
  #fit metab model
  mm <- metab(bayes_specs, data=dat)
  writeRDS(mm, paste0("../data/metabolism/modeled/",sitename,"_metabolism.rds"))
}

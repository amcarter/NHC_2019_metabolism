############################################
# Estimate gas exchange from Raymond et al 2012
# Code from M Vlah

# edited A Carter
# 07-06-2020

library(dplyr)
library(readr)
library(lubridate)

site_deets <- read_csv(paste0(metab_projdir, '/data/siteData/NHCsite_metadata.csv')) 
dat <- read_csv(paste0(metab_projdir,'/data/metabolism/processed/CBP.csv'))
site_deets <- site_deets[1:7,] %>% select(sitecode, S=slope, W=width_m)
dat <- format_data_for_lake_metabolizer(dat)

calc_raymond_K600_eq7 <- function(dat, site_deets, site){
    dat$S <- site_deets$S[site_deets$sitecode==site]
    dat$V = dat$discharge/(dat$z.mix*site_deets$W[site_deets$sitecode==site])
    dat$D = dat$z.mix
    dat$Q = dat$discharge
    
    #get min and max for all terms in raymond eqn 7 (from table 2)
    dat$t1a = (dat$V * dat$S)^(0.86 + 0.016)
    dat$t1b = (dat$V * dat$S)^(0.86 - 0.016)
    dat$t2a = dat$Q^(-0.14 + 0.012)
    dat$t2b = dat$Q^(-0.14 - 0.012)
    dat$t3a = dat$D^(0.66 + 0.029)
    dat$t3b = dat$D^(0.66 - 0.029)
    dat$t1min = apply(dat[,c("t1a", "t1b")], 1, min)
    dat$t1max = apply(dat[,c("t1a", "t1b")], 1, max)
    dat$t2min = apply(dat[,c("t2a", "t2b")], 1, min)
    dat$t2max = apply(dat[,c("t2a", "t2b")], 1, max)
    dat$t3min = apply(dat[,c("t3a", "t3b")], 1, min)
    dat$t3max = apply(dat[,c("t3a", "t3b")], 1, max)
    
    #calculate k600 (m/d)
    dat$raymondk600min = (4725 - 445) * dat$t1min * dat$t2min * dat$t3min
    dat$raymondk600max = (4725 + 445) * dat$t1max * dat$t2max * dat$t3max
    dat$raymondk600mean = 4725 * (dat$V * dat$S)^(0.86) * dat$Q^(-0.14) * dat$D^(0.66)
    
    #convert to K600 (1/d)
    K600s = data.frame(
        apply(select(dat, starts_with('raymond')), 2, function(x) x / dat$D))
    colnames(K600s) = sub('k', 'K', colnames(K600s))
    dat = cbind(dat$datetime, K600s)
}

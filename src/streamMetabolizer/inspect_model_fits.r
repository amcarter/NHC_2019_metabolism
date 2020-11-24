# assess the fit and K values for SM metabolism model runs
# 11/17/2020
library(ks)
library(zoo)
library(tidyverse)
#setwd('C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism/data')


plot_rhats <- function(fit){
  rh <- get_fit(fit)$daily %>%
    select(date, ends_with('Rhat'))
  ylim = range(c(rh$K600_daily_Rhat, rh$GPP_daily_Rhat, rh$ER_daily_Rhat, 1.05), 
               na.rm = T)
  plot(x = rh$date, y = rh$K600_daily_Rhat, 
       ylab = "rhat (convergence metric)", xlab = "date",
       type = "l", lwd = 1.5, ylim = ylim)        
  lines(rh$date, rh$GPP_daily_Rhat, col = "forestgreen", lwd = 1.5)
  lines(rh$date, rh$ER_daily_Rhat, col = "sienna", lwd = 1.5)
  abline(h = 1.05, lty = 2, col = "red")
  mtext("Rhat below 1.05 is good", 3, 0, adj = 0, cex = .8)
  legend("topleft", 
         legend = c("K600", "GPP", "ER"),
         col = c(1, "forestgreen", "sienna"),
         lty = 1, bty = "n", lwd = 1.5)
}
plot_binning <- function(fit){
  mm_fit <- get_fit(fit)
  
  SM_output <- mm_fit$daily
  SM_KQbin <-  mm_fit$KQ_binned
  SM_day <- get_data_daily(fit)
  SM_specs <- get_specs(fit)
  
  day <- data.frame(SM_day$discharge.daily, 
                    SM_output$K600_daily_50pct, 
                    SM_output$GPP_50pct,
                    SM_output$K600_daily_Rhat,
                    rep('daily', dim(SM_output)[1]))
  
  colnames(day)<-c('Q', 'K600', 'GPP','Rhat', 'Group')
  
  # gg<-ggplot(day, aes(x=log(Q), y = GPP, col=Rhat))+
  #   geom_point() +
  #   geom_hline(yintercept = 0.2)
  
  nodes<-data.frame(exp(SM_specs$K600_lnQ_nodes_centers), 
                    exp(SM_KQbin$lnK600_lnQ_nodes_50pct),
                    exp(SM_specs$K600_lnQ_nodes_meanlog))
  colnames(nodes)<-c('Q', 'K600', 'K600_prior')
  pm <- par()$mar
  plot(log(day$Q), day$K600,
       xlab = "logQ", ylab = "K600", 
       col = "grey", pch = 19, cex = 1.5)
  points(log(nodes$Q), nodes$K600, col = "brown3", pch = 19, cex = 1.5)
  points(log(nodes$Q), nodes$K600_prior, col = "brown3", cex = 1.5)
  par(new = T, mar = c(0,0,0,0))
  plot(1,1, type = 'n', axes = FALSE, xlab = '', ylab = '')
  legend("top", legend = c("daily", "prior", "posterior"),
         col = c("grey", "brown3", "brown3"),
         pch = c(19, 1, 19), cex = 1,  bty = 'n', ncol = 3)
  par(mar = pm)
}
plot_KvER <- function(fit){
  KvER <- get_fit(fit)
  pcor <- round(cor(KvER$daily$K600_daily_mean, 
                    KvER$daily$ER_daily_mean, 
                    use = "na.or.complete"),2)
  plot(-KvER$daily$ER_daily_mean, KvER$daily$K600_daily_mean,
       xlab = "ER (gO2/m2/d)", ylab = "K600 (/d)")
  mtext(paste0("pearson's correlation = ", pcor),
        side = 3, line = 0, adj = 1, cex = .8)
}
plot_metab <- function(met, ylim = NULL, doy = F){

  yrange = range(c(met$GPP.upper, met$ER.lower), na.rm = T)
  
  if(!is.null(ylim)){yrange = ylim}
  if(doy == T){
    met <- met %>%
      mutate(date = doy) %>%
      arrange(date)
    }
  plot(met$date, met$GPP, 
       type = "l", lwd = 2, col = "forestgreen",
       ylim = yrange, xlab = "date", ylab = "gO2/m2/d")  
  lines(met$date, met$ER, 
        lwd = 2, col = "sienna")
  polygon(na.approx(c(met$date, rev(met$date)), na.rm = F), 
          na.approx(c(met$GPP.lower, rev(met$GPP.upper)), na.rm = F),
          col = alpha("forestgreen", 0.4), border = NA)
  polygon(na.approx(c(met$date, rev(met$date)), na.rm = F), 
          na.approx(c(met$ER.lower, rev(met$ER.upper)), na.rm = F),
          col = alpha("sienna", 0.4), border = NA)
  abline(h = 0)
  
}
plot_hall_metab <- function(met, ylim = NULL, 
                            site = c("CBP", "WB", "UNHC"), doy = F) {

  ss <- data.frame(hall_site = c("Concrete", "Wood Bridge", "Blackwood"),
                   site = c("CBP", "WB", "UNHC"))
  if(!site %in% ss$site){
    print(paste(site, "not in Hall data"))
    break
  }
  
  hall <- read_csv("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/hall_50yl/data/hall/hall_table_15.csv") %>%
    rename(hall_site = site) %>%
    mutate(ER_gO2m2d = -ER_gO2m2d) %>%
    left_join(ss) %>%
    select(-hall_site, -newdate) %>%
    filter(site %in% !!site) %>%
    mutate(doy = format(as.Date(date, format = "%m/%d/%Y"), "%j")) %>%
    select(-date)
  
  if(!"doy" %in% colnames(met)){
    met <- met %>%
      mutate(doy = format(date, "%j"))
  }
  
  met <- met %>%
    full_join(hall, by = "doy")
  
  if(is.null(ylim)){
    ylim <- range(c(met$GPP.upper, met$ER.lower, 
                    met$GPP_gO2m2d, met$ER_gO2m2d), na.rm = T)
  }
  
  plot_metab(met, ylim, doy = doy)
  
  if(doy == T){
    met <- met %>%
      mutate(date = doy) %>%
      arrange(date)
  }
  
  points(met$date, met$GPP_gO2m2d, pch = 19, col = "forestgreen")
  points(met$date, met$ER_gO2m2d, pch = 19, col = "sienna")
}
plot_kde_metab <- function(met, lim = NULL){
  
  kernel <- kde(na.omit(met[,c("GPP","ER")]))
  if(is.null(lim)) {
    lim <- quantile(c(met$GPP, -met$ER), .99, na.rm = T) 
  } 
  
  plot(kernel, xlab = "GPP (gO2m2d)", ylab = "ER (gO2m2d)", 
       ylim = c(-lim, 0), xlim = c(0, lim), 
       display = "filled.contour",
       cont=c(30,60,90), #lwd = 1,
       col = c(NA, 
               alpha("grey25", .25), 
               alpha("grey25", .5), 
               alpha("grey25", .75)))
       
  abline(0,-1)
}
plot_kde_hall_metab <- function(met, lim = NULL, 
                                site = c("CBP", "WB", "UNHC")){
  
  ss <- data.frame(hall_site = c("Concrete", "Wood Bridge", "Blackwood"),
                   site = c("CBP", "WB", "UNHC"))
  if(!site %in% ss$site){
    print(paste(site, "not in Hall data"))
    break
  }
  
  hall <- read_csv("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/hall_50yl/data/hall/hall_table_15.csv") %>%
    rename(hall_site = site) %>%
    mutate(ER_gO2m2d = -ER_gO2m2d) %>%
    left_join(ss) %>%
    select(-hall_site, -newdate) %>%
    filter(site %in% !!site)
    
  kernel_hall <- kde(na.omit(hall[,c("GPP_gO2m2d","ER_gO2m2d")]))
  
  if(is.null(lim)) {
    lim <- quantile(c(hall$GPP_gO2m2d, -hall$ER_gO2m2d, 
                   met$GPP, -met$ER), .99, na.rm = T) 
  }
  
  plot_kde_metab(met, lim)
  par(new=T)
  plot(kernel_hall, 
       display = "filled.contour",
       xlab = "", ylab = "", 
       ylim = c(-lim, 0), xlim = c(0, lim),
       cont = c(30,60,90), 
       #lwd = 1,
       col = c(NA, 
               alpha("darkred", .25),                 
               alpha("darkred", .5), 
               alpha("darkred", .75)))
  abline(0,-1)
  legend("topright", cex = 1.4,
         c(paste0("Today  n = ",nrow(na.omit(met))),
           paste0("Hall  n = ", nrow(na.omit(hall)))),
         fill = c(alpha("grey25", .75), alpha("darkred", .75)), 
         border = NA, bty = "n")
   
}
plot_diagnostics <- function(fit, site, ylim = NULL, lim = NULL){
  m <- rbind(c(1,1,1,1,2,2),
             c(3,3,4,4,5,5))
  layout(m)
  par(mar = c(4,4,2,2))
  met = predict_metab(fit)
  
  plot_hall_metab(met, ylim = ylim)
  plot_kde_hall_metab(met, lim = lim)
  plot_rhats(fit)
  plot_binning(fit)
  plot_KvER(fit)
  mtext(site, outer = T, line = -1.5)
}

#site <- "CBP"

# Look at the metabolism predictions and fits for each year of data
# fit <- readRDS(paste0("metabolism/modeled/", site, "_nreg_v1.rds"))
# plot_metab_preds(predict_metab(fit))
# plot_binning(fit)
# plot_rhats(fit)
# plot_KvER(fit)
# 
# plot_DO_preds(fit, style = "dygraphs", y_var = "conc")
# 
# get_fit(fit)$overall %>%
#   select(ends_with('Rhat'))


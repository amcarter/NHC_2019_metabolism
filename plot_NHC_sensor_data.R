###########################################################
# Make plots of NHC sensor data from 2019-2020
# AMC 20200907
#
#  1. Read in processed NHC sensor timeseries
#  2. Plot individual site variables
#  3. Plot 3d plots with z axis as distance

library(tidyverse)
library(lubridate)
library(plotly)
library(zoo)
library(ggplot2)
setwd(metab_projdir)

# site metadata:
mdat <- read_csv("data/siteData/NHCsite_metadata.csv")
sites <- mdat[c(1:5,7),] %>%
  select(site=sitecode, distance_m)

sites$distance_m<- c(0,1550,3450,5950,6100,8450)

# make a df with all NHC sensor data ####
ss <- "NHC"
dat <- read_csv(paste0("data/metabolism/processed/",ss, ".csv"), guess_max = 100000)
dat$datetime <- with_tz(dat$DateTime_UTC, "EST")

# calculate wrt and velocity using coarse crossection data
#   cross section is level * measured width in march
#   velocity is discharge/cross section
#   MRT for 100m reach is 100/velocity/60 s/ 60 m

alldat <- dat %>%
  mutate(spc_uScm = 1000*SpecCond_mScm) %>%
  select(datetime, DO_mgL, DO_sat = DO.sat, spc_uScm, temp_C = WaterTemp_C, discharge_m3s = Discharge_m3s, depth, level_m) %>%
  mutate(
    site = ss, 
    distance_m = sites$distance_m[sites$site==ss],
    width_m = mdat$width_mar_m[mdat$sitecode==ss],
    length_m = mdat$habitat_length_m[mdat$sitecode==ss],
    velocity_ms = discharge_m3s/(width_m*level_m),
    WRT_hr = 100/velocity_ms/60/60)


for(ss in sites$site[-1]){
  dat <- read_csv(paste0("data/metabolism/processed/",ss, ".csv"), guess_max = 100000)
  dat$datetime <- with_tz(dat$DateTime_UTC, "EST")
  if(ss=="UNHC"){
    dat$SpecCond_uScm <- NA
  }
  dat <- dat %>%
      select(datetime, DO_mgL, DO_sat = DO.sat, spc_uScm=SpecCond_uScm, 
             temp_C = WaterTemp_C, discharge_m3s = Discharge_m3s, 
             depth, level_m) %>%
    mutate(
      site = ss, 
      distance_m = sites$distance_m[sites$site==ss],
      width_m = mdat$width_mar_m[mdat$sitecode==ss],
      length_m = mdat$habitat_length_m[mdat$sitecode==ss],
      velocity_ms = discharge_m3s/(width_m*level_m),
      WRT_hr = 100/velocity_ms/60/60)
  alldat <- bind_rows(alldat, dat)
}

alldat<- alldat[-which(is.na(alldat$datetime)),]
# summarize daily values of dataset ####
ddat <- alldat %>%
  mutate(
    date = as.Date(substr(datetime, 1, 10)),
    DO_psat = DO_mgL/DO_sat) %>%
  select(-datetime) %>%
  group_by(date, site) %>%
  summarize_all(mean, na.rm=TRUE) %>%
  mutate(DOY = format(date, "%j")) %>%
  ungroup() 
ddat <- left_join(ddat,sites[,c(1,3)])
ddat$site <- factor(ddat$site, levels=sites$site)
ddat$DOY <- as.numeric(ddat$DOY)
# plot in 2D panels to see if the data are all there

# DO plot at all of the sites
ggplot(data=ddat, aes(x=date, y=DO_psat)) +
  geom_path(aes(color=site), lwd=1) +
  scale_color_grey() +
  theme_classic()

# DO "contors" with y as distance and x variable

# Mikes code ####
# 
zz <- t(apply(z_mat, 1, na.approx, na.rm = FALSE))
colnames(zz) <- colnames(z_mat)
# # zz <- z_var
# # zz[is.na(zz)] <- 0
# # zz <- zz %>%
# #   as_tibble() %>%
# #   group_by(distance_m, DOY) %>%
# #   summarize(DO_psat = mean(DO_psat, na.rm = TRUE)) %>%
# #   ungroup()
# 
library(viridis)
ncolors <- 12

facmat <- cut(as.vector(zz),
              breaks = ncolors,
              labels = viridis(ncolors)) %>%
  as.character() %>%
  matrix(nrow = nrow(zz),
         ncol = ncol(zz))

# facmat <- matrix(cut(as.vector(zz),
#                      breaks = ncolors,
#                      labels = viridis(ncolors)),
#                  nrow = nrow(zz),
#                  ncol = ncol(zz))
zz<-zz[,-171]
xlims <- range(as.numeric(colnames(zz)))
ylims <- c(1, nrow(zz))
zrow <- nrow(zz)
plot(rep(0, zrow), seq(1:365), pch=15, col=facmat[, 1], xlim = xlims,
     ylim = ylims)
for(i in 1:ncol(zz)){
  dep <- colnames(zz)[i]
  points(rep(as.numeric(dep), zrow), seq(1:365), pch=15, col=facmat[, i])
}
# 
# pivot_longer(zz, )
#   # stat_contour(bins = 100, na.rm = TRUE)
# # sum(apply(z_var, 1, function(x) any(is.na(x))))
# plot(zz$distance_m, zz$DOY, pch=15, col=as.character(cutfac))
# 
# plot(ddat$distance_m, ddat$discharge_m3s, log="y", type="l")


# make contour plots of daily data ####

# make table of just the x,y,z you want
ddat$ln_discharge_m3s <- log(ddat$discharge_m3s)

generate_z_mat <- function(dat, x="DOY", y="distance_m", z, 
                           x_vals=NA, y_vals=NA, smooth=F, n=3){
  z_var <- dat %>%
    select(!!x, !!y, !!z) %>%
    group_by(!!sym(x), !!sym(y)) %>%
    summarize_all(mean, na.rm=T) %>%
    ungroup() %>%
    pivot_wider(names_from = !!y, values_from = !!z)
  z_var <- z_var[, c(x, as.character(sort(as.numeric(colnames(z_var[-1])))))]
  
  # set x and y values  
  if(x=="DOY"){x_vals <- 1:365}
  if(y=="distance_m"){ y_vals <- seq(0,8450,by=50)}
  
  if(is.na(x_vals)){
    xrange <- range(dat[,x], na.rm=T)
    x_vals <- seq(xrange[1], xrange[2], length.out=100)
    x_breaks <- seq((x_vals[1]-diff(x_vals)[1]/2),
                   (x_vals[100]+diff(x_vals)[1]/2),
                   by=diff(x_vals)[1])
    zz <-data.frame(ln_discharge_m3s = x_vals,
                    fac = cut(x_vals,x_breaks))
    breaks <- unique(as.numeric(c(sub("\\((.+),.*", "\\1", zz$fac), 
                                  sub("[^,]*,([^]]*)\\]", "\\1", zz$fac))))
    for(col in colnames(z_var[,-1])){
      tmp<-z_var[,c(x,col)]
      tmp$fac<- cut(pull(tmp[,1]),breaks)
      tt <- tmp[,-1] %>% group_by(fac) %>%
        summarize_all(mean, na.rm=T) %>%
        ungroup()
      zz <- left_join(zz,tt)  
    }
    z_var <- zz[,-2] 
  }
  

  if(smooth==T){
    z_filt <- bind_rows(z_var[(nrow(z_var)-n+1):nrow(z_var),], z_var, z_var[1:n,] )
    for(i in 2:ncol(z_filt)){
      z_filt[,i] <- na.approx(z_filt[,i], na.rm=F) 
      z_filt[,i] <- stats::filter(z_filt[,i],rep(1/n,n), sides=2)
    }
    z_var <- z_filt[(n+1):(nrow(z_filt)-n),]
  }
  
  z_mat <- matrix(data=NA, nrow=length(x_vals), ncol = length(y_vals), 
                  dimnames=list(NULL,as.character(y_vals)))

  for(col in colnames(z_var[,-1])){
    z_mat[,col] <- (z_var[,col]) 
  }

  z_mat <- as.data.frame(t(apply(z_mat, MARGIN = 1, na.approx, na.rm = F)))
  
  colnames(z_mat) <- y_vals
  z_mat$X <- x_vals
  return(z_mat)
}

z_mat <- generate_z_mat(ddat, x="temp_C",z = "spc_uScm")
z_mat<- rename(z_mat, temp_C = X)
zz <- pivot_longer(z_mat,cols = -temp_C,
                    names_to = "distance_m", values_to = "Sp_Cond")
zz$distance_m <- as.numeric(zz$distance_m)
 

g_spcT<- ggplot(data=zz, aes(y=distance_m, x=temp_C, fill=Sp_Cond)) +
  geom_raster()+
     scale_fill_gradientn(colors =viridis(5))

figure <-  ggarrange(g_DOQ, g_tempQ, g_spcQ,
                      ncol = 1, nrow = 3)
ggsave("spaceQ_contours.pdf",figure, path="figures", width=8.5,height=11)

# DO vs DOY ####
z_var <- ddat %>%
  select(DOY, distance_m, DO_psat) %>%
  group_by(DOY, distance_m) %>%
  summarize_all(mean, na.rm=T) %>%
  ungroup() %>%
  pivot_wider(names_from = distance_m, values_from = DO_psat)
z_var <- z_var[,c("DOY",as.character(sort(as.numeric(colnames(z_var[-1])))))]

#z_var$DO_psat[is.nan(z_var$DO_psat)] = NA

z_filt <- bind_rows(z_var[(nrow(z_var)-2):nrow(z_var),], z_var, z_var[1:3,] )
for(i in 2:ncol(z_filt)){
  z_filt[,i] <- na.approx(z_filt[,i], na.rm=F) 
  z_filt[,i] <- stats::filter(z_filt[,i],rep(1/3,3), sides=2)
}

z_filt <- z_filt[4:(nrow(z_filt)-3),-1]
dist_vec <- seq(0,8450,by=50)

z_mat <- matrix(data=NA, nrow=nrow(z_filt), ncol = length(dist_vec), 
                dimnames=list(NULL,as.character(dist_vec)))

for(col in colnames(z_filt)){
  z_mat[,col] <- pull(z_filt[,col]) 
}

z_mat <- as.data.frame(t(apply(z_mat, MARGIN = 1, na.approx, na.rm = F)))
colnames(z_mat) <- dist_vec
z_mat$DOY <- seq(1,365)
z <- pivot_longer(z_mat, cols = -DOY, 
                   names_to="distance_m", values_to = "DO_psat")
z$distance_m <- as.numeric(z$distance_m)


g_DO <- ggplot(data=z, aes(y=distance_m, x=DOY, color=DO_psat)) +
  geom_point() +
  scale_color_gradientn(colors=viridis(5))

library(ggpubr)
figure <- ggarrange(g_DO, g_temp, g_spc,
                    #labels = c("A", "B", "C"),
                    ncol = 1, nrow = 3)

ggsave("spaceTime_contours.pdf",figure, path="figures", width=8.5,height=11)
 
# 
# DO_DOY_plt <- plot_ly(z=z, type="surface") %>%
#   layout(
#     title = 'NHC Dissolved Oxygen',
#     scene = list(
#       xaxis = list(title = 'position',
#                    tickmode = 'array',
#                    tickvals = 0:6,
#                    ticktext = colnames(y)),
#       yaxis = list(title = 'DOY'),
#       zaxis = list(title = 'DO (% sat)')
#     ))
#     
# DO_DOY_plt
# 
# 
# 

# 
# # DO vs MRT ####
# y_var <- ddat %>%
#   select(WRT_hr, distance_m, DO_psat) %>%
#   group_by(WRT_hr, distance_m) %>%
#   summarize_all(mean, na.rm=T) %>%
#   ungroup() %>%
#   pivot_wider(names_from = distance_m, values_from = DO_psat)
# y_var <- y_var[,c("WRT_hr",as.character(sort(as.numeric(colnames(y_var[-1])))))]
# 
# for(i in 2:ncol(y_var)){
#   y_var[,i]<- na.approx(y_var[,i], na.rm=F) 
# }
# 
# 
# x <- as.numeric(y_var$WRT_hr)
# z <- as.numeric(colnames(y_var[,-1]))
# y <- as.matrix(y_var[,-1])
# 
# 
# 
# DO_WRT_plt <- plot_ly(z=y, type="surface") %>%
#   layout(
#     title = 'NHC Dissolved Oxygen',
#     scene = list(
#       xaxis = list(title = 'position',
#                    tickmode = 'array',
#                    tickvals = 0:6,
#                    ticktext = colnames(y)),
#       yaxis = list(title = 'WRT (hr)',
#                    tickmode = 'array',
#                    tickvals = 0:2,
#                    ticktext =c(0,y[nrow(y),1])),
#       zaxis = list(title = 'DO (% sat)')
#     ))
# 
# DO_WRT_plt
# 


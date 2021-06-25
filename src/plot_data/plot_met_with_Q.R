library(pracma)
# setwd("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism/data")
source("NHC_2019_metabolism/src/streamMetabolizer/inspect_model_fits.r")


dat <- readRDS("NHC_2019_metabolism/data/metabolism/compiled/met_preds_stream_metabolizer.rds")

preds_19 <- dat$preds %>% 
  filter(year == 2019)
preds_nhc <- dat$preds %>%
  filter(site %in% c("NHC", "UNHC"))
preds <- dat$preds %>% filter(date > as.Date('2017-01-01'))
qq <- read_csv('NHC_2019_metabolism/data/rating_curves/interpolatedQ_allsites.csv')
plot(qq$DateTime_UTC, qq$NHC.Q, log = 'y', type = 'l')
qq <- qq %>%
  mutate(date = as.Date(with_tz(DateTime_UTC, tz = "EST")))%>%
  group_by(date) %>%
  select(-DateTime_UTC, -notes, -PWC.Q) %>%
  summarize_all(mean, na.rm = T)

# plot annual met with discharge


sites$sscode <- c("NHC_8.5", "NHC_6.9", "NHC_5", "NHC_2.5", "NHC_2.3", "NHC_0")

  ylim = c(-10,4)
  Qlim = c(.01, max(qq$NHC.Q, na.rm = T) * 5000)
  xlim = c(as.Date(c("2019-03-06", "2020-03-20")))
png("figures/metQ_across_sites_SM_2019.png",
    width = 10, height = 6, units = 'in', res = 300)
  par(oma = c(4,4,4,4),
      mfrow = c(3,2),
      mar = c(0,0,0,0))
  for(s in c(6,3,5,2,4,1)){
    ss <- sites$sitecode[s]
    sn <- sites$sscode[s]
    tmp <- preds_19 %>%
      filter(site == ss)
    tmpq <- qq[,c(1, s+1)] %>%
      filter(date >= xlim[1], 
             date <= xlim[2])
    colnames(tmpq) <- c('date', 'discharge')
    plot_metab(tmp, ylim = ylim, xlim = xlim, xaxt = 'n', yaxt = 'n')
    if(s %in% c(4,5,6)) {
      axis(2, at = c(-6, -3, 0, 3))
    } 
    if(s %in% c(1,4)){
      axis(1, at = seq(as.Date('2019-03-01'), by = 'month', length.out = 13),
           labels = month.abb[c(3:12, 1:3)])
      # axis(1, at = seq(as.Date('2019-03-01'), by = '2 months', length.out = 7),
      #      labels = month.abb[c(3,5,7,9,11,1,3)])
    } 
    par(new = T)
    plot(tmpq$date, tmpq$discharge, log = "y", type = "l", ylim = Qlim, xlim = xlim,
         axes = FALSE, xlab = '', ylab = '')
    mtext(sn, cex = 1, line = -1.5, adj = 0.02)
    if(s %in% c(1,2,3)){
      axis(4, at = c(0.01, 1, 100), labels = c(0.01, 1, 100))
    } 
  }
  par(new = T, mfrow = c(1,1), oma = c(4,4,0,4))
  plot(1,1, axes = F, ann = F, type = 'n')
  mtext("gC/m2/d", 2, 2.8)
  mtext("Q (m3/s)", 4, 2.8)
  # mtext("Date", 1, 2.5)
  mtext("Metabolism across Sites in 2019", 3, -1)
dev.off()

png("figures/metQ_across_years_SM_2019.png",
    width = 10, height = 6, units = 'in', res = 300)
  par(oma = c(4,4,4,4),
      mfrow = c(3,2),
      mar = c(0,0,0,0))
  florence <- as.Date("2018-09-14")
  xlim = c(as.Date(c("2017-03-01", "2018-03-01")))
  for(y in c(2017,2018,2019)){
    tmp <- preds_nhc %>%
      filter(site == "UNHC",
             year == y)
    tmpq <- qq[,c(1, 7)] %>%
      filter(date >= xlim[1], 
             date <= xlim[2])
    colnames(tmpq) <- c('date', 'discharge')
    plot_metab(tmp, ylim = ylim, xlim = xlim, xaxt = 'n', yaxt = 'n')
    axis(2, at = c(-6, -3, 0, 3))
    abline(v = florence, lty = 2)
    if(y == 2017){ mtext("NHC_0",3, line = .5) }
    if(y == 2019){
      axis(1, at = seq(as.Date('2019-03-01'), by = 'month', length.out = 13),
           labels = month.abb[c(3:12, 1:3)])
    } 
    par(new = T)
    plot(tmpq$date, tmpq$discharge, log = "y", type = "l", ylim = Qlim, xlim = xlim,
         axes = FALSE, xlab = '', ylab = '')
    mtext(y, cex = 1, line = -1.5, adj = 0.02)

    tmp <- preds_nhc %>%
      filter(site == "NHC",
             year == y)
    tmpq <- qq[,c(1, 2)] %>%
      filter(date >= xlim[1], 
             date <= xlim[2])
    colnames(tmpq) <- c('date', 'discharge')
    plot_metab(tmp, ylim = ylim, xlim = xlim, xaxt = 'n', yaxt = 'n')
    abline(v = florence, lty = 2, )
    if(y == 2017){ mtext("NHC_8.5",3, line =.5) }
    if(y == 2019){
      axis(1, at = seq(as.Date('2019-03-01'), by = 'month', length.out = 13),
           labels = month.abb[c(3:12, 1:3)])
    } 
    par(new = T)
    plot(tmpq$date, tmpq$discharge, log = "y", type = "l", ylim = Qlim, xlim = xlim,
         axes = FALSE, xlab = '', ylab = '')
    mtext(y, cex = 1, line = -1.5, adj = 0.02)
    axis(4, at = c(0.01, 1, 100), labels = c(0.01, 1, 100))
    xlim = xlim + 365
  }
  par(new = T, mfrow = c(1,1), oma = c(4,4,0,4))
  plot(1,1, axes = F, ann = F, type = 'n')
  mtext("gC/m2/d", 2, 2.8)
  mtext("Q (m3/s)", 4, 2.8)
  # mtext("Date", 1, 2.5)
  mtext("Metabolism across years", 3, -1)
dev.off()


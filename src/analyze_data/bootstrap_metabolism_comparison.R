# Bootstrap Metabolism comparison NHC SP data and Hall data #####
library(viridis)
library(beanplot)
library(scales)

#setup ####
#historic data 
dat <- readRDS("NHC_2019_metabolism/data/metabolism/compiled/met_preds_direct_calc.rds")
nhc_68_70 = dat$preds %>%
    filter(era == "then", 
           site == "CBP") %>%
    select(-era)

nhc_68 <- nhc_68_70 %>%
    filter(year == 1968)
nhc_69 <- nhc_68_70 %>%
    filter(year == 1969)

gpp_68_70 = nhc_68_70$GPP
gpp_68 = nhc_68$GPP
gpp_69 = nhc_69$GPP

er_68_70 = nhc_68_70$ER
er_68 = nhc_68$ER
er_69 = nhc_69$ER

nep_68_70 = gpp_68_70 + er_68_70
nep_68 = gpp_68 + er_68
nep_69 = gpp_69 + er_69

# contemporary data
nhc_new <- dat$preds %>%
    filter(era == "now",
           site == "CBP")

gpp_new = nhc_new$GPP
er_new = nhc_new$ER
nep_new = gpp_new + er_new

dates_new = nhc_new$date

#time-series comparison of means then and now? ####
plot(gpp_new, type='l')
acf(gpp_new, na.action=na.pass)
pacf(gpp_new, na.action=na.pass)
#strong autocorrelation and partial autocorr;

#will have to model error as an autoregressive process
#not stationary; can't use pure AR
qqnorm(er_new); abline(1, 1, col='red', lty=2)
#normalish; no need to go bayesian

#nhc_17 is irregular, so can't use arima; only option would be GAM

#distribution plots ####
png(width=9, height=6, units='in', type='cairo', res=300,
    filename='figures/metab_distributions_v2.png')

defpar = par(mfrow=c(2,3))

#plot GPP dists, then and now
plot(density(gpp_68_70, na.rm=TRUE), bty='l', col='sienna3',
     main='GPP 1968-70 vs. 2019', xlab='GPP')
lines(density(gpp_new, na.rm=TRUE), col='blue')
legend('topright', 
       legend=c('68-70; n=76', 
                paste0('19; n=', length(which(!is.na(gpp_new))))),
       col = c('sienna3','blue'), lty = 1, bty = 'n',
       seg.len = 1, cex = 0.9, lwd = 2)

#plot ER dists, then and now
plot(density(er_68_70, na.rm=TRUE), bty='l', col='sienna3',
     main='ER 1968-70 vs. 2019', xlab='ER')
lines(density(er_new, na.rm=TRUE), col='blue')
legend('topleft', 
       legend=c('68-70; n=76', 
                paste0('19; n=', length(which(!is.na(er_new))))),
       col = c('sienna3','blue'), 
       lty = 1, bty = 'n', seg.len = 1, cex = 0.9, lwd = 2)

#plot NEP dists, then and now
plot(density(nep_68_70, na.rm=TRUE), bty='l', col='sienna3',
    main='NEP 1968-70 vs. 2019', xlab='NEP', ylim = c(0,4.5))
lines(density(nep_new, na.rm=TRUE), col='blue')
legend('topleft', 
       legend=c('68-70; n=76', 
                paste0('19; n=', length(which(!is.na(nep_new))))),
       col = c('sienna3','blue'), 
       lty = 1, bty = 'n', seg.len = 1, cex = 0.9, lwd = 2)

#plot GPP dists by year
# cols = viridis(6)
cols = c(rep('sienna3', 2), 'blue')
plot(density(gpp_68, na.rm=TRUE),bty='l', col=cols[1],
    main='GPP by year', xlab='GPP', ylim = c(0,3), xlim = c(-.1, 1))
lines(density(gpp_69, na.rm=TRUE), col=cols[2])
lines(density(gpp_new, na.rm=TRUE), col=cols[3])
legend('topright',
    legend=c(paste0('68; n=', length(gpp_68)),
             paste0('69; n=', length(gpp_69)),
             paste0('19; n=', length(which(!is.na(gpp_new))))),
    col = cols, lty = 1, bty = 'n', 
    seg.len = 1, cex = 0.9, lwd = 2)

#plot ER dists by year
plot(density(er_68, na.rm=TRUE), bty='l', col=cols[1],
    main='ER by year', xlab='ER', xlim = c(-1.2, .2))
lines(density(er_69, na.rm=TRUE), col=cols[2])
lines(density(er_new, na.rm=TRUE), col=cols[3])
legend('topleft',
    legend=c(paste0('68; n=', length(gpp_68)),
             paste0('69; n=', length(gpp_69)),
             paste0('19; n=', length(which(!is.na(er_new))))), 
    col = cols, lty = 1, bty = 'n', 
    seg.len = 1, cex = 0.9, lwd = 2)

#plot NEP dists by year
plot(density(nep_68, na.rm=TRUE), bty='l', col=cols[1],
    main='NEP by year', xlab='NEP', ylim = c(0, 4.5), xlim = c(-.8, .2))
lines(density(nep_69, na.rm=TRUE), col=cols[2])
lines(density(nep_new, na.rm=TRUE), col=cols[3])
legend('topleft',
    legend = c(paste0('68; n=', length(gpp_68)),
               paste0('69; n=', length(gpp_69)),
               paste0('19; n=', length(which(!is.na(nep_new))))),
    col=cols, lty=1, bty='n', seg.len=1, cex=0.9, lwd=2)

dev.off()
               

#plot temporal coverage for historic data ####

historic_dates = nhc_68_70$date[! is.na(nhc_68_70$GPP)]
historic_year_agg = as.character(historic_dates)
substr(historic_year_agg, 1, 4) = '2020'
historic_year_agg = as.Date(historic_year_agg)
hy_num = as.numeric(historic_year_agg)

new_dates = nhc_new$date[!is.na(nhc_new$GPP)]
new_year_agg = as.character(new_dates)
substr(new_year_agg, 1, 4) = '2020'
new_year_agg = as.Date(new_year_agg)
ny_num = as.numeric(new_year_agg)


png(width=7, height=5, units='in', type='cairo', res=300,
    filename='figures/seasonalcoverage.png')
lims = c(min(ny_num), max(ny_num))

par(mfrow=c(2,1), mar=c(0,0,0,0), oma=c(3,4,3,0))

beanplot(ny_num, horizontal = TRUE, col = 'gray', xaxt = 'n', 
         frame.plot = FALSE, ylim = lims)
mtext('2019', 2, cex = 1.5)
mtext('CBP Annual Coverage', 3)
legend('topright', legend = paste('n =', length(! is.na(ny_num))),
       bty = 'n', cex = 1.3, text.font = 2)
beanplot(hy_num, horizontal = TRUE, col='brown3', xaxt = 'n',
       frame.plot = FALSE, ylim = lims)
mtext('1969', 2, cex = 1.5)
legend('topright', legend = paste('n =', length(!is.na(hy_num))),
       bty = 'n', cex = 1.3, text.font = 2)
axis(1, at=seq(as.Date('2020-01-01'), as.Date('2020-12-31'), length.out=13)[1:12],
    labels=month.abb)
dev.off()

#compare overall distributions then and now ####

#first assess normality
# png(width=7, height=6, units='in', type='cairo', res=300,
#     filename='~/Dropbox/streampulse/figs/NHC_comparison/normality_assessment.png')

par(mfrow=c(2,2), mar=c(0, 0, 0, 0), oma=rep(4, 4))
qqnorm(gpp_new, lty=2, xlab='', ylab='', main='', xaxt='n', yaxt='n', bty='o')
abline(0, 1, col='red', lty=2)
legend('bottom', legend='GPP 2019', bty='n', cex=1.3)
qqnorm(er_new, lty=2, xlab='', ylab='', main='', xaxt='n', yaxt='n', bty='o')
abline(0, 1, col='red', lty=2)
legend('bottom', legend='ER 2019', bty='n', cex=1.3)
qqnorm(gpp_68_70, lty=2, xlab='', ylab='', main='', xaxt='n', yaxt='n', bty='o')
abline(0, 1, col='red', lty=2)
legend('top', legend='GPP 1968-70', bty='n', cex=1.3)
qqnorm(er_68_70, lty=2, xlab='', ylab='', main='', xaxt='n', yaxt='n', bty='o')
abline(0, 1, col='red', lty=2)
legend('top', legend='ER 1968-70', bty='n', cex=1.3)
mtext('Theoretical Quantiles', 1, outer=TRUE, line=1.5)
mtext('Sample Quantiles', 2, outer=TRUE, line=1.5)
mtext('Normal Q-Q Plots (red line = 1:1)', 3, outer=TRUE, line=1.5)

# dev.off()

#nonnormal, but CLT probably applies.
#let's assess equality of variance with an F-test
var.test(gpp_68_70, gpp_new) #equal: p = 0.1783
var.test(er_68_70, er_new) #equal: p = 0.7358

# equal variance, so could do a 2-sample t-test, but not gonna right now.
#not i.i.d., so welch's t-test is out (could transform)
#can't do Mann-Whitney-Wilcoxon Test either because of autocorr

#bootstrap 2-samp t-test for GPP (and we'll go with Welch's) ####

#get observed t-statistic
t_obs_gpp = t.test(gpp_68_70, gpp_new, var.equal=FALSE)$statistic

#artificially make both sample means identical (satisfy the null)
gpp_68_70_mod = gpp_68_70 - mean(gpp_68_70, na.rm=TRUE) +
    mean(c(gpp_68_70, gpp_new), na.rm=TRUE)
gpp_new_mod = gpp_new - mean(gpp_new, na.rm=TRUE) +
    mean(c(gpp_68_70, gpp_new), na.rm=TRUE)

#verify
round(mean(gpp_68_70_mod, na.rm=TRUE), 7) ==
    round(mean(gpp_new_mod, na.rm=TRUE), 7)

#get historic monthly data coverage to use for sample weights
nhc_68_70 <- filter(nhc_68_70, !is.na(GPP))
month_counts_68_70 = tapply(rep(1, nrow(nhc_68_70)),
    substr(nhc_68_70$date, 6, 7), sum)
month_proportions = month_counts_68_70 / nrow(nhc_68_70)
#split gpp vector by month for each dataset
#commented portions are for uniformly distributing monthly draw weights
gpp_68_70_bymo <- split(gpp_68_70_mod[!is.na(gpp_68_70_mod)],
                        factor(substr(nhc_68_70$date, 6, 7)))
# gpp_68_70_bymo = split(gpp_68_70_mod,
#     factor(rep(c('01','02','03','04','05','06','07','08','10','11','12'),
#     length.out=length(gpp_68_70_mod))))
gpp_new_bymo = split(gpp_new_mod[!is.na(gpp_new_mod)],
                     factor(substr(dates_new[!is.na(gpp_new_mod)], 6, 7)))
# gpp_new_bymo = split(gpp_new_mod,
#     factor(rep(c('01','02','03','04','05','06','07','08','10','11','12'),
#     length.out=length(gpp_new_mod))))
nsamp_new = sum(sapply(gpp_new_bymo, length))
nsamp_68_70 = length(which(!is.na(gpp_68_70_mod)))


#determine monthly sample sizes for modern dataset; deal with remainders
month_samp_new = month_proportions * nsamp_new
extra_sample_probs = month_samp_new %% 1
month_samp_new = floor(month_samp_new)

#get bootstrap estimate of sampling distribution of the t-stat if H0 is true;
#i.e. bootstrap the null distribution (weight draws by historic monthly coverage)
nsamp = 20000

t_vect_gpp = vector(length=nsamp)
for(i in 1:nsamp){
    samp_68_70_gpp = samp_new_gpp = c()
    remainder_months = sample(c(1:8, 10:12), size=sum(extra_sample_probs),
        prob=extra_sample_probs)
    for(j in c(1:8, 10:12)){
        extra_samp = ifelse(j %in% remainder_months, 1, 0)
        j = sprintf('%02d', j)
        samp_68_70_gpp = append(samp_68_70_gpp, sample(gpp_68_70_bymo[[j]],
            size=month_counts_68_70[j], replace=TRUE))
        samp_new_gpp = append(samp_new_gpp, sample(gpp_new_bymo[[j]],
            size=month_samp_new[j] + extra_samp, replace=TRUE))
    }
    # samp_68_70_gpp = sample(gpp_68_70_mod, size=nsamp_68_70, replace=TRUE)
    # samp_new_gpp = sample(gpp_new_mod, size=nsamp_new, replace=TRUE)
    t_vect_gpp[i] = t.test(samp_68_70_gpp, samp_new_gpp,
        var.equal=FALSE)$statistic
}

#p-val is proportion of times observed t-statistic >= bootstrap t-statistic
pval_gpp = (sum(t_vect_gpp <= t_obs_gpp) + 1) / (nsamp + 1)
if(pval_gpp == 1){
    pval_gpp = (sum(t_vect_gpp >= t_obs_gpp) + 1) / (nsamp + 1)
}

#bootstrap Welch's t-test for ER ####

#get observed t-statistic
t_obs_er = t.test(er_68_70, er_new, var.equal=FALSE)$statistic

#artificially make both sample means identical (satisfy the null)
er_68_70_mod = er_68_70 - mean(er_68_70, na.rm=TRUE) +
    mean(c(er_68_70, er_new), na.rm=TRUE)
er_new_mod = er_new - mean(er_new, na.rm=TRUE) +
    mean(c(er_68_70, er_new), na.rm=TRUE)

#verify
mean(er_68_70_mod, na.rm=TRUE) == mean(er_new_mod, na.rm=TRUE)

#split er vector by month for each dataset
er_68_70_bymo = split(er_68_70_mod[!is.na(er_68_70_mod)],
                      factor(substr(nhc_68_70$date, 6, 7)))
er_new_bymo = split(er_new_mod, factor(substr(dates_new, 6, 7)))
er_new_bymo = lapply(er_new_bymo, na.omit)
nsamp_new = sum(sapply(er_new_bymo, length))
nsamp_68_70 = length(er_68_70_mod[!is.na(er_68_70_mod)])

#determine monthly sample sizes for modern dataset; deal with remainders
month_samp_new = month_proportions * nsamp_new
extra_sample_probs = month_samp_new %% 1
month_samp_new = floor(month_samp_new)

#get bootstrap estimate of sampling distribution of the t-stat if H0 is true;
#i.e. bootstrap the null distribution (weight draws by historic monthly coverage)
nsamp = 20000
t_vect_er = vector(length=nsamp)
for(i in 1:nsamp){
    samp_68_70_er = samp_new_er = c()
    remainder_months = sample(c(1:8, 10:12), size=sum(extra_sample_probs),
        prob=extra_sample_probs)
    for(j in c(1:8, 10:12)){
        extra_samp = ifelse(j %in% remainder_months, 1, 0)
        j = sprintf('%02d', j)
        samp_68_70_er = append(samp_68_70_er, sample(er_68_70_bymo[[j]],
            size=month_counts_68_70[j], replace=TRUE))
        samp_new_er = append(samp_new_er, sample(er_new_bymo[[j]],
            size=month_samp_new[j] + extra_samp, replace=TRUE))
    }
    # samp_68_70_er = sample(er_68_70_mod, size=nsamp_68_70, replace=TRUE)
    # samp_new_er = sample(er_new_mod, size=nsamp_new, replace=TRUE)
    t_vect_er[i] = t.test(samp_68_70_er, samp_new_er,
        var.equal=FALSE)$statistic
}

#p-val is proportion of times observed t-statistic >= bootstrap t-statistic
pval_er = (sum(t_vect_er >= t_obs_er) + 1) / (nsamp + 1)
if(pval_er == 1){
    pval_er = (sum(t_vect_er >= t_obs_er) + 1) / (nsamp + 1)
}
#visualize GPP hypothesis test ####

# png(width=7, height=6, units='in', type='cairo', res=300,
#     filename='~/Dropbox/streampulse/figs/NHC_comparison/bootstrap_welch_t_weighted_filtered.png')

# par(mfrow=c(2,1), mar=c(4,4,1,2), oma=c(0,0,3,0))
# 
# plot(density(t_vect_gpp), xlab='t-value', main='', xlim = c(-6,2))
# qs = quantile(t_vect_gpp, probs=c(0.025, 0.975))
# dd = density(t_vect_gpp)
# ddo = order(dd$x)
# xdens = dd$x[ddo]
# ydens = dd$y[ddo]
# xdens_lt = xdens[xdens <= qs[1]]
# ydens_lt = ydens[xdens <= qs[1]]
# polygon(c(xdens_lt, rev(xdens_lt)), c(ydens_lt, rep(0, length(ydens_lt))),
#     col='lightgreen', border='lightgreen')
# xdens_ut = xdens[xdens >= qs[2]]
# ydens_ut = ydens[xdens >= qs[2]]
# polygon(c(xdens_ut, rev(xdens_ut)), c(ydens_ut, rep(0, length(ydens_ut))),
#     col='lightgreen', border='lightgreen')
# abline(v=t_obs_gpp, lty=2, col='red', lwd=2)
# legend('topleft', legend='GPP', bty='n', text.font=2, cex=1)
# legend('topleft', legend=paste('\np =', round(pval_gpp, 3)), bty='n',
#     text.font=1, cex=1)

#why bimodality in the null dist?
#1. not from skewed historic draw weights; artificially uniformified them to test.
#2. is from skewed modern draw weights; artificially uniformified them to test.
#3. there is multimodality in some of the modern monthly GPP dists.
#if most draws come from one GPP peak or another,
#the t-val may land in one H0 peak or another

#visualize ER hypothesis test ####

# plot(density(t_vect_er), xlim=c(-2, 5), xlab='t-value', main='')
# qs = quantile(t_vect_er, probs=c(0.025, 0.975))
# dd = density(t_vect_er)
# ddo = order(dd$x)
# xdens = dd$x[ddo]
# ydens = dd$y[ddo]
# xdens_lt = xdens[xdens <= qs[1]]
# ydens_lt = ydens[xdens <= qs[1]]
# polygon(c(xdens_lt, rev(xdens_lt)), c(ydens_lt, rep(0, length(ydens_lt))),
#     col='lightgreen', border='lightgreen')
# xdens_ut = xdens[xdens >= qs[2]]
# ydens_ut = ydens[xdens >= qs[2]]
# polygon(c(xdens_ut, rev(xdens_ut)), c(ydens_ut, rep(0, length(ydens_ut))),
#     col='lightgreen', border='lightgreen')
# abline(v=t_obs_er, lty=2, col='red', lwd=2)
# legend('topleft', legend='ER', bty='n', text.font=2, cex=1)
# legend('topleft', legend=paste('\np =', round(pval_er, 3)), bty='n',
#     text.font=1, cex=1)
# 
# mtext("Observed Welch's t-values (red) relative to bootstrapped null dists", 3,
#     outer=TRUE, line=1, font=2, cex=1.3)

# dev.off()

#verify with Mann-Whitney-Wilcoxon Test? ####

#visualize dists again with non-bootstrapped means ####

# png(width=7, height=6, units='in', type='cairo', res=300,
#     filename='~/Dropbox/streampulse/figs/NHC_comparison/means_raw_boxplot.png')
# par(mfrow = c(1,1))
# gppHmean = paste('mean =', round(mean(gpp_68_70, na.rm=TRUE), 2))
# gppCmean = paste('mean =', round(mean(gpp_new, na.rm=TRUE), 2))
# erHmean = paste('mean =', round(mean(er_68_70, na.rm=TRUE), 2) * -1)
# erCmean = paste('mean =', round(mean(er_new, na.rm=TRUE), 2) * -1)
# boxplot(gpp_68_70, gpp_new,  -1* er_68_70, -1*er_new,
#     ylab='', col='gray',
#     names=c('GPP 1968-70', 'GPP 2017-19', 'ER 1968-70', 'ER 2017-19'))
# axis(1, at=1:4, labels=c(gppHmean, gppCmean, erHmean, erCmean),
#     line=1.5, col='transparent', tcl=0, font=2)
# mtext(expression(paste("gm"^"-2" * " d"^"-1")), 2, line=2)
# mtext('Another look at distributions, then and now (not bootstrapped)', 3,
#     cex=1, font=2)

# dev.off()

#bootstrap some confidence bounds ####
nsamp = 20000
mean_vect_er_68_70 = mean_vect_er_17_19 = mean_vect_gpp_68_70 =
    mean_vect_gpp_17_19 = vector(length=nsamp)
for(i in 1:nsamp){
    samp_68_70_er = sample(er_68_70, size=length(er_68_70), replace=TRUE)
    samp_17_19_er = sample(er_new, size=length(er_new), replace=TRUE)
    samp_68_70_gpp = sample(gpp_68_70, size=length(er_68_70), replace=TRUE)
    samp_17_19_gpp = sample(gpp_new, size=length(gpp_new), replace=TRUE)
    mean_vect_er_68_70[i] = mean(samp_68_70_er, na.rm=TRUE)
    mean_vect_er_17_19[i] = mean(samp_17_19_er, na.rm=TRUE)
    mean_vect_gpp_68_70[i] = mean(samp_68_70_gpp, na.rm=TRUE)
    mean_vect_gpp_17_19[i] = mean(samp_17_19_gpp, na.rm=TRUE)
}

# plot(density(mean_vect_er_68_70 * -1))
CI = data.frame('CI95_lower'=numeric(4), 'median'=numeric(4),
                'CI95_upper'=numeric(4),
                row.names=c('GPP_then', 'GPP_now', 'ER_then', 'ER_now'))
CI[1,] = quantile(sort(mean_vect_gpp_68_70), probs=c(0.025, 0.5, 0.975))
CI[2,] = quantile(sort(mean_vect_gpp_17_19), probs=c(0.025, 0.5, 0.975))
CI[3,] = -quantile(sort(mean_vect_er_68_70) , probs=c(0.025, 0.5, 0.975))
CI[4,] = -quantile(sort(mean_vect_er_17_19), probs=c(0.025, 0.5, 0.975))


# weighted by month CI, proportional to Hall sampling
gpp_68_70_bymo = split(gpp_68_70[!is.na(gpp_68_70)],
                      factor(substr(nhc_68_70$date, 6, 7)))
er_68_70_bymo = split(er_68_70[!is.na(er_68_70)],
                      factor(substr(nhc_68_70$date, 6, 7)))
gpp_new_bymo = split(gpp_new[!is.na(gpp_new)],
                     factor(substr(dates_new[!is.na(gpp_new)], 6, 7)))
er_new_bymo = split(er_new[!is.na(er_new)],
                     factor(substr(dates_new[!is.na(er_new)], 6, 7)))
nsamp_new = length(gpp_new[!is.na(gpp_new)])/length(gpp_68_70[!is.na(gpp_68_70)])
nsamp = 20000
mean_vect_er_68_70 = mean_vect_er_new = mean_vect_gpp_68_70 =
    mean_vect_gpp_new = c()
for(i in 1:nsamp){
    samp_68_70_er = samp_new_er = 
        samp_68_70_gpp = samp_new_gpp = c()
    
    for(m in names(month_counts_68_70)){
        mn = month_counts_68_70[m]
        t_er_68_70 <- er_68_70_bymo[[m]]
        t_gpp_68_70 <- gpp_68_70_bymo[[m]]
        t_er_new <- er_new_bymo[[m]]
        t_gpp_new <- gpp_new_bymo[[m]]
        
        samp_68_70_er = c(samp_68_70_er,
                          sample(t_er_68_70, 
                                 size = mn, 
                                 replace=TRUE))
        samp_new_er = c(samp_new_er,
                        sample(t_er_new, 
                               size = round(nsamp_new * mn, 0), 
                               replace=TRUE))
        samp_68_70_gpp = c(samp_68_70_gpp, 
                           sample(t_gpp_68_70, 
                                  size = mn, 
                                  replace=TRUE))
        samp_new_gpp = c(samp_new_gpp,
                             sample(t_gpp_new, 
                                    size = round(nsamp_new * mn, 0),     
                                    replace=TRUE))
        }
    mean_vect_er_68_70[i] = mean(samp_68_70_er)
    mean_vect_er_new[i] = mean(samp_new_er)
    mean_vect_gpp_68_70[i] = mean(samp_68_70_gpp)
    mean_vect_gpp_new[i] = mean(samp_new_gpp)
    if(i %% 1000 == 0){ print(i) }
}

 # plot(density(mean_vect_er_68_70 * -1))
CI_prop = data.frame('CI95_lower'=numeric(4), 'median'=numeric(4),
                  'CI95_upper'=numeric(4), 
                  'met' = c(rep("GPP", 2), rep("ER", 2)),
                  'era' = rep(c('then', 'now'), 2), 
                  'prop' = "by hall sampling")
CI_prop[1,1:3] = quantile(sort(mean_vect_gpp_68_70), probs=c(0.025, 0.5, 0.975))
CI_prop[2,1:3] = quantile(sort(mean_vect_gpp_new), probs=c(0.025, 0.5, 0.975))
CI_prop[3,1:3] = -quantile(sort(mean_vect_er_68_70), probs=c(0.025, 0.5, 0.975))
CI_prop[4,1:3] = -quantile(sort(mean_vect_er_new), probs=c(0.025, 0.5, 0.975))
# bootstrapped proportional to year, excluding september
    
n_68 = length(gpp_68_70[!is.na(gpp_68_70)])
n_new = length(er_new[!is.na(er_new)])
month_props <- c(31, 28, 31, 30, 31, 30, 31, 31, 0, 31, 30, 31)
names(month_props) <- c('01', '03', '02', '04', '05', '06',
                        '07', '08', '09', '10', '11', '12')
mean_vect_er_68_70 = mean_vect_er_new = mean_vect_gpp_68_70 =
    mean_vect_gpp_new = c()
turboset <- list(c(1:12), c(10,11), c(3,4), c(7,8), c(1,2))
names <- c("year", "oct_nov", "mar_apr", "jul_aug", "jan_feb")
for(ss in 1:5){
    set <- turboset[[ss]]
    for(i in 1:nsamp){
        samp_68_70_er = samp_new_er = 
            samp_68_70_gpp = samp_new_gpp = c()
        
        for(m in names(month_props[set])){
            nn <- month_props[m]/sum(month_props[set])
            t_er_68_70 <- er_68_70_bymo[[m]]
            t_gpp_68_70 <- gpp_68_70_bymo[[m]]
            t_er_new <- er_new_bymo[[m]]
            t_gpp_new <- gpp_new_bymo[[m]]
            
            samp_68_70_er = c(samp_68_70_er,
                              sample(t_er_68_70, 
                                     size = round(nn * n_68, 0), 
                                     replace=TRUE))
            samp_new_er = c(samp_new_er,
                            sample(t_er_new, 
                                   size = round(nn * n_new, 0), 
                                   replace=TRUE))
            samp_68_70_gpp = c(samp_68_70_gpp, 
                               sample(t_gpp_68_70, 
                                      size = round(nn * n_68, 0), 
                                      replace=TRUE))
            samp_new_gpp = c(samp_new_gpp,
                                 sample(t_gpp_new, 
                                        size = round(nn * n_new, 0),     
                                        replace=TRUE))
            }
        mean_vect_er_68_70[i] = mean(samp_68_70_er)
        mean_vect_er_new[i] = mean(samp_new_er)
        mean_vect_gpp_68_70[i] = mean(samp_68_70_gpp)
        mean_vect_gpp_new[i] = mean(samp_new_gpp)
        if(i %% 1000 == 0){ print(i) }
    }

CI_p = data.frame('CI95_lower'=numeric(4), 'median'=numeric(4),
    'CI95_upper'=numeric(4), 'met' = c(rep("GPP", 2), rep("ER", 2)),
    'era' = rep(c('then', 'now'), 2), 'prop' = names[ss])
CI_p[1,1:3] = quantile(sort(mean_vect_gpp_68_70), probs=c(0.025, 0.5, 0.975))
CI_p[2,1:3] = quantile(sort(mean_vect_gpp_new), probs=c(0.025, 0.5, 0.975))
CI_p[3,1:3] = -quantile(sort(mean_vect_er_68_70), probs=c(0.025, 0.5, 0.975))
CI_p[4,1:3] = -quantile(sort(mean_vect_er_new), probs=c(0.025, 0.5, 0.975))

CI_prop <- bind_rows(CI_prop, CI_p)
}
write.csv(CI_prop, 'NHC_2019_metabolism/data/metabolism/compiled/bootstrapped/met_means_bootstrapped_by_month_proportions.csv')


png('figures/bootstrapped_CI_daily_mean_direct_calc_byhall.png', 
    height = 5, width = 7, res = 300, units = 'in')
    par(mfrow = c(1,2), mar = c(3,2,4,1), oma = c(0,3,1,0))
    tmp <- CI_prop %>%
        filter(prop == "by hall sampling") 
    boxplot(t(tmp[,1:3]), col = c(alpha("brown3",.8),"grey35" ), 
            ylim = c(0.24, 0.6),
            xaxt = "n", xlab = "")
    axis(1, at=c(1.5, 3.5), labels=c("GPP", "ER"), line=0)
    mtext("Proportional to Hall 1972 sampling")
    mtext("CI around mean (g C/m2/d)", side = 2, line = 3)
    
    tmp <- CI_prop %>%
        filter(prop == "year") 
    boxplot(t(tmp[,1:3]), col = c(alpha("brown3",.8),"grey35" ), 
            ylim = c(0.24, 0.6),
            xaxt = "n", xlab = "")
    axis(1, at=c(1.5, 3.5), labels=c("GPP", "ER"), line=0)
    mtext("Proportional to days per month")
    legend("bottomright",cex=1, bty = "n",
           c("then", "now"),
           fill = c(alpha("brown3",.75),"grey35" ))
    par(new = T, mfrow = c(1,1))
    mtext('Bootstrapped daily metabolism estimates at Concrete Bridge', 
          line = 2, cex = 1.1)
 
dev.off()

png('figures/bootstrapped_CI_daily_mean_direct_calc_byseason.png', 
    height = 4, width = 8, res = 300, units = 'in')
    par(mfrow = c(1,4), mar = c(3,2,4,1), oma = c(0,3,1,0))
    tmp <- CI_prop %>%
        filter(prop == "jan_feb") 
    boxplot(t(tmp[,1:3]), col = c(alpha("brown3",.8),"grey35" ), 
            ylim = c(0.15, 0.7),
            xaxt = "n", xlab = "")
    axis(1, at=c(1.5, 3.5), labels=c("GPP", "ER"), line=0)
    mtext("Jan - Feb")
    mtext("CI around mean (g C/m2/d)", side = 2, line = 3)
    # legend("bottom",cex=1, bty = "n",
    #        c("then", "now"),
    #        fill = c(alpha("brown3",.75),"grey35" ))
    
    tmp <- CI_prop %>%
        filter(prop == "mar_apr") 
    boxplot(t(tmp[,1:3]), col = c(alpha("brown3",.8),"grey35" ), 
            ylim = c(0.15, 0.7),
            xaxt = "n", xlab = "")
    axis(1, at=c(1.5, 3.5), labels=c("GPP", "ER"), line=0)
    mtext("Mar - Apr")
    
    tmp <- CI_prop %>%
        filter(prop == "jul_aug") 
    boxplot(t(tmp[,1:3]), col = c(alpha("brown3",.8),"grey35" ), 
            ylim = c(0.15, 0.7),
            xaxt = "n", xlab = "")
    axis(1, at=c(1.5, 3.5), labels=c("GPP", "ER"), line=0)
    mtext("Jul - Aug")
    
    tmp <- CI_prop %>%
        filter(prop == "oct_nov") 
    boxplot(t(tmp[,1:3]), col = c(alpha("brown3",.8),"grey35" ), 
            ylim = c(0.15, 0.7),
            xaxt = "n", xlab = "")
    axis(1, at=c(1.5, 3.5), labels=c("GPP", "ER"), line=0)
    mtext("Oct - Nov")

    par(new = T, mfrow = c(1,1))
    mtext('Bootstrapped daily metabolism estimates at Concrete Bridge', 
          line = 3.4, cex = 1.1)
 
dev.off()
         
# CI by months
CI_months = data.frame()

for(m in names(month_counts_68_70)){
    mean_vect_er_68_70 = mean_vect_er_new = mean_vect_gpp_68_70 =
        mean_vect_gpp_new = c()
    for(i in 1:nsamp){
        samp_68_70_er = samp_new_er = 
            samp_68_70_gpp = samp_new_gpp = c()
    
        t_er_68_70 <- er_68_70_bymo[[m]]
        t_gpp_68_70 <- gpp_68_70_bymo[[m]]
        t_er_new <- er_new_bymo[[m]]
        t_gpp_new <- gpp_new_bymo[[m]]
        
        nm <- month_counts_68_70[m]
        samp_68_70_er = c(samp_68_70_er,
                          sample(t_er_68_70, 
                                 size = nm, 
                                 replace=TRUE))
        samp_new_er = c(samp_new_er,
                        sample(t_er_new, 
                               size = nsamp_new * nm, 
                               replace=TRUE))
        samp_68_70_gpp = c(samp_68_70_gpp, 
                           sample(t_gpp_68_70, 
                                  size = nm, 
                                  replace=TRUE))
        samp_new_gpp = c(samp_new_gpp,
                             sample(t_gpp_new, 
                                    size = nsamp_new * nm, 
                                    replace=TRUE))

        mean_vect_er_68_70[i] = mean(samp_68_70_er)
        mean_vect_er_new[i] = mean(samp_new_er)
        mean_vect_gpp_68_70[i] = mean(samp_68_70_gpp)
        mean_vect_gpp_new[i] = mean(samp_new_gpp)
    }
    print(m) 
    nn <- length(t_gpp_new)
    CI = data.frame('month' = rep(m, 4), 
                    'met' = rep(c('gpp', 'er'),2),
                    'era' = c(rep('then', 2), rep('now',2)),
                    'CI95_lower' = numeric(4), 
                    'median' = numeric(4), 
                    'CI95_upper' = numeric(4),
                    'n' = c(nm, nm, nn, nn))
    CI[1, 4:6] = quantile(sort(mean_vect_gpp_68_70), probs=c(0.025, 0.5, 0.975))
    CI[2, 4:6] = -quantile(sort(mean_vect_er_68_70), probs=c(0.025, 0.5, 0.975))
    CI[3, 4:6] = quantile(sort(mean_vect_gpp_new), probs=c(0.025, 0.5, 0.975))
    CI[4, 4:6] = -quantile(sort(mean_vect_er_new), probs=c(0.025, 0.5, 0.975))
    if(nm == 1) { 
        CI[1:2,c(4,6)] <- NA
        }
    CI_months = bind_rows(CI_months, CI)    
}
CI = data.frame('month' = rep('09', 4), 
                'met' = rep(c('gpp', 'er'),2),
                'era' = c(rep('then', 2), rep('now',2)))

CI_months <- bind_rows(CI_months, CI) %>%
    as_tibble() %>%
    mutate(month = as.character(month)) %>%
    arrange(month)

CI_months$met <- factor(CI_months$met, levels = c("gpp", "er"))
write.csv(CI_months, 'NHC_2019_metabolism/data/metabolism/compiled/bootstrapped/monthly_bootstrapped_met_CIs_direct_calc.csv')


then <- nhc_68_70 %>%
    arrange(doy) %>%
    group_by(doy) %>%
    summarize(ER = mean(ER, na.rm = T),
              GPP = mean(GPP, na.rm = T)) %>%
    mutate(across(-doy, ~na.approx(., x = as.numeric(doy))))
now <- nhc_new %>%
    arrange(doy) %>%
    group_by(doy) %>%
    summarize(ER = mean(ER, na.rm = T),
              GPP = mean(GPP, na.rm = T)) %>%
    mutate(across(-doy, ~na.approx(., x = as.numeric(doy))))
par(mfrow = c(2,1), mar = c(2,4,2,2))
plot(now$doy, now$GPP, type = "l", ylim = c(-2, 1.2))
lines(now$doy, now$ER)
lines(then$doy, then$GPP, lty = 2)
lines(then$doy, then$ER, lty = 2)
abline(h = 0)
plot(now$doy, -now$GPP/now$ER, type = 'l')
lines(then$doy, -then$GPP/then$ER, lty = 2)



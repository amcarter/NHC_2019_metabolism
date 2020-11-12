## simulations with variation all at sub-daily scale
# prepare input data (DO used only to pick first DO of each day)
dat <- streamMetabolizer::data_metab('3', res='15')
dat_daily <- data.frame(date=as.Date(paste0("2012-09-", 18:20)),
                        GPP.daily=2, ER.daily=-3, K600.daily=21, stringsAsFactors=FALSE)

# define simulation parameters
mm <- metab_sim(
  specs(mm_name('sim'), err_obs_sigma=0.1, err_proc_sigma=2,
        GPP_daily=NULL, ER_daily=NULL, K600_daily=NULL),
  data=dat, data_daily=dat_daily)
# actual simulation happens during prediction - different each time
get_params(mm)

#>         date K600.daily GPP.daily ER.daily err.obs.sigma err.obs.phi
#> 1 2012-09-18         21         2       -3           0.1           0
#> 2 2012-09-19         21         2       -3           0.1           0
#> 3 2012-09-20         21         2       -3           0.1           0
#>   err.proc.sigma err.proc.phi discharge.daily
#> 1              2            0        19.13472
#> 2              2            0        22.90962
#> 3              2            0        20.21814
predict_metab(mm)
#>         date GPP GPP.lower GPP.upper ER ER.lower ER.upper msgs.fit warnings
#> 1 2012-09-18   2        NA        NA -3       NA       NA       NA         
#> 2 2012-09-19   2        NA        NA -3       NA       NA       NA         
#> 3 2012-09-20   2        NA        NA -3       NA       NA       NA         
#>   errors
#> 1       
#> 2       
#> 3       

predict_DO(mm)[seq(1,50,by=10),]
#>          date          solar.time   DO.sat depth temp.water     light  DO.pure
#> 1  2012-09-18 2012-09-18 04:05:58 9.083329  0.16       3.60    0.0000 8.410000
#> 11 2012-09-18 2012-09-18 06:35:58 9.174032  0.16       3.23  361.1082 8.074219
#> 21 2012-09-18 2012-09-18 09:05:58 8.638056  0.16       5.51 1339.0007 8.873394
#> 31 2012-09-18 2012-09-18 11:35:58 7.784199  0.16       9.68 1778.8622 9.108354
#> 41 2012-09-18 2012-09-18 14:05:58 7.309720  0.16      12.35 1498.8112 8.494834
#>      DO.mod   DO.obs
#> 1  8.410000 8.288769
#> 11 8.120033 8.108756
#> 21 9.088803 8.966921
#> 31 9.373993 9.310303
#> 41 8.645758 8.474404

predict_DO(mm)[seq(1,50,by=10),]
#>          date          solar.time   DO.sat depth temp.water     light  DO.pure
#> 1  2012-09-18 2012-09-18 04:05:58 9.083329  0.16       3.60    0.0000 8.410000
#> 11 2012-09-18 2012-09-18 06:35:58 9.174032  0.16       3.23  361.1082 8.074219
#> 21 2012-09-18 2012-09-18 09:05:58 8.638056  0.16       5.51 1339.0007 8.873394
#> 31 2012-09-18 2012-09-18 11:35:58 7.784199  0.16       9.68 1778.8622 9.108354
#> 41 2012-09-18 2012-09-18 14:05:58 7.309720  0.16      12.35 1498.8112 8.494834
#>      DO.mod   DO.obs
#> 1  8.410000 8.249887
#> 11 7.992529 8.004143
#> 21 8.414859 8.334371
#> 31 8.808968 8.773604
#> 41 8.486494 8.380405
# or same each time if seed is set
mm@specs$sim_seed <- 236
predict_DO(mm)$DO.obs[seq(1,50,by=10)]
#> [1] 8.514605 8.226884 9.242439 9.682286 8.660459
predict_DO(mm)$DO.obs[seq(1,50,by=10)]
#> [1] 8.514605 8.226884 9.242439 9.682286 8.660459
# fancy GPP equation
dat_daily <- data.frame(date=as.Date(paste0("2012-09-", 18:20)),
                        Pmax=8, alpha=0.01, ER.daily=-3, K600.daily=21, stringsAsFactors=FALSE)
mm <- metab_sim(
  specs(mm_name('sim', GPP_fun='satlight'), err_obs_sigma=0.1, err_proc_sigma=2,
        Pmax=NULL, alpha=NULL, ER_daily=NULL, K600_daily=NULL),
  data=dat, data_daily=dat_daily)
get_params(mm)
#>         date K600.daily Pmax alpha ER.daily err.obs.sigma err.obs.phi
#> 1 2012-09-18         21    8  0.01       -3           0.1           0
#> 2 2012-09-19         21    8  0.01       -3           0.1           0
#> 3 2012-09-20         21    8  0.01       -3           0.1           0
#>   err.proc.sigma err.proc.phi discharge.daily
#> 1              2            0        19.36654
#> 2              2            0        19.95024
#> 3              2            0        21.22849

predict_metab(mm) # metab estimates are for data without errors
#>         date      GPP GPP.lower GPP.upper ER ER.lower ER.upper msgs.fit
#> 1 2012-09-18 3.182805        NA        NA -3       NA       NA       NA
#> 2 2012-09-19 3.165193        NA        NA -3       NA       NA       NA
#> 3 2012-09-20 3.147224        NA        NA -3       NA       NA       NA
#>   warnings errors
#> 1                
#> 2                
#> 3                
predict_DO(mm)[seq(1,50,by=10),]
#>          date          solar.time   DO.sat depth temp.water     light  DO.pure
#> 1  2012-09-18 2012-09-18 04:05:58 9.083329  0.16       3.60    0.0000 8.410000
#> 11 2012-09-18 2012-09-18 06:35:58 9.174032  0.16       3.23  361.1082 8.273214
#> 21 2012-09-18 2012-09-18 09:05:58 8.638056  0.16       5.51 1339.0007 9.939735
#> 31 2012-09-18 2012-09-18 11:35:58 7.784199  0.16       9.68 1778.8622 9.910489
#> 41 2012-09-18 2012-09-18 14:05:58 7.309720  0.16      12.35 1498.8112 9.203392
#>       DO.mod    DO.obs
#> 1   8.410000  8.461259
#> 11  8.212799  8.004816
#> 21 10.180817 10.165990
#> 31 10.134104 10.111424
#> 41  9.394355  9.253832
## simulations with variation at both sub-daily and multi-day scales
sp <- specs(mm_name('sim', pool_K600='none'),
            K600_daily = function(n, ...) pmax(0, rnorm(n, 10, 3))) # n is available within sim models
mm <- metab(sp, dat)
get_params(mm)
#>         date K600.daily GPP.daily  ER.daily err.obs.sigma err.obs.phi
#> 1 2012-09-18  10.973189  0.872383 -7.136814          0.01           0
#> 2 2012-09-19   8.874616  6.162460 -6.991966          0.01           0
#> 3 2012-09-20   9.813784  6.693289  0.000000          0.01           0
#>   err.proc.sigma err.proc.phi discharge.daily
#> 1            0.2            0        19.37273
#> 2            0.2            0        20.94919
#> 3            0.2            0        20.99108

predict_metab(mm)
#>         date       GPP GPP.lower GPP.upper         ER ER.lower ER.upper
#> 1 2012-09-18  5.755216        NA        NA  -1.801561       NA       NA
#> 2 2012-09-19  5.479744        NA        NA -10.943692       NA       NA
#> 3 2012-09-20 17.948548        NA        NA  -4.696068       NA       NA
#>   msgs.fit warnings errors
#> 1       NA                
#> 2       NA                
#> 3       NA 

## K~Q model
dat <- data_metab('10','15')
sp <- specs(mm_name('sim', pool_K600='binned'))
mm <- metab(sp, dat)
pars <- get_params(mm)
attr(pars, 'K600_eqn')
#> $K600_lnQ_nodes_centers
#> [1] 2.7 2.9 3.1 3.3
#> 
#> $K600_lnQ_cnode_meanlog
#>  [1] 1.791759 1.791759 1.791759 1.791759 1.791759 1.791759 1.791759 1.791759
#>  [9] 1.791759 1.791759
#> 
#> $K600_lnQ_cnode_sdlog
#>  [1] 1 1 1 1 1 1 1 1 1 1
#> 
#> $K600_lnQ_nodediffs_meanlog
#>  [1] 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2
#> 
#> $K600_lnQ_nodediffs_sdlog
#>  [1] 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5
#> 
#> $lnK600_lnQ_nodes
#> [1] 1.8055839 1.7903077 1.2187666 0.6656684
#> 
#> $K600_daily_predlog
#>  [1] 1.2645783 1.6185911 0.9744718 1.7094699 1.3963858 1.7913175 0.9430280
#>  [8] 1.4406586 0.9215569 1.4184119
#> 
if (FALSE) {
  plot_DO_preds(predict_DO(mm))
  plot_DO_preds(mm)
  library(ggplot2)
}

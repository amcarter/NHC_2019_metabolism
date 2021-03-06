# Initial GHG analyses
# NHC gas data from 11/2019 - 3/2020

library(tidyverse)
library(lubridate)
library(ggpubr)

setwd("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism")


dat <- read_csv("data/gas_data/compiled_ghg_dataframe_allvariables_filled.csv") %>%
  mutate(distance_m = 8450 - distance_upstream_m,
         datetime = force_tz(datetime, tz = "EST"), 
         date = as.Date(datetime))
d2 <- read_csv("data/gas_data/compiled_ghg_dataframe_short.csv") %>%
  mutate(force_tz(datetime, tz = "EST"))

# Mass Balance Calculations ####
# F_out = F_in + F_instream_prod + F_gw - F_evasion 

#filter dataframe
gas <- dat %>%
  select(date, site, CH4.ugL, CO2.ugL, N2O.ugL) %>% 
  pivot_longer(ends_with("ugL"), names_to = "gas", 
               values_to = "concentration_ugl") %>%
  mutate(gas = substr(gas, 1,3)) 
sat <- dat %>%
  select(date, site, ends_with("sat")) %>% 
  pivot_longer(ends_with("sat"), names_to = "gas", 
               values_to = "sat_ugl") %>%
  mutate(gas = substr(gas, 1,3)) 
k <- dat %>%
  select(date, site, starts_with("K_")) %>% 
  pivot_longer(starts_with("K_"), names_to = "gas", 
               values_to = "K_gas_d") %>%
  mutate(gas = substr(gas, 3,5)) 

mb <- dat %>%
  select(date, site, distance_m, width_m = width_march_m, 
         ws_area_km2, discharge_m3s, depth_m = depth,
         watertemp_C, AirPres_kPa, AirTemp_C, GPP, ER,
         ends_with("ugld")) %>%
  pivot_longer(ends_with("ugld"), names_to = "gas", values_to = "flux_ugld") %>%
  mutate(gas = substr(gas, 1,3)) %>%
  left_join(gas, by = c("date", "site", "gas")) %>%
  left_join(k, by = c("date", "site", "gas")) %>%
  left_join(sat, by = c("date", "site", "gas")) %>%
  arrange(date, distance_m)

mb$site <- factor(mb$site, levels = c("NHC", "PM", "CBP", "WB", "WBP","UNHC"))
sites <- dat %>%
  select(site, distance_m, width_march_m) %>%
  group_by(site) %>%
  summarize_all(mean) %>%
  ungroup() %>%
  arrange(distance_m)
  
sites <- sites %>%
  mutate(sa_m2 = c(NA, diff(sites$distance_m)) * width_march_m) %>%
  select(site, sa_m2)
  

# F_in is the mass flux across the input point:
# F_in = pGas (mg/m3) * discharge (m3/s)
# first, CO2

# 
# dd <- mb %>%
#   filter(date == dates[1]) %>%
#   # arrange(distance_m) %>%
#   mutate(discharge_m3s = na.approx(discharge_m3s, ws_area_km2, na.rm = T),
#          depth = na.approx(depth, distance_m, na.rm = F),
#          GPP = na.approx(GPP, distance_m, na.rm = F),
#          ER = na.approx(ER, distance_m, na.rm = F),
#          CH4.ugL = na.approx(CH4.ugL, distance_m, na.rm = F),
#          CO2.ugL = na.approx(CO2.ugL, distance_m, na.rm = F),
#          N2O.ugL = na.approx(N2O.ugL, distance_m, na.rm = F),
#          CH4.flux_ugld = na.approx(CH4.flux_ugld, distance_m, na.rm = F),
#          CO2.flux_ugld = na.approx(CO2.flux_ugld, distance_m, na.rm = F),
#          N2O.flux_ugld = na.approx(N2O.flux_ugld, distance_m, na.rm = F)) 
#   


dates <- unique(mb$date)
flux = data.frame()
# Calculate for each date:
for(i in 1:length(dates)){
  dd <- mb %>%
    filter(date == dates[i]) %>%
    left_join(sites, by = "site") 
  n <- nrow(dd)-3
  tmp <- dd %>%
    slice(c(-1, -2, -3)) %>%
    mutate(F_in_mgm2d = dd$discharge_m3s[1:n] * 
             dd$concentration_ugl[1:n]/sa_m2*60*60*24,
           F_out_mgm2d = discharge_m3s * concentration_ugl/sa_m2*60*60*24,
           F_evas_mgm2d = (flux_ugld + dd$flux_ugld[1:n]) * depth_m / 2,
           F_instr_mgm2d = ifelse(gas == "CO2", (ER - GPP) * (44 / 32) *1000, 0),
           F_ot_mgm2d =  F_out_mgm2d - F_in_mgm2d - F_instr_mgm2d + 
             F_evas_mgm2d,
           gw_in_m3s = discharge_m3s - dd$discharge_m3s[1:n]) %>%
    select(c(site, date, gas, gw_in_m3s, ends_with("mgm2d")))
  
  dd <- dd %>% left_join(tmp, by = c("site", "date", "gas"))
  flux = bind_rows(flux, dd)
}

flux$site <- factor(flux$site, levels = c("UNHC", "WBP", "WB","CBP", "PM", "NHC"))

ggplot(flux, aes(x = site, F_ot_mgm2d)) +
  geom_boxplot() +
  facet_wrap(gas~., scales = "free_y", ncol = 1) +
  geom_hline(yintercept = 0, col = "grey30")

longf <- flux %>%
  as.tibble() %>%
  mutate(F_evas_mgm2d = -F_evas_mgm2d) %>%
  pivot_longer(starts_with("F_"), names_to = "flux_category", values_to = "mgm2d")
png("figures/gas/flux_components_by_site.png",
    width = 8, height = 6, res = 300, units = "in")
  longf %>%
    filter(site != "UNHC", 
           !(flux_category %in% c("F_in_mgm2d", "F_out_mgm2d"))) %>%
    ggplot(aes(site, (mgm2d), fill = flux_category)) +
      geom_boxplot(position = "dodge") +
      facet_wrap(gas~., scales = "free_y", ncol = 1) +
      geom_hline(yintercept = 0, col = "gray30")
dev.off()

ggplot(longf, aes(site, gw_in_m3s)) +#, col = factor(date))) +
  geom_boxplot()

longf %>%
  filter(site != "UNHC",
         flux_category == "F_ot_mgm2d") %>%
  ggplot(aes(gw_in_m3s/sa_m2 *60*60*24, mgm2d, col = site)) +
    geom_point() +
    facet_wrap(gas~., scales = "free_y") +
    xlab("ground water input (m/d)") +
    ylab("missing flux (mg/m2/d)") +
    geom_hline(yintercept = 0, col = "grey30")

write_csv(longf, "data/gas_data/mass_balance_nhc_ghg_data.csv")
chem <- dat %>% 
  mutate(NH4NO3 = nh4n_mgl/no3n_mgl,
         CH4CO2 = (CH4.ugL/16)/(CO2.ugL/44),
         CH4O2 = (CH4.ugL/16)/(DO.obs/32)) %>%
  select(site, date, no3n_mgl, NH4NO3, CH4CO2, CH4O2, DO.obs, doc_mgl) 

flux %>%
  left_join(chem, by = c("site", "date")) %>%
  ggplot(aes(CH4CO2, F_ot_mgm2d, col = site)) +
    geom_point() +
    facet_wrap(gas~., scales = "free_y", ncol = 1) +
    xlab("CH4:CO2") +
    xlim(0,.0087) +
    # xlim(0,1.8) +
    ylab("missing flux (mg/m2/d)") +
    geom_hline(yintercept = 0, col = "grey30")

png("figures/gas/missing_flux_by_gw_anaerobic_metric.png",
    width = 8, height = 7, res = 300, units = "in")
  flux %>%
    left_join(chem, by = c("site", "date")) %>%
    mutate('accumulated gw (m3/d)' = gw_in_m3s/60*60*24) %>%
    filter(NH4NO3 < 0.5, 
           CH4CO2 < 0.0087) %>%
    pivot_longer(c("NH4NO3","CH4CO2", "accumulated gw (m3/d)"), 
                 names_to = "element_ratio", values_to = "ratio_value") %>%
    ggplot(aes(ratio_value, (F_ot_mgm2d), col = site)) +
      geom_point() +
      facet_grid(gas~element_ratio, scales = "free") +
      xlab("") +
      # xlim(0,.5) +
      ylab("missing flux (mg/m2/d)") +
      geom_hline(yintercept = 0, col = "grey30")
dev.off() 


ggplot(chem, aes(site, CH4O2)) +
  geom_boxplot()
-CO2_evas_mgm2d),col = "steelblue", lwd = 2) +
  geom_line(aes(y = CO2_gw_mgm2d), col = "sienna3", lwd = 2)# +
  geom_line(aes(y = CO2_ot_mgm2d), lty = 2, lwd = 2)# +
  legend()

# dc <- bind_rows(dd[1,], dc)
par(mfrow = c(3,1), mar = c(1,4,0,2), oma = c(4,0,4,0))
plot(dc$distance_m, dc$CO2_out_mgm2d - dc$CO2_in_mgm2d, lwd = 2, 
     type = "l", ylim = c(-6000, 12000), ylab = "CO2 mg/m2/d",  
      xaxt = 'n', col.bor = "grey30") 
lines(dc$distance_m, dc$CO2_instr_mgm2d, col = "forestgreen", lwd = 2) 
lines(dc$distance_m, dc$CO2_evas_mgm2d, col = "steelblue", lwd = 2) 
lines(dc$distance_m, dc$CO2_gw_mgm2d, col = "sienna3", lwd = 2)
abline(h = 0, col = "grey30")
legend("topleft",
       legend = c("total flux", "In stream production", "evasion", 
                  "missing flux (gw, anaerobic)"),
       lty = 1, lwd = 2, bty = 'n', 
       col = c("black", "forestgreen", "steelblue","sienna3"))
mtext(dc$date[1], outer = T, cex = 1.2)

plot(dc$distance_m, dc$CH4_out_mgm2d - dc$CH4_in_mgm2d, lwd = 2, 
     type = "l", ylim = c(-60, 30), ylab = "CH4 mg/m2/d", 
     xaxt = 'n', col.bor = "grey30") 
lines(dc$distance_m, dc$CH4_evas_mgm2d, col = "steelblue", lwd = 2) 
lines(dc$distance_m, dc$CH4_ot_mgm2d, col = "sienna3", lwd = 2)
abline(h = 0, col = "grey30")
# legend("topleft",
#        legend = c("total flux", "evasion", "missing flux (gw, in stream)"),
#        lty = 1, lwd = 2, bty = 'n', 
#        col = c("black", "steelblue","sienna3"))

plot(dc$distance_m, dc$N2O_out_mgm2d - dc$N2O_in_mgm2d, lwd = 2, 
     type = "l", ylim = c(-.1, 1.5), ylab = "N2O mg/m2/d", xlab = "distance", bor.col = "grey30") 
lines(dc$distance_m, dc$N2O_evas_mgm2d, col = "steelblue", lwd = 2) 
lines(dc$distance_m, dc$N2O_ot_mgm2d, col = "sienna3", lwd = 2)
abline(h = 0, col = "grey30")
# legend("topleft",
#        legend = c("total flux", "evasion", "missing flux (gw, in stream)"),
#        lty = 1, lwd = 2, bty = 'n', 
#        col = c("black", "steelblue","sienna3"))





ds <- data.frame(distance_m = seq(0, dc$distance_m[6], by = 10)) %>%
  left_join(dc, by = "distance_m") %>%
  mutate(across(c(-datetime, -date, -site), ~na.approx(.,na.rm = F)))

plot(ds$distance_m, ds$CO2_pred_ugL, type = "l", ylim = c(1, 7500), log = "y")

for(r in 2:nrow(ds)){
  lambda = (ds$K_CO2[r-1]/60/60/24)/
    (ds$discharge_m3s[r-1]/ds$depth_m[r-1]/ds$width_m[r-1])
  F_instr = (ds$ER[r] - ds$GPP[r])*(44/32) * ds$width_m[r] * 10
  F_out = ds$CO2_pred_ugL[r-1] * exp(-10* lambda) + F_instr/ds$discharge_m3s[r]/60/60/24
  ds$CO2_pred_ugL[r] <- F_out
}

lines(ds$distance_m, ds$CO2_pred_ugL, col = 2)

  g)
red_ugL = (CO2_in_mg + CO2_instr_mg - CO2_evas_mg)/
            discharge_m3s/60/60/24/ttime_d,
         CH4_pred_ugL = (CH4_in_mg - CH4_evas_mg)/
            discharge_m3s/60/60/24/ttime_d,
         N2O_pred_ugL = (N2O_in_mg - N2O_evas_mg)/
           discharge_m3s/60/60/24/ttime_d) %>%

ggplot(flux, aes(distance_m, CO2.ugL, color = factor(date))) + 
  geom_point() + geom_line() + ylim(0,7500) +
  geom_point(aes(y = CO2_pred_ugL), color = 1)

plot(flux$distance_m, flux$CO2.ugL, type = "b")
points(flux$distance_m, flux$CO2_pred_ugL, col = "brown3")
t, distance_upstream_m, watertemp_C, GPP, ER, 
         discharge_m3s, DO.obs, CH4.flux_ugld, CO2.flux_ugld, N2O.flux_ugld) %>%
  pivot_longer(ends_with("ugld"), names_to = "gas", values_to = "flux_ugld") %>%
  mutate(gas = substr(gas, 1,3)) %>%
  left_join(gas, by = c("datetime", "site", "gas")) %>%
  mutate(date = as.Date(datetime),
         date = case_when(date == as.Date("2019-12-04") ~ 
                            as.Date("2019-12-03"),
                          date == as.Date("2020-01-30") ~ 
                            as.Date("2020-01-29"), 
                          TRUE ~ date)) %>%
  arrange(date, distance_upstream_m)

flux$site <- factor(flux$site, levels = c("NHC", "PM", "CBP", "WB", "WBP","UNHC"))

# png("figures/gas/ghg_boxplots_bydate.png", height = 5, width = 7, units = "in", res = 300)
  ggplot(flux, aes(x = habitat, y = flux_ugld, 
                   fill = habitat)) +
    geom_boxplot() +
    facet_wrap(gas~., scales = "free_y", ncol = 1)
# dev.off()
png("figures/gas/ghg_longitudinal_bydate.png", height = 5, width = 7, units = "in", res = 300)
ggplot(flux, aes(x = -distance_upstream_m, y = concentration_ugl, 
                 color = factor(date))) +
  geom_line() +
  geom_point() +
  xlab("upstream --> downstream (m)")+
  facet_wrap(gas~., scales = "free_y", ncol = 1)
dev.off()
png("figures/gas/ghgflux_longitudinal_boxplots.png", height = 5, width = 4, units = "in", res = 300)
ggplot(flux, aes(x = site, y = flux_ugld))+ 
  geom_boxplot() +
  scale_x_discrete(limits = rev(levels(flux$site)))+
  xlab("upstream --> downstream (site)")+
  facet_wrap(gas~., scales = "free_y", ncol = 1)
dev.off()

ggplot(flux, aes(x = watertemp_C, concentration_ugl, color = date)) +
  geom_point(size = 2) +
  # geom_smooth(method = lm, color = "grey30") +
  # theme(legend.position = "bottom") +
  # scale_color_gradient(name = "log(Q)") +
  facet_wrap(gas~., scales = "free", ncol = 1)


# ccc <- ggplot(flux, aes(x = GPP+ER, y = concentration_ugl, 
#                  color = watertemp_C)) +
#   geom_point(size = 2) +
#   geom_smooth(method = lm, color = "grey30") +
#   theme(legend.position = "bottom") +
#   scale_color_gradient(name = "temp C") +
#   facet_wrap(gas~., scales = "free_y", ncol = 1)
# fff <- ggplot(flux, aes(x = GPP+ER, y = flux_ugld, 
#                  color = watertemp_C)) +
#   geom_point(size = 2) +
#   geom_smooth(method = lm, color = "grey30") +
#   theme(legend.position = "bottom") +
#   scale_color_gradient(name = "temp C") +
#   facet_wrap(gas~., scales = "free_y", ncol = 1)
# png("figures/gas/linear_ghg_by_NEPtemp.png", height = 7, width = 6, units = "in", res = 300)
# ggarrange(cc, ff)
# ggarrange(ccc, fff)
# dev.off()


# PCA ####
dat.pca <- d2 %>%
  select(-site, -habitat, -date, -K600, -ends_with("ugld")) %>%
  prcomp()

autoplot(dat.pca, data = d2, colour = "habitat")



# linear mixed effects models ####

# rescale covariates, normalize each to the mean
scaled <- dat %>%
  mutate(logQ = log(discharge_m3s),
         no3n_mgl = ifelse(no3n_mgl > 0.6, NA, no3n_mgl), 
         doc_mgl = ifelse(doc_mgl < 0, 0, doc_mgl)) %>%
  select(site, habitat, logQ, watertemp_C, ER, GPP, DO.obs, no3n_mgl, doc_mgl, 
         ends_with("ugld"), ends_with("ugL")) %>%
  mutate(across(-c(site, habitat), ~scale(.)[,1, drop = T]), 
         across(c(site, habitat), ~factor(.)), 
         ER = -ER, 
         NEP = ER - GPP)
apply()
mm <- lmer(N2O.flux_ugld ~ CH4.flux_ugld + no3n_mgl + (1|site), data = scaled)

# steps_r <- mm_step$random %>%
#   as.tibble() %>%
#   mutate(predictor = rownames(mm_step$random), 
#          effect = "random")
# 
# steps_lme <- mm_step$fixed %>%
#   as.tibble() %>%
#   mutate(predictor = rownames(mm_step$fixed),
#          effect = "fixed") %>%
#   bind_rows(steps_r)

mm <- lmer(CH4.flux_ugld ~ logQ + watertemp_C + GPP + ER + DO.obs + 
            CO2.ugL + (1|site), data = scaled)
mm_step <- lmerTest::step(mm)
mGPPER <- get_model(mm_step)
confint(mGPPER)
summary(mGPPER) 
# r.squaredGLMM(mGPPER)
steps_r <- mm_step$random %>%
  as.tibble() %>%
  mutate(predictor = rownames(mm_step$random),
         effect = "random")

CH4_steps <- mm_step$fixed %>%
  as.tibble() %>%
  mutate(predictor = rownames(mm_step$fixed),
         effect = "fixed") %>%
  bind_rows(steps_r) %>%
  mutate(gas = "CH4")

mm <- lmer(CO2.flux_ugld ~ logQ + watertemp_C + GPP + ER + DO.obs + 
             (1|site), data = scaled)
mm_step <- lmerTest::step(mm)
mGPPER <- get_model(mm_step)
confint(mGPPER)
summary(mGPPER)
steps_r <- mm_step$random %>%
  as.tibble() %>%
  mutate(predictor = rownames(mm_step$random),
         effect = "random")

CO2_steps <- mm_step$fixed %>%
  as.tibble() %>%
  mutate(predictor = rownames(mm_step$fixed),
         effect = "fixed") %>%
  bind_rows(steps_r) %>%
  mutate(gas = "CO2") 

mm <- lmer(N2O.flux_ugld ~ logQ + watertemp_C + GPP + ER + DO.obs + 
              (1|site), data = scaled)
mm_step <- lmerTest::step(mm)
mGPPER <- get_model(mm_step)
confint(mGPPER)
summary(mGPPER) 
steps_r <- mm_step$random %>%
  as.tibble() %>%
  mutate(predictor = rownames(mm_step$random),
         effect = "random")

model_steps <- mm_step$fixed %>%
  as.tibble() %>%
  mutate(predictor = rownames(mm_step$fixed),
         effect = "fixed") %>%
  bind_rows(steps_r) %>%
  mutate(gas = "N2O") %>%
  bind_rows(CO2_steps,
            CH4_steps) %>%
  slice(c(-6, -7, -13, -14, -21, -22)) %>%
  select(gas, predictor, Eliminated, sum_sq = 'Sum Sq', 
         NumDF, DenDF, F_value = 'F value', 'Pr(>F)')


write_csv(model_steps, "data/gas_data/Satterthwaite_DFmethod_lme_steps.csv")

# Analyze water chem data ####
spchem <- read_csv("data/water_chemistry/all_grab_data.csv") %>%
  filter(siteID %in% c("NHC", "UNHC")) %>%
  select(-flagID, -flagComment, -methodDetail, -writeInMethod, -regionID, -method) %>% 
  group_by(siteID, dateTimeUTC, variable) %>%
  summarize_all(mean, na.rm = T) %>%
  ungroup() %>%#data.frame()
  pivot_wider(names_from = "variable", values_from = "value") %>%
  select(-TDN, -TN, -NH4, -PO4) %>%
  filter(across(-c(siteID, dateTimeUTC, TOC), ~!is.na(.)))

# pair with discharge and temperature
raw_dat <- data.frame()
for(site in c("NHC", "UNHC")){
  dat <- read_csv(paste0("data/metabolism/processed/", site, ".csv"), 
                  guess_max = 10000)
  dat$site <- site
  raw_dat = bind_rows(raw_dat, dat)
}
         
chem <- raw_dat %>% 
  as.tibble() %>%
  select(dateTimeUTC = DateTime_UTC, siteID = site, discharge, temp.water) %>%
  mutate(across(c(discharge, temp.water), na.approx, na.rm = T)) %>%
  right_join(spchem, by = c("dateTimeUTC", "siteID")) %>%
  filter(across(c(discharge, temp.water), ~!is.na(.))) %>%
  mutate(logQ = log(discharge),
         across(-c(siteID, dateTimeUTC, discharge), ~scale(.)[,1, drop = T]))

apply(chem, 2, function(x) sum(is.na(x)))
  
chem.pca <- chem %>%
  select(-siteID, -dateTimeUTC, -discharge, -TOC, -temp.water) %>%
  prcomp()

png("figures/gas/waterchem_pca.png", height = 5, width = 5, units = "in", res = 300)
autoplot(chem.pca, data = chem, colour = "siteID", size = 2, 
         loadings = TRUE, loadings.label.colour = 'black',
         loadings.label = TRUE, loadings.label.size = 5) +
  # scale_color_gradient(low = "black", high = "red") +
  labs(title = "Water Chem by Discharge @ NHC, UNHC")
dev.off()
# biplot(chem.pca, pch = 20)
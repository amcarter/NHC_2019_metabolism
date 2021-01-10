# Initial GHG analyses
# NHC gas data from 11/2019 - 3/2020

library(tidyverse)
library(lubridate)
library(ggpubr)

setwd("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism")


dat <- read_csv("data/gas_data/compiled_ghg_dataframe_allvariables.csv")
d2 <- read_csv("data/gas_data/compiled_ghg_dataframe_short.csv")

# Mass Balance Calculations ####
# F_out = F_in + F_instream_prod + F_gw - F_evasion 

#filter dataframe
mb <- dat %>%
  filter(!is.na(datetime),
         site != "MC751") %>%
  mutate(distance_m = 8450 - distance_upstream_m) %>%
  select(datetime, site, distance_m, width_march_m, ws_area_km2,
         AirPres_kPa, AirTemp_C, GPP, ER, discharge_m3s, watertemp_C,
         ends_with(c("ugL", "sat","ugld")), starts_with("K")) 

# F_in is the mass flux across the input point:
# F_in = pGas (mg/m3) * discharge (m3/s)
# first, CO2

mb %>% 
  mutate(across(ends_with("ugL"), ~discharge_m3s * ., .names = "{col}_in_mgs"))

onger(ends_with("ugL"), names_to = "gas", 
               values_to = "concentration_ugl") %>%
  mutate(gas = substr(gas, 1,3)) 

flux <- dat %>%
  filter(!is.na(datetime),
         site != "MC751") %>%
  select(datetime, site,  habitat, distance_upstream_m, watertemp_C, GPP, ER, 
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
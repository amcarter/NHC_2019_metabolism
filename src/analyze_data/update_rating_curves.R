# Update rating curve file with new discharge measurement

# A Carter
# 12/15/2020

setwd("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism/")
source("src/analyze_data/calc_discharge_from_crosssection.r")

# 1. load and munge data ####
qdat <- read_csv("data/rating_curves/discharge_measurements.csv")
sp_qdat <- read_csv("data/rating_curves/discharge_measurements_NHC_UNHC.csv")
sdat <- read_csv("data/siteData/NHCsite_metadata.csv")
                 
# 2. initial calcs to create rc_sheet ####
# setting negative V equal to zero to be consistent with E Moore's data
qdat$velocity_ms[qdat$velocity_ms<0] <- 0

qq <- data.frame()
for(site in unique(qdat$site)){
  sdat <- qdat %>%
  filter(site == !!site)
  
  for(date in unique(sdat$date)){
    dat <- sdat %>%
      filter(date == !!date)
  
    plot_xc(dat)
    out <- calc_xc_discharge(dat$distance_m, dat$depth_m, dat$velocity_ms)
    out <- dat %>%
      slice(1) %>%
      select(date, time, site, stage_m, notes) %>%
      bind_cols(out)
    qq <- bind_rows(qq, out)
  }
}

qq[,-5]

# 3. run on old SP data for NHC and UNHC ####
#qq <- data.frame()
qdat <- sp_qdat
for(site in unique(qdat$site)){
  sdat <- qdat %>%
  filter(site == !!site)
  
  for(date in unique(sdat$date)){
    dat <- sdat %>%
      filter(date == !!date)
  
    plot_xc(dat)
    out <- calc_xc_discharge(dat$distance_m, dat$depth_m, dat$velocity_ms)
    out <- dat %>%
      slice(1) %>%
      select(date, time, site, stage_m, notes) %>%
      bind_cols(out)
    qq <- bind_rows(qq, out)
  }
}

# qq[,-5]
dev.off()
qq <- sdat %>% 
  select(site = sitecode, location_m = distance_m) %>%
  right_join(qq, by = "site")
write_csv(qq, "data/rating_curves/calculated_discharge.csv")  



# 3. Add a new point to the Q file ####
# this step is likely to be custom for each site/date, so we'll do one at a time



write_csv(qq, "data/rating_curves/calculated_discharge.csv")  

# 4. add stages ####
# for the points that don't have it get from LL
# for the points that do, compare with LL

# NHC
nhc <- read_csv("data/metabolism/processed/NHC.csv") %>%
  group_by(date = as.Date(DateTime_EST)) %>%
  select(date, level_m, NHC_Q = discharge) %>%
  summarize_all(mean, na.rm = T) %>%
  mutate(site = "NHC") %>%
  ungroup() 
nhcq <- read_csv("data/metabolism/processed/UNHC.csv") %>%
  group_by(date = as.Date(DateTime_EST)) %>%
  select(date, UNHC_Q = discharge) %>%
  summarize_all(mean, na.rm = T) %>%
  full_join(nhc, by = "date") %>%
  select(-site, -level_m) %>%
  ungroup()
pm <- read_csv("data/metabolism/processed/PM_lvl.csv") %>%
  group_by(date = as.Date(DateTime_EST)) %>%
  select(date, level_m) %>%
  summarize_all(mean, na.rm = T) %>%
  mutate(site = "PM") %>%
  ungroup() 
cbp <- read_csv("data/metabolism/processed/CBP_lvl.csv") %>%
  group_by(date = as.Date(DateTime_EST)) %>%
  select(date, level_m) %>%
  summarize_all(mean, na.rm = T) %>%
  mutate(site = "CBP") %>%
  ungroup() 
wb <- read_csv("data/metabolism/processed/WB_lvl.csv") %>%
  group_by(date = as.Date(DateTime_EST)) %>%
  select(date, level_m) %>%
  summarize_all(mean, na.rm = T) %>%
  mutate(site = "WB") %>%
  ungroup() 
wbp <- read_csv("data/metabolism/processed/WBP_lvl.csv") %>%
  group_by(date = as.Date(DateTime_EST)) %>%
  select(date, level_m) %>%
  summarize_all(mean, na.rm = T) %>%
  mutate(site = "WBP") %>%
  ungroup() 
unhc <- read_csv("data/metabolism/processed/UNHC.csv") %>%
  group_by(date = as.Date(DateTime_EST)) %>%
  select(date, level_m) %>%
  summarize_all(mean, na.rm = T) %>%
  mutate(site = "UNHC") %>%
  ungroup() 

all <- nhc %>%
  select(-NHC_Q) %>%
  bind_rows(pm, cbp, wb, wbp, unhc) %>%
  pivot_wider(names_from = "site", values_from = "level_m") %>%
  arrange(NHC)

plot(all$NHC, all$UNHC,pch = 20, xlim = c(.4, 1.5), ylim = c(.2, 1.5))
points(all$NHC, all$PM,pch = 20, col = 2)
points(all$NHC, all$CBP, pch = 20,col = 3)
points(all$NHC, all$WB, pch = 20,col = 4)
points(all$NHC, all$WBP, pch = 20,col = 5)

all %>%
  filter(!is.na(CBP)) %>%
  ggplot( aes(NHC, CBP, color = date)) +
    geom_point()
# there is an approximate linear relationship between level at NHC and
# at the other sites, I will use this to fill in the level when it is missing
qq <- read_csv("data/rating_curves/calculated_discharge.csv")

m <- lm(PM ~ NHC, all)
l.pm <- summary(m)$coefficients[,1]
m <- lm(CBP ~ NHC, all)
l.cbp <- summary(m)$coefficients[,1]
m <- lm(WB ~ NHC, all)
l.wb <- summary(m)$coefficients[,1]
m <- lm(WBP ~ NHC, all)
l.wbp <- summary(m)$coefficients[,1]
m <- lm(UNHC ~ NHC, all)
l.unhc <- summary(m)$coefficients[,1]

tmp <- all %>%
  rename(NHC_level_m = NHC) %>%
  right_join(qq, by = c("date")) %>%
  mutate(level_m = case_when(grepl("PM", site) ~ PM,
                             grepl("CBP", site) ~ CBP,
                             grepl("WBP", site) ~ WBP,
                             grepl("WB", site) ~ WB, 
                             grepl("UNHC", site) ~ UNHC,
                             grepl("NHC", site) ~ NHC_level_m)) %>%
  select(-PM, -CBP, -WB, -WBP, -UNHC) %>%
  mutate(level_m = case_when(!is.na(level_m) ~ level_m,
                             grepl("PM", site) ~ l.pm[1] + l.pm[2] * NHC_level_m,
                             grepl("CBP", site) ~ l.cbp[1] + l.cbp[2] * NHC_level_m,
                             grepl("WBP", site) ~ l.wbp[1] + l.wbp[2] * NHC_level_m, 
                             grepl("WB", site) ~ l.wb[1] + l.wb[2] * NHC_level_m,
                             grepl("UNHC", site) ~ l.unhc[1] + l.unhc[2] * NHC_level_m))
  


# the 9/22/2016 date is problematic
w <- which(tmp$date == as.Date("2016-09-22"))
tmp$level_m[w] <- NA
plot(tmp$stage_m, tmp$level_m)
abline(0,1)

m <- lm(level_m ~ stage_m, tmp)
ll <- summary(m)$coefficients[,1]
qq_l <- tmp %>%
  mutate(level_m = case_when(is.na(level_m) ~ ll[1] + ll[2]* stage_m, 
                             !is.na(level_m) ~ level_m)) %>%
  left_join(nhcq, by = "date")

# inspect data ####
qq_l %>%
  filter(site == "NHC") %>%
  ggplot(aes(level_m, discharge)) +
  geom_point(aes(color = date), size = 2)

# Assign proper discharges to the pool sections on days when a pool and 
# neighboring riffle were measured

# PM: that q is just wrong, but we didn't measure a riffle. Will pair with 
# interpolated Q to see what we get

# CBP
qq_l %>%
  filter(grepl("CBP", site))
# these two Q's look comprable, I will average them.
w <- which(grepl("CBP", qq_l$site) & qq_l$date == as.Date("2020-06-19"))
qq_l$discharge[w] <- mean(qq_l$discharge[w], na.rm = T)
w <- which(grepl("CBP", qq_l$site) & qq_l$date == as.Date("2020-08-23"))
qq_l$discharge[w] <- mean(qq_l$discharge[w], na.rm = T)

# WBP
qq_l %>%
  filter(grepl("WBP", site))
# here I am going to trust the riffle velocity
w <- which(grepl("WBP", qq_l$site) & qq_l$date == as.Date("2020-06-19"))
qq_l$discharge[w[1]] <- qq_l$discharge[w[2]]

# now correct the average velocites based on new Qs 
# and get rid of the riffle points
qq_l$velocity_avg <- qq_l$discharge/qq_l$xc_area
qq_l <- qq_l %>%
  filter(!grepl("riffle", site))

# look across sites on the two dates where multiple were measured.

a <- qq_l %>%
  filter(date == as.Date("2020-06-19")) 
data.frame(location_m = c(sdat$distance_m[sdat$sitecode =="NHC"],
                          sdat$distance_m[sdat$sitecode =="UNHC"]),
           discharge = c(a$NHC_Q[1], a$UNHC_Q[1])) %>%
  bind_rows(a) %>%
  ggplot(aes(location_m, discharge)) +
  geom_point()

# 6/19/20 is on the falling limb of a storm, while 8/23/20 is just low.
# none of this makes any sense.
write_csv(qq_l, "data/rating_curves/calculated_discharge_modified.csv")

# Build Rating Curves ####
nhc <- qq_l %>%
  filter(site =="NHC")
unhc <- qq_l %>%
  filter(site =="UNHC")
m <- lm(log(discharge)~log(level_m), data = nhc)
m_coef <- summary(m)$coefficients[,1]
plot(nhc$level_m, nhc$discharge)
lines(seq(0,2, by = .01), exp(m_coef[1]) * seq(0,2, by = .01) ^ m_coef[2])
points(nhc$stage_m, nhc$discharge, col = 3)
m <- lm(log(discharge)~log(stage_m), data = nhc)
m_coef <- summary(m)$coefficients[,1]
lines(seq(0,2, by = .01), exp(m_coef[1]) * seq(0,2, by = .01) ^ m_coef[2], lty = 2)


build_powerlaw_rc <- function(l, q){
  m <- lm(log(q)~log(l))
  coef <- summary(m)$coefficients[,1]
  return(coef)
}

rc <- data.frame(site = c("NHC", "UNHC"),
           a = NA,
           b = NA,
           formula = "log(q) = a + b * log(level)")
rc[rc$site =="NHC",2:3] <- build_powerlaw_rc(nhc$level_m, nhc$discharge)
rc[rc$site =="UNHC",2:3] <- build_powerlaw_rc(unhc$level_m, unhc$discharge)

write_csv(rc, "data/rating_curves/modified_ZQ_curves.csv")

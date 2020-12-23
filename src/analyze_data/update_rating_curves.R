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

qq <- read_csv("data/rating_curves/calculated_discharge.csv") 
qq <- qq%>%
  mutate(DateTime_EST = round_date(ymd_hms(paste(date, time), tz = "EST"),
                                   unit = "15 minutes"))

# 4. add stages ####
# for the points that don't have it get from LL
# for the points that do, compare with LL

# NHC
levels <- read_csv("data/rating_curves/NHC_UNHC_corrected_level.csv") %>%
  mutate(DateTime_EST = with_tz(DateTime_UTC, tz = "EST")) %>%
  select(DateTime_EST, level_unhc, level_nhc) 
levels <- read_csv("data/metabolism/processed/PM_lvl.csv") %>%
  mutate(DateTime_EST = with_tz(DateTime_EST, tz = "EST")) %>%
  select(DateTime_EST, level_pm = level_m) %>%
  full_join(levels)
levels <- read_csv("data/metabolism/processed/CBP_lvl.csv") %>%
  mutate(DateTime_EST = with_tz(DateTime_EST, tz = "EST")) %>%
  select(DateTime_EST, level_cbp = level_m) %>%
  full_join(levels)
levels <- read_csv("data/metabolism/processed/WB_lvl.csv") %>%
  mutate(DateTime_EST = with_tz(DateTime_EST, tz = "EST")) %>%
  select(DateTime_EST, level_wb = level_m) %>%
  full_join(levels)
levels <- read_csv("data/metabolism/processed/WBP_lvl.csv") %>%
  mutate(DateTime_EST = with_tz(DateTime_EST, tz = "EST")) %>%
  select(DateTime_EST, level_wbp = level_m) %>%
  full_join(levels)
levels <- read_csv("data/metabolism/processed/PWC.csv") %>%
  mutate(DateTime_EST = with_tz(DateTime_EST, tz = "EST")) %>%
  select(DateTime_EST, level_pwc = level_m) %>%
  full_join(levels) %>%
  arrange(DateTime_EST)
par(mfrow = c(2,3))
plot((levels$level_nhc), (levels$level_unhc), pch = 20, xlim = c(.5, .7), ylim = c(0.4, 0.6))
unhcab <- lm(log(level_unhc) ~ log(level_nhc), levels)$coefficients
abline(unhcab)
plot(log(levels$level_nhc), log(levels$level_pm), pch = 20, col = 2) 
pmab <- lm(log(level_pm) ~log(level_nhc), levels)$coefficients
abline(pmab)
plot(log(levels$level_nhc), log(levels$level_cbp), pch = 20, col = 3)
cbpab <- lm(log(level_cbp) ~log(level_nhc), levels)$coefficients
abline(cbpab)
plot(log(levels$level_nhc), log(levels$level_wb), pch = 20, col = 4)
wbab <- lm(log(level_wb) ~log(level_nhc), levels)$coefficients
abline(wbab)
plot(log(levels$level_nhc), log(levels$level_wbp), pch = 20, col = 5)
wbpab <- lm(log(level_wbp) ~log(level_nhc), levels)$coefficients
abline(wbpab)
plot(log(levels$level_nhc), log(levels$level_pwc), pch = 20, col = 6)
pwcab <- lm(log(level_pwc) ~log(level_nhc), levels)$coefficients
abline(pwcab)

# grab level for the day/time of Q measurement
qql <- qq %>% 
  left_join(levels) %>%
  mutate(level = case_when(site == "NHC" ~ level_nhc,
                           site == "PM" ~ level_pm,
                           site == "CBP" ~ level_cbp,
                           site == "WB" ~ level_wb,
                           site == "WBP" ~ level_wbp,
                           site == "PWC" ~ level_pwc,
                           site == "UNHC" ~ level_unhc),
         level_mod = case_when(site == "PM" ~ 
                                 exp(pmab[1] + pmab[2] * log(level_nhc)),
                               site == "CBP" ~
                                 exp(cbpab[1] + cbpab[2] * log(level_nhc)),
                               site =="WB" ~ 
                                 exp(wbab[1] + wbab[2] * log(level_nhc)),
                               site == "WBP" ~
                                 exp(wbpab[1] + wbpab[2] * log(level_nhc)),
                               site == "UNHC" ~
                                 exp(unhcab[1] + unhcab[2] * log(level_nhc)))) %>%
  select(-level_pm, -level_cbp, -level_wb,
         -level_wbp, -level_pwc, -notes)

write_csv(qql, "data/rating_curves/calculated_discharge_with_levels.csv")
qql <- read_csv("data/rating_curves/calculated_discharge_with_levels.csv")

par(mfrow = c(1,1))
plot(qql$stage_m, qql$level)
abline(0,1)


# # the 9/22/2016 date is problematic
# w <- which(tmp$date == as.Date("2016-09-22"))
# tmp$level_m[w] <- NA
ggplot(qql, aes(stage_m, level, color = site)) +
  geom_point()+
  geom_smooth(method = "lm")

qq_l <- qql %>%
  mutate(level = case_when(is.na(level) ~  stage_m,
           #                level > 1.6 ~ stage_m,
                           !is.na(level) ~ level))

# inspect data ####
qq_l %>%
  filter(site == "UNHC") %>%
  ggplot(aes(level, discharge)) +
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
# qq_l %>%
#   rename(DateTime_UTC = DateTime_EST) %>%
# write_csv(qq_l, "data/rating_curves/calculated_discharge_modified.csv")
# Build Rating Curves ####
qq_l <- read_csv("data/rating_curves/calculated_discharge_modified.csv")
nhc <- qq_l %>%
  filter(site =="NHC")
nhc <- bind_rows(nhc, data.frame(level = .53, discharge = .05)) # this is the min stage and the min flow
unhc <- qq_l %>%
  filter(site =="UNHC")%>%
  slice(c(-2,-3))
  unhc$level[1]<- c(.38)
m <- lm(log(discharge)~log(level), data = nhc)
m_coef <- summary(m)$coefficients[,1]
plot(nhc$level, nhc$discharge)
lines(seq(0,2, by = .01), exp(m_coef[1]) * seq(0,2, by = .01) ^ m_coef[2])
m <- lm(log(discharge)~log(stage_m), data = nhc)
m_coef <- summary(m)$coefficients[,1]
points(nhc$level, nhc$discharge,pch = 20, col = 5)
lines(seq(0,2, by = .01), exp(m_coef[1]) * seq(0,2, by = .01) ^ m_coef[2], col = 5)

m <- lm(log(discharge)~log(level), data = unhc)
m_coef <- summary(m)$coefficients[,1]
plot(unhc$level, unhc$discharge)
lines(seq(0,2, by = .01), exp(m_coef[1]) * seq(0,2, by = .01) ^ m_coef[2])
lines(seq(0,2, by = .01), exp(m_coef[1]) * seq(0,2, by = .01) ^ m_coef[2], col = 2)
points(unhc$level, unhc$discharge,pch = 20, col = 2)

build_powerlaw_rc <- function(l, q, site){
  m <- lm(log(q)~log(l))
  coefs <- summary(m)$coefficients[,1]
  out<- data.frame(site = site,
                   a = coefs[1, drop = T],
                   b = coefs[2, drop = T],
                   formula = "log(q) = a + b * log(level)",
                   min_Q = min(q, na.rm = T),
                   max_Q = max(q, na.rm = T),
                   min_l = min(l, na.rm = T),
                   max_l = max(l, na.rm = T),
                   n = length(!is.na(q)),
                   row.names = NULL)
  return(out)
}

rc <- data.frame()
rc <- bind_rows(rc, build_powerlaw_rc(nhc$level, nhc$discharge, "NHC"))
#rc <- bind_rows(rc, build_powerlaw_rc(nhc$stage_m, nhc$discharge, "NHC"))
rc <- bind_rows(rc, build_powerlaw_rc(unhc$level, unhc$discharge, "UNHC"))

rc$notes <- c("with low flow min stage point included (based on UNHC)",
              "with two high flow low level points excluded",
              "all Q measurements")

write_csv(rc, "data/rating_curves/modified_ZQ_curves.csv")
# look across sites on the two dates where multiple were measured. ####

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
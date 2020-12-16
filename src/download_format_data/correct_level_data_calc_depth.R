# Correct level data from NHC sites and 
# calculate average depth as best as we can for now


# this process is idosincratic for each site and will be based
# on times that the sensor was moved or when there is bad data.

library(tidyverse)
library(xts)
library(dygraphs)
library(lubridate)
library(zoo)
setwd('C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism/data')

norm_depths = data.frame(site = c("NHC","PM","CBP","WB","WBP","PWC","UNHC"),
                         depth_cm = c(50, 65, 54, 68, 70, 70, 30),
                         level_cm = c(80, 90, 75, 60, 90, 90, 30))

# 1. Load Data ####
# load field notes
fnotes <- read_csv("water_chemistry/nhc_ysi_data.csv") %>%
  mutate(DateTime_EST = round_date(as.POSIXct(paste(date, as.character(time)), 
                                   format = "%m/%d/%Y %H:%M:%S", tz = "EST"),
                                   unit = '15 minutes'))
spdepths <- read_csv("water_chemistry/streampulse_depth_measurements_2017.csv") %>%
  filter(!is.na(Date)) %>%
  mutate(DateTime_EST = round_date(as.POSIXct(paste(Date, as.character(Time)), 
                                   format = "%m.%d.%Y %H:%M", tz = "EST"),
                                   unit = '15 minutes'),
         waterdepth_cm = stage_at_sensor *100) %>%
  select(site = Site, DateTime_EST, waterdepth_cm) %>%
  filter(site %in% c("NHC", "UNHC"))


fnotes <- bind_rows(fnotes, spdepths)
depths <- fnotes %>%
  mutate(date = as.Date(DateTime_EST)) %>%
  select(site, date, waterdepth_cm) %>%
  filter(date > as.Date("2018-01-01") & !(site %in% c("PWC","MC751"))) %>%
  group_by(date, site) %>%
  summarize_all(mean, na.rm = T) %>%
  ungroup() %>%
  pivot_wider(names_from = site, values_from = waterdepth_cm) %>%
  mutate_all(na.approx, na.rm = F) %>%
  pivot_longer(cols = -date, names_to = "site", values_to = "depth_cm")

ggplot(depths, aes(x = date, y = depth_cm, color = site)) +
  geom_line()
ggplot(spdepths, aes(x = DateTime_EST, y = waterdepth_cm, color = site)) +
  geom_point()
# CBP ####
dat <- read_csv("metabolism/processed/CBP.csv") %>%
  mutate(DateTime_EST = force_tz(DateTime_EST, tz = "EST"))
fn <- fnotes %>% filter(site == "CBP")

# explore data
lvl <- xts(x = dat$level_m, order.by = dat$DateTime_EST)
dygraph(lvl)

plot(dat$DateTime_EST, dat$level_m, type = "l")
points(fn$DateTime_EST, fn$waterdepth_cm/100, col = 2)
abline(h = .27) # this looks like the level where points below are just sensor out of water
dat$level_m[which(dat$level_m < 0.27)] <- NA

# examine field notes for times the sensor was moved:
select(fn, DateTime_EST, sensoroffset_cm, waterdepth_cm, notes) %>%
  data.frame()

# first, move the section before May 8 up to water depths
s1 <- fn %>%
  select(DateTime_EST, waterdepth_cm, sensoroffset_cm, depthtosensor_cm) %>%
  full_join(dat) %>% 
  filter(DateTime_EST <= as.POSIXct("2019-06-03 16:00:00", tz = "EST")) %>%
  mutate(delta = waterdepth_cm/100 - level_m,
         level_m = level_m + mean(delta, na.rm = T)) %>%
  select(DateTime_EST, level_m) %>%
  arrange(DateTime_EST)

s1$level_m[s1$DateTime_EST >= as.POSIXct("2019-05-08 15:30:00", tz = "EST")] <- .095 + 
  s1$level_m[s1$DateTime_EST >= as.POSIXct("2019-05-08 15:30:00", tz = "EST")]


# The data in Sept look shifted down
# fix this by a linear shift
d1 <- as.POSIXct("2019-09-10 21:15:00", tz = "EST")
d2 <- as.POSIXct("2019-10-05 19:00:00", tz = "EST")

c1 <- as.POSIXct("2019-09-13 00:15:00", tz = "EST")
c2 <- as.POSIXct("2019-10-05 01:45:00", tz = "EST")
df <- data.frame(d = c(d1,c1,c2,d2),
                 l = c(dat$level_m[which(dat$DateTime_EST == d1)], NA, NA, 
                       dat$level_m[which(dat$DateTime_EST == d2)])) %>%
  transform(l = na.approx(l, d))

df1 <- data.frame(d = seq(from = c1, to = c2, by = "15 min"))
delta = c(dat$level_m[which(dat$DateTime_EST == c1)] - df$l[2], 
          .3 - df$l[3])
df1$delta <- rep(NA, nrow(df1))
df1$delta[1] <- delta[1]
df1$delta[nrow(df1)] <- delta[2]
s2 <- df1 <- transform(df1, delta = na.approx(delta, d)) %>%
  rename(DateTime_EST = d) %>%
  left_join(dat) %>%
  mutate(level_m = level_m - delta) %>% 
  select(DateTime_EST, level_m)

plot(dat$DateTime_EST, dat$level_m, type = "l")
points(fn$DateTime_EST, fn$waterdepth_cm/100, col = 2)
lines(s1$DateTime_EST, s1$level_m, col = 3)  
lines(s2$DateTime_EST, s2$level_m, col = 3)  
ss <- bind_rows(s1, s2) %>%
  rename(l = level_m) %>%
  right_join(dat)
ss$level_m[which(!is.na(ss$l))] <- ss$l[which(!is.na(ss$l))]
dat <- ss %>% select(-l)
# last section is the one at the end, I will do another linear shift.

d1 <- as.POSIXct("2020-02-12 17:45:00", tz = "EST")
delta <- dat$level_m[which(dat$DateTime_EST == d1 - 60*60)] -
  dat$level_m[which(dat$DateTime_EST == d1)]
w <- which(dat$DateTime_EST >= d1)
dat$level_m[w] <- delta + dat$level_m[w]
dat$level_m[(w-3):(w-1)] <- NA
dat$level_m[(w-4):w] <- na.approx(dat$level_m[(w-4):w])

ss <- dat %>% 
  filter(DateTime_EST >= as.POSIXct("2019-06-01 00:00:00", tz = "EST")) %>%
  left_join(fn, by = "DateTime_EST") %>%
  mutate(delta = level_m - waterdepth_cm/100,
         l = level_m - mean(delta, na.rm = T)) %>%
  select(DateTime_EST, l)
ss <- right_join(ss, dat)
ss$level_m[which(!is.na(ss$l))] <- ss$l[which(!is.na(ss$l))]
dat <- ss %>% select(-l) %>%
  arrange(DateTime_EST)

lines(dat$DateTime_EST, dat$level_m, col = 4)

# get rid of a few more impossible points:
dat <- dat %>% 
  mutate(depth = ifelse(dat$discharge < 0.01 & dat$level_m > 1, NA, depth),
         discharge = ifelse(dat$discharge < 0.01 & dat$level_m > 1, NA, discharge))

# # last thing, the final section needs to shift down so that the overbank flow 
# # hits at the same place:
# 
# ggplot(dat, aes(level_m, depth)) +
#   geom_point(aes(col = DateTime_EST)) +
#   ylim(0,2) +
#   geom_vline(xintercept = c(.95, 1.12))
# 
# 

# CBP depth ####
# calc depth function in stream metabolizer can be calibrated by providing a
# depth at unit discharge, this is the discharge = 1 m3s
# unit discharge is at a stage of ~ 1 m in this stream, 
# zero flow corresponds to a stage of ~ 0.25 m
# Hall declared standard water level to be 50 cm above zero flow,
# which was the normal spring flow and corresponded to a depth of 0.54 m at CBP

ggplot(dat, aes(level_m, depth)) +
  geom_point(aes(col = DateTime_EST)) +
  ylim(0,2) + 
  geom_vline(xintercept = c(.71,1.05)) +
  geom_hline(yintercept = .88)
# Given all that, I am going to say level 75 should correspond to depth 54. 
# in my data, level 75 corresponds to Q of 0.25

# using the default exponent from Leopold and Maddock 1953, f = 0.294, 
# I can calculate the average depth at unit Q from the above observation:
# d = c Q^f   0.54 = c * 0.25^0.294

c <- .54/(.25 ^.294)
tmp <- dat %>% 
  mutate(depth2 = calc_depth(discharge, c)) 
ggplot(tmp, aes(level_m, depth2)) +
  geom_point() +
  geom_point(aes(y = depth), col = "blue") +
  ylim(0,2) 

dat <- tmp %>% 
  select(-depth) %>%
  rename(depth = depth2)
 
# once we reach overbank flow, I will assume that changes in level continue
# to cause a linear change in depth with the same slope as before.
# This relationship shifts part way through the year, so I'll do it twice.
# The date where this divide happens is in Oct during low flow, I'll use
# the end of the break at Oct 6

dl <- dat %>% 
  filter(DateTime_EST > as.POSIXct("2019-06-01 00:00:00", tz = "EST")) 

ggplot(dl, aes(level_m, depth)) +
  geom_point(aes(col = DateTime_EST)) +
  ylim(0,2) +
  geom_vline(xintercept = 1.12) +
  geom_hline(yintercept = .9)

dl <- dat %>% 
  filter(DateTime_EST > as.POSIXct("2019-06-01 00:00:00", tz = "EST") &
           level_m < 1.12 & depth < 0.9) 

logEstimate <- lm(depth~log(level_m),data=dl)
x <- seq(0.1, 2, by = .1)
plot(dl$level_m, dl$depth, ylim = c(0,3), xlim = c(0,2))
coef <- logEstimate$coefficients
curve(coef[1] + coef[2] * log(x), add=TRUE, col=4, lwd = 3)

# linEstimate <- lm(depth~level_m,data=dl)
# x <- seq(0.1, 2, by = .1)
# plot(dl$level_m, dl$depth, ylim = c(0,1.5), xlim = c(0,2))
# coef <- linEstimate$coefficients
# curve(coef[1] + coef[2] * (x), add=TRUE, col=4, lwd = 3)
# 
dat$depth2 <- coef[1] + coef[2] * log(dat$level_m)
# dat$depth3 <- coef[1] + coef[2] * (dat$level_m)

plot(dat$DateTime_EST, dat$depth, type = "l", ylim = c(0,2))
lines(dat$DateTime_EST, dat$depth2, col = 3)
#lines(dat$DateTime_EST, dat$depth3, col = 4)
lines(dat$DateTime_EST, dat$level_m, col = 2)

dat_lvl <- dat %>%
  mutate(depth = ifelse(depth2 < 0.15, 0.15, depth2)) %>%
  select(-depth2)

write_csv(dat_lvl, "metabolism/processed/CBP_lvl.csv")

# PM  ####

dat <- read_csv("metabolism/processed/PM.csv") %>%
  mutate(DateTime_EST = force_tz(DateTime_EST, tz = "EST"))
fn <- fnotes %>% filter(site == "PM")

# explore data
lvl <- xts(x = dat$level_m, order.by = dat$DateTime_EST)
dygraph(lvl)

plot(dat$DateTime_EST, dat$level_m, type = "l")
points(fn$DateTime_EST, fn$waterdepth_cm/100, col = 2)
abline(h = .55) # this looks like the level where points below are just sensor out of water
dat$level_m[which(dat$level_m < 0.55)] <- NA

# examine field notes for times the sensor was moved:
select(fn, DateTime_EST, sensoroffset_cm, waterdepth_cm, notes) %>%
  data.frame()

# first, move the section before May 8 up to water depths
s1 <- dat %>% 
  filter(DateTime_EST <= as.POSIXct("2019-05-08 13:45:00")) %>%
  mutate(level_m = level_m + .22) %>%
  select(DateTime_EST, level_m ) %>%
  arrange(DateTime_EST)


# Chunks that need to be snapped with the rest of the data: 
d1 <- as.POSIXct("2019-11-20 16:45:00", tz = "EST") 
d2 <- as.POSIXct("2019-11-26 19:00:00", tz = "EST") 

df <- data.frame(d = c(d1, (d1 + 14*15*60), (d2 - 3*15*60), d2),
                 l = c(dat$level_m[which(dat$DateTime_EST == d1)], NA, NA, 
                       dat$level_m[which(dat$DateTime_EST == d2)])) %>%
  transform(l = na.approx(l, d))

df1 <- data.frame(d = seq(from = d1 + 14*15*60, to = d2 - 3*15*60, by = "15 min"))
delta = c(dat$level_m[which(dat$DateTime_EST == d1 + 14*15*60)] - df$l[2], 
          dat$level_m[which(dat$DateTime_EST == d2 - 3*15*60)] - df$l[3])
df1$delta <- rep(NA, nrow(df1))
df1$delta[1] <- delta[1]
df1$delta[nrow(df1)] <- delta[2]

s2 <- df1 <- transform(df1, delta = na.approx(delta, d)) %>%
  rename(DateTime_EST = d) %>%
  left_join(dat) %>%
  mutate(level_m = level_m - delta) %>% 
  select(DateTime_EST, level_m)

d1 <- as.POSIXct("2020-02-27 16:00:00", tz = "EST") 
d2 <- as.POSIXct("2020-02-28 04:30:00", tz = "EST") 

df <- data.frame(d = c(d1, (d1 + 2*15*60), (d2 - 1*15*60), d2),
                 l = c(dat$level_m[which(dat$DateTime_EST == d1)], NA, NA, 
                       dat$level_m[which(dat$DateTime_EST == d2)])) %>%
  transform(l = na.approx(l, d))

df1 <- data.frame(d = seq(from = d1 + 2*15*60, to = d2 - 1*15*60, by = "15 min"))
delta = c(dat$level_m[which(dat$DateTime_EST == d1 + 2*15*60)] - df$l[2], 
          dat$level_m[which(dat$DateTime_EST == d2 - 1*15*60)] - df$l[3])
df1$delta <- rep(NA, nrow(df1))
df1$delta[1] <- delta[1]
df1$delta[nrow(df1)] <- delta[2]

s3 <- df1 <- transform(df1, delta = na.approx(delta, d)) %>%
  rename(DateTime_EST = d) %>%
  left_join(dat) %>%
  mutate(level_m = level_m - delta) %>% 
  select(DateTime_EST, level_m)

plot(dat$DateTime_EST, dat$level_m, type = "l")
points(fn$DateTime_EST, fn$waterdepth_cm/100, col = 2)
lines(s1$DateTime_EST, s1$level_m, col = 3)  
lines(s2$DateTime_EST, s2$level_m, col = 3)  
lines(s3$DateTime_EST, s3$level_m, col = 3)  

ss <- bind_rows(s1, s2, s3) %>%
  rename(l = level_m)

# move data up tp match point measures

dat <- dat %>% left_join(ss) %>%
  mutate(level_m = ifelse(!is.na(l), l, level_m)) %>%
  select(-l) 

dat <- fn %>% select(DateTime_EST, waterdepth_cm) %>%
  right_join(dat) %>%
  mutate(delta = waterdepth_cm/100 - level_m,
         level_m = level_m + mean(delta, na.rm = T)) %>%
  select(-waterdepth_cm, -delta) %>%
  arrange(DateTime_EST)

lines(dat$DateTime_EST, dat$level_m, col = 4)

plot(dat$level_m, dat$depth, ylim = c(0,1))
abline(h = .4, v = .9)

dat <- dat %>%
  mutate(discharge = ifelse(depth > .4 & level_m < .9, NA, discharge),
         depth = ifelse(depth >.4 & level_m <.9, NA, depth))

# # compare to CBP q to get coefficient for unit discharge
# cbp <- read_csv("metabolism/processed/CBP_lvl.csv") %>%
#   select(DateTime_EST, cbp.depth = depth, cbp.q = discharge) %>%
#   left_join(dat)
# plot(cbp$cbp.q, cbp$discharge, log = "xy")
# dat$depth2 <- calc_depth(dat$discharge, c)
# plot(cbp$cbp, cbp$depth)
# abline = c(0,1)
# 
# plot(dat$depth, dat$depth2)
# plot(dat$level_m, dat$depth, ylim = c(0,2.5))

# PM depth ####
# calc depth function in stream metabolizer can be calibrated by providing a
# depth at unit discharge, this is the discharge = 1 m3s
# unit discharge is at a stage of ~ 1 m in this stream, 
# zero flow corresponds to a stage of ~ 0.25 m
# Hall declared standard water level to be 50 cm above zero flow,
# which was the normal spring flow and corresponded to a depth of 0.54 m at CBP
# .35 at UNHC, .7 at WBP, .68 at WB, .65 at PM, .5 at NHC
# normal depths are CBP = PM = NHC = WBP = 75, WB = 50, UNHC = 25
nd <- norm_depths %>%
  filter(site =="PM")

ggplot(dat, aes(level_m, log(discharge))) +
  geom_point(aes(col = DateTime_EST)) +
 # ylim(0,2) + 
  geom_vline(xintercept =  1) + 
  geom_hline(yintercept = log(.8))

# Given all that, I am going to say level 1 should correspond to depth 65. 
# in my data, level 100 corresponds to Q of 0.8

# using the default exponent from Leopold and Maddock 1953, f = 0.294, 
# I can calculate the average depth at unit Q from the above observation:
# d = c Q^f   0.65 = c * 0.8^0.294

c <- .65/(.8 ^.294)
tmp <- dat %>% 
  mutate(depth2 = calc_depth(discharge, c)) 
ggplot(tmp, aes(level_m, depth2)) +
  geom_point() +
  geom_point(aes(y = depth), col = "blue") +
  ylim(0,2) 

dat <- tmp %>% 
  select(-depth) %>%
  rename(depth = depth2)

# once we reach overbank flow, I will assume that changes in level continue
# to cause a linear change in depth with the same slope as before.
# This relationship shifts part way through the year, so I'll do it twice.
# The date where this divide happens is in Oct during low flow, I'll use
# the end of the break at Oct 6

dl <- dat %>% 
  filter(DateTime_EST > as.POSIXct("2019-06-01 00:00:00", tz = "EST")) 

ggplot(dat, aes(level_m, depth)) +
  geom_point(aes(col = DateTime_EST)) +
  ylim(0,2) +
  geom_vline(xintercept = 1.12) +
  geom_hline(yintercept = .9)

dl <- dat %>% 
  filter(DateTime_EST > as.POSIXct("2019-06-01 00:00:00", tz = "EST") &
           level_m < 1.12 & depth < 0.9) 

logEstimate <- lm(depth~log(level_m),data=dl)
x <- seq(0.1, 2, by = .1)
plot(dl$level_m, dl$depth, ylim = c(0,3), xlim = c(0,2))
coef <- logEstimate$coefficients
curve(coef[1] + coef[2] * log(x), add=TRUE, col=4, lwd = 3)

# linEstimate <- lm(depth~level_m,data=dl)
# x <- seq(0.1, 2, by = .1)
# plot(dl$level_m, dl$depth, ylim = c(0,1.5), xlim = c(0,2))
# coef <- linEstimate$coefficients
# curve(coef[1] + coef[2] * (x), add=TRUE, col=4, lwd = 3)
# 
dat$depth2 <- coef[1] + coef[2] * log(dat$level_m)
# dat$depth3 <- coef[1] + coef[2] * (dat$level_m)

plot(dat$DateTime_EST, dat$depth, type = "l", ylim = c(0,2))
lines(dat$DateTime_EST, dat$depth2, col = 3)
#lines(dat$DateTime_EST, dat$depth3, col = 4)
lines(dat$DateTime_EST, dat$level_m, col = 2)

dat_lvl <- dat %>%
  mutate(depth = ifelse(depth2 < 0.15, 0.15, depth2)) %>%
  select(-depth2)



write_csv(dat, "metabolism/processed/PM_lvl.csv")
# WB  ####

dat <- read_csv("metabolism/processed/WB.csv") %>%
  mutate(DateTime_EST = force_tz(DateTime_EST, tz = "EST"))
fn <- fnotes %>% filter(site == "WB")

# explore data
lvl <- xts(x = dat$level_m, order.by = dat$DateTime_EST)
dygraph(lvl)

plot(dat$DateTime_EST, dat$level_m, type = "l")
points(fn$DateTime_EST, fn$waterdepth_cm/100, col = 2)
abline(h = .29) # this looks like the level where points below are just sensor out of water
dat$level_m[which(dat$level_m < 0.29)] <- NA

# examine field notes for times the sensor was moved:
select(fn, DateTime_EST, sensoroffset_cm, waterdepth_cm, notes) %>%
  data.frame()

# first, move the section before May 15 up to water depths
s1 <- dat %>% 
  filter(DateTime_EST <= as.POSIXct("2019-05-15 13:45:00", tz = "EST")) %>%
  mutate(level_m = 
           ifelse(DateTime_EST >= as.POSIXct("2019-03-28 00:00:00", tz = "EST") &
                    DateTime_EST <= as.POSIXct("2019-04-16 00:00:00", tz = "EST"),
                          level_m + .6,
                          level_m + .4)) %>%
  select(DateTime_EST, level_m ) %>%
  arrange(DateTime_EST)


# Chunks that need to be snapped with the rest of the data: 
d1 <- as.POSIXct("2019-06-02 16:45:00", tz = "EST") 
d2 <- as.POSIXct("2019-06-18 20:00:00", tz = "EST") 

s2 <- dat %>% 
  filter(DateTime_EST <= d2 & DateTime_EST >= d1) %>%
  mutate(level_m = level_m + .06) %>%
  select(DateTime_EST, level_m ) %>%
  arrange(DateTime_EST)

plot(dat$DateTime_EST, dat$level_m, type = "l")
points(fn$DateTime_EST, fn$waterdepth_cm/100, col = 2)
lines(s1$DateTime_EST, s1$level_m, col = 3)  
lines(s2$DateTime_EST, s2$level_m, col = 3)  

ss <- bind_rows(s1, s2) %>%
  rename(l = level_m)

# move data up tp match point measures

dat <- dat %>% left_join(ss) %>%
  mutate(level_m = ifelse(!is.na(l), l, level_m)) %>%
  select(-l) 
dd <- as.POSIXct("2019-06-02 00:00:00", tz = "EST")
dat <- fn %>% 
  select(DateTime_EST, waterdepth_cm) %>%
  filter(DateTime_EST > dd)  %>%
  right_join(dat) %>%
  mutate(delta = waterdepth_cm/100 - level_m,
         level_m = level_m + mean(delta, na.rm = T)) %>%
  select(-waterdepth_cm, -delta) %>%
  arrange(DateTime_EST)

lines(dat$DateTime_EST, dat$level_m, col = 4)
ggplot(dat, aes(level_m, depth)) +
  geom_point(aes(col = DateTime_EST)) +
  ylim(0,4) 

write_csv(dat, "metabolism/processed/WB_lvl.csv")

# WBP  ####

dat <- read_csv("metabolism/processed/WBP.csv") %>%
  mutate(DateTime_EST = force_tz(DateTime_EST, tz = "EST"))
fn <- fnotes %>% filter(site == "WBP")

# explore data
lvl <- xts(x = dat$level_m, order.by = dat$DateTime_EST)
dygraph(lvl)

plot(dat$DateTime_EST, dat$level_m, type = "l")
points(fn$DateTime_EST, fn$waterdepth_cm/100, col = 2)
abline(h = .67) # this looks like the level where points below are just sensor out of water
dat$level_m[which(dat$level_m < 0.67)] <- NA

# examine field notes for times the sensor was moved:
select(fn, DateTime_EST, sensoroffset_cm, waterdepth_cm, notes) %>%
  data.frame()

# first, move the section before May 15 up to water depths
s1 <- dat %>% 
  filter(DateTime_EST <= as.POSIXct("2019-05-08 22:45:00", tz = "EST")) %>%
  left_join(fn) %>%
  mutate(delta = waterdepth_cm/100 - level_m,
         level_m = level_m + mean(delta, na.rm = T)) %>%
  select(DateTime_EST, level_m ) %>%
  arrange(DateTime_EST)


# Chunks that need to be snapped with the rest of the data: 
d1 <- as.POSIXct("2019-11-20 18:00:00", tz = "EST") 
d2 <- as.POSIXct("2019-12-04 01:45:00", tz = "EST") 

df <- data.frame(d = c(d1, (d1 + 9*15*60), (d2 - 8*60*60), d2),
                 l = c(dat$level_m[which(dat$DateTime_EST == d1)], NA, NA, 
                       dat$level_m[which(dat$DateTime_EST == d2)])) %>%
  transform(l = na.approx(l, d))

df1 <- data.frame(d = seq(from = d1 + 9*15*60, to = d2 - 8*60*60, by = "15 min"))
delta = c(dat$level_m[which(dat$DateTime_EST == d1 + 14*15*60)] - df$l[2], 
          dat$level_m[which(dat$DateTime_EST == d2 - 3*15*60)] - df$l[3])
df1$delta <- rep(NA, nrow(df1))
df1$delta[1] <- delta[1]
df1$delta[nrow(df1)] <- delta[2]

s2 <- df1 <- transform(df1, delta = na.approx(delta, d)) %>%
  rename(DateTime_EST = d) %>%
  left_join(dat) %>%
  mutate(level_m = level_m - delta) %>% 
  select(DateTime_EST, level_m)

s3 <- dat %>% 
  filter(DateTime_EST >= as.POSIXct("2020-02-12 15:45:00", tz = "EST")) %>%
  mutate(level_m = level_m + .07) %>%
  select(DateTime_EST, level_m ) %>%
  arrange(DateTime_EST)

plot(dat$DateTime_EST, dat$level_m, type = "l")
points(fn$DateTime_EST, fn$waterdepth_cm/100, col = 2)
lines(s1$DateTime_EST, s1$level_m, col = 3)  
lines(s2$DateTime_EST, s2$level_m, col = 3)  
lines(s3$DateTime_EST, s3$level_m, col = 3)  

ss <- bind_rows(s1, s2, s3) %>%
  rename(l = level_m)

# move data up tp match point measures

dat <- dat %>% left_join(ss) %>%
  mutate(level_m = ifelse(!is.na(l), l, level_m)) %>%
  select(-l) 
dd <- as.POSIXct("2019-05-08 22:45:00", tz = "EST")
dat <- fn %>% 
  select(DateTime_EST, waterdepth_cm) %>%
  filter(DateTime_EST > dd)  %>%
  right_join(dat) %>%
  mutate(delta = waterdepth_cm/100 - level_m,
         level_m = ifelse(DateTime_EST > dd,
           level_m + mean(delta, na.rm = T),
           level_m)) %>%
  select(-waterdepth_cm, -delta) %>%
  arrange(DateTime_EST)

lines(dat$DateTime_EST, dat$level_m, col = 4)

ggplot(dat, aes(level_m, log(discharge))) +
  geom_point(aes(col = DateTime_EST))# +
  ylim(0,4) 

write_csv(dat, "metabolism/processed/WBP_lvl.csv")

# PWC  ####

dat <- read_csv("metabolism/processed/PWC.csv") %>%
  mutate(DateTime_EST = force_tz(DateTime_EST, tz = "EST"))
fn <- fnotes %>% filter(site == "PWC")

# explore data
lvl <- xts(x = dat$level_m, order.by = dat$DateTime_EST)
dygraph(lvl)

plot(dat$DateTime_EST, dat$level_m, type = "l")
points(fn$DateTime_EST, fn$waterdepth_cm/100, col = 2)
abline(h = .67) # this looks like the level where points below are just sensor out of water
dat$level_m[which(dat$level_m < 0.67)] <- NA

# examine field notes for times the sensor was moved:
select(fn, DateTime_EST, sensoroffset_cm, waterdepth_cm, notes) %>%
  data.frame()

# first, move the section before May 15 up to water depths
s1 <- dat %>% 
  filter(DateTime_EST <= as.POSIXct("2019-05-08 22:24:00", tz = "EST")) %>%
  mutate(level_m = ifelse(DateTime_EST <= as.POSIXct("2019-03-27 18:45:00", tz = "EST"),
                          level_m + .22,
                          level_m + .14)) %>%
  select(DateTime_EST, l = level_m ) %>%
  arrange(DateTime_EST)


plot(dat$DateTime_EST, dat$level_m, type = "l")
points(fn$DateTime_EST, fn$waterdepth_cm/100, col = 2)
lines(s1$DateTime_EST, s1$l, col = 3)  

dat <- dat %>% left_join(s1) %>%
  mutate(level_m = ifelse(!is.na(l), l, level_m)) %>%
  select(-l) 

ggplot(dat, aes(level_m, depth)) +
  geom_point(aes(col = DateTime_EST))# +
ylim(0,4) 

write_csv(dat, "metabolism/processed/PWC_lvl.csv")

# NHC ####
dat <- read_csv("metabolism/processed/NHC.csv") %>%
  mutate(DateTime_EST = force_tz(DateTime_EST, tz = "EST"))
fn <- fnotes %>% 
  filter(site == "NHC",
         !is.na(waterdepth_cm))%>%
  select(DateTime_EST, waterdepth_cm) 

fn <- left_join(fn, dat, by = "DateTime_EST") %>%
  mutate(delta = waterdepth_cm/100 - level_m)
# explore data
lvl <- xts(x = dat$level_m, order.by = dat$DateTime_EST)
dygraph(lvl)
abline(h = .4) # this looks like the level where points below are just sensor out of water
dat$level_m[which(dat$level_m < 0.4)] <- NA

plot(fn$DateTime_EST, fn$waterdepth_cm/100, col = 2)

plot(fn$DateTime_EST, fn$waterdepth_cm/100, col = 2, pch = 19, ylim = c(.5, 1.5))
lines(dat$DateTime_EST, dat$level_m)

par(new = T)
plot(fn$DateTime_EST, fn$delta, type = "l", col = 3)
axis(4)

hist(fn$delta)
plot(fn$level_m, fn$delta)
abline(h=0)

plot(fn$waterdepth_cm/100, fn$level_m)
abline(0,1)

dat$level_m <- dat$level_m + mean(fn$delta, na.rm = T)

# 12/04/2020 I am not saving this file yet, I am waiting to see if there is 
# updated field notes first.
# write_csv(dat, "metabolism/processed/NHC.csv")

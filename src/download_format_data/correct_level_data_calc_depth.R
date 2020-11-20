# Correct level data from NHC sites and 
# calculate average depth as best as we can for now


# this process is idosincratic for each site and will be based
# on times that the sensor was moved or when there is bad data.

library(tidyverse)
library(xts)
library(dygraphs)
library(zoo)
setwd('C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism/data')

# load field notes
fnotes <- read_csv("water_chemistry/nhc_ysi_data.csv") %>%
  mutate(DateTime_EST = round_date(as.POSIXct(paste(date, as.character(time)), 
                                   format = "%m/%d/%Y %H:%M:%S", tz = "EST"),
                                   unit = '15 minutes'))

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
  filter(DateTime_EST <= as.POSIXct("2019-06-03 16:00:00")) %>%
  mutate(delta = waterdepth_cm/100 - level_m,
         level_m = level_m + mean(delta, na.rm = T)) %>%
  select(DateTime_EST, level_m) %>%
  arrange(DateTime_EST)

s1$level_m[s1$DateTime_EST >= as.POSIXct("2019-05-08 15:30:00")] <- .095 + 
  s1$level_m[s1$DateTime_EST >= as.POSIXct("2019-05-08 15:30:00")]


# The data in Sept look shifted down
# fix this by a linear shift
d1 <- as.POSIXct("2019-09-10 21:15:00", tz = "EST")
d2 <- as.POSIXct("2019-10-05 19:00:00", tz = "EST")

c1 <- as.POSIXct("2019-09-13 00:15:00", tz = "EST")
c2 <- as.POSIXct("2019-10-04 01:45:00", tz = "EST")
df <- data.frame(d = c(d1,c1,c2,d2),
                 l = c(dat$level_m[which(dat$DateTime_EST == d1)], NA, NA, 
                       dat$level_m[which(dat$DateTime_EST == d2)])) %>%
  transform(l = na.approx(l, d))

df1 <- data.frame(d = seq(from = c1, to = c2, by = "15 min"))
delta = c(dat$level_m[which(dat$DateTime_EST == c1)] - df$l[2], 
          dat$level_m[which(dat$DateTime_EST == c2)] - df$l[3])
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

# 
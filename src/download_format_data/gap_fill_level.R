# gap fill level data for unhc and nhc
# plot_pres <- function(NHC, waterpres = "level_m", airpres = "waterdepth_m", extra = NA){
#   if(is.na(extra)){
#   NHC %>% select(all_of(waterpres), all_of(airpres)) %>%
#     xts(order.by = NHC$DateTime_UTC) %>%
#     dygraph() %>%
#     dyRangeSelector()
#   } else {
#   NHC %>% select(all_of(waterpres), all_of(airpres), all_of(extra)) %>%
#     xts(order.by = NHC$DateTime_UTC) %>%
#     dygraph() %>%
#     dyRangeSelector()
#     
#   }
# }

# read in corrected level data:
lvl_long <- read_csv("NHC_2019_metabolism/data/rating_curves/all_sites_level_corrected.csv") %>%
  group_by(site, DateTime_UTC) %>%
  summarize_all(mean, na.rm = T) %>%
  ungroup()

# ggplot(lvl_long, aes(DateTime_UTC, level_m, color = site)) +
#   geom_line()
lvl <- lvl_long %>%
  pivot_wider(names_from = "site", values_from = "level_m") %>%
  arrange(DateTime_UTC) %>% 
  select(DateTime_UTC, NHC, UNHC) %>%
  mutate(month = as.numeric(format(DateTime_UTC, "%m")))

lvl_s <- lvl
n <- nrow(lvl)
lvl_s$UNHC <- c( rep(NA_real_, 8),lvl_s$UNHC[1:(n-8)])

mnhc <- lm(NHC ~ UNHC, data = lvl_s)
munhc <- lm(UNHC ~ NHC, data = lvl_s)
lvl_s$nhc_mod <- predict(mnhc, newdata = lvl_s)
lvl_s$unhc_mod <- predict(munhc, newdata = lvl_s)


# qm <- lvl_s %>%
#   mutate(UNHC = c(lvl_s$UNHC[9:n], rep(NA_real_, 8)),
#          unhc_mod = c(lvl_s$unhc_mod[9:n], rep(NA_real_, 8)),
#          UNHC = case_when((DateTime_UTC > ymd_hms("2017-08-01 04:00:00") &
#                              DateTime_UTC < ymd_hms("2018-01-24 23:15:00")) ~
#                             NA_real_,
#                           TRUE ~ UNHC)) %>%
#   rename(level_nhc = NHC, level_unhc = UNHC) %>%
#   mutate(nhc_mod = ifelse(is.na(level_nhc), nhc_mod, level_nhc),
#          unhc_mod = ifelse(is.na(level_unhc), unhc_mod, level_unhc),
#          notes = case_when(is.na(level_nhc) ~ "nhc modeled",
#                            is.na(level_unhc) ~ "unhc modeled"))
qm <- lvl_s %>%
  mutate(UNHC = c(lvl_s$UNHC[9:n], rep(NA_real_, 8)),
         unhc_mod = c(lvl_s$unhc_mod[9:n], rep(NA_real_, 8)),
         NHC = case_when(((DateTime_UTC > ymd_hms('2017-10-15 00:00:00') &
                            DateTime_UTC < ymd_hms("2018-01-24 23:15:00"))) ~
                            NA_real_,
                          TRUE ~ NHC)) %>%
  rename(level_nhc = NHC, level_unhc = UNHC) %>%
  mutate(nhc_mod = ifelse(is.na(level_nhc), nhc_mod, level_nhc),
         unhc_mod = ifelse(is.na(level_unhc), unhc_mod, level_unhc),
         notes = case_when(is.na(level_nhc) ~ "nhc modeled",
                           is.na(level_unhc) ~ "unhc modeled"))

# snap interpolated levels  to their neighbors ####

nhc_gaps <- rle_custom(is.na(qm$level_nhc))
unhc_gaps <- rle_custom(is.na(qm$level_unhc))

# fill in and plot modeled data
# par(mfrow = c(1,1))
# plot(qm$DateTime_UTC, qm$nhc_mod, log="y", main = "NHC", pch = 20)
# points(qm$DateTime_UTC[qm$notes == "nhc modeled"],
#        qm$nhc_mod[qm$notes=="nhc modeled"], pch = 20, col = "red")
# plot(qm$DateTime_UTC, qm$unhc_mod, log="y", main = "UNHC", pch = 20)
# points(qm$DateTime_UTC[qm$notes == "unhc modeled"],
#        qm$unhc_mod[qm$notes=="unhc modeled"], pch = 20, col = "red")

# find endpoints of measured and modeled data in the gaps
nhc_gaps <- nhc_gaps[nhc_gaps$values==1,]
unhc_gaps <- unhc_gaps[unhc_gaps$values==1,]

# Don't allow interpolated discharge to be lower than NHC min flow
# NHCmin <- min(qm$level_nhc, na.rm=T)
# m <- min(qm$nhc_mod, na.rm=T)
# it isnt.

for(i in 1:nrow(nhc_gaps)){
  a <- nhc_gaps[i,]$starts
  b <- nhc_gaps[i,]$stops
  
  if(a==1) next
  startdiff <- qm$level_nhc[a-1] - qm$nhc_mod[a]
  if(b==nrow(qm)){
    enddiff <- startdiff
  } else{
    enddiff <- qm$level_nhc[b+1] - qm$nhc_mod[b]
  }
  if(is.na(startdiff)||is.na(enddiff)) next
  
  diffQ <- seq(startdiff, enddiff, length.out=nhc_gaps[i,]$lengths)
  
  qm$nhc_mod[a:b]<- qm$nhc_mod[a:b] + diffQ
}

# snap UNHC gaps
# Don't allow interpolated discharge to be lower than NHC min flow
UNHCmin <- min(qm$unhc_mod, na.rm=T)
w <- which.min(qm$unhc_mod)

for(i in 1:nrow(unhc_gaps)){
  a <- unhc_gaps[i,]$starts
  b <- unhc_gaps[i,]$stops
  
  if(a==1) next
  startdiff <- qm$level_unhc[a-1] - qm$unhc_mod[a]
  if(w > a & w < b) {
    enddiff <- UNHCmin - qm$unhc_mod[w]
    diffQ <- seq(startdiff, enddiff, length.out = w - a + 1)
    qm$unhc_mod[a:w] <- qm$unhc_mod[a:w] + diffQ
    a <- w + 1  
    startdiff <- enddiff
  }
  if(b==nrow(qm)){
    enddiff <- startdiff
  } else{
    enddiff <- qm$level_unhc[b+1] - qm$unhc_mod[b]
  }
  if(is.na(startdiff)||is.na(enddiff)) next
  
  diffQ <- seq(startdiff, enddiff, length.out = b - a + 1)
  
  qm$unhc_mod[a:b] <- qm$unhc_mod[a:b] + diffQ
}

# double check that everything looks okay
# par(mfrow = c(1,1))
# plot(qm$DateTime_UTC, qm$nhc_mod, log="y", main = "NHC", pch = 20)
# points(qm$DateTime_UTC[qm$notes == "nhc modeled"],
#        qm$nhc_mod[qm$notes=="nhc modeled"], pch = 20, col = "red")
# plot(qm$DateTime_UTC, qm$unhc_mod, log="y", main = "UNHC", pch = 20)
# points(qm$DateTime_UTC[qm$notes == "unhc modeled"],
#        qm$unhc_mod[qm$notes=="unhc modeled"], pch = 20, col = "red")

qm_l <- qm %>%
  select(DateTime_UTC, NHC = nhc_mod, UNHC = unhc_mod, notes) %>%
  pivot_longer(cols = -c(DateTime_UTC, notes), names_to = "site", values_to = "level_m")
lvl_long<- lvl_long %>%
  filter(!(site %in% c("NHC", "UNHC"))) %>%
  bind_rows(qm_l)
write_csv(lvl_long, "NHC_2019_metabolism/data/rating_curves/all_sites_level_corrected.csv")

nhc <- read_csv("NHC_2019_metabolism/data/metabolism/corrected_level/NHC_lvl.csv") %>%
  select(-level_m)
nhc <- qm %>%
  select(DateTime_UTC, level_m = nhc_mod) %>%
  right_join(nhc) %>%
  arrange(DateTime_UTC)
write_csv(nhc, "NHC_2019_metabolism/data/metabolism/corrected_level/NHC_lvl.csv")

unhc <- read_csv("NHC_2019_metabolism/data/metabolism/corrected_level/UNHC_lvl.csv") %>%
  select(-level_m)
unhc <- qm %>%
  select(DateTime_UTC, level_m = unhc_mod) %>%
  right_join(unhc) %>%
  arrange(DateTime_UTC)
write_csv(unhc, "NHC_2019_metabolism/data/metabolism/corrected_level/UNHC_lvl.csv")

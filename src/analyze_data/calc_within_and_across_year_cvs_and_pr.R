library(agricolae)


dat <- read_csv("NHC_2019_metabolism/data/metabolism/compiled/metabolism_summary_table_gC.csv") %>%
  filter(method == 'uninformed_raymond'| site %in% c("CBP", 'all') & is.na(year)) %>%
  mutate(year = ifelse(method == 'uninformed_raymond', year, 1969)) %>%
  select(-c(3:11))
# build summary table####
sum <- dat %>% filter(year != 2020, site != 'PWC') %>%
  mutate(across(starts_with('peak'), ~as.numeric(format(., '%j'))+4, 
                .names = '{col}_doy'))


# Get CV's of annual metabolism: ####
sum <- dat %>% filter(year != 1969, year != 2020, site != 'PWC') %>%
  mutate(across(starts_with('peak'), ~as.numeric(format(., '%j'))+4, 
                .names = '{col}_doy'))
sum %>% filter(year == 2019) %>% #summary()
  mutate(nep = er_mean - gpp_mean) %>% arrange(nep) %>% select(site, nep)
  summarize(across(ends_with('cum'), 
                   .fns = list(mean = ~mean(.),
                               sd = ~sd(.),
                               cv = ~sd(.)/mean(.)),
                   .names = '{col}_{fn}'))
sum %>% 
  filter(site %in% c('NHC','UNHC')) %>%
  group_by(site) %>%
  summarize(across(ends_with('cum'), 
                   .fns = list(mean = ~mean(.),
                               sd = ~sd(.),
                               cv = ~sd(.)/mean(.)),
                   .names = '{col}_{fn}')) 

yy <- sum %>% 
  filter(site %in% c('NHC','UNHC')) 
aov(er_mean~year, data=yy) %>%
  HSD.test("year", group=TRUE)
# Get CV's of within year metabolism: ####
dat <- readRDS("NHC_2019_metabolism/data/metabolism/compiled/met_preds_stream_metabolizer.rds")$preds %>%
  filter(era == 'now', site !='PWC', year != 2020)

cvs <- dat %>% 
  group_by(site, year) %>%
  summarize(across(c(GPP, ER), 
                   .fns = list(cv = ~sd(., na.rm = T)/mean(., na.rm = T) *12/32,
                               min = ~min(.,na.rm = T),
                               max = ~max(., na.rm = T),
                               mean = ~mean(., na.rm = T)),
                   .names = '{col}_{fn}'))
cvs %>% filter(year == 2019) %>% summary()

read_csv('NHC_2019_metabolism/data/metabolism/compiled/metabolism_summary_table_gC.csv') %>%
  select(site, year, gpp_cv, gpp_min, er_max, gpp_mean)
# get within year P:Rs ####
dat %>% 
  # filter(ER !=0) %>%
  mutate(PR = -GPP/ER) %>%
  # group_by(site, year) %>%
  summarize(PR = median(PR, na.rm = T))

cbp <- dat %>%
  filter(site == 'CBP' )%>%
  select(date, GPP, ER)%>%
  mutate(PR = -GPP/ER,
         PR = ifelse(is.infinite(PR), NA, PR)) 

plot(cbp$date, cbp$PR)

# compare magnitude of OCt resp to rest of year ####
dat %>%
  select(site, year, date, GPP, ER)%>%
  mutate(month = format(date, '%b'),
         across(c(GPP, ER), 
                .fns = list(all = ~ifelse(month %in% c('Sep', 'Oct','Nov'), NA, .),
                            oct = ~ifelse(month %in% c('Sep','Oct','Nov'), ., NA)),
                .names = '{col}_{fn}')) %>%
  group_by(site, year) %>%
  summarize(across(starts_with(c('GPP', 'ER')), 
                   .fns = list(mean = ~mean(., na.rm = T),
                               max = ~min(., na.rm = T)))) %>%
  mutate(mult = ER_oct_max/ER_mean) %>%
  select(site, year, mult) # mult is the max fall resp divided by mean annual resp

# find the end of the fall Resp Peak:

dat %>%
  mutate(month = format(date, '%b')) %>%
  filter(year == 2019, 
         month %in% c('Aug','Sep','Oct')) %>%
  select(site, date, ER, discharge) %>%
  mutate(discharge = ifelse(site == 'NHC', discharge, NA)) %>%
  group_by(site) %>%
  mutate(ER = na.approx(ER, na.rm = F, maxgap = 3))%>%
  ggplot(aes(date, ER, col = site)) +
  geom_line() + geom_smooth() 
  geom_line(aes(y = log(discharge)))
  pivot_wider(names_from = site, values_from = ER) %>%

######################################################
# pull out storms from NHC data (already flagged on portal, this doesn't find them for you)

library(tidyr)
library(dplyr)
library(lubridate)

setwd("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism")

filelist <- list.files("data/metabolism/processed")
allsites_storms <- data.frame()
for(i in 1:length(filelist)){
    # 1. get data for your river reach
    sitename <- substr(filelist[i], 1, nchar(filelist[i])-4)
    dat <- read.csv(paste0("data/metabolism/processed/",sitename,".csv"), 
                    header=TRUE, stringsAsFactors = FALSE)
    dat$DateTime_UTC <- ymd_hms(dat$DateTime_UTC)
    
    #remove duplicate rows
    dat <- dat[! duplicated(dat),]
    
    #get indices of all points that have been flagged with the word "storm".
    #i just flagged all likely storms for Eno 2017 and NHC 2017
    storm_rows = grep('storm', dat$flagcomment, ignore.case=TRUE)
    
    #get datetimes associated with those i just flagged (comment == 'storm')
    storm_dt = dat$DateTime_UTC[storm_rows]
    
    # create boolean column for storm or not storm
    dat$storm = as.numeric(dat$DateTime_UTC %in% storm_dt)
    
    #bring in metabolism estimates; deal with \\N characters from MySQL
    modres = readRDS(paste0("data/metabolism/condensed/condensed_",sitename,".rds"))
    modres <- modres$metab%>%
        select(-msgs.fit, -warnings, -errors)
    
    dat$waterpres_kPa <- (dat$WaterPres_kPa-dat$AirPres_kPa)
    dat <- dat %>% select(solar.time, DO_mgL, watertemp_C=WaterTemp_C, 
                          waterpres_kPa,storm)
    #average storm data by day and merge with metab estimates
    stormdaily = group_by(dat, 'date'=as.Date(substr(solar.time, 1, 10))) %>%
        summarize_if(is.numeric, mean, na.rm=TRUE) %>%
        left_join(modres, by='date') %>%
        as.data.frame()
    
    # Select days with >20% flagged as storm
    stormdaily$storm[stormdaily$storm>=.2]<-1
    stormdaily$storm[stormdaily$storm<.2]<-0
    
    # Only flag storm day with peak flow
    #get storm indices
    stormbool = as.logical(stormdaily$storm)
    
    stormrle = rle(stormbool)
    stormends = cumsum(stormrle$lengths)[c(FALSE, TRUE)]
    stormstarts = stormends - stormrle$lengths[c(FALSE, TRUE)] + 1
    
    stormchunks = mapply(`:`, stormstarts, stormends)
    stormmax = c()
    for(i in 1:length(stormchunks)){
        maxQ <- max(stormdaily$waterpres_kPa[stormchunks[[i]]], na.rm= TRUE)
        max_storm <- stormchunks[[i]][which(stormdaily$waterpres_kPa[stormchunks[[i]]]==maxQ)]
        stormmax = append(stormmax,max_storm)
    }
    stormdaily$storm[stormdaily$storm==1]<- .5
    stormdaily$storm[stormmax]<-1
    
    # find the % increase in flow for each storm
    Qincrease <- c()
    for(i in 1:length(stormmax)){
        stormQ <- stormdaily$discharge.m3s[stormmax[i]]
        preQ <- stormdaily$discharge.m3s[stormmax[i]-1]
        prestorm <- stormdaily$storm[stormmax[i]-1]
        if(prestorm == 0.5){
            preQ <- stormdaily$discharge.m3s[stormmax[i]-2] 
        }
        dQ <- stormQ/preQ
        Qincrease <- append(Qincrease, dQ)
    }
    
    # remove storms where stormflow isn't greater than pre stormflow
    stormyness <- data.frame(stormmax, Qincrease)
    stormyness <- stormyness[rev(order(Qincrease)),]
    # select the 10 greatest flow changes:
    stormdaily$storm<-0
    stormdaily$storm[stormyness$stormmax[1:10]]<-1
    stormdaily$site <- sitename
    
    allsites_storms <- rbind(allsites_storms, stormdaily)
}

#save csv
write.csv(allsites_storms, 'data/stormmet.csv', row.names = F)


#####################################################################################
# 2. some basic exploration ####
#read in daily storm data
allsites_storms<- read_csv('data/stormmet.csv')
sites <- c("UNHC","WBP", "WB", "CBP", "PM", "NHC")
# how many days do you want to look at?
predays <- 4
postdays <- 6

storms = which(allsites_storms$storm==1)
starts <- storms-predays
ends <- storms+postdays
allsites_storms$stormday <- NA
for(i in 1:length(starts)){
    n <-starts[i]
    allsites_storms$stormday[n:(n+predays+postdays)]<- seq(-2,5, by=1)
}
storms<- allsites_storms[!is.na(allsites_storms$stormday),]

#plot storm occurance relative to metab series
ylims = range(storms$GPP, storms$ER, na.rm=TRUE)

par(mfrow = c(2,3))
for(i in sites){
    subset <- storms[storms$site==i,]
    stormstarts <- which(subset$storm==1)-2
    plot( seq(-predays, postdays),subset$GPP[stormstarts[1]:(stormstarts[1]+predays+postdays)],ylim=c(-10,5), 
          main=i, ylab="gO2/m2/d", xlab="days since storm", type="n")
    abline(v=0, col="gray", lty=2, lwd=2)
    abline(h=0)
    for(j in 1:length(stormstarts)){
        n <- stormstarts[j]
        lines(seq(-predays,postdays), subset$GPP[n:(n+predays+postdays)], col = alpha("forestgreen",.5), lwd=2)
        lines(seq(-predays,postdays), subset$ER[n:(n+predays+postdays)], col = alpha("brown3",.5), lwd=2)
    }
}



plot(stormdaily$date, stormdaily$GPP, type='l', col='red', ylim=ylims,
    xlab='', ylab='g O2 m^-2 d^-1', main='CBP storms')
lines(stormdaily$date, stormdaily$ER, col='blue')
abline(v=stormdaily$date[stormdaily$storm > 0], col='gray', lty=2)


#get metab averages during storm days
duringGPP = tapply(stormdaily$GPP[stormbool], stormfac,
    mean, na.rm=TRUE)
duringER = tapply(stormdaily$ER[stormbool], stormfac,
    mean, na.rm=TRUE) * -1

#get pre- and post-storm index collections
during = which(stormbool)
day1pre = stormstarts - 1
day2pre = stormstarts - 2
day3pre = stormstarts - 3
day1post = stormends + 1
day2post = stormends + 2
day3post = stormends + 3

#get storm averages and pre - during, post - during differences
#pre - post would also be good to know
storm_summary = matrix(NA, 10, 7,
    dimnames=list(c('GPPmn', 'ERmn', 'GPP:ERmn', 'GPPse', 'ERse',
        'difGPPmn', 'difERmn', 'difGPP:ERmn', 'difGPPse', 'difERse'),
    c('day3pre', 'day2pre', 'day1pre', 'during', 'day1post', 'day2post', 'day3post')))
duringratio = mean(duringGPP, na.rm=TRUE) / mean(duringER, na.rm=TRUE)
for(i in 1:ncol(storm_summary)){
    day = colnames(storm_summary)[i]
    gppset = stormdaily$GPP[get(day)]
    erset = stormdaily$ER[get(day)] * -1

    storm_summary[1, i] = mean(gppset, na.rm=TRUE)
    storm_summary[2, i] = mean(erset, na.rm=TRUE)
    storm_summary[3, i] = storm_summary[1, i] / storm_summary[2, i]
    storm_summary[4, i] = sd(gppset, na.rm=TRUE) / sqrt(length(stormstarts))
    storm_summary[5, i] = sd(erset, na.rm=TRUE) / sqrt(length(stormstarts))

    if(day != 'during'){
        # gppdif = gppset - duringGPP
        # erdif = erset - duringER

        # storm_summary[6, i] = mean(gppdif, na.rm=TRUE)
        # storm_summary[7, i] = mean(erdif, na.rm=TRUE)
        storm_summary[6, i] = storm_summary[1,i] - mean(duringGPP, na.rm=TRUE)
        storm_summary[7, i] = storm_summary[2,i] - mean(duringER, na.rm=TRUE)
        storm_summary[8, i] = storm_summary[3, i] - duringratio
        storm_summary[9, i] = storm_summary[1,i] - sd(duringGPP, na.rm=TRUE) #wrong
        storm_summary[10, i] = storm_summary[2,i] - sd(duringER, na.rm=TRUE) #wrong
        # storm_summary[9, i] = sd(gppdif, na.rm=TRUE)
        # storm_summary[10, i] = sd(erdif, na.rm=TRUE)
    } else {
        storm_summary[6:10, i] = NA
    }
}

#plot
error.bars <- function(x, y, upper, lower=upper, cap.length=0.1, horiz=F,...){
    if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
        stop("One or more vectors is not the same length")

    if(horiz==F) {
        arrows(x,y+upper, x, y-lower, angle=90, code=3, length=cap.length, ...)
    } else if (horiz==T) {
        arrows(x+upper,y, x-lower, y, angle=90, code=3, length=cap.length, ...)
    }
}

# png(width=14, height=10, units='in', res=300,
#     filename='~/Desktop/NHC_stormMetabDif.png', type='cairo')

par(mfrow=c(2, 1), mar=c(3, 5, 3, 1))
gppymax = colSums(storm_summary[c(1,4),])
erymax = colSums(storm_summary[c(2,5),])
barplot(storm_summary[1:3,], beside=TRUE, names.arg=rep('', 7), main='Day mean',
    ylim=c(0, max(c(gppymax, erymax))), col=c('red3', 'blue2', 'purple3'),
    border=c('red3', 'blue2', 'purple3'), ylab='g O2 m^-2 d^-1',
    cex.names=1.5, cex.lab=1.5, cex.main=1.5)
legend('topright', legend=c('GPP', 'ER', 'GPP:ER'),
    fill=c('red3', 'blue2', 'purple3'), bty='n')
barx = seq(1.5, 26.5, 4)
error.bars(barx, storm_summary[1,], upper=storm_summary[4,],
    lower=storm_summary[4,], cap.length=0.01)
error.bars(barx + 1, storm_summary[2,], upper=storm_summary[5,],
    lower=storm_summary[5,], cap.length=0.01)
barplot(storm_summary[6:8,], beside=TRUE, main='Day mean minus storm mean',
    col=c('red3', 'blue2', 'purple3'), ylab='g O2 m^-2 d^-1',
    border=c('red3', 'blue2', 'purple3'), cex.names=1.5, cex.lab=1.5, cex.main=1.5)
#dev.off()

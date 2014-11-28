# This script takes the processed point count, guild, and land cover files and
# makes the data frames necessary for gdist analyses.
#--------------------------------------------------------------------------------*
# ---- SET-UP ----
#================================================================================*

# Load packages:

library(reshape2)
library(plyr)

# Source scripts:

source('analysis-scripts/RFuns.R')

# Get data:

pc = read.csv('derived-data/pc_10_14.csv')
g = read.csv('derived-data/guilds.csv')
lc = read.csv('derived-data/pts_lc100.csv')

#-------------------------------------------------------------------------*
# ---- Days since May 1 for each year ----
#-------------------------------------------------------------------------*

# 2006 is a typo (most likely 2009)! 

pc[pc$year == 2006,'year'] <- 2009

# Changing my decimal date to days of the year

month = as.numeric(substr(as.character(pc$ddate),0,1))

dday = pc$ddate - month

df = data.frame(year = pc$year, month, dday)

# Convert the decimal day to the day of the month:

df$day = ifelse(df$month == 6, df$dday*30,df$dday*31)

df$day = df$day + 1

# Create date object:

df$date = as.Date(as.character(paste(df$year, df$month, df$day, sep = '/')))

# Determine the day of year for a sampling event:

pc$ddate = sapply(df$date, dayOfYear)

# And change the column name:

names(pc)[3] <- 'day'

#-------------------------------------------------------------------------*
# ---- Dealing with bad records and removing exotic species ----
#=========================================================================*

# Remove Bert Drake (all zero counts?):

pc = pc[pc$site!='DRAKBERMD1',]
lc = lc[lc$site!='DRAKBERMD1',]

# Reduce site and count frames to just the sites with point counts:

lc = merge(lc, data.frame(site = unique(pc$site)), all = F)

pc = merge(pc, data.frame(site = unique(lc$site)), all = F)

# I can already write the lc csv to a file, first remove unnecessary 
# columns:

lc = lc[,-c(2:3)]

# NE2 doesn't have any count data ... remove:

lc = lc[lc$site!='NE2',]

write.csv(lc, 'gdist_data/site_covs.csv', row.names = F)

#-------------------------------------------------------------------------*
# ---- Remove exotic species ----
#-------------------------------------------------------------------------*

# Remove eust, ropi, and hosp:

pc = pc[pc$species!='eust' & pc$species!='hosp' & pc$species!='ropi',]

#-------------------------------------------------------------------------*
# ---- Remove fly-overs ----
#-------------------------------------------------------------------------*

# Coded in many ways!

pc = pc[pc$detection != 'flyover',]

pc = pc[pc$detection != 'flyover-v',]

pc = pc[pc$detection != 'v, flyover',]

pc = pc[pc$detection != 'flyover (a)',]

pc = pc[pc$detection != 'flyover (b)',]

# Remove detection type:

pc = pc[,-12]

#-------------------------------------------------------------------------*
# ---- Sky data ----
#-------------------------------------------------------------------------*

summary(pc$sky)

pc[is.na(pc$sky) & pc$site == 'CAMPLIBMD1',13]<-0

pc[is.na(pc$sky) & pc$site == 'ROHRSALMD1',13]<-2

#-------------------------------------------------------------------------*
# ---- Drop wind ----
#-------------------------------------------------------------------------*
# Note: No way to fix and maintain counts, which are more important.

pc = pc[,-14]

#-------------------------------------------------------------------------*
# ---- Fix temperature data ----
#-------------------------------------------------------------------------*
# Remove missing records:

pc = pc[pc$temp != '<NA>',]

# Find offending records:

table(pc$temp)

# First, n:

offending.recs = subset(pc,pc$temp == 'n')

good.recs = df.matcher(pc, offending.recs, 'no')

match.recs = good.recs[good.recs$day == unique(offending.recs$day) &
                 good.recs$year == unique(offending.recs$year),c(1,12)]

match.recs = na.omit(match.recs)

match.temp = mean(as.numericF(unique(match.recs)[,2]))

new.temps = rep(match.temp,dim(offending.recs)[1])

offending.recs$temp<-as.factor(new.temps)

pc = rbind(good.recs, offending.recs)

# -9999's ... Ali:

offending.recs = subset(pc,pc$temp == '-9999' & pc$observer == 'ajr')

good.recs = df.matcher(pc, offending.recs, 'no')

match.recs = good.recs[good.recs$day == unique(offending.recs$day) &
                         good.recs$year == unique(offending.recs$year),c(1,12)]

match.recs = na.omit(match.recs)

match.temp = mean(as.numericF(unique(match.recs)[,2]))

new.temps = rep(match.temp,dim(offending.recs)[1])

offending.recs$temp<-as.factor(new.temps)

pc = rbind(good.recs, offending.recs)

# -9999's ... Nora:

offending.recs = subset(pc,pc$temp == '-9999' & pc$observer == 'ned')

good.recs = df.matcher(pc, offending.recs, 'no')

match.recs = good.recs[good.recs$day == unique(offending.recs$day) &
                         good.recs$year == unique(offending.recs$year),c(1,12)]

match.recs = na.omit(match.recs)

match.temp = mean(as.numericF(unique(match.recs)[,2]))

new.temps = rep(match.temp,dim(offending.recs)[1])

offending.recs$temp<-as.factor(new.temps)

pc = rbind(good.recs, offending.recs)

# And then there's blanks:

offending.recs = subset(pc,pc$temp == '')

good.recs = df.matcher(pc, offending.recs, 'no')

match.recs = good.recs[good.recs$day == unique(offending.recs$day) &
                         good.recs$year == unique(offending.recs$year),c(1,12)]

match.recs = na.omit(match.recs)

match.temp = mean(as.numericF(unique(match.recs)[,2]))

new.temps = rep(match.temp,dim(offending.recs)[1])

offending.recs$temp<-as.factor(new.temps)

pc = rbind(good.recs, offending.recs)

# I've added NA's along the way:

pc = pc[!is.na(pc$site),]

# Convert temperature to numeric:

pc$temp = as.numericF(pc$temp)

#-------------------------------------------------------------------------*
# ---- Fix observer data ----
#-------------------------------------------------------------------------*

table(pc$observer)

pc[pc$observer == 'cw','observer']<-'ckw'

#-------------------------------------------------------------------------*
# ---- Multiple counts for some sites  ----
#-------------------------------------------------------------------------*

# A few still needed fixing ... multiple temp listings for CAMPLIBMD1:

pc[pc$site == 'CAMPLIBMD1' & pc$year == 2010,]

pc[pc$site == 'CAMPLIBMD1' & pc$year == 2010,'temp']<-66

# One Nora observation for MCCOCHMD1?

pc[pc$site == 'MCCOCHRMD1' & pc$year == 2009,]

row.names(pc) = 1:nrow(pc)

pc = pc[-3162,]

row.names(pc) = 1:nrow(pc)

# ODELVINMD2, temperature:

pc[pc$site == 'ODELVINMD2' & pc$year == 2009,'temp']<-68

# ROHRSALMD1, temperature:

pc[pc$site == 'ROHRSALMD1' & pc$year == 2010,]

pc[pc$site == 'ROHRSALMD1' & pc$year == 2010,'temp']<-71

#-------------------------------------------------------------------------*
# ---- Write observation-level covariates  ----
#-------------------------------------------------------------------------*

obs = pc[,c(1:4,12:13)]

obs = unique(obs)

obs$year = factor(obs$year)
obs$sky = factor(obs$sky)

write.csv(obs, 'gdist_data/obs_covs.csv', row.names = F)

#-------------------------------------------------------------------------*
# ---- Make and write observation frame  ----
#-------------------------------------------------------------------------*

y0 = pc[c(1:2,6:11)]

y0 =   ddply(y0, c('site','year','species'), numcolwise(sum))

# The remainder must be done on the fly ... depends on species and lumping:

write.csv(y0, 'gdist_data/y0.csv', row.names = F)



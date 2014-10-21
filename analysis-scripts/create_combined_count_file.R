# Libraries

library(plyr)
library(lubridate)

# Avian prep include distance and time:

setwd('C:/Users/Brian/Documents/community_analysis/raw_data')
list.files()

pc09 = read.csv('pc_2009.csv')
pc10 = read.csv('pc_2010.csv')
pc10s = read.csv('pc_supplement_2010.csv')
pc12 = read.csv('2012 Nestwatch Data.csv')  # Records 2012

setwd('C:/Users/Brian/Documents/community_analysis/derived_data')
lc = read.csv('pts_lc.csv')

#---------------------------------------------------------------------------------------
# Add impervious surface data to counts 
#---------------------------------------------------------------------------------------

# Fix the site  code naming conventions in the 2009 point count by merging with lc data:

pc09 = merge(pc09,lc,by = c('X','Y'), all = F)

# Basic merges to include imp data:

pc10 = merge(pc10, lc)
pc10s = merge(pc10s, lc)
pc12 = merge(pc12, lc, by.x = 'Site', by.y = 'site')

#---------------------------------------------------------------------------------------
# Simplify frames, removing unwanted columns
#---------------------------------------------------------------------------------------

# colnames = site, x, y, imp, date, observer, time.interval, species, d10,d20,d30,d40,d50,detect

pc09 = pc09[,c(23,1,2,24,7:16)]

pc10 = pc10[,c(22,1,2,23,4,5,7:14)]

pc10s = pc10s[,c(21,1,2,22,4,5,7:14)]

pc12 = pc12[,c(1,16,17,23,2,3,5:12)]

# Change column names for pc12:

colnames(pc12) = colnames(pc10s)

# Make a list of the point count data frames:

pc.list = list(pc09,pc10,pc10s,pc12)

# Bind the rows across data frames:

pc = ldply(pc.list, data.frame)

#---------------------------------------------------------------------------------------
# Change date formats and add year and day of year columns
#---------------------------------------------------------------------------------------

pc$Date = as.Date(pc$Date, format = '%m/%d/%Y')

pc$year = year(pc$Date)

# Some of the years were mis-labeled as 2006, making the change:

pc$year = ifelse(pc$year == 2006,2009,pc$year)

pc$jdate = yday(pc$Date)

#---------------------------------------------------------------------------------------
# Change the order of the columns and rename columns
#---------------------------------------------------------------------------------------

names(pc)

pc = pc[,c(1:5,15:16,6:14)]

colnames(pc) = c('site','x','y','imp','date','year','jdate','observer','time.int',
                  'spec','d10','d20','d30','d40','d50','obs.meth')

#---------------------------------------------------------------------------------------
# Fix case formats
#---------------------------------------------------------------------------------------

pc$observer = tolower(pc$observer)
pc$spec = tolower(pc$spec)
pc$obs.meth = tolower(pc$obs.meth)


#---------------------------------------------------------------------------------------
# Fix names (mispellings and misID)
#---------------------------------------------------------------------------------------

# Source the file with the fix names function:

source('C:/Users/Brian/Documents/community_analysis/r_scripts/fix_names_and_subset.R')

# Fix the names

pc = fix.sp.names(pc)

unique(pc[order(pc$spec),'spec'])


#---------------------------------------------------------------------------------------
# Change -9999 to NA
#---------------------------------------------------------------------------------------

make.na = function(x) ifelse(x == -9999, NA, x)

pc$d10 = make.na(pc$d10)
pc$d20 = make.na(pc$d20)
pc$d30 = make.na(pc$d30)
pc$d40 = make.na(pc$d40)
pc$d50 = make.na(pc$d50)

#---------------------------------------------------------------------------------------
# Write to file
#---------------------------------------------------------------------------------------

setwd('C:/Users/Brian/Documents/community_analysis/derived_data')
write.csv(pc, 'pc.csv', row.names = F)

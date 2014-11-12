#================================================================================*
# ---- PREPARE DATA FOR ANALYSIS ----
#================================================================================*
# This script takes the processed point count, guild, and land cover files and
# creates a wide frame of log-tranformed count data, linked with land cover data.

#--------------------------------------------------------------------------------*
# ---- SET-UP ----
#================================================================================*

# Load packages:

library(vegan)
library(reshape)
library(plyr)

# Get data:

setwd('~/Community_analysis')

pc = read.csv('derived-data/pc_10_14.csv')
g = read.csv('derived-data/guilds.csv')
lc = read.csv('derived-data/pts_lc100.csv')

#--------------------------------------------------------------------------------*
# ---- ENSURE LC AND PC FILES HAVE THE SAME SITES LISTED ----
#================================================================================*
# Note: Bert Drake removed due to no counts.

pc.site = data.frame(site = sort(unique(pc[pc$site!='DRAKBERMD1','site'])))
lc.site = data.frame(site = sort(unique(lc[lc$site!='DRAKBERMD1','site'])))

pc = merge(pc, lc.site, all = F)
lc = merge(lc, pc.site, all = F)

#--------------------------------------------------------------------------------*
# ---- PREPARE SPECIES POINT COUNT DATA FOR ANALYSIS ----
#================================================================================*

# For this run, I will use the mean count per location across distance classes:

pc$counts = rowSums(pc[,7:11])

# Prepare count frame

counts1 =  pc[,-c(3:5,7:15)]

t1 = melt(counts1, id = c('site','year','species'))

t2 = data.frame(cast(t1, site ~ species ~ variable, mean))

t3 = as.matrix(t2)

t3[is.na(t3)]<-0

t4 = data.frame(t3)

colnames(t4) = gsub('.counts','',colnames(t4))

counts = t4

#--------------------------------------------------------------------------------*
# ---- Log-transformed count data ----
#--------------------------------------------------------------------------------*

counts.lt = counts

for(i in 1:dim(counts)[2]){
  counts.lt[,i] = log(1+counts.lt[,i])
}

#--------------------------------------------------------------------------------*
# ---- Counts by guild ----
#================================================================================*

# Will use t1 (molten data frame) as the starting point, merging with
# trophic, nest, and foraging-trophic guilds. First, I will need to sum the
# counts across species within a guild for a given year, then I will need to
# take the average values across years.

# Create a new row of the guild frame htat combines foraging and trophic
# guilds:

  g$foraging_trophic = paste(g$foraging, g$trophic, sep ='_')

# Merge the frames (also removes species not in both lists):

  g$species = tolower(g$species)
  
  pcg = merge(t1,g, all = F)

# Function to create the count frame for a given guild:

  guild.shape = function(guild){
    t1a = aggregate(pcg$value, by = list(pcg$site,pcg$year,pcg[,guild]),sum)
      names(t1a) = c('site','year','guild','count')
    t1 = melt(t1a, id = c('site','year','guild'))
    t2 = data.frame(cast(t1, site ~ guild ~ variable, mean))
    t3 = as.matrix(t2)
    t3[is.na(t3)]<-0
    t4 = data.frame(t3)
    colnames(t4) = gsub('.count','',colnames(t4))
    return(t4)
    }

nest.counts = guild.shape('nest')
trophic.counts = guild.shape('trophic')
foraging_trophic.counts = guild.shape('foraging_trophic')

###############################################################################################################################################
# WORKS TO THIS POINT
###############################################################################################################################################

  guild.log.transform = function(life.history){
    # Get data
     pc2 = prep.guild(life.history)
    # Convert to log-transformed counts:
      for(i in 1:dim(pc2)[2]){
        pc2[,i] = log(1+pc2[,i])
      }
    return(pc2)
    }

# Return raw point count abundances by guild:

  pc.trophic = prep.guild(trophic)
  pc.nest = prep.guild(nest)

# Return log-transformed abundances by guild:

  trophic.logtransformed = guild.log.transform(trophic)
  nest.logtransformed = guild.log.transform(nest)

#--------------------------------------------------------------------------------*
# ---- PREPARE SPECIES POINT COUNT DATA FOR ABUNDANCE ANALYSIS ----
#================================================================================*

pc = read.csv('derived-data/pc_10_14.csv')
g = read.csv('derived-data/guilds.csv')
lc = read.csv('derived-data/pts_lc100.csv')

#--------------------------------------------------------------------------------*
# ---- ABUNDANCE ACROSS SPECIES ----
#================================================================================*

# Subset pc frame into site, year, species, and sum of detections across 
# distance classes:

pc2 = pc[,c(1,2,6)]

pc2$counts = rowSums(pc[7:11])

# Because multiple time intervals are present (and removal method is used), sum
# by year:

nyrs = ddply(pc2, .(site), summarize, c = length(unique(year)))
t.count = ddply(pc2, .(site), summarize, c = sum(counts))

pc3 = merge(nyrs,t.count, by = 'site')
pc3$t.a = pc3[,3]/pc3[,2]

pc3 = pc3[,-3]

# merge with lc frame:

pc.t.abund = merge(pc3, lc, all = F)

pc.t.abund = pc.t.abund[,-c(4,5)]

pc.abund = pc.t.abund[,c(1:2,4:5,3)]

names(pc.abund)[2] = 'nyrs'

#--------------------------------------------------------------------------------*
# ---- Relative abundance by guild ----
#================================================================================*

# Get data:

g$species = tolower(g$species)

# Merge count and guild data:

pcg = merge(pc2,g, all = F)

# Remove Bert Drake:

pcg = pcg[pcg$site!='DRAKBERMD1',]

# Combine feeding and trophic guilds:

pcg$trophic = paste(pcg$foraging, pcg$trophic, sep='-')

# Calculate the sum by trophic and nest guild

trophic.count = ddply(pcg, .(site, trophic), summarize, c = sum(counts))

nest.count = ddply(pcg, .(site, nest), summarize, c = sum(counts))

# Convert to wide format:

troph.melt = melt(trophic.count, id = c('site','trophic'))
nest.melt = melt(nest.count, id = c('site','nest'))

troph.cast = cast(troph.melt, site ~ trophic ~ variable,sum)
nest.cast = cast(nest.melt, site ~ nest ~ variable,sum)

troph.df = data.frame(troph.cast)
nest.df = data.frame(nest.cast)

troph.df$site = row.names(troph.df)
nest.df$site = row.names(nest.df)

# Merge with the count frame:

pc.abund2 = merge(pc.abund, troph.df)
pc.abund2 = merge(pc.abund2, nest.df)

pc.abund = pc.abund2

# Calculate the abundance scaled by years
# 
# for (i in 6:length(names(pc.abund))){
#   pc.abund[,i] = pc.abund[,i]/pc.abund[,2]
# }

# Make separate frames for trophic and nesting guilds:

troph.ab = pc.abund[,c(1,6:18)]

troph.ab$t = rowSums(troph.ab[,2:14])

nest.ab = pc.abund[,c(1,19:24)]

nest.ab$t = rowSums(nest.ab[,2:7])


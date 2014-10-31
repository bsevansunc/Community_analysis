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
# ---- PREPARE SPECIES POINT COUNT DATA FOR ANALYSIS ----
#================================================================================*

# For this run, I will use the max count per location across distance classes:

pc$counts = pc[,7]+pc[,8]+pc[,9]+pc[,10]+pc[,11]

# Prepare count frame

pc = pc[,-c(3:5,7:15)]

t = melt(pc, id = c('site','year','species'))

t2 = cast(t, site ~ species ~ variable,max)

t3 = data.frame(t2)

t3m = as.matrix(t3)

t3m[t3m == -Inf]<-0

t4 = data.frame(t3m)

colnames(t4) = gsub('.counts','',colnames(t4))

# One of the sites has never had a count, remove it:

t4$sums = rowSums(t4)

t4 = t4[t4$sums!=0,-56]

# Ensure all lc and pc files are covered (and vice-versa):

lc = unique(lc)

lc.sites = data.frame(lc[,1])
colnames(lc.sites) = 'site'

t4$site = row.names(t4)

t5 = merge(t4, lc.sites, all = F)

pc.sites = data.frame(t5$site)
names(pc.sites) = 'site'

pc = t5[,-1]

lc1 = merge(lc, pc.sites, all = F, incomparables = NA)

lc = lc1

#--------------------------------------------------------------------------------*
# ---- Log-transform count data ----
#--------------------------------------------------------------------------------*

pc2 = pc[,-53]

for(i in 1:dim(pc2)[2]){
  pc2[,i] = log(1+pc2[,i])
}

pc.logtransformed = pc2

pc = pc[,-53]

#--------------------------------------------------------------------------------*
# ---- Counts by guild ----
#================================================================================*

# Get data:

  pc1 = read.csv('derived-data/pc_10_14.csv')
  g$species = tolower(g$species)
  
  pcg = merge(pc1,g, all = F)

# For this run, I will use the max count per location across distance classes:

  pcg$counts = pcg[,7]+pcg[,8]+pcg[,9]+pcg[,10]+pcg[,11]

# Remove Bert Drake:

  pcg = pcg[pcg$site!='DRAKBERMD1',]

# Separate by life history trait

  pcg.troph = pcg
  pcg.nest = pcg

# Combine feeding and trophic guilds:

  pcg$trophic = paste(pcg.troph$foraging, pcg.troph$trophic, sep='-')

# Aggregate count data by guild:

  trophic = aggregate(pcg$counts,by = list(pcg$site,pcg$year,pcg$trophic),sum)
    names(trophic) = c('site','year','guild','count')
  nest = aggregate(pcg$counts,by = list(pcg$site,pcg$year,pcg$nest),sum)
    names(nest) = c('site','year','guild','count')

# Functions to prepare count frames for raw and log-transformed abundances:

  prep.guild = function(life.history){
    t = melt(life.history, id = c('site','year','guild'))
    t2 = cast(t, site ~ guild,max)
    t3 = data.frame(t2)[,-1]
    t3m = as.matrix(t3)
    t3m[t3m == -Inf]<-0
    t4 = data.frame(t3m)
    t4$site = t2$site
    lc = unique(lc)
    lc.sites = data.frame(lc[,1])
    colnames(lc.sites) = 'site'
    t5 = merge(t4, lc.sites, all = F)
    pc.sites = data.frame(t2$site)
    names(pc.sites) = 'site'
    pc = t5[,-1]
    return(pc)
  }

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


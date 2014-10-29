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

head(pc)

pc2 = pc[,-53]

for(i in 1:dim(pc2)[2]){
  pc2[,i] = log(1+pc2[,i])
}

pc = pc2

#--------------------------------------------------------------------------------*
# ---- Counts by guild ----
#================================================================================*

pc1 = read.csv('derived-data/pc_10_14.csv')
g$species = tolower(g$species)

pcg = merge(pc1,g, all = F)


#### Prepare count frame ----

# For this run, I will use the max count per location across distance classes:

pcg$counts = pcg[,7]+pcg[,8]+pcg[,9]+pcg[,10]+pcg[,11]

# Remove Bert Drake:

pcg = pcg[pcg$site!='DRAKBERMD1',]

# Separate by life history trait

pcg.troph = pcg
pcg.nest = pcg
pcg.migr = pcg

pcg$trophic = paste(pcg.troph$foraging, pcg.troph$trophic, sep='-')


trophic = aggregate(pcg$counts,by = list(pcg$site,pcg$year,pcg$trophic),sum)
names(trophic) = c('site','year','guild','count')
nest = aggregate(pcg$counts,by = list(pcg$site,pcg$year,pcg$nest),sum)
names(nest) = c('site','year','guild','count')
mig = aggregate(pcg$counts,by = list(pcg$site,pcg$year,pcg$migratory),sum)
names(mig) = c('site','year','guild','count')


# Function to prepare count frames

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
  t5 = t5[,-1]
  # Convert to log-transformed counts:
  pc2 = t5  
  for(i in 1:dim(pc2)[2]){
    pc2[,i] = log(1+pc2[,i])
  }
  pc2
}

trophic = prep.guild(trophic)
nest = prep.guild(nest)
mig = prep.guild(mig)
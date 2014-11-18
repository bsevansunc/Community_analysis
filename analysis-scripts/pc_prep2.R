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

counts.sp = t4

#--------------------------------------------------------------------------------*
# ---- Total abundance (mean counts across species) ----
#--------------------------------------------------------------------------------*

counts.total = data.frame(count = rowSums(counts.sp))

#--------------------------------------------------------------------------------*
# ---- Log-transformed count data ----
#--------------------------------------------------------------------------------*

counts.lt = counts.total

for(i in 1:dim(counts.lt)[2]){
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

#--------------------------------------------------------------------------------*
# Log-transformed count data
#--------------------------------------------------------------------------------*

  guild.log.transform = function(life.history){
    # Get data
     pc2 = life.history
    # Convert to log-transformed counts:
      for(i in 1:dim(pc2)[2]){
        pc2[,i] = log(1+pc2[,i])
      }
    return(pc2)
    }

nest.counts.lt = guild.log.transform(nest.counts)
trophic.counts.lt = guild.log.transform(trophic.counts)
foraging_trophic.counts.lt = guild.log.transform(foraging_trophic.counts)

#--------------------------------------------------------------------------------*
# ---- Relative abundance by guild ----
#--------------------------------------------------------------------------------*

guild.rel.abund = function(guild.df){
  rel.g.abund = list()
    for (i in 1:ncol(guild.df)){
      rel.g.abund[[i]] = guild.df[i]/(rowSums(guild.df))
    }
    df = do.call('cbind', rel.g.abund)
    return(df)
}

nest.ra = guild.rel.abund(nest.counts)
trophic.ra = guild.rel.abund(trophic.counts)
foraging_trophic.ra = guild.rel.abund(foraging_trophic.counts)

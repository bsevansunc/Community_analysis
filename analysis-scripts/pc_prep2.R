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
g = read.csv('guild_data/niche.csv')
nw = read.csv('guild_data/niche_width.csv')
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

# Remove species not part of this analysis:

pc = pc[pc$species!='rsha'&pc$species!='hosp'&pc$species!='eust'&pc$species!='ropi'&pc$species!='chsw'&pc$species!='eaki',]

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

# Remerge with lc

pc.names = data.frame(site = row.names(counts.sp))

lc2 = merge(lc, pc.names)

lc = lc2[,c(1,5)]

#--------------------------------------------------------------------------------*
# ---- Total abundance (mean counts across species) ----
#--------------------------------------------------------------------------------*

counts.total = data.frame(count = rowSums(counts.sp))

#########################################################################################

# niche width by site and niche width dimension

nwOutFun = function(niche.dimension, main.name){
j.out = numeric()
t.count = numeric()
count.m = matrix(nrow = length(counts.sp[,1]),ncol = length(nw$sp))
scale.count.m = matrix(nrow = length(counts.sp[,1]),ncol = length(nw$sp))

for (j in 1:length(counts.sp[,1])){
  for(i in 1:length(nw$sp)){
    count.m[j,i] = counts.sp[j,i]
    scale.count.m[j,i] = count.m[j,i]*nw[i,niche.dimension]
}}
nw2 = rowSums(scale.count.m)/rowSums(count.m)
nw_lc = data.frame(site = lc$site, imp = lc$imp, nw = nw2)

mod1 = lm(nw~imp, data = nw_lc)
mod2 = lm(nw~imp+I(imp^2), data = nw_lc)

lc_imp = data.frame(imp = seq(min(lc$imp),max(lc$imp), by = .01))
mod.preds = predict(mod2, newdata = lc_imp, se.fit = T, interval = 'confidence')$fit
mod.preds = data.frame(mod.preds)

p = plot(nw~imp, data = nw_lc,
     pch = 19, col = 'gray70',
     #type = 'l', lwd = 2,
     ylim = c(min(nw_lc$nw),max(nw_lc$nw)),
     main = main.name,
     xlab = '% Impervious',
     ylab = 'Niche Width (Rank percentile)',
     bty = 'l')
p1 = p + lines(fit~lc_imp$imp, data = mod.preds, lwd = 2) +
  lines(lwr~lc_imp$imp, data = mod.preds, lwd = 1.5, lty = 2) +
  lines(upr~lc_imp$imp, data = mod.preds, lwd = 1.5, lty = 2)
return(list(summary(mod1), summary(mod2),p1))
}

pdf('guild_data/guild_out/niche_width_scatterplots.pdf', 
    width = 6.5, height = 5.5, onefile = T)

nwOutFun('nw', 'Combined Eltonian and Grinnelian niche width')

nwOutFun('enw', 'Eltonian niche width')

nwOutFun('dietEven', 'Diet breadth')

nwOutFun('fsEven', 'Foraging strata breadth')

nwOutFun('varHeight', 'Nest height variability')

nwOutFun('offspring_max', 'Maximum number of offspring')

nwOutFun('gnw', 'Grinnelian niche width')

nwOutFun('Brange', 'Breeding range')

nwOutFun('LNB', 'Local niche breadth')

nwOutFun('RNB', 'Regional niche breadth')

nwOutFun('BiomeH', 'Biome diversity')

cor(nw$BiomeH, nw$Brange)
cor(nw$BiomeH, nw$RNB)
cor(nw$BiomeH, nw$LNB)

nw$new_gnw = pRank(nw$BiomeH+nw$RNB+nw$LNB)
nw$new_enw = pRank(nw$dietEven+nw$varHeight+nw$offspring_max)
nw$new_nw = pRank(nw$new_gnw+nw$new_enw)

nwOutFun('new_gnw','Grinnelian niche width 2')
nwOutFun('new_enw','Eltonian niche width 2')
nwOutFun('new_nw','Combined niche width 2')

nw$new_nw2 = pRank(nw$new_enw + nw$gnw)

nwOutFun('new_nw2', 'Combined niche width 3')

dev.off()

#########################################################################################

# niche by site and niche width dimension
# niche width by site and niche width dimension

nwOutFun = function(niche.dimension, main.name, ylab){
  j.out = numeric()
  t.count = numeric()
  count.m = matrix(nrow = length(counts.sp[,1]),ncol = length(nw$sp))
  scale.count.m = matrix(nrow = length(counts.sp[,1]),ncol = length(nw$sp))
  
  for (j in 1:length(counts.sp[,1])){
    for(i in 1:length(nw$sp)){
      count.m[j,i] = counts.sp[j,i]
      scale.count.m[j,i] = count.m[j,i]*nw[i,niche.dimension]
    }}
  nw2 = rowSums(scale.count.m)/rowSums(count.m)
  nw_lc = data.frame(site = lc$site, imp = lc$imp, nw = nw2)
  
  mod1 = lm(nw~imp, data = nw_lc)
  mod2 = lm(nw~imp+I(imp^2), data = nw_lc)
  
  lc_imp = data.frame(imp = seq(min(lc$imp),max(lc$imp), by = .01))
  mod.preds = predict(mod2, newdata = lc_imp, se.fit = T, interval = 'confidence')$fit
  mod.preds = data.frame(mod.preds)
  
  p = plot(nw~imp, data = nw_lc,
           pch = 19, col = 'gray70',
           #type = 'l', lwd = 2,
           ylim = c(min(nw_lc$nw),max(nw_lc$nw)),
           main = main.name,
           xlab = '% Impervious',
           ylab = ylab,
           bty = 'l')
  p1 = p + lines(fit~lc_imp$imp, data = mod.preds, lwd = 2) +
    lines(lwr~lc_imp$imp, data = mod.preds, lwd = 1.5, lty = 2) +
    lines(upr~lc_imp$imp, data = mod.preds, lwd = 1.5, lty = 2)
    abline(h = 0, col ='red', lty = 3, lwd = 3)
  return(list(summary(mod1), summary(mod2),p1))
}

pdf('guild_data/guild_out/population_trajectory.pdf', 
    width = 6.5, height = 5.5, onefile = T)

nwOutFun('trend', 'Population trend (1966-2012) by site', 'Population trajectory scaled by abundance')
dev.off()

pdf('guild_data/guild_out/diet_plots.pdf', 
    width = 6.5, height = 5.5, onefile = T)

nw$inv_per = t1$Diet.Inv

nwOutFun('inv_per','Proportion of invertebrates by impervious','% Invertebrates')

nw$fruit_per = t1$Diet.Fruit

nwOutFun('fruit_per','Proportion of fruit by impervious','% Fruit')

nw$seed_per = t1$Diet.Seed

nwOutFun('seed_per','Proportion of seed by impervious','% Seed')




nwOutFun('enw', 'Eltonian niche width')

nwOutFun('dietEven', 'Diet breadth')

nwOutFun('fsEven', 'Foraging strata breadth')

nwOutFun('varHeight', 'Nest height variability')

nwOutFun('offspring_max', 'Maximum number of offspring')

nwOutFun('gnw', 'Grinnelian niche width')

nwOutFun('Brange', 'Breeding range')

nwOutFun('LNB', 'Local niche breadth')

nwOutFun('RNB', 'Regional niche breadth')

nwOutFun('BiomeH', 'Biome diversity')

cor(nw$BiomeH, nw$Brange)
cor(nw$BiomeH, nw$RNB)
cor(nw$BiomeH, nw$LNB)

nw$new_gnw = pRank(nw$BiomeH+nw$RNB+nw$LNB)
nw$new_enw = pRank(nw$dietEven+nw$varHeight+nw$offspring_max)
nw$new_nw = pRank(nw$new_gnw+nw$new_enw)

nwOutFun('new_gnw','Grinnelian niche width 2')
nwOutFun('new_enw','Eltonian niche width 2')
nwOutFun('new_nw','Combined niche width 2')

nw$new_nw2 = pRank(nw$new_enw + nw$gnw)

nwOutFun('new_nw2', 'Combined niche width 3')

dev.off()
#########################################################################################

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

  pcg$nft = paste(pcg$nest, pcg$foraging_trophic)

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
nft.counts = guild.shape('nft')

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
nft.lt = guild.log.transform(nft.counts)


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

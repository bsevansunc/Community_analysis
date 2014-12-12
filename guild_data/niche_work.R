#---------------------------------------------------------------------------------*
# ---- Combined niche width files ----
#=================================================================================*

# Get data:

gnw = read.csv('guild_data/grinnelian_niche_width.csv')
enw = read.csv('guild_data/eltonian_niche_width.csv')

# Merge files:

mg1 = merge(enw, gnw, by.x  = 'sp', by.y = 'alpha')

mg1$nw = pRank(mg1$enw + mg1$gnw)

mg1[order(mg1$nw),]

# Write to file:

write.csv(mg1, 'guild_data/niche_width.csv', row.names = F)

#---------------------------------------------------------------------------------*
# ---- Combined guild files ----
#=================================================================================*

gg = read.csv('guild_data/guild_grinnelian_raw.csv')

ge = read.csv('guild_data/guild_eltonian.csv')

# Dealing with naming problems!

name.test = data.frame(common = gg[,3],dgg = rep(1, length(gg[,1])))

ge_names = data.frame(ge[,1:3],dge = rep(1, length(ge[,1])))

test = merge(ge_names, name.test, all.x = T, all.y = F)

gg[,3]<-as.character(gg[,3])

gg[gg$AOU == 7310,3]<-'Eastern Tufted Titmouse'

gg[gg$AOU == 4050,3]<-'Piliated Woodpecker'

gg[gg$AOU == 4120,3]<-'Yellow-shafted Flicker'

# Now merge gg and ge_names, only keeping matching records

gg1 = merge(ge_names, gg, by.x = 'common', by.y = 'CommonName', all = F)

mg2 = mg1[,c(2,14,15,17,20,25,28)]

gg2 = gg1[,c(2,12,14:15,22,32:33)]

# Merge gg and ge

niche1 = merge(gg2, ge)

niche2 = niche1[,c(1:7,13,14:16,41)]

# Change cup-cavity (carw) to cavity:

niche2[niche2$nest_type == 'Cup-cavity','nest_type']<-'cavity'

niche2[niche2$nest_loc == 'shrub-tree','nest_loc']<-'tree-shrub'

# Make combined nest category

niche2$nest = factor(paste(niche2$nest_type, niche2$nest_loc, sep = '-'))

niche2[niche2$nest == 'cup-tree-shrub','nest']<-'cup-generalist'

niche2[niche2$nest == 'cup-ground-shrub','nest']<-'cup-ground'

niche2$nest = as.character(niche2$nest)

niche2[niche2$nest == 'cavity-shrub-building'|niche2$nest == 'cup-building','nest']<-'building'

niche2$nest = factor(niche2$nest)

# Simplify foraging category

niche2$Foraging = as.character(niche2$Foraging)

niche2[niche2$Foraging == 'aerial foraging'|niche2$Foraging == 'hawks','Foraging']<-'aerial'

niche2[niche2$Foraging == 'hover/glean','Foraging']<-'foliage glean'

niche2[niche2$Foraging == 'foliage glean','Foraging']<-'foliage'

niche2[niche2$Foraging == 'ground glean','Foraging']<-'ground'

niche2[niche2$Foraging == 'bark glean','Foraging']<-'bark'

niche2$Foraging = factor(niche2$Foraging)

table(niche2$Foraging)
summary(niche2)

write.csv(niche2, 'niche.csv', row.names = F)

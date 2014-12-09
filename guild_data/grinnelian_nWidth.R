gg = read.csv('guild_grinnelian_raw.csv')

ge = read.csv('guild_eltonian.csv')

# Dealing with naming problems!

name.test = data.frame(common = gg[,3],dgg = rep(1, length(gg[,1])))

ge_names = data.frame(ge[,1:3],dge = rep(1, length(ge[,1])))

test = merge(ge_names, name.test, all.x = T, all.y = F)

gg[,3]<-as.character(gg[,3])

gg[gg$AOU == 7310,3]<-'Eastern Tufted Titmouse'

gg[gg$AOU == 4050,3]<-'Piliated Woodpecker'

gg[gg$AOU == 4120,3]<-'Yellow-shafted Flicker'

# Now merge gg and ge, only keeping matching records

mg1 = merge(ge_names, gg, by.x = 'common', by.y = 'CommonName', all = F)

mg2 = mg1[,c(2,14,15,17,20,25,28)]

# Function for percentile rank

pRank <- function(x) trunc(rank(x))/length(x)

mg2$Brange = 1-pRank(mg2$Brange_Area_km2)

mg2$LNB = 1-pRank(mg2$new_Im)

mg2$RNB = 1-pRank(mg2$Tol)

mg2$BiomeH = 1-pRank(mg2$BiomeH_noMex)

# Remove unnecessary columns:

mg2 = mg2[,-c(3:7)]

# Change the 1966-2004 trend record:

names(mg2)[2] <- 'trend'

# Grinnelian Niche width:

mg2$gnw = rowSums(mg2[,3:6])

mg2$gnw = pRank(mg2$gnw)

# Check out the results!

mg2[order(mg2$gnw),]

# Take away the trend data:

mg2 = mg2[,-2]

# Write the file:

write.csv(mg2, 'guild_grinnelian.csv', row.names = F)

# Get eltonian niche width file, merge them:

enw = read.csv('eltonian_niche_width.csv')

mg3 = merge(enw, mg2, by.x  = 'sp', by.y = 'alpha')

mg3$nw = pRank(mg3$enw + mg3$gnw)

mg3[order(mg3$nw),]

# Write to file:

write.csv(mg3, 'niche_width.csv', row.names = F)

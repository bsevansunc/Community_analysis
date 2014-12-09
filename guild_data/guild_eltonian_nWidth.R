library(vegan)

ge = read.csv('guild_eltonian.csv')

# Add min and max nest height data for the Carolina Wren (report from All About Birds):

ge[ge$alpha == 'carw',4:5]<-c(0,2)

# Remove the red-shouldered hawk:

ge = ge[ge$alpha!='rsha',]

# Diet diversity

diet = data.frame(row.names = ge[,2], ge[,20:28])

diet.div = diversity(diet)

diet.even = diet.div/log(length(diet))

diet$div = diet.div

diet$even = diet.even

# Foraging strata diversity:

fs = data.frame(row.names = ge[,2], ge[,30:34])

fs.div = diversity(fs)

fs.even = fs.div/log(length(fs))

fs$div = fs.div

fs$even = fs.even

# Nest height variation (alas, diversity or evenness not possible with the available data):

varHeight = ge$nest_max - ge$nest_min

# Variation in incubation time:

inc_var = ge$Inc_period_max - ge$Inc_period_min

# Variation in nestling time:

nestling_var = ge$Nestling_period_max - ge$Nestling_period_min

# Variation in clutch size:

clutch_var = ge$Clutch_max - ge$Clutch_min

# Variation in broods per year:

brood_var = ge$Brood_max

# Summarizing Eltonian niche width:

enw.df = data.frame(sp = ge$alpha, dietEven = diet.even, 
                    fsEven = fs.even, varHeight, 
                    inc_var, nestling_var, clutch_var, brood_var)

# Because bhco is a nest parasite (with up to 40 broods per year, I will give an nHeight and breeding var of 0):

enw.df[enw.df$sp == 'bhco',c(4,8)]<- c(max(na.omit(enw.df$varHeight)),max(na.omit(enw.df$brood_var)))

# Dobkon et al. 1995 says that House Wren exhibit 11.5 m range in nest height:

enw.df[enw.df$sp == 'howr',4]<-11.5

# Convert to percentiles:

pRank <- function(x) trunc(rank(x))/length(x)

enwP.df = enw.df
for (i in 2:length(enw.df)){
  enwP.df[,i] = 1-pRank(enw.df[,i])
}

enwP.df$breeding_var = pRank(rowSums(enwP.df[,5:8]))

# Remove specific breeding info:

enwP.df = enwP.df[,-c(5:8)]

# Come up with a single index of niche width:

enwP.df$enw1 = pRank(rowSums(enwP.df[,2:5]))

enwP.df[order(enwP.df$enw1),]

# Change to "1 - ..." so 0 represents the most generalist:

for(i in 2:length(enwP.df)){
  enwP.df[,i] = 1 - enwP.df[,i]
}

# Come up with a single index of niche width:

enw.df$enw1 = rowSums(enw.df[,2:5])

# Take a look!

enw.df[order(enw.df$enw1),]

# Scaled diet and foraging strata evenness (so 0's are generalist, recorded as: 1 -)

diet.even = 1-pRank(diet.even)
fs.even = 1-pRank(fs.even)

nHeight_var = 1- pRank(varHeight)



# Write file:

write.csv(enw.df, 'eltonian_niche_width.csv', row.names = F)


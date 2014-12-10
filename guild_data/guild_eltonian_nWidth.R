library(vegan)

ge = read.csv('guild_eltonian.csv')

# Take away the red-shouldered hawk:

ge = ge[ge$alpha!='rsha',]

# Add min and max nest height data for the Carolina Wren (report from All About Birds):

ge[ge$alpha == 'carw',4:5]<-c(0,2)

# Remove the red-shouldered hawk:

ge = ge[ge$alpha!='rsha',]

# Diet diversity

diet = data.frame(row.names = ge[,2], ge[,20:28])

dietEven = diversity(diet)/log(length(diet))

# Foraging strata diversity:

fs = data.frame(row.names = ge[,2], ge[,30:34])

fsEven = diversity(fs)/log(length(fs))

# Nest height variation (alas, diversity or evenness not possible with the available data):

varHeight = ge$nest_max - ge$nest_min

# Minimum incubation time:

inc_min = ge$Inc_period_min # Lower values represent generalist

# Minimum nestling time:

nestling_min = ge$Nestling_period_min # Lower values should represent generalist

# Maximum clutch size:

clutch_max = ge$Clutch_max # Higher values represent generalist

# Maximum number of broods per year:

brood_max = ge$Brood_max # Higher values represent generalist 

# Maximum number of offspring per year

offspring_max = clutch_max*brood_max # Higher values represent generalist

# Summarizing Eltonian niche width:

enw.df = data.frame(sp = ge$alpha, dietEven,fsEven, varHeight, offspring_max)

# Because bhco is a nest parasite (with up to 40 broods per year, I will give an nHeight and breeding var of 0):

enw.df[enw.df$sp == 'bhco',c(4,5)]<- c(max(na.omit(enw.df$varHeight)),40)

# Dobkon et al. 1995 says that House Wren exhibit 11.5 m range in nest height:

enw.df[enw.df$sp == 'howr',4]<-11.5

# Function to convert to percentiles:

pRank <- function(x) trunc(rank(x))/length(x)

enwP.df = enw.df
for (i in 2:length(enw.df)){
  enwP.df[,i] = 1-pRank(enw.df[,i])
}

# Come up with a single index of Eltonian niche width:

enwP.df$enw = pRank(rowSums(enwP.df[,2:5]))

# Take a look!

enwP.df[order(enwP.df$enw),]

# Write file:

write.csv(enwP.df, 'eltonian_niche_width.csv', row.names = F)


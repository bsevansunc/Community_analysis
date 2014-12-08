library(vegan)

ge = read.csv('guild_eltonian.csv')

scale01Fun = function(x){
  x2 = na.omit(x)
  (x-min(x2))/(max(x2)-min(x2))
}

# Add min and max nest height data for the Carolina (report from All About Birds):

ge[ge$alpha == 'carw',4:5]<-c(0,2)

# Diet diversity

diet = data.frame(row.names = ge[,2], ge[,20:28] )

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

# Scaled diet and foraging strata evenness

diet.even = scale01Fun(diet.even)
fs.even = scale01Fun(fs.even)

# Nest height variation (alas, diversity or evenness not possible with the available data):

varHeight = ge$nest_max - ge$nest_min

nHeight_var = varHeight/max(na.omit(varHeight))

nHeight_var = scale01Fun(varHeight)

# Variation in incubation time:

inc_var = ge$Inc_period_max - ge$Inc_period_min

inc_var = scale01Fun(inc_var)

# Variation in nestling time:

nestling_var = ge$Nestling_period_max - ge$Nestling_period_min

nestling_var = scale01Fun(nestling_var)

# Variation in clutch size:

clutch_var = ge$Clutch_max - ge$Clutch_min

clutch_var = scale01Fun(clutch_var)

# Variation in broods per year:

brood_var = ge$Brood_max - ge$Brood_min

brood_var = scale01Fun(brood_var)

# Summarizing nesting:

breeding_var = numeric()
for(i in 1:length(ge[,1])){
 breeding_var[i] = sum(inc_var[i], nestling_var[i],
                       clutch_var[i], brood_var[i]) 
}

breeding_var = scale01Fun(breeding_var)

# Summarizing Eltonian niche width:

enw.df = data.frame(sp = ge$alpha, diet.even, fs.even, nHeight_var, 
           breeding_var)

# Because bhco is a nest parasite (with up to 40 broods per year, I will give an nHeight and breeding var of 1):

enw.df[enw.df$sp == 'bhco',4:5]<- c(1,1)

# Come up with a single index of niche width:

enw = rowMeans(enw.df[,2:5], na.rm = T)

enw = scale01Fun(enw)

enw.df$enw = enw

# Change enw so that 0 represents generalists:

enw.df$enw = 1-enw

# Take a look!

enw.df[order(enw.df$enw),]


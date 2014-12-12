# Prep guild file

setwd('~/Community_analysis/guild_data')

list.files()

diet = read.table('diet_foraging', header=T, sep = '\t', quote = "\"", fill=T)

nest = read.csv('nest_territory.csv')

# Missing info on the Northern Flicker!

ysfl = nest[nest$alpha == 'ysfl',]

ysfl[,c(15:22)]<-c(11,13,24,27,5,8,1,1)

nest = nest[-45,]

nest = rbind(nest, ysfl)

# Subset frames to necessary columns:

names(diet)

d = diet[,c(8:20,24:30,36)]

names(nest)

nest = nest[,c(1:6,10:11,13:22)]

# Try merge, fix merge:

g = merge(nest, d, by.x = 'common',by.y = 'English', all = F)

g1 = data.frame(sp = g$alpha, count = rep(1, length(g$alpha)))
n1 = data.frame(sp = nest$alpha, c1 = rep(1, length(nest$alpha)))

merge(g1, n1, all = T)

d$English = as.character(d$English)

nest[nest$alpha == 'bggn','common']
d[d$English == 'Blue-grey Gnatcatcher','English']<-'Blue-gray Gnatcatcher'

nest[nest$alpha == 'eawp','common']
nest[nest$alpha == 'eawp','latin']
d[d$Scientific == 'Contopus virens','English']<-'Eastern Wood-Pewee'

d[d$English == 'Grey Catbird','English']<-'Gray Catbird'

nest[nest$alpha == 'piwo','common']
nest[nest$alpha == 'piwo','latin']
d[d$Scientific == 'Dryocopus pileatus','English']<-'Piliated Woodpecker' # Yes, I know this is not the right spelling!


nest[nest$alpha == 'tuti','common']
nest[nest$alpha == 'tuti','latin']
d[d$Scientific == 'Baeolophus bicolor','English']<-'Eastern Tufted Titmouse' 

nest[nest$alpha == 'ysfl','common']
nest[nest$alpha == 'ysfl','latin']
d[d$Scientific == 'Colaptes auratus','English']<-'Yellow-shafted Flicker' 

d$English = factor(d$English)

# Merge (data checked, names are now equivalent):

g = merge(nest, d, by.x = 'common',by.y = 'English', all = F)

# Change Nest location categories:
g$nest_loc = as.character(g$nest_loc)
g$nest_loc[1] <- 'Tree'


gsub(g$nest_loc[i], 'Tree-branch', 'Tree')

for(i in 1:length(g[,1])){
  if (g$nest_loc[i] == 'Tree-branch') g$nest_loc[i] <-  'Tree'
  if (g$nest_loc[i] == 'Ground-tree-branch-shrub') g$nest_loc[i] <-  'Generalist'
  if (g$nest_loc[i] == 'Tree-trunk') g$nest_loc[i] <-  'Tree'
  if (g$nest_loc[i] == 'Shrub-tree-branch') g$nest_loc[i] <- 'Shrub-tree'
}

g$nest_loc = tolower(g$nest_loc)
g$nest_loc = factor(g$nest_loc)

# There's one blank nest type:

g[g$nest_type =='','nest_type']<-'cup'

# Some of the diet and foraging columns have no values, remove fields:

summary(g)
names(g)

g = g[,-c(24,31,32)]

# Reset the factor levels for 5-category dietary guilds:

g$Diet.5Cat = factor(g$Diet.5Cat)

# Missing info on the Eastern Phoebe!

eaph = data.frame(common = 'Eastern Phoebe', alpha = 'eaph', latin = 'Sayornis phoebe',
                  nest_min = 1, nest_max = 5, nest_mean = 3, nest_type = 'cup', 
                  nest_loc = 'building', summer_flock = 'N',
                  broods = 'double', 15, 16, 16, 20, 1, 4, 1, 2, 'Sayornis phoebe',
                  90, 0,0,0,0,10,0,0,0,'Invertebrate',
                  0,50,50,0,0,19.7)

names(eaph) = names(g)

g = rbind(g, eaph)

# Missing info on the Red-bellied Woodpecker!

rbwo = data.frame(common = 'Red-bellied Woodpecker', alpha = 'rbwo', latin = 'Melanerpes carolinus',
                  nest_min = 2, nest_max = 18, nest_mean = 7.6, nest_type = 'cavity', 
                  nest_loc = 'tree', summer_flock = 'N',
                  broods = 'double', 12, 12, 24, 27, 2, 6, 1, 3, 'Melanerpes carolinus',
                  30, 10,10,0,0,0,20,10,20,'Omnivore',
                  0,0,20,60,20,69.5)

names(rbwo) = names(g)

g = rbind(g, rbwo)

g = g[order(as.character(g$alpha)),]

# Write cleaned guild data to file:

write.csv(g, 'guild_eltonian.csv', row.names = F)

#################################################################

nest_timeMin = g$Inc_period_min +  g$Nestling_period_min

nest_timeMax = g$Inc_period_max +  g$Nestling_period_max

varInc = g$Inc_period_max - g$Inc_period_min

varNestling = g$Nestling_period_max - g$Nestling_period_min

varNest_time = nest_timeMax - nest_timeMin

hist(varNestling)

hist(varNest_time)

plot(nest_timeMin~g$BodyMass.Value)
nestTimeMod = lm(nest_timeMin~g$BodyMass.Value)

names(g)

nType = g$nest_type
nType[nType == '']<-'cup'
nType[nType == 'Cup-cavity']<-'cavity'
nType = factor(nType)

nTimeByMass = nest_timeMin/g$BodyMass.Value


df1 = data.frame(nType, nTimeByMass)

plot(df1$nType,df1$nTimeByMass)
summary(lm(nTimeByMass~nType, data = df1))

plot(varNest_time~g$BodyMass.Value)
summary(lm(varNest_time~g$BodyMass.Value))
summary(lm(varNest_time~nType))
summary(lm(varNest_time~nest_timeMin))

var_nTimeByMass = varNest_time/g$BodyMass.Value

plot(varNest_time~nType)
plot(varNestling~nType)
plot(varInc~nType)



plot(var_nTimeByMass~nType)

varHeight = g$nest_max - g$nest_min

plot(varHeight~g$nest_mean)
summary(lm(varHeight~g$nest_mean))

cor.test(varHeight,g$nest_mean)

varHeightbyMean = varHeight/g$nest_mean

varHeight = na.omit(varHeight)

varHeight = g$nest_max - g$nest_min
varHeightScaled = numeric()
for (i in 1:length(varHeight)){
  if (is.na(varHeight[i])) varHeightScaled[i] = NA 
  else  varHeightScaled[i] = varHeight[i]/max(na.omit(varHeight))
  }

hist(varHeightScaled)
hist(varHeightScaled)

maxOff = numeric()
minOff = numeric()
varOff = numeric()
varOffScaled = numeric()
for(i in 1:length(nest[,1])){
  if 
}
n2$Clutch_max*n2$Brood_max
minOff = n2$Clutch_min*n2$Brood_min

varOff = (maxOff-minOff)/maxOff

varOff = 



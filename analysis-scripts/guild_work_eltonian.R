# Prep guild file

setwd('/Users/bsevans/gits/Community_analysis/guild_data')

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

# Remove a few unnecessary columns:

names(g)

# Change Nest location categories:
g$nest_loc[1] <- 'Tree'
g$nest_loc = as.character(g$nest_loc)

gsub(g$nest_loc[i], 'Tree-branch', 'Tree')

for(i in 1:length(g[,1])){
  if (g$nest_loc[i] == 'Tree-branch') g$nest_loc[i] <-  'Tree'
  if (g$nest_loc[i] == 'Ground-tree-branch-shrub') g$nest_loc[i] <-  'Generalist'
  if (g$nest_loc[i] == 'Tree-trunk') g$nest_loc[i] <-  'Tree'
  if (g$nest_loc[i] == 'Shrub-tree-branch') g$nest_loc[i] <- 'Shrub-tree'
}

g$nest_loc = tolower(g$nest_loc)
g$nest_loc = factor

for(i in 1:length(g[,1])){
  if (g$nest_loc[i] == 'Tree-branch') gsub('Tree-branch', 'Tree',g$nest_loc[i])
}

# Remove a few unnecessary columns:


hist(g$Inc_period_min)

g = g[,-c(31:32)]

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



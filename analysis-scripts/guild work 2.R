setwd('C:/Users/Brian/Documents/community_analysis/derived_data')

g = read.csv('guilds.csv')

fix(g)
g

g$habitat = as.character(rep(NA,length(g[1])))
g$migratory = as.character(rep(NA,length(g[1])))

fix(g)

table(g$foraging)



g1 = g

g = g1

table(g$nest)

g$nest = ifelse(g$nest == 'cup-tree ','cup-tree',g$nest)
g$nest = ifelse(g$nest == 'cup-gen ','cup-shrub',g$nest)

table(g$trophic)

#g$trophic = ifelse(g$trophic == 'insct/om', 'omnivore',g$trophic)

write.csv(g,'guilds2.csv')

# Add population trend data:



g$trend = as.character(rep(NA,length(g[1])))
g$trend.sig = as.character(rep(NA,length(g[1])))

g1 = g
g = g1

g$trend = as.numeric(g$trend)
g$trend
fix(g)

write.csv(g, 'pop.trends.csv', row.names = F)

g1 = g
g = g[g$trend.sig == 'y',]

table(g[g$trend>0,]$migratory)

table(g[g$trend<0,]$migratory)


table(g[g$trend>0,]$habitat)

table(g[g$trend<0,]$habitat)

table(g1[g1$trend.sig=='n',]$habitat)

g1


table(g[g$trend>0,]$trophic)

table(g[g$trend<0,]$trophic)

# ###
se = function(x) sd(x)/sqrt(length(x))

g2= g1[g1$trend.sig!= 'NA',]

g2 = g2[g2$trend.sig == 'y',]

tapply(g2$trend,g2$nest, mean)

tapply(g2$trend,g2$nest, se)

tapply(g2$trend,g2$foraging, mean)

tapply(g2$trend,g2$foraging, se)

head(g2)


tapply(g2$trend,g2$trophic, mean)

tapply(g2$trend,g2$trophic, se)



tapply(g2$trend,g2$habitat, mean)

tapply(g2$trend,g2$habitat, se)

g2$migratory = factor(g2$migratory)
summary(g2)
tapply(g2$trend,g2$migratory, mean)

#

mig.means = tapply(g2$trend,g2$migratory, mean)
mig.se = tapply(g2$trend,g2$migratory, se)
y = 1:3

plot(mig.means,y, pch = 19, xlim = c(-2,2), cex = 1.5)
means = tapply(g2$trend,g2$migratory, mean)
cis = mig.se*1.96
lcl = means - cis
ucl = means + cis
lse = means - mig.se
use = means + mig.se

arrows(x0 = lcl, y0 = y,x1 = ucl, angle = 0)
arrows(x0 = lse, y0 = y,x1 = use, angle = 0, lwd = 3)

abline(v = 0, lty =2)

#
g2 = g2[!is.na(g2$trend),]
g2$trophic = factor(g2$trophic)


means = tapply(g2$trend,g2$trophic, mean)
ses = tapply(g2$trend,g2$trophic, se)
y = 1:5

plot(means,y, pch = 19, xlim = c(-5,5), cex = 1.5)
cis = ses*1.96
lcl = means - cis
ucl = means + cis
lse = means - ses
use = means + ses

arrows(x0 = lcl, y0 = y,x1 = ucl, angle = 0)
arrows(x0 = lse, y0 = y,x1 = use, angle = 0, lwd = 3)

abline(v = 0, lty =2)

#
g2.all = g2
g2$trend.sig = factor(g2$trend.sig)
g2 = g2[g2$trend.sig =='y',]
g2 = g2[!is.na(g2$trend),]

g2$trophic = factor(g2$trophic)

# Nest

means = tapply(g2$trend,g2$nest, mean)
ses = tapply(g2$trend,g2$nest, se)
y = 1:5

par(mar = c(5,5,4,3)+.1)
plot(means,y, pch = 19, xlim = c(-5,5), cex = 1.25, ylim = c(.5,5.5),
     yaxt = 'n', ylab = '',xlab='Mean population change 1966-2012 (%)', cex.lab = 1.5)
cis = ses*1.96
lcl = means - cis
ucl = means + cis
lse = means - ses
use = means + ses


arrows(x0 = lcl, y0 = y,x1 = ucl, angle = 0, lwd = 1.5)
arrows(x0 = lse, y0 = y,x1 = use, angle = 0, lwd = 3)
axis(2, c(1,2,3,4,5),c('Cavity','Cup-S','Cup-T','Edif.','Ground'),cex.axis = 1.4, las =1,hadj =1)

abline(v = 0, lty =2)

# Trophic


means = tapply(g2$trend,g2$trophic, mean)
means = means[c(3,5)]
ses = tapply(g2$trend,g2$trophic, se)
ses = ses[c(3,5)]


mmeans = tapply(g2$trend,g2$migratory, mean)
mses = tapply(g2$trend,g2$migratory, se)

nmeans = tapply(g2$trend,g2$nest, mean)
nses = tapply(g2$trend,g2$nest, se)

means = c(means, mmeans,nmeans)
ses = c(ses, mses,nses)
cis = ses*1.96
lcl = means - cis
ucl = means + cis
lse = means - ses
use = means + ses

y = 1:10


par(mar = c(5,8,4,3)+.1)
plot(means,y, pch = 19, xlim = c(-8,8), cex = 1.25, ylim = c(.5,10.5),
     yaxt = 'n', ylab = '',xlab='Mean population change 1966-2012 (%)', cex.lab = 1.5)


arrows(x0 = lcl, y0 = y,x1 = ucl, angle = 0, lwd = 1.5)
arrows(x0 = lse, y0 = y,x1 = use, angle = 0, lwd = 3)
axis(2, y,c('InsectiV','OmniV','NeotM','Res','S-DM','Cavity', 'Cup-Shb','Cup-Tr','Edif','Gnd'),cex.axis = 1.2, las =1,hadj =1)
axis(2, c(1.5,4,8.5),c('Troph','Migrat.','Nesting'),cex.axis = 1.5, line = 4.5, tick =F)
abline(v = 0, lty =2)
abline(h=2.5)
abline(h=5.5)



############
# Habitat
#############

g2$habitat = factor(g2$habitat)

means = tapply(g2$trend,g2$habitat, mean)
ses = tapply(g2$trend,g2$habitat, se)
y = 1:4

plot(means,y, pch = 19, xlim = c(-5,5), cex = 1.25, ylim = c(.5,4.5),
     yaxt = 'n', ylab = '',xlab='Mean population change 1966-2012 (%)', cex.lab = 1.5)
cis = ses*1.96
lcl = means - cis
ucl = means + cis
lse = means - ses
use = means + ses

arrows(x0 = lcl, y0 = y,x1 = ucl, angle = 0, lwd = 1.5)
arrows(x0 = lse, y0 = y,x1 = use, angle = 0, lwd = 3)
axis(2, c(1,2,3,4),c('Forest','Open','Scrub','Water'),cex.axis = 1.4, las =1,hadj =1)

abline(v = 0, lty =2)



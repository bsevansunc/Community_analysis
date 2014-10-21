# Community analysis

# Load data:

setwd('C:/Users/Brian/Documents/community_analysis/derived_data')

pc = read.csv('point_counts.csv')
lc = read.csv('pts_lc.csv')


names(lc) = c('site','x','y','imp')

library(outliers)

#====================================================================
# Data prep
#====================================================================

# Merge count and land cover data:

p = merge(pc, lc,  by = 'site', all.x = F, all.y =F)

# Remove unused columns:
names(p)

p = p[,c(1,3,4,6:8)]

colnames(p) = c('site','species','abund','x','y','imp')

# Remove duplicate columns

p = unique(p)

# Add a presence/absence column for richness

p$pa = ifelse(p$abund>0,1,0)

# Re-order the columns:

p = p[,c(1:3,7,4:6)]

# Remove unused factor levels:

p$site = factor(p$site)

#====================================================================
# Total abundance
#====================================================================

t1 = data.frame(unique(p$site),tapply(p$abund,p$site, sum))

t1.imp = unique(data.frame(p$site,p$imp))

t1$imp = t1.imp$p.imp

t1 = t1[,-1]

colnames(t1) = c('abund','imp')

head(t1)

tab = t1

hist(tab$abund)

summary(tab$abund)
plot(tab$abund~tab$imp)

grubbs.test(tab$abund)

tab1 = tab[tab$abund<110,]

grubbs.test(tab1$abund)

tab1 = tab1[tab1$abund<90,]

grubbs.test(tab1$abund)

tab1 = tab1[tab1$abund<70,]

grubbs.test(tab1$abund)

tab1 = tab1[tab1$abund<60,]

grubbs.test(tab1$abund)

hist(tab1$abund)

plot(tab1$abund~tab1$imp)

summary(lm(tab1$abund~tab1$imp+(I(tab1$imp^2))))

summary(tab1$abund)

# abline(lm(tab1$abund~tab1$imp+(I(tab1$imp^2))))

# Ack! Richness?
p1 = p[p$pa>0,]

p1 = p1[p1$species!='UNID',]
p1$species = factor(p1$species)

p2 = data.frame(unique(p1$site),tapply(p1$pa,p1$site, sum))

head(p2)

colnames(p2) = c('site','rich')

r1 = merge(p2,lc, all = F)

mod = glm(rich~imp+I(imp^2),data = r1,family = 'poisson')

tframe = data.frame(predict(mod, type = 'response', se.fit = T))
tframe$imp = r1$imp

tframe = tframe[order(tframe$imp),]

tf2 = tframe
tf2$lcl = tf2$fit-1.96*tf2$se.fit
tf2$ucl = tf2$fit+1.96*tf2$se.fit

summary(lm(rich~imp+I(imp^2), data = r1))
plot(rich~imp, data = r1, 
     xlab = '',ylab = '',
     bty = 'l',
     pch = 19, col ='gray80', cex = .75)
mtext('Impervious surface within 100 m (%)', 1,3, cex = 1.3)
mtext('Species richness', 2,2.5, cex = 1.3)

lines(tf2$fit~tf2$imp, lwd = 2)
lines(tf2$lcl~tf2$imp, lwd = 1, lty =2)
lines(tf2$ucl~tf2$imp, lwd = 1, lty =2)

summary(lm(rich~imp+I(imp^2), data = r1))

mod.imp.imp2 = glm(rich~imp+I(imp^2),data = r1,family = 'poisson')
mod.imp = glm(rich~imp,data = r1,family = 'poisson')
mod.null = glm(rich~1,data = r1,family = 'poisson')

anova(mod.imp,mod.null, test ='Chisq')
anova(mod.imp,mod.imp.imp2, test ='Chisq')



# R-squared?

y = r1$rich
y.mod = predict(mod,type='response')

r2 = 1-(sum((y-y.mod)^2)/sum((y-mean(y))^2))

#===============================================================================================
head(r1)
head(p)
head(g2)

g2 = g1

t1 = merge(p,g2, all =F)

head(t1)
table(t1$trophic)
table(t1$trophic)/sum(table(t1$trophic))*100

table(t1$nest)/sum(table(t1$nest))*100

table(t1$habitat)/sum(table(t1$habitat))*100

table(t1$nest)/sum(table(t1$nest))*100

table(t1$migratory)/sum(table(t1$migratory))*100


troph = aggregate(t1$abund, by = list(t1$site, t1$trophic), sum)
  colnames(troph) = c('site','troph', 'count')
site.sum = aggregate(t1$abund, by = list(t1$site), sum)
  colnames(site.sum) = c('site', 't.abund')
imps = unique(data.frame(t1$site,t1$imp))
imps = imps[order(imps$t1.site),]
  colnames(imps) = c('site','imp')

troph1 = merge(site.sum, troph)
troph1 = merge(troph1, imps)
troph1$r.abund = troph1$count/troph1$t.abund

###############################################################
mod = glm(rich~imp+I(imp^2),data = r1,family = 'poisson')

tframe = data.frame(predict(mod, type = 'response', se.fit = T))
tframe$imp = r1$imp

tframe = tframe[order(tframe$imp),]

tf2 = tframe
tf2$lcl = tf2$fit-1.96*tf2$se.fit
tf2$ucl = tf2$fit+1.96*tf2$se.fit

plot(rich~imp, data = r1, 
     xlab = '',ylab = '',
     bty = 'l',
     pch = 19, col ='gray80', cex = .75)
mtext('Impervious surface within 100 m (%)', 1,3, cex = 1.3)
mtext('Species richness', 2,2.5, cex = 1.3)

lines(tf2$fit~tf2$imp, lwd = 2)
lines(tf2$lcl~tf2$imp, lwd = 1, lty =2)
lines(tf2$ucl~tf2$imp, lwd = 1, lty =2)
###############################################################
# OMNIVORES
###############################################################

ins_om = troph1[troph1$troph == 'omnivore',]

mod.omni = lm(I(log(1+count))~imp+I(imp^2),data = ins_om)

omni.frame = data.frame(predict(mod.omni, type = 'response', se.fit = T))
omni.frame$imp = ins_om$imp

omni.frame = omni.frame[order(omni.frame$imp),]

tf2 = omni.frame
tf2$lcl = tf2$fit-1.96*tf2$se.fit
tf2$ucl = tf2$fit+1.96*tf2$se.fit

plot(I(log(1+count))~imp, data = ins_om, 
     xlab = '',ylab = '',
     bty = 'l',
     pch = 19, col ='gray80', cex = .75,
     main = 'Omnivore', cex.main = 2)
mtext('Impervious surface within 100 m (%)', 1,3, cex = 1.3)
mtext('Log abundance', 2,2.5, cex = 1.3)

lines(tf2$fit~tf2$imp, lwd = 2)
lines(tf2$lcl~tf2$imp, lwd = 1, lty =2)
lines(tf2$ucl~tf2$imp, lwd = 1, lty =2)

summary(mod.omni)

legend('topright', c('p < 0.001',expression(R^2*' = 0.07')), border = 'white', box.col = 'white', cex = 1.5)

###############################################################
# insectivoreS
###############################################################

ins = troph1[troph1$troph == 'insectivore',]

mod.ins = lm(I(log(1+count))~imp,data = ins)

ins.frame = data.frame(predict(mod.ins, type = 'response', se.fit = T))
ins.frame$imp = ins$imp

ins.frame = ins.frame[order(ins.frame$imp),]

tf2 = ins.frame
tf2$lcl = tf2$fit-1.96*tf2$se.fit
tf2$ucl = tf2$fit+1.96*tf2$se.fit

plot(I(log(1+count))~imp, data = ins, 
     xlab = '',ylab = '',
     bty = 'l',
     pch = 19, col ='gray80', cex = .75,
     main = 'Insectivore', cex.main = 2)
mtext('Impervious surface within 100 m (%)', 1,3, cex = 1.3)
mtext('Log abundance', 2,2.5, cex = 1.3)

lines(tf2$fit~tf2$imp, lwd = 2)
lines(tf2$lcl~tf2$imp, lwd = 1, lty =2)
lines(tf2$ucl~tf2$imp, lwd = 1, lty =2)

summary(mod.ins)

legend('topright', c('p < 0.001',expression(R^2*' = 0.23')), border = 'white', box.col = 'white', cex = 1.5)

#
###############################################################
# granivoreS
###############################################################

gran = troph1[troph1$troph == 'granivore',]

mod.gran = lm(I(log(1+count))~imp,data = gran)

gran.frame = data.frame(predict(mod.gran, type = 'response', se.fit = T))
gran.frame$imp = gran$imp

gran.frame = gran.frame[order(gran.frame$imp),]

tf2 = gran.frame
tf2$lcl = tf2$fit-1.96*tf2$se.fit
tf2$ucl = tf2$fit+1.96*tf2$se.fit

plot(I(log(1+count))~imp, data = gran, 
     xlab = '',ylab = '',
     bty = 'l',
     pch = 19, col ='gray80', cex = .75,
     main = 'Granivore', cex.main = 2)
mtext('Impervious surface within 100 m (%)', 1,3, cex = 1.3)
mtext('Log abundance', 2,2.5, cex = 1.3)

lines(tf2$fit~tf2$imp, lwd = 2)
lines(tf2$lcl~tf2$imp, lwd = 1, lty =2)
lines(tf2$ucl~tf2$imp, lwd = 1, lty =2)

summary(mod.gran)

legend('topright', c('p < 0.001',expression(R^2*' = 0.23')), border = 'white', box.col = 'white', cex = 1.5)

#

#

#

#################################################################
# NESTING GUILDS
#################################################################

nest = aggregate(t1$abund, by = list(t1$site, t1$nest), sum)
colnames(nest) = c('site','nest', 'count')
site.sum = aggregate(t1$abund, by = list(t1$site), sum)
colnames(site.sum) = c('site', 't.abund')
imps = unique(data.frame(t1$site,t1$imp))
imps = imps[order(imps$t1.site),]
colnames(imps) = c('site','imp')

nest1 = merge(site.sum, nest)
nest1 = merge(nest1, imps)

nest1$r.abund = nest1$count/nest1$t.abund

nest1$labund = log(nest1$count+1)

ground = nest1[nest1$nest == 'ground',]
plot(ground$labund~ground$imp)
abline(lm(ground$r.abund~ground$imp))
summary(lm(ground$r.abund~ground$imp))

cavity = nest1[nest1$nest == 'cavity',]
plot(cavity$labund~cavity$imp)
abline(lm(cavity$r.abund~cavity$imp))
summary(lm(cavity$r.abund~cavity$imp))

mod.cav = lm(labund~imp, data = cavity)
cav.frame = data.frame(predict(mod.cav, type = 'response', se.fit = T))
cav.frame$imp = cavity$imp

cav.frame = cav.frame[order(cav.frame$imp),]

tf2 = cav.frame
tf2$lcl = tf2$fit-1.96*tf2$se.fit
tf2$ucl = tf2$fit+1.96*tf2$se.fit

plot(labund~imp, data = cavity, 
     xlab = '',ylab = '',
     bty = 'l',
     pch = 19, col ='gray80', cex = .75,
     main = 'Cavity nesters', cex.main = 2)
mtext('Impervious surface within 100 m (%)', 1,3, cex = 1.3)
mtext('Log abundance', 2,2.5, cex = 1.3)

lines(tf2$fit~tf2$imp, lwd = 2)
lines(tf2$lcl~tf2$imp, lwd = 1, lty =2)
lines(tf2$ucl~tf2$imp, lwd = 1, lty =2)

summary(mod.cav)

legend('topright', c('p < 0.001',expression(R^2*' = 0.14')), border = 'white', box.col = 'white', cex = 1.5)


# Edificarians

edif = nest1[nest1$nest == 'edif',]
plot(edif$labund~edif$imp)
abline(lm(edif$r.abund~edif$imp))
summary(lm(edif$r.abund~edif$imp))

mod.edif = lm(labund~imp, data = edif)
edif.frame = data.frame(predict(mod.edif, type = 'response', se.fit = T))
edif.frame$imp = edif$imp

edif.frame = edif.frame[order(edif.frame$imp),]

tf2 = edif.frame
tf2$lcl = tf2$fit-1.96*tf2$se.fit
tf2$ucl = tf2$fit+1.96*tf2$se.fit

plot(labund~imp, data = edif, 
     xlab = '',ylab = '',
     bty = 'l',
     pch = 19, col ='gray80', cex = .75,
     main = 'Edificarians', cex.main = 2)
mtext('Impervious surface within 100 m (%)', 1,3, cex = 1.3)
mtext('Log abundance', 2,2.5, cex = 1.3)

lines(tf2$fit~tf2$imp, lwd = 2)
lines(tf2$lcl~tf2$imp, lwd = 1, lty =2)
lines(tf2$ucl~tf2$imp, lwd = 1, lty =2)

summary(mod.edif)

legend('topright', c('p < 0.001',expression(R^2*' = 0.16')), border = 'white', box.col = 'white', cex = 1.5)


# Cup shrub


cup_shrub = nest1[nest1$nest == 'cup-shrub',]
plot(labund~imp)
abline(lm(cup_shrub$r.abund~cup_shrub$imp))
summary(lm(labund~imp+I(imp^2),data = cup_shrub))


cup_tree = nest1[nest1$nest == 'cup-tree',]
plot(cup_tree$r.abund~cup_tree$imp)
abline(lm(cup_tree$r.abund~cup_tree$imp))
summary(lm(labund~imp+I(imp^2),data = cup_tree))

#################################################################

habitat = aggregate(t1$abund, by = list(t1$site, t1$habitat), sum)
colnames(habitat) = c('site','habitat', 'count')
site.sum = aggregate(t1$abund, by = list(t1$site), sum)
colnames(site.sum) = c('site', 't.abund')
imps = unique(data.frame(t1$site,t1$imp))
imps = imps[order(imps$t1.site),]
colnames(imps) = c('site','imp')

habitat1 = merge(site.sum, habitat)
habitat1 = merge(habitat1, imps)
habitat1$labund = log(habitat1$count+1)

# Forest


forest = habitat1[habitat1$habitat == 'forest',]


mod.forest = lm(labund~imp, data = forest)
forest.frame = data.frame(predict(mod.forest, type = 'response', se.fit = T))
forest.frame$imp = forest$imp

forest.frame = forest.frame[order(forest.frame$imp),]

tf2 = forest.frame
tf2$lcl = tf2$fit-1.96*tf2$se.fit
tf2$ucl = tf2$fit+1.96*tf2$se.fit

plot(labund~imp, data = forest, 
     xlab = '',ylab = '',
     bty = 'l',
     pch = 19, col ='gray80', cex = .75,
     main = 'Forest habitat', cex.main = 2)
mtext('Impervious surface within 100 m (%)', 1,3, cex = 1.3)
mtext('Log abundance', 2,2.5, cex = 1.3)

lines(tf2$fit~tf2$imp, lwd = 2)
lines(tf2$lcl~tf2$imp, lwd = 1, lty =2)
lines(tf2$ucl~tf2$imp, lwd = 1, lty =2)

summary(mod.forest)

legend('topright', c('p < 0.001',expression(R^2*' = 0.12')), border = 'white', box.col = 'white', cex = 1.5)


# Open habitat

# Open


open = habitat1[habitat1$habitat == 'open',]


mod.open = lm(labund~imp+I(imp^2), data = open)
open.frame = data.frame(predict(mod.open, type = 'response', se.fit = T))
open.frame$imp = open$imp

open.frame = open.frame[order(open.frame$imp),]

tf2 = open.frame
tf2$lcl = tf2$fit-1.96*tf2$se.fit
tf2$ucl = tf2$fit+1.96*tf2$se.fit

plot(labund~imp, data = open, 
     xlab = '',ylab = '',
     bty = 'l',
     pch = 19, col ='gray80', cex = .75,
     main = 'open habitat', cex.main = 2)
mtext('Impervious surface within 100 m (%)', 1,3, cex = 1.3)
mtext('Log abundance', 2,2.5, cex = 1.3)

lines(tf2$fit~tf2$imp, lwd = 2)
lines(tf2$lcl~tf2$imp, lwd = 1, lty =2)
lines(tf2$ucl~tf2$imp, lwd = 1, lty =2)

summary(mod.open)

legend('topright', c('p < 0.001',expression(R^2*' = 0.12')), border = 'white', box.col = 'white', cex = 1.5)


#################################################################

habitat = aggregate(t1$abund, by = list(t1$site, t1$habitat), sum)
colnames(habitat) = c('site','habitat', 'count')
site.sum = aggregate(t1$abund, by = list(t1$site), sum)
colnames(site.sum) = c('site', 't.abund')
imps = unique(data.frame(t1$site,t1$imp))
imps = imps[order(imps$t1.site),]
colnames(imps) = c('site','imp')

habitat1 = merge(site.sum, habitat)
habitat1 = merge(habitat1, imps)

habitat1$r.abund = habitat1$count/habitat1$t.abund
habitat1$habitat = factor(habitat1$habitat)
habitat1$labund = log(habitat1$count + 1)

forest = habitat1[habitat1$habitat == 'forest',]
plot(labund~imp, data = forest)
abline(lm(forest$r.abund~forest$imp))
summary(lm(forest$r.abund~forest$imp))

mod = lm(labund~imp, data = forest)


open = habitat1[habitat1$habitat == 'open',]
plot(labund~imp, data = open)
abline(lm(open$r.abund~open$imp))
summary(lm(labund~imp+I(imp^2), data = open))

#

migratory = aggregate(t1$abund, by = list(t1$site, t1$migratory), sum)
colnames(migratory) = c('site','migratory', 'count')
site.sum = aggregate(t1$abund, by = list(t1$site), sum)
colnames(site.sum) = c('site', 't.abund')
imps = unique(data.frame(t1$site,t1$imp))
imps = imps[order(imps$t1.site),]
colnames(imps) = c('site','imp')

migratory1 = merge(site.sum, migratory)
migratory1 = merge(migratory1, imps)

migratory1$r.abund = migratory1$count/migratory1$t.abund
migratory1$migratory = factor(migratory1$migratory)
migratory1$labund = log(migratory1$count+1)

nt = migratory1[migratory1$migratory == 'nt',]
plot(labund~imp, data = nt)
abline(lm(nt$r.abund~nt$imp))
summary(lm(labund~imp, data = nt))


res = migratory1[migratory1$migratory == 'r',]
plot(res$r.abund~res$imp)
abline(lm(res$r.abund~res$imp))
summary(lm(labund~imp, data = res))



tmig = migratory1[migratory1$migratory == 't',]
plot(tmig$r.abund~tmig$imp)
abline(lm(tmig$r.abund~tmig$imp))
summary(lm(tmig$r.abund~tmig$imp))

summary(lm(trend~trophic+nest+habitat+migratory, data = g))

g$trend.sig = factor(g$trend.sig)

summary(g)

summary(lm(trend~nest, data = g))

summary(lm(trend~migratory, data = g))

summary(lm(trend~trophic, data = g))


summary(lm(abund~imp+I(imp^2), data = t1))
plot(abund~imp, data = t1)

summary(t1)
t1$habitat = factor(t1$habitat)
t1$migratory = factor(t1$migratory)
t1a = t1
t1$trophic = ifelse(t1$trophic == 'insct/om','omnivore',t1$trophic)
summary(t1)
t1$trophic = factor(t1$trophic)
summary(t1)
t1= t1a
summary(t1)

mod = lm(abund~imp+I(imp^2)+imp*trophic+imp*nest+imp*habitat+imp*migratory, data = t1)
mod = lm(abund~imp*nest, data = t1)

# 
t1
head(t1)
summary(t1)


t.insectivore = t1[t1$trophic == 'insectivore',]
d.insectivore = density(t.insectivore$imp, from = 0, to = 100)
plot(d.insectivore, col = 'red', lwd = 2, main = '', bty = 'l', 
     xlab = 'Impervious surface within 100 m (%)', cex.lab = 1.5)
    legend('topright','Insectivore',bty = 'n', cex = 2)

t.omnivore = t1[t1$trophic == 'insct/om'|t1$trophic == 'omnivore',]
d.omnivore = density(t.omnivore$imp)
plot(d.omnivore, col = 'blue', lwd = 2, main = '', bty = 'l',
  xlab = 'Impervious surface within 100 m (%)', cex.lab = 1.5)
  legend('topright','Omnivore',bty = 'n', cex = 2)
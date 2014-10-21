#--------------------------------------------------------------------------------*
# ---- SET-UP ----
#================================================================================*

# Load packages:

library(vegan)
library(reshape)

# Get data:

setwd('~/Community_analysis')

pc = read.csv('derived-data/pc_10_14.csv')
g = read.csv('derived-data/guilds.csv')
lc = read.csv('derived-data/pts_lc100.csv')

#--------------------------------------------------------------------------------*
# ---- PREPARE DATA FOR ANALYSIS ----
#================================================================================*

# For this run, I will use the max count per location across distance classes:

pc$counts = pc[,7]+pc[,8]+pc[,9]+pc[,10]+pc[,11]

# Prepare count frame

pc = pc[,-c(3:5,7:15)]

t = melt(pc, id = c('site','year','species'))

t2 = cast(t, site ~ species ~ variable,max)

t3 = data.frame(t2)

t3m = as.matrix(t3)

t3m[t3m == -Inf]<-0

t4 = data.frame(t3m)

colnames(t4) = gsub('.counts','',colnames(t4))

# One of the sites has never had a count, remove it:

t4$sums = rowSums(t4)

t4 = t4[t4$sums!=0,-56]

# Ensure all lc and pc files are covered (and vice-versa):

lc = unique(lc)

lc.sites = data.frame(lc[,1])
  colnames(lc.sites) = 'site'

t4$site = row.names(t4)

t5 = merge(t4, lc.sites, all = F)

pc.sites = data.frame(t5$site)
  names(pc.sites) = 'site'

pc = t5[,-1]

lc1 = merge(lc, pc.sites, all = F, incomparables = NA)

lc = lc1

#--------------------------------------------------------------------------------*
# ---- ORDINATION ----
#================================================================================*

mod1 = metaMDS(pc)
e1 = envfit(mod1~imp, data=lc)

plot(mod1, type ='n')
points(mod1, display = "sites", cex = 0.5, pch=19, col="red")
plot(e1)


lc$lc.li = lc$lc21+lc$lc22
lc$lc.hi = lc$lc23+lc$lc24


ord <- cca(pc ~ lc.li+lc.hi+lc.forest+lc.ag, data=lc)

ord2 <- cca(pc ~ imp+can, data=lc)

# Exploring output

plot(ord, type = "n")
points(ord, display = "sites", cex = 0.5, pch=19, col="red")
text(ord, display = "spec", cex=0.7, col="darkgreen")
ord.fit <- envfit(ord ~ lc.li+lc.hi+lc.forest+lc.ag, data=lc, perm=1000)
plot(ord.fit)

# Plot again, zooming on center region:
plot(ord, type = "n", xlim = c(-2,2),ylim = c(-2,2))
points(ord, display = "sites", cex = 0.5, pch=19, col="red")
text(ord, display = "spec", cex=0.7, col="darkgreen")
ord.fit <- envfit(ord ~ lc.li+lc.hi+lc.forest+lc.ag, data=lc, perm=1000)
plot(ord.fit)

# High intensity:
plot(ord, type = "n", xlim = c(0,2),ylim = c(0,2))
points(ord, display = "sites", cex = 0.5, pch=19, col="red")
text(ord, display = "spec", cex=0.7, col="darkgreen")

# Low intensity:
plot(ord, type = "n", xlim = c(0,2),ylim = c(-2,0))
points(ord, display = "sites", cex = 0.5, pch=19, col="red")
text(ord, display = "spec", cex=0.7, col="darkgreen")

# AG:
plot(ord, type = "n", xlim = c(-2,0),ylim = c(0,2))
points(ord, display = "sites", cex = 0.5, pch=19, col="red")
text(ord, display = "spec", cex=0.7, col="darkgreen")

# Forest:
plot(ord, type = "n", xlim = c(-2,0),ylim = c(-2,0))
points(ord, display = "sites", cex = 0.5, pch=19, col="red")
text(ord, display = "spec", cex=0.7, col="darkgreen")

# Output
anova(ord)

# Explained variance:

ord

#Ugh! Only 8%

# Partitioning:

ord.no.lc.hi = cca(pc ~ lc.li+lc.forest+lc.ag, data=lc)
ord.no.lc.hi
ord.only.lc.hi = cca(pc ~ lc.hi, data=lc)
ord.only.lc.hi

ord.no.lc.li = cca(pc ~ lc.hi+lc.forest+lc.ag, data=lc)
ord.no.lc.li
ord.only.lc.li = cca(pc ~ lc.li, data=lc)
ord.only.lc.li

ord.no.for = cca(pc ~ lc.hi+lc.li+lc.ag, data=lc)
ord.no.for
ord.only.for = cca(pc ~ lc.forest, data=lc)
ord.only.for

ord.no.ag = cca(pc ~ lc.hi+lc.li+lc.forest, data=lc)
ord.no.ag
ord.only.ag = cca(pc ~ lc.ag, data=lc)
ord.only.ag

# Ordination of imp and can only

ord2

plot(ord2, type = "n")
points(ord2, display = "sites", cex = 0.5, pch=19, col="red")
text(ord2, display = "spec", cex=0.7, col="darkgreen")
ord2.fit <- envfit(ord2 ~ imp+can, data=lc, perm=1000)
plot(ord2.fit)

# Upper left quadrant, high can, low imp:

plot(ord2, type = "n", xlim = c(-2,0),ylim = c(0,1.5))
points(ord2, display = "sites", cex = 0.5, pch=19, col="red")
text(ord2, display = "spec", cex=0.7, col="darkgreen")

# Lower left quadrant, low can, low imp:

plot(ord2, type = "n", xlim = c(-2,0),ylim = c(-1,0))
points(ord2, display = "sites", cex = 0.5, pch=19, col="red")
text(ord2, display = "spec", cex=0.7, col="darkgreen")

# Upper right quadrant, "high" can, high imp:

plot(ord2, type = "n", xlim = c(0,1),ylim = c(0,.5))
points(ord2, display = "sites", cex = 0.5, pch=19, col="red")
text(ord2, display = "spec", cex=0.7, col="darkgreen")

# Lower right quadrant, "high" can, high imp:

plot(ord2, type = "n", xlim = c(0,3),ylim = c(-1,0))
points(ord2, display = "sites", cex = 0.5, pch=19, col="red")
text(ord2, display = "spec", cex=0.7, col="darkgreen")

anova(ord2)

# Explained variance

ord2 # Ack! 0.0559!


ord2.imp <- cca(pc ~ imp, data=lc)

ord2.can <- cca(pc ~ can, data=lc)

#--------------------------------------------------------------------------------*
# ---- ORDINATION, GUILDS ----
#================================================================================*

pc = read.csv('derived-data/pc_10_14.csv')
g$species = tolower(g$species)

pcg = merge(pc,g, all = F)


#### Prepare count frame ----

# For this run, I will use the max count per location across distance classes:

pcg$counts = pcg[,7]+pcg[,8]+pcg[,9]+pcg[,10]+pcg[,11]

# Remove Bert Drake:

pcg = pcg[pcg$site!='DRAKBERMD1',]

# Separate by life history trait

pcg.troph = pcg
pcg.nest = pcg
pcg.migr = pcg

pcg$trophic = paste(pcg.troph$foraging, pcg.troph$trophic, sep='-')


trophic = aggregate(pcg$counts,by = list(pcg$site,pcg$year,pcg$trophic),sum)
  names(trophic) = c('site','year','guild','count')
nest = aggregate(pcg$counts,by = list(pcg$site,pcg$year,pcg$nest),sum)
names(nest) = c('site','year','guild','count')
mig = aggregate(pcg$counts,by = list(pcg$site,pcg$year,pcg$migratory),sum)
names(mig) = c('site','year','guild','count')


# Function to prepare count frames

prep.guild = function(life.history){
  t = melt(life.history, id = c('site','year','guild'))
  t2 = cast(t, site ~ guild,max)
  t3 = data.frame(t2)[,-1]
  t3m = as.matrix(t3)
  t3m[t3m == -Inf]<-0
  t4 = data.frame(t3m)
  t4$site = t2$site
  lc = unique(lc)
  lc.sites = data.frame(lc[,1])
  colnames(lc.sites) = 'site'
  t5 = merge(t4, lc.sites, all = F)
  pc.sites = data.frame(t2$site)
  names(pc.sites) = 'site'
  t5[,-1]
}

trophic = prep.guild(trophic)
nest = prep.guild(nest)
mig = prep.guild(mig)

ord.troph <- cca(trophic ~ lc.li+lc.hi+lc.forest+lc.ag, data=lc)
ord.nest <- cca(nest ~ lc.li+lc.hi+lc.forest+lc.ag, data=lc)
ord.mig <- cca(mig ~ lc.li+lc.hi+lc.forest+lc.ag, data=lc)

### 

# Exploring output

plot(ord.troph, type = "n",xlim = c(-1,1),ylim = c(-1,1))
ord.fit = envfit(ord ~ lc.li+lc.hi+lc.forest+lc.ag, data=lc, perm=1000)
plot(ord.fit, col =1)
text(ord.troph, display = "spec", cex=0.7, col="darkgreen")
text(ord.nest, display = "spec", cex=0.7, col="red")
text(ord.mig, display = "spec", cex=0.7, col="blue")

ord2.troph <- cca(trophic ~ imp + can, data=lc)
ord2.nest <- cca(nest ~ imp + can, data=lc)
ord2.mig <- cca(mig ~ imp + can, data=lc)

w = 2 # w stands for window

ord = cca(trophic)
plot(ord, type = "n",xlim = c(-w,w),ylim = c(-w,w))
ord.env = envfit(ord2 ~ imp + can, data=lc, perm=1000)
plot(ord.env, col =1)


head(lc)

urb = factor(ifelse(lc$imp < 5,'R',ifelse(lc$imp>60,'U','S')))
ordispider(ord, urb, col = 'blue')
points(ord, display = "sites", cex = 0.5, pch=19, col="red")

text(ord, display = "spec", cex=0.7, col="darkgreen")

text(ord2, display = "spec", cex=0.7, col="blue")
text(ord2.nest, display = "spec", cex=0.7, col="red")
text(ord2.mig, display = "spec", cex=0.7, col="blue")

# Explained variance:

ord.troph # .113
ord.nest # .1371
ord.mig # 0.039 --> very low, dropping

# Partitioning:

ord.troph.no.lc.hi = cca(trophic~ lc.li+lc.forest+lc.ag, data=lc)
ord.troph.no.lc.hi
ord.troph.no.lc.li = cca(trophic~ lc.hi+lc.forest+lc.ag, data=lc)
ord.troph.no.lc.li
ord.troph.no.forest = cca(trophic~ lc.hi+lc.li+lc.ag, data=lc)
ord.troph.no.forest
ord.troph.no.ag = cca(trophic~ lc.hi+lc.li+lc.forest, data=lc)
ord.troph.no.ag

ord.nest.no.lc.hi = cca(nest~ lc.li+lc.forest+lc.ag, data=lc)
ord.nest.no.lc.hi
ord.nest.no.lc.li = cca(nest~ lc.hi+lc.forest+lc.ag, data=lc)
ord.nest.no.lc.li
ord.nest.no.forest = cca(nest~ lc.hi+lc.li+lc.ag, data=lc)
ord.nest.no.forest
ord.nest.no.ag = cca(nest~ lc.hi+lc.li+lc.forest, data=lc)
ord.nest.no.ag

ord2.troph.no.imp = cca(trophic~ can, data=lc)
ord2.troph.no.can = cca(trophic~ imp, data=lc)

ord2.troph
ord2.troph.no.imp
ord2.troph.no.can

ord2.nest.no.imp = cca(nest~ can, data=lc)
ord2.nest.no.can = cca(nest~ imp, data=lc)
ord2.nest.no.imp
ord2.nest.no.can


plot(ord2.troph, type = "n",xlim = c(-.75,.75),ylim = c(-.75,.75))
ord2.fit = envfit(ord2 ~ imp + can, data=lc, perm=1000)
plot(ord2.fit, col =1)
text(ord2.troph, display = "spec", cex=0.7, col="darkgreen")
text(ord2.nest, display = "spec", cex=0.7, col="red")


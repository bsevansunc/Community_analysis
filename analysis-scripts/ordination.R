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
# ---- Log-transform count data ----
#================================================================================*

head(pc)

pc2 = pc[,-53]

for(i in 1:dim(pc2)[2]){
  pc2[,i] = log(1+pc2[,i])
}

#--------------------------------------------------------------------------------*
# ---- ORDINATION BY SPECIES ----
#================================================================================*

ord2 <- cca(pc ~ imp+can, data=lc)

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

# Lower right quadrant, low can, high imp:

plot(ord2, type = "n", xlim = c(0,5),ylim = c(-1,0))
points(ord2, display = "sites", cex = 0.5, pch=19, col="red")
text(ord2, display = "spec", cex=0.7, col="darkgreen")

anova(ord2)

# Explained variance

ord2 # Ack! 0.0719!

ord2.imp <- cca(pc ~ imp, data=lc) # .049 ... can explains .0229 that imp does not

ord2.can <- cca(pc ~ can, data=lc) # .0451 ... imp explains .0268 that can does not

#--------------------------------------------------------------------------------*
# ---- ORDINATION BY GUILD ----
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
  t5 = t5[,-1]
  # Convert to log-transformed counts:
  pc2 = t5  
  for(i in 1:dim(pc2)[2]){
    pc2[,i] = log(1+pc2[,i])
  }
  pc2
}

trophic = prep.guild(trophic)
nest = prep.guild(nest)
mig = prep.guild(mig)

ord.troph <- cca(trophic ~ can+imp, data=lc)
ord.nest <- cca(nest ~ can+imp, data=lc)
ord.mig <- cca(mig ~ can+imp, data=lc)

# Explore explained variance

ord.troph # .1034

cca(trophic ~ imp, data=lc) # .0861, .0173 explained by can not explained by imp

cca(trophic ~ can, data=lc) # .0652, .0382 explained by imp not explained by can

#

ord.nest # .2005

cca(nest ~ imp, data=lc) # .1826, .0179 explained by can not explained by imp

cca(nest ~ can, data=lc) # .1447, .0558 explained by imp not explained by can

#

ord.mig # .05710

cca(mig ~ imp, data=lc) # .0526, .0045 explained by can not explained by imp

cca(mig ~ can, data=lc) # .0196, .0375 explained by imp not explained by can

plot(ord.troph, type = "n",xlim = c(-.75,.75),ylim = c(-.75,.75))
ord.fit = envfit(ord.troph ~ imp + can, data=lc, perm=1000)
plot(ord.fit, col =1)
text(ord.troph, display = "spec", cex=0.75, col="darkgreen")

plot(ord.nest, type = "n",xlim = c(-.75,.75),ylim = c(-.75,.75))
ord.fit = envfit(ord.nest ~ imp + can, data=lc, perm=1000)
plot(ord.fit, col =1)
text(ord.nest, display = "spec", cex=1, col="red")

plot(ord.mig, type = "n",xlim = c(-.15,.15),ylim = c(-.15,.15))
ord.fit = envfit(ord.mig ~ imp + can, data=lc, perm=1000)
plot(ord.fit, col =1)
text(ord.mig, display = "spec", cex=1.5, col="blue")


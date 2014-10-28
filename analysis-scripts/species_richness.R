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

pc = pc2

#--------------------------------------------------------------------------------*
# ---- Counts by guild ----
#================================================================================*

pc1 = read.csv('derived-data/pc_10_14.csv')
g$species = tolower(g$species)

pcg = merge(pc1,g, all = F)


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

#--------------------------------------------------------------------------------*
# ---- Species richness ----
#================================================================================*

sr = specnumber(pc)

plot(sr~imp, data = lc, xlab = 'Impervious surface (%)', ylab = 'Species Richness')

mod.i2c2.full = glm(sr~imp + can + imp:can + I(imp^2) + I(can^2), data = lc)
mod.i2.inter = glm(sr~imp + can + imp:can + I(imp^2), data = lc)
mod.c2.inter = glm(sr~imp + can + imp:can + I(can^2), data = lc)
mod.ic.inter = glm(sr~imp + can + imp:can, data = lc)
mod.i2c2 = glm(sr~imp + can + I(imp^2) + I(can^2), data = lc)
mod.ic2 = glm(sr~imp + can + I(can^2), data = lc)
mod.c2 = glm(sr~can + I(can^2), data = lc)
mod.ic = glm(sr~can + imp, data = lc)
mod.c = glm(sr~can, data = lc)
mod.i2c = glm(sr~imp + can + I(imp^2), data = lc)
mod.i2 = glm(sr~imp + I(imp^2), data = lc)
mod.i = glm(sr~imp, data = lc)
mod.null = glm(sr~1, data = lc)


mod.list = list(mod.i2c2.full, mod.i2.inter,mod.c2.inter, mod.ic.inter, 
                mod.i2c2, mod.ic2, mod.c2, mod.ic, mod.c, mod.i2c, mod.i2,
                mod.i, mod.null)



aic_out = numeric()
formulas = character()
deviance = numeric()

for(i in 1:length(mod.list)){
  formulas[i] = as.character(mod.list[[i]]$formula[3])
  aic_out[i] = AIC(mod.list[[i]])
  deviance[i] = mod.list[[i]]$deviance
}

out.df = data.frame(formulas, aic_out, deviance)

out.df = out.df[order(out.df$aic_out),]

out.df$dAIC = out.df$aic_out-min(out.df$aic_out)

out.df$w = round(exp(-0.5*out.df$dAIC)/sum(exp(-0.5*out.df$dAIC)),3)

out.df = out.df[,c(1:2, 4:5, 3)]

# Output:

out.df

summary(mod.list[[2]])

summary(lm(sr~imp + can + imp:can + I(imp^2), data = lc))

imp = seq(0,100, 1)
can = seq(0,100, 1)

mod.predictions = predict(mod.i2.inter, type = 'response')

plot(mod.predictions~lc$imp)
plot(mod.predictions~lc$can)

x = 1:100
y = 1:100

grid = mesh(x,y)
g1 = grid$x + grid$y
g1 = ifelse(g1>100, NA, 1)

z =with(grid, 14.817241-0.287151*x+0.0023908*x^2-0.0317364*y+0.0060869*x*y)

z = z*g1

t = melt(z)
t = na.omit(t)
names(t) = c('imp','can','richness')

t2 = ggplot(t, aes(imp,can,z = richness))
t2 + geom_tile(aes(fill = richness))+
  scale_fill_continuous(name = 'Richness',low = "yellow", high = "red")+
  stat_contour(bins = 20,size = .25)+
  xlab('% Impervious') + ylab('% Canopy') + 
  # Add themes:
  theme(axis.text = element_text(size=14, color = 1),
        axis.title.x = element_text(vjust = -.5),
        axis.title.y = element_text(vjust = .5),
        axis.title = element_text(size=18, vjust = -1),
        axis.line = element_line(colour = "black"),
        panel.background = element_blank())


#--------------------------------------------------------------------------------*
# ---- Guild richness ----
#================================================================================*

gr.troph = specnumber(trophic)
gr.nest = specnumber(nest)
gr.mig = specnumber(mig)

#--------------------------------------------------------------------------------*
# ---- Trophic ----
#--------------------------------------------------------------------------------*

plot(gr.troph~imp, data = lc, xlab = 'Impervious surface (%)', ylab = 'Species Richness')

mod.i2c2.full = glm(gr.troph~imp + can + imp:can + I(imp^2) + I(can^2), data = lc)
mod.i2.inter = glm(gr.troph~imp + can + imp:can + I(imp^2), data = lc)
mod.c2.inter = glm(gr.troph~imp + can + imp:can + I(can^2), data = lc)
mod.ic.inter = glm(gr.troph~imp + can + imp:can, data = lc)
mod.i2c2 = glm(gr.troph~imp + can + I(imp^2) + I(can^2), data = lc)
mod.ic2 = glm(gr.troph~imp + can + I(can^2), data = lc)
mod.c2 = glm(gr.troph~can + I(can^2), data = lc)
mod.ic = glm(gr.troph~can + imp, data = lc)
mod.c = glm(gr.troph~can, data = lc)
mod.i2c = glm(gr.troph~imp + can + I(imp^2), data = lc)
mod.i2 = glm(gr.troph~imp + I(imp^2), data = lc)
mod.i = glm(gr.troph~imp, data = lc)
mod.null = glm(gr.troph~1, data = lc)


mod.list = list(mod.i2c2.full, mod.i2.inter,mod.c2.inter, mod.ic.inter, 
                mod.i2c2, mod.ic2, mod.c2, mod.ic, mod.c, mod.i2c, mod.i2,
                mod.i, mod.null)


aic_out = numeric()
formulas = character()
deviance = numeric()

for(i in 1:length(mod.list)){
  formulas[i] = as.character(mod.list[[i]]$formula[3])
  aic_out[i] = AIC(mod.list[[i]])
  deviance[i] = mod.list[[i]]$deviance
}

out.df = data.frame(formulas, aic_out, deviance)

out.df = out.df[order(out.df$aic_out),]

out.df$dAIC = out.df$aic_out-min(out.df$aic_out)

out.df$w = round(exp(-0.5*out.df$dAIC)/sum(exp(-0.5*out.df$dAIC)),3)

out.df = out.df[,c(1:2, 4:5, 3)]

# Output:

out.df

summary(mod.list[[2]])

mod.predictions = predict(mod.i2.inter, type = 'response')

plot(mod.predictions~lc$imp)
plot(mod.predictions~lc$can)

x = seq(min(lc$imp), max(lc$imp), by = .5)
y = seq(min(lc$can), max(lc$can), by = .5)

grid<-mesh(x,y)

g1 = grid$x+grid$y

g1 = ifelse(g1>100, NA, 1)

grid$x = g1 * grid$x
grid$y = g1 * grid$y


z =with(grid, 6.868165184-0.116676485*x-0.005092819*y+0.000872133*x^2+0.001824743*x*y)

grid[[3]] = z

persp3D(z = z, x = x, y = y, theta = 25, phi = 5, 
        xlab = 'Impervious surface (%)',
        ylab = 'Canopy cover (%)',
        zlab = 'Trophic richness',ticktype ='detailed',cex.axis = .65,
        cex.lab = .8)

# Hmmmm ... this makes little sense because you can't be over 100% land cover (summed). This 
# plot, however, sugests highest richness with highest canopy cover and impervious surface.
# Need to do this with real data!

grid<-mesh(x,y)

g1 = grid$x+grid$y

g1 = ifelse(g1>100, NA, 1)


x = 1:100
y = 1:100

grid = mesh(x,y)
g1 = grid$x + grid$y
g1 = ifelse(g1>100, NA, 1)

z =with(grid, 6.868165184-0.116676485*x-0.005092819*y+0.000872133*x^2+0.001824743*x*y)

z = z*g1

t = melt(z)
t = na.omit(t)
  names(t) = c('imp','can','richness')

t2 = ggplot(t, aes(imp,can,z = richness))
  t2 + geom_tile(aes(fill = richness))+
    scale_fill_continuous(name = 'Richness',low = "yellow", high = "red")+
    stat_contour(bins = 20,size = .5)+
    xlab('% Impervious') + ylab('% Canopy') + 
    # Add themes:
  theme(axis.text = element_text(size=14, color = 1),
        axis.title.x = element_text(vjust = -.5),
        axis.title.y = element_text(vjust = .5),
        axis.title = element_text(size=18, vjust = -1),
        axis.line = element_line(colour = "black"),
        panel.background = element_blank())


#--------------------------------------------------------------------------------*
# ---- Nest ----
#--------------------------------------------------------------------------------*

plot(gr.nest~imp, data = lc, xlab = 'Impervious surface (%)', ylab = 'Species Richness')

mod.i2c2.full = glm(gr.nest~imp + can + imp:can + I(imp^2) + I(can^2), data = lc)
mod.i2.inter = glm(gr.nest~imp + can + imp:can + I(imp^2), data = lc)
mod.c2.inter = glm(gr.nest~imp + can + imp:can + I(can^2), data = lc)
mod.ic.inter = glm(gr.nest~imp + can + imp:can, data = lc)
mod.i2c2 = glm(gr.nest~imp + can + I(imp^2) + I(can^2), data = lc)
mod.ic2 = glm(gr.nest~imp + can + I(can^2), data = lc)
mod.c2 = glm(gr.nest~can + I(can^2), data = lc)
mod.ic = glm(gr.nest~can + imp, data = lc)
mod.c = glm(gr.nest~can, data = lc)
mod.i2c = glm(gr.nest~imp + can + I(imp^2), data = lc)
mod.i2 = glm(gr.nest~imp + I(imp^2), data = lc)
mod.i = glm(gr.nest~imp, data = lc)
mod.null = glm(gr.nest~1, data = lc)


mod.list = list(mod.i2c2.full, mod.i2.inter,mod.c2.inter, mod.ic.inter, 
                mod.i2c2, mod.ic2, mod.c2, mod.ic, mod.c, mod.i2c, mod.i2,
                mod.i, mod.null)


aic_out = numeric()
formulas = character()
deviance = numeric()

for(i in 1:length(mod.list)){
  formulas[i] = as.character(mod.list[[i]]$formula[3])
  aic_out[i] = AIC(mod.list[[i]])
  deviance[i] = mod.list[[i]]$deviance
}

out.df = data.frame(formulas, aic_out, deviance)

out.df = out.df[order(out.df$aic_out),]

out.df$dAIC = out.df$aic_out-min(out.df$aic_out)

out.df$w = round(exp(-0.5*out.df$dAIC)/sum(exp(-0.5*out.df$dAIC)),3)

out.df = out.df[,c(1:2, 4:5, 3)]

# Output:

out.df

summary(mod.list[[4]])

mod.predictions = predict(mod.list[[4]], type = 'response')

plot(mod.predictions~lc$imp)
plot(mod.predictions~lc$can)

x = seq(0, round(max(lc$imp),0), by = 1)
y = seq(0, round(max(lc$can),0), by = 1)

x = 1:100
y = 1:100

grid = mesh(x,y)
g1 = grid$x + grid$y
g1 = ifelse(g1>100, NA, 1)

z =with(grid, 4.3346876-0.0206286*x-0.0103436*y+0.0008807*x*y)

z = z*g1

t = melt(z)
t = na.omit(t)
names(t) = c('imp','can','richness')

t2 = ggplot(t, aes(imp,can,z = richness))
t2 + geom_tile(aes(fill = richness))+
  scale_fill_continuous(name = 'Richness',low = "yellow", high = "red")+
  stat_contour(bins = 20,size = .25)+
  xlab('% Impervious') + ylab('% Canopy') + 
  # Add themes:
  theme(axis.text = element_text(size=14, color = 1),
        axis.title.x = element_text(vjust = -.5),
        axis.title.y = element_text(vjust = .5),
        axis.title = element_text(size=18, vjust = -1),
        axis.line = element_line(colour = "black"),
        panel.background = element_blank())

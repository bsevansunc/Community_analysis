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

summary(lm(sr~imp + can + imp:can + I(imp^2) + I(can^2), data = lc))


imp = seq(0,100, 1)
can = seq(0,100, 1)

mod.predictions = predict(mod.i2.inter, type = 'response')

plot(mod.predictions~lc$imp)
plot(mod.predictions~lc$can)

?cluster
pca1 = princomp(~imp + can, data = lc)
plot(lc$imp~pca1$scores[,2])
plot(lc$can~pca1$scores[,2])
plot(lc$imp~pca1$scores[,1])
plot(lc$can~pca1$scores[,1])

plot(sr~pca1$scores[,1])
plot(sr~pca1$scores[,2])





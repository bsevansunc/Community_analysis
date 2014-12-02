# This script takes the processed point count, guild, and land cover files and
# creates a wide frame of log-tranformed count data, linked with land cover data.
#--------------------------------------------------------------------------------*
# ---- SET-UP ----
#================================================================================*

# Load packages:

library(vegan)
library(reshape2)
library(plyr)
library(unmarked)
library(AICcmodavg)
library(ggplot2)

# Get data:

pc = read.csv('derived-data/pc_10_14.csv')
g = read.csv('derived-data/guilds.csv')
lc = read.csv('derived-data/pts_lc100.csv')

# Remove Bert Drake (all zero counts?):

pc = pc[pc$site!='DRAKBERMD1',]
lc = lc[lc$site!='DRAKBERMD1',]

# Reduce site and count frames to just the sites with point counts:

lc = merge(lc, data.frame(site = unique(pc$site)), all = F)

pc = merge(pc, data.frame(site = unique(lc$site)), all = F)

# 2006 is a typo! A little digging reveals that the year is most likely 2009:

pc[pc$year == 2006,'year'] <- 2009

# Remove eust, ropi, and hosp:

pc = pc[pc$species!='eust' & pc$species!='hosp' & pc$species!='ropi',]

# Maintain only necessary columns:

pc1 = pc[,c(1,2,6:11)]

# Function to subset to a given species by species:

pc.sub = function(sp, year){
  # Subset to year and species
  df = pc1[pc1$species == sp & pc1$year == year,]
  # Remove year and species columns:
  df = df[,-c(2,3)]
  # Combine counts at a given site and distance class:
  df = ddply(df, c('site'), numcolwise(sum))
  # Merge with lc to get 0's (and convert NA's to 0):
  df = merge(df, lc, all =T)
  df[is.na(df)] <-0
  # Remove spatial location columns:
  df = df[,-c(7:8)]
  return(df)
  } 

#--------------------------------------------------------------------------------*
# ---- CALCULATE TOTAL SPECIES ABUNDANCE BY IMPERVIOUS AND CANOPY COVER ----
#================================================================================*

#---------------------------------------------------------------------*
# ---- CREATE UNMARKED FRAME ----
#---------------------------------------------------------------------*

# Function to create unmarked frame:

umf.fun = function(sp, yr) {
  pc = pc.sub(sp, yr)
  lc = pc[,c(1,7:8)]
  row.names(pc) = pc[,1]  # Sets sites as row names
  pc = pc[,-c(1, 7:8)]       # Removes sites and lc from the count frame
  # Make a new covariate frame:
  covs = data.frame(can = scale(lc$can), 
                    imp = scale(lc$imp), row.names = lc$site)
  # Create unmarked frame:
  umf <- unmarkedFrameDS(y=as.matrix(pc), 
                       siteCovs=covs, survey="point",
                       dist.breaks=c(0, 10, 20, 30, 40, 50), 
                       unitsIn="m")
  return(umf)
}

umf.fun('eato', 2009)

#---------------------------------------------------------------------*
# ---- RUN MODELS ----
#---------------------------------------------------------------------*

# Vector of candidate model formulas:

models = c('~can ~1',
           '~can ~ imp',
           '~can ~ imp + I(imp^2)',
           '~1 ~1')
           
           
# Run models for a given year and species:

mod.run = function(sp, yr){
mod.outs = list()
for(i in 1:length(models)){
  mod.outs[[i]] = distsamp(formula(models[i]),umf.fun(sp, yr), output = 'abund')
}
names(mod.outs) = models
# Model table:
fl = fitList(fits = mod.outs)
mod.tab = modSel(fl, nullmod = '~1 ~1')
return(list(fl, mod.tab))
}


# Predict model at sites:
mod.plot = function(sp, yr){
  mods = mod.run(sp, yr)
  mod.preds.site = predict(mods[[1]], type="state")
  test = cbind(lc, mod.preds.site)
  t1 = test[ order(test[,5]), ]
  plot(Predicted~imp, data = t1, type = 'l', lwd = 2, 
     xlab = '% Impervious', ylab = 'Predicted abundance',
     bty = 'l', ylim = c(0,max(upper)))
  lines(lower~imp, data = t1, type = 'l', lwd = 1, lty = 2)
  lines(upper~imp, data = t1, type = 'l', lwd = 1, lty = 2)
  }

mod.run('sosp',2009)
mod.plot('sosp',2009)


#---------------------------------------------------------------------*
# ---- Model-averaged parameter estimates ----
#---------------------------------------------------------------------*

#---------------------------------------------------------------------*
# ---- Predict model at sites ----
#---------------------------------------------------------------------*





# Backtransform function to turn z values back to original scale:

bt.cov = function(x, x.scaled) x.scaled * sd(x) + mean(x)

# Convert can and imp to original scales and add to pred frame:

mod.preds.site$can = bt.cov(lc$can, covs$can)
mod.preds.site$imp = bt.cov(lc$imp, covs$imp)

df1 = mod.preds.site
z= df1$Predicted

p = ggplot(df1, aes(imp, can))
p + geom_point(aes(color = z, size = z))+
  scale_color_gradient(low = 'yellow', high = 'red')+
  scale_shape(solid = FALSE)+
  xlab('% Impervious') + ylab('% Canopy') + 
  # Add themes:
  theme(axis.text = element_text(size=14, color = 1),
        axis.title.x = element_text(vjust = -.5),
        axis.title.y = element_text(vjust = .5),
        axis.title = element_text(size=18, vjust = -1),
        axis.line = element_line(colour = "black"),
        legend.position="none",
        panel.background = element_blank())

#---------------------------------------------------------------------*
# ---- Predict model for canopy cover (imp at mean value) ----
#---------------------------------------------------------------------*

lc.can = data.frame(can = seq(min(covs$can),max(covs$can), by = 1), imp = mean(covs$imp))
lc.can = data.frame(can = lc.can, imp = rep(0,length(lc.can)))

mod.preds.can = predict(fl, newdata = lc.can, type = 'state')

mod.preds.can$can = bt.cov(lc$can, lc.can$can)

plot(Predicted~can, data = mod.preds.can, 
     type = 'l', lwd = 2,
     ylim = c(min(mod.preds.can$lower),max(mod.preds.can$upper)),
     main = 'Across species: Canopy',
     xlab = '% Canopy',
     ylab = 'Predicted abundance',
     bty = 'l')
lines(lower~can, data = mod.preds.can,
      type = 'l', lty = 2,lwd = 2)
lines(upper~can, data = mod.preds.can,
      type = 'l', lty = 2,lwd = 2)

#---------------------------------------------------------------------*
# ---- Predict model for impervious cover (can at mean value) ----
#---------------------------------------------------------------------*

lc.imp = seq(min(covs$imp),max(covs$imp), by = .01)
lc.imp = data.frame(imp = lc.imp, can = rep(0,length(lc.imp)))

mod.preds.imp = predict(fl, newdata = lc.imp, type = 'state')

mod.preds.imp$imp = bt.cov(lc$imp, lc.imp$imp)

plot(Predicted~imp, data = mod.preds.imp, 
     type = 'l', lwd = 2,
     ylim = c(min(mod.preds.imp$lower),max(mod.preds.imp$upper)),
     main = 'Across species: Impervious',
     xlab = '% Impervious',
     ylab = 'Predicted abundance',
     bty = 'l')
lines(lower~imp, data = mod.preds.imp,
      type = 'l', lty = 2,lwd = 2)
lines(upper~imp, data = mod.preds.imp,
      type = 'l', lty = 2,lwd = 2)


################################################
backTransform(best.mod['det'])

# backTransform(linearComb(best.mod['det'], c(1,10)))

site.level.density <- predict(m1, type="state")

sld = site.level.density
sld$imp = lc$imp
sld$can = lc$can

plot(Predicted~imp, data = sld[order(sld$imp),], 
     type = 'l', lwd = 2,
     ylim = c(min(sld$lower),max(sld$upper)),
     main = 'Across species',
     xlab = '% Impervious',
     ylab = 'Predicted abundance',
     bty = 'l')
lines(lower~imp, data = sld[order(sld$imp),],
      type = 'l', lty = 2,lwd = 2)
lines(upper~imp, data = sld[order(sld$imp),],
      type = 'l', lty = 2,lwd = 2)


names(mod.outs) = models

mod.outs[[1]]

null <- distsamp(~1~1, umf, output = 'density')

can_imp <- distsamp(~can~imp), umf, output = 'density')

can_null <- distsamp(~1~1, umf, output = 'density')

hn_can <- distsamp(~can~1), umf, output = 'density')

hn_null
hn_can

backTransform(hn_null, type="state")

backTransform(hn_null, type="det")

backTransform(hn_can, type="state")

backTransform(linearComb(hn_can['det'], c(1,5)))




site.level.density <- predict(hn_null, type="state")$Predicted

plot(site.level.density~lc$imp)







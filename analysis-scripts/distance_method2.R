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

# For each site, get the number of surveys and add to the lc frame:

visits = ddply(pc, .(site), summarize, length(unique(year)))[,2]
visits = as.matrix(visits, ncol = 1)

#--------------------------------------------------------------------------------*
# ---- CALCULATE TOTAL SPECIES ABUNDANCE BY IMPERVIOUS AND CANOPY COVER ----
#================================================================================*

# For total species abundance, make a frame summarizing distance by site:

pc.ta = ddply(pc, .(site), summarize, d10 =  mean(d10))
pc.ta$d20 = ddply(pc, .(site), summarize, mean(d20))[,2]
pc.ta$d30 = ddply(pc, .(site), summarize, mean(d30))[,2]
pc.ta$d40 = ddply(pc, .(site), summarize, mean(d40))[,2]
pc.ta$d50 = ddply(pc, .(site), summarize, mean(d50))[,2]

#---------------------------------------------------------------------*
# ---- CREATE UNMARKED FRAME ----
#---------------------------------------------------------------------*

row.names(pc.ta) = pc.ta[,1]  # Sets sites as row names
pc.ta = pc.ta[,-1]            # Removes sites from the count frame

names(pc.ta) = c('[0,10]', '(10,20]', '(20,30]', '(30,40]', '(40,50]')

covs = data.frame(can = scale(lc$can), imp = scale(lc$imp), row.names = lc$site)

# covs = data.frame(can = lc$can, imp = lc$imp, row.names = lc$site)

umf <- unmarkedFrameDS(y=as.matrix(pc.ta), 
                       siteCovs=covs, survey="point",
                       dist.breaks=c(0, 10, 20, 30, 40, 50), 
                       unitsIn="m")

hist(umf, xlab="distance (m)", main="", cex.lab=0.8, cex.axis=0.8)

#---------------------------------------------------------------------*
# ---- RUN MODELS ----
#---------------------------------------------------------------------*

# Vector of candidate model formulas:

# models = c('~1 ~1',
#            '~can  ~1',
#            '~imp ~1',
#            '~can + imp ~1',
#            '~can + imp + imp:can ~1')

models = c('~1 ~1',
           '~can + imp + imp:can ~1',
           '~can + imp + imp:can  ~imp',
           '~can + imp + imp:can ~can',
           '~can + imp + imp:can ~imp + can',
           '~can + imp + imp:can ~imp + can + imp:can')
           
           
# Run models:

mod.outs = list()
for(i in 1:length(models)){
  mod.outs[[i]] = distsamp(formula(models[i]),umf, output = 'abund')
}

names(mod.outs) = models

#---------------------------------------------------------------------*
# ---- Model table ----
#---------------------------------------------------------------------*

fl = fitList(fits = mod.outs)

aictab(mod.outs, modnames = models)

modSel(fl, nullmod = '~1 ~1')

#---------------------------------------------------------------------*
# ---- Model-averaged parameter estimates ----
#---------------------------------------------------------------------*

#---------------------------------------------------------------------*
# ---- Predict model at sites ----
#---------------------------------------------------------------------*

mod.preds.site = predict(fl, type="state")

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







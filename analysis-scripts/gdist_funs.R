# Load packages:

library(reshape2)
library(plyr)
library(unmarked)
library(AICcmodavg)
library(ggplot2)

# Get data:

y0 = read.csv('gdist_data/y0.csv')
site_covs = read.csv('gdist_data/site_covs.csv')
obs_covs = read.csv('gdist_data/obs_covs.csv')

#--------------------------------------------------------------------------------*
# ---- DATA PREP ----
#================================================================================*

#--------------------------------------------------------------------------------*
# ---- Prep observation (y0) frame ----
#--------------------------------------------------------------------------------*
# Function to pull out year and make d matrix:

y.SpYr = function(species.of.interest, yr){
  # Vector of unique sites:
  sites = unique(y0[1])
  # Vector of sites for a given year:
  site.yr = unique(subset(y0, year == yr)[1])
  # Subset observation frame to year of interest:
  y.yr = subset(y0, year == yr)[,-2]
  # Subset observation frame to species of interest. If statement allows for SOI
  # to be across all species:
  if (species.of.interest == 'total'){
    y0a = y.yr[,-2]
    y.sub =   ddply(y0a, 'site', numcolwise(sum))
  } else {
    y.sub = subset(y.yr, species == species.of.interest)[,-2]
  }
  # Merge with the sampled sites for a given year:
  y.sub2 = merge(y.sub, site.yr, all = T)
  # Turn NA's into 0's for sites sampled with no observations:
  y.sub2[is.na(y.sub2)]<-0
  # Merge with samples frame to create NA's for sites not sampled that year:
  return(merge(sites, y.sub2, all = T))
}

# Create observation matrix of counts (input)

y.Sp = function(species.of.interest){
  y09 = y.SpYr(species.of.interest, 2009)
  y10 = y.SpYr(species.of.interest, 2010)
  y12 = y.SpYr(species.of.interest, 2012)
  df = cbind(y09,y10[,-1],y12[,-1])
  row.names(df) = df[,1]
  df = df[,-1]
  return(as.matrix(df))
}

#--------------------------------------------------------------------------------*
# ---- Prep site covariates frame ----
#--------------------------------------------------------------------------------*

site_covsFrame = function(){
  data.frame(can = scale(site_covs$can), 
  imp = scale(site_covs$imp), 
  row.names = site_covs$site)
}

#--------------------------------------------------------------------------------*
# ---- Prep yearly site covariates frame ----
#--------------------------------------------------------------------------------*

# Function to convert a obs_cov to a scaled variable in wide format:

wideObs = function(variable){
  # Subset frame to variable of interest and grouping variables:
  df1 = data.frame(site = obs_covs$site, year = obs_covs$year, 
                   var = obs_covs[,variable])
  # Scale variable if numeric or integer:
  if (variable == 'day'|variable == 'temp') {
    df1$var = scale(df1$var)
  }
  # Reshape to wide format, with year as columns:
  melt.df = melt(df1, id.vars = c('site','year'), measure.vars = 'var')
  cast.df = cast(melt.df, site ~ year)
  cast.df = cast.df[,-1]
  return(as.data.frame(cast.df))
}

yearlyCovFun = function(){
  obsCov1 = c('day','observer','temp','sky')
  out.list = list()
  for (i in 1:length(obsCov1)){
    out.list[[i]] = wideObs(obsCov1[i])
  }
  names(out.list) = obsCov1
  return(out.list)
  }

t = yearlyCovFun()

t2  = t[[4]]

t2[is.na(t2)]<-NA


#--------------------------------------------------------------------------------*
# ---- CREATE UNMARKED FRAME ----
#--------------------------------------------------------------------------------*

umfMaker = function(species.of.interest){
  y1 = y.Sp(species.of.interest)
  s_covs = site_covsFrame()
  y_covs = yearlyCovFun()
  distances = seq(0,50,10)
  # Testing if NA's are the problem:
  y1 = na.omit(y1)
  y_covs = na.omit(y_covs)
  s_covs = merge(y0, s_covs, by.x = row.names(y1), by.y = row.names(s_covs), all = F)
  unmarkedFrameGDS(y=y0, siteCovs=s_covs, yearlySiteCovs=y_covs, numPrimary=3, 
                   dist.breaks=distances, survey="point", unitsIn="m")
  }

umf.cach = umfMaker('cach')
t = na.omit(umf.cach)

#--------------------------------------------------------------------------------*
# ---- CALCULATE TOTAL SPECIES ABUNDANCE BY IMPERVIOUS AND CANOPY COVER ----
#================================================================================*

umf.cach = umfMaker('cach')

mod$Null <- gdistsamp(lambdaformula=~1, phiformula=~1, pformula=~1, 
                       umf.cach, output="abund", K=20)

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







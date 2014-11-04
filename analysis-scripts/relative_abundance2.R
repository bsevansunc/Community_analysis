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

# Get data:

setwd('~/Community_analysis')

pc = read.csv('derived-data/pc_10_14.csv')
g = read.csv('derived-data/guilds.csv')
lc = read.csv('derived-data/pts_lc100.csv')

# Remove Bert Drake (all zero counts?):

pc = pc[pc$site!='DRAKBERMD1',]
lc = lc[lc$site!='DRAKBERMD1',]

# Convert guild species to lower case

g$species = tolower(g$species)

# Combine foraging and trophic guilds:

g$ft = paste(g$foraging, g$trophic, sep = '_')

# Reduce site and count frames to just the sites with point counts:

lc = merge(lc, data.frame(site = unique(pc$site)), all = F)

pc = merge(pc, data.frame(site = unique(lc$site)), all = F)

# For each site, get the number of surveys and add to the lc frame:

visits = ddply(pc, .(site), summarize, length(unique(year)))[,2]
visits = as.matrix(visits, ncol = 1)

#--------------------------------------------------------------------------------*
# ---- CALCULATE TOTAL SPECIES ABUNDANCE BY IMPERVIOUS AND CANOPY COVER ----
#================================================================================*

# Make a frame summarizing distance by site, year, and species (sum of counts
# across removal time intervals):

pc.sa = ddply(pc, .(site, year, species), summarize, d10 = sum(d10))
pc.sa$d20 = ddply(pc, .(site, year, species), summarize, sum(d20))[,4]
pc.sa$d30 = ddply(pc, .(site, year, species), summarize, sum(d30))[,4]
pc.sa$d40 = ddply(pc, .(site, year, species), summarize, sum(d40))[,4]
pc.sa$d50 = ddply(pc, .(site, year, species), summarize, sum(d50))[,4]

guilds = merge(pc.sa, g, all = F)

guilds.possible.ft = sort(unique(guilds$ft))
guilds.possible.nest = sort(unique(guilds$nest))

# Function to create unmarked frame for a given guild:

guild.umf.fun = function(guild.class, guild){
  # Sum the counts for each guild, year, and distance class:
    guilds$gc = guilds[,guild.class]
    pc.g = ddply(guilds, .(site, year, gc), summarize, d10 = sum(d10))
    pc.g$d20 = ddply(guilds, .(site, year, gc), summarize, sum(d20))[,4]
    pc.g$d30 = ddply(guilds, .(site, year, gc), summarize, sum(d30))[,4]
    pc.g$d40 = ddply(guilds, .(site, year, gc), summarize, sum(d40))[,4]
    pc.g$d50 = ddply(guilds, .(site, year, gc), summarize, sum(d50))[,4]
  # Calculate the max count for each guild across years at a site:
    pc.g2 = ddply(pc.g, .(site, gc), summarize, d10 = max(d10))
    pc.g2$d20 = ddply(pc.g, .(site, gc), summarize, max(d20))[,3]
    pc.g2$d30 = ddply(pc.g, .(site, gc), summarize, max(d30))[,3]
    pc.g2$d40 = ddply(pc.g, .(site, gc), summarize, max(d40))[,3]
    pc.g2$d50 = ddply(pc.g, .(site, gc), summarize, max(d50))[,3] 
  # Extract the guild of interest:
    goi = guild
  # Subset the dataframe to just that guild:
    pc.g3 = pc.g2[pc.g2[,'gc'] == goi,-2]
  # Create a dataframe of just sites:
    site.df = data.frame(site = sort(unique(pc$site)))
  # Merge site and guild frames:
    pc.g4 = merge(pc.g3, site.df, all = T)
  # Change NA's to 0's
    pc.g4[is.na(pc.g4)] = 0
  # Set sites as rownames:
    row.names(pc.g4) = pc.g4[,1]
  # Remove sites from the count frame:
    pc.g4 = pc.g4[,-1]
  # Name the columns:
    names(pc.ta) = c('[0,10]', '(10,20]', '(20,30]', '(30,40]', '(40,50]')
  # Create scaled covariate dataframe:
    covs = data.frame(can = scale(lc$can), imp = scale(lc$imp),
                      row.names = lc$site)
  # Create unmarked frame:
    umf = unmarkedFrameDS(y=as.matrix(pc.g4), 
            siteCovs=covs, survey="point",
            dist.breaks=c(0, 10, 20, 30, 40, 50), 
            unitsIn="m")
    return(umf)
  }

# Create lists of unmarked frames across guilds:

umf.ft = list()
umf.nest = list()

  for (i in 1:length(guilds.possible.ft)){
    umf.ft[[i]] = guild.umf.fun('ft',guilds.possible.ft[i])
  }

  for (i in 1:length(guilds.possible.nest)){
    umf.nest[[i]] = guild.umf.fun('nest',guilds.possible.nest[i])
  }

names(umf.ft) = guilds.possible.ft
names(umf.nest) = guilds.possible.nest

#---------------------------------------------------------------------*
# ---- RUN MODELS ----
#---------------------------------------------------------------------*

run.guild.mods = function(guild.umf){
# Vector of candidate model formulas:
  models = c('~1 ~1',
             '~1 ~imp',
             '~1 ~imp + can',
             '~1 ~imp + I(imp^2)',
             '~1 ~imp + I(imp^2) + can',
             '~1 ~imp + I(imp^2) + can + I(can^2)',
             '~1 ~can',
             '~1 ~can + I(can^2)',
             '~1 ~imp + can + I(can^2)',
             '~1 ~imp + can + imp:can',
             '~1 ~imp + can + imp:can + I(imp^2)',
             '~1 ~imp + can + imp:can +I(can^2)',
             '~1 ~imp + can + imp:can +I(imp^2)+ I(can^2)')
  # Run models:
    umf = guild.umf
    mod.outs = list()
    for(j in 1:length(models)){
      mod.outs[[j]] = distsamp(formula(models[j]),umf, output = 'density')
    }
    # Create fitlist
      fl = fitList(fits = mod.outs)
    # Create model selection tables
      aic.tab = aictab(mod.outs, modnames = models)
      mod.sel.tab = modSel(fl, nullmod = '1')
    # Predict models at sites:
      mod.preds.site = predict(fl, type="state")
    # Backtransform function to turn z values back to original scale:
      bt.cov = function(x, x.scaled) x.scaled * sd(x) + mean(x)
    # Convert can and imp to original scales and add to pred frame:
      mod.preds.site$can = bt.cov(lc$can, covs$can)
      mod.preds.site$imp = bt.cov(lc$imp, covs$imp)
    # Return list of models, fitlist, aic.tab, model selection tab, 
    # and model predictions:
      list.out = list(models, mod.outs, fl, aic.tab, mod.sel.tab, mod.preds.site)
      names(list.out) = c('models','model.outs','fitList','AICtable',
                     'ModelSelectionTable','Model.prediction.by.site')
      return(list.out)
  }

guild.mods.troph = list()
guild.mods.nest = list()
for (i in 1:length(umf.ft)){
  guild.mods.troph[[i]] = run.guild.mods(umf.ft[[i]])
}
for (i in 1:length(umf.nest)){
  guild.mods.nest[[i]] = run.guild.mods(umf.nest[[i]])
}

test = run.guild.mods(umf.nest[[6]])

#---------------------------------------------------------------------*
# ---- CREATE DATA FRAME OF PREDICTIONS BY GUILD AND SITE ----
#---------------------------------------------------------------------*

troph.df = data.frame(guild.mods.troph[[1]][[6]])
troph.df = troph.df[,c(5:6,1)]
for(i in 2:length(guild.mods.troph)){
  troph.df[,i+2] = guild.mods.troph[[i]][[6]][,1]
}
names(troph.df)[3:15] = guilds.possible.ft[1:length(guilds.possible.ft)]

nest.df = data.frame(guild.mods.nest[[1]][[6]])
nest.df = nest.df[,c(5:6,1)]
for(i in 2:length(guild.mods.nest)){
  nest.df[,i+2] = guild.mods.nest[[i]][[6]][,1]
}
names(nest.df)[3:8] = c(as.character(guilds.possible.nest[1:6]))
names(nest.df)[4:6] = c('cup_gen','cup_shrub','cup_tree')

# Create proportions:

troph.prop.df = troph.df
for (i in 3:length(troph.df)){
  troph.prop.df[,i] = troph.df[,i]/rowSums(troph.df[,3:15])
}

nest.prop.df = nest.df
for (i in 3:length(nest.df)){
  nest.prop.df[,i] = nest.df[,i]/rowSums(nest.df[,3:8])
}

# Create lcl and ucl frames:

troph.lcl.df = data.frame(guild.mods.troph[[1]][[6]])
troph.lcl.df = troph.lcl.df[,c(5:6,3)]
for(i in 2:length(guild.mods.troph)){
  troph.lcl.df[,i+2] = guild.mods.troph[[i]][[6]][,3]
}
names(troph.lcl.df)[3:15] = guilds.possible.ft[1:length(guilds.possible.ft)]

troph.ucl.df = data.frame(guild.mods.troph[[1]][[6]])
troph.ucl.df = troph.ucl.df[,c(5:6,4)]
for(i in 2:length(guild.mods.troph)){
  troph.ucl.df[,i+2] = guild.mods.troph[[i]][[6]][,4]
}
names(troph.ucl.df)[3:15] = guilds.possible.ft[1:length(guilds.possible.ft)]

nest.lcl.df = data.frame(guild.mods.nest[[1]][[6]])
nest.lcl.df = nest.lcl.df[,c(5:6,3)]
for(i in 2:length(guild.mods.nest)){
  nest.lcl.df[,i+2] = guild.mods.nest[[i]][[6]][,3]
}
names(nest.lcl.df)[3:8] = c(as.character(guilds.possible.nest[1:6]))
names(nest.lcl.df)[4:6] = c('cup_gen','cup_shrub','cup_tree')

nest.ucl.df = data.frame(guild.mods.nest[[1]][[6]])
nest.ucl.df = nest.ucl.df[,c(5:6,4)]
for(i in 2:length(guild.mods.nest)){
  nest.ucl.df[,i+2] = guild.mods.nest[[i]][[6]][,4]
}
names(nest.ucl.df)[3:8] = c(as.character(guilds.possible.nest[1:6]))
names(nest.ucl.df)[4:6] = c('cup_gen','cup_shrub','cup_tree')


# Function to extract lcl and ucl from frame:

prop.fun.cl = function(lcl.df, ucl.df,guild){
  high = ucl.df[,guild]
  low = lcl.df[,guild]  
  lcl1 = lcl.df[ , -which(names(lcl.df) %in% c(guild))]
  ucl1 = lcl.df[ , -which(names(ucl.df) %in% c(guild))]
  lcl = low/rowSums(ucl1[3:length(ucl1)])
  ucl = high/rowSums(lcl1[3:length(lcl1)])
  data.frame(can = lcl.df[,1], imp = lcl.df[,2], lcl = lcl, ucl = ucl)
}

###########

# Graphical output (function):

scatterout = function(df, guild){
  df1 = data.frame(df$imp, df$can, df[,guild])
  names(df1) = c('imp','can','z')
  p = ggplot(df1, aes(imp, can))
  p + geom_point(aes(color = z, size = z))+
    scale_color_gradient(low = 'blue', high = 'red')+
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
}

scatterout.lc(nest.prop.df,'cavity','imp')

scatterout.lc = function(df.in, response, predictor){
  # Make data frame:
  df = data.frame(df.in[,response], pc.abund[,predictor])
#   df.cl = prop.fun.cl(lcl.df, ucl.df,response)
#   df.ucl = data.frame(ucl = df.cl$ucl, predictor = df.cl[,predictor])
#   df.lcl = data.frame(lcl = df.cl$lcl, predictor = df.cl[,predictor])
  names(df) = c('response','predictor')
  # Get modeled data and construct grid to plot:
  model = lm(response~predictor+I(predictor^2), data = df)
#   model.ucl = lm(ucl~predictor+I(predictor^2), data = df.ucl)
#   model.lcl = lm(lcl~predictor+I(predictor^2), data = df.lcl)
  grid <- with(df, expand.grid(predictor = seq(min(predictor), max(predictor), length = 50)))
  grid$response <- stats::predict(model, newdata=grid)
#    grid$ucl <-stats::predict(model.ucl, newdata=grid)
#    grid$lcl <-stats::predict(model.lcl, newdata=grid)
  # Plot label for x-axis:
  xlabel = if (predictor == 'imp') '% Impervious' else '% Canopy'
  # Plot data:
  ggplot(df, aes(predictor, response)) + geom_point(color = 1, size = 1) + 
#     geom_smooth(aes(ymin = lcl, ymax = ucl), data = grid, stat='identity', col = 1) +
    #geom_line(data = grid, stat = 'identity', col = 1, size = 1.5)+
    ylab('Proportional abundance') + 
    xlab(xlabel) + 
    # Add themes:
    theme(axis.text = element_text(size=10, color = 1),
          axis.title.x = element_text(vjust = -.5),
          axis.title.y = element_text(vjust = 1),
          axis.title = element_text(size=10, vjust = -1),
          axis.line = element_line(colour = "black"),
          legend.position="none",
          panel.background = element_blank())
}

# Write a function that combines the above into a single plot:

plot.outs = function(df.in, response){
  plot.list = list(scatterout(df.in, response), 
                   scatterout.lc(df.in,response,'imp'), 
                   scatterout.lc(df.in,response,'can'))
  plot.layout = matrix(c(1,1,3,1,1,2), nrow = 2, byrow = T)
  multiplot(plotlist = plot.list, layout = plot.layout)
}

plot.outs(nest.prop.df, 'cavity')

plot.fun = function(df.in, guild.name){
  out_dir = 'output/plots/relative_abundance/'
  out_name = paste('rel_abund_',names(df.in[i]),'.pdf', sep ='')
  out = paste(out_dir, out_name, sep = '')
  pdf(out, width = 7.68, height = 4.8)
  plot.outs(df.in, guild.name)
  dev.off()
}

# Plot output, nest guilds:

for (i in 3:length(names(nest.prop.df))){
  plot.fun(nest.prop.df, names(nest.prop.df[i]))
}

# Plot output, feeding-trophic guilds:

names(troph.prop.df)[c(9,12)] = c('foliage_insct_om','ground_insct_om')

for (i in 3:length(names(troph.prop.df))){
  plot.fun(troph.prop.df, names(troph.prop.df[i]))
}


# guild relative abundance

#--------------------------------------------------------------------------------*
# ---- SETUP ----
#================================================================================*

# Load packages:

library(ggplot2)
library(car)

# Run source script to prepare data for analysis:

source('analysis-scripts/pc_prep2.R')
source('analysis-scripts/multiplot_function.R')

#--------------------------------------------------------------------------------*
# ---- FUNCTIONS FOR ANALYSES ----
#================================================================================*

# Function to scale proportional value to (0,1):

props.scaled = function(df, response){
  x = df[,response]
  n = length(x)
  s = .5
  (x*(n-1)+s)/n
}

# Create list of potential models:

models = list('response~imp + can + imp:can + I(imp^2) + I(can^2)',
              'response~imp + can + imp:can + I(imp^2)',
              'response~imp + can + imp:can + I(can^2)',
              'response~imp + can + imp:can',
              'response~imp + can + I(imp^2) + I(can^2)',
              'response~imp + can + I(can^2)',
              'response~can + I(can^2)',
              'response~can + imp',
              'response~can',
              'response~imp + can + I(imp^2)',
              'response~imp + I(imp^2)',
              'response~imp',
              'response~1')  

# Run model list function:

run.mods = function(df, response){
#   response = cbind(df[,response], rowSums(df) - df[,response])
  response = props.scaled(df, response)
  df = cbind(response, lc[-c(1:3)])
#   response = logit(df[,response])
  out.list = list()
  for(i in 1:length(models)){
#     out.list[[i]] = glm(as.formula(models[[i]]), data = df, family = binomial(link = 'logit'))
#     out.list[[i]] = glm(as.formula(models[[i]]), data = pc.abund)
    out.list[[i]] = betareg(as.formula(models[[i]]), data = df)
  }
  names(out.list) = models
  return(out.list)
}

#--------------------------------------------------------------------------------*
# Nest relative abundance
#--------------------------------------------------------------------------------*

cav.edif.mods = run.mods(nest.ra,'cavity.edificarian')
cav.tree.mods = run.mods(nest.ra,'cavity.tree')
cup.ground.mods = run.mods(nest.ra,'cup.ground')
cup.shrub.mods = run.mods(nest.ra,'cup.shrub')
cup.tree.mods = run.mods(nest.ra,'cup.tree')

nest.mod.list = list(cav.edif.mods, cav.tree.mods, cup.ground.mods, cup.shrub.mods, cup.tree.mods)
names(nest.mod.list) = c('cavity.edificarian','cavity.tree', 'cup.ground','cup.shrub', 'cup.tree')

aic.nest.mods = list()
for (i in 1:length(nest.mod.list)){
  aic.nest.mods[[i]] = aictab(nest.mod.list[[i]], modnames = names(nest.mod.list[[i]]))
}

names(aic.nest.mods) = names(nest.mod.list)

aic.nest.mods

coef.mods = function(guild.mod.list){
  modavg(guild.mod.list, 'imp', exclude = c('imp:can','I(imp^2)'))->modavg.imp
    mod.avg.imp = data.frame(modavg.imp[c(1, 3:7)])
  modavg(guild.mod.list, 'can', exclude = c('imp:can','I(can^2)'))->modavg.can
    mod.avg.can = data.frame(modavg.can[c(1, 3:7)])
  modavg(guild.mod.list, 'imp:can')->modavg.imp.can
    mod.avg.impcan = data.frame(modavg.imp.can[c(1, 3:7)])
  modavg(guild.mod.list, 'I(imp^2)')->modavg.imp2
    mod.avg.imp2 = data.frame(modavg.imp2[c(1, 3:7)])
  modavg(guild.mod.list, 'I(can^2)')->modavg.can2
    mod.avg.can2 = data.frame(modavg.can2[c(1, 3:7)])
  out.df = rbind(mod.avg.imp, mod.avg.can, mod.avg.impcan, mod.avg.imp2, mod.avg.can2)
  return(out.df)
  }

coef.nest.mods = list()
for (i in 1:length(nest.mod.list)){
  coef.nest.mods[[i]] = coef.mods(nest.mod.list[[i]])
}

names(coef.nest.mods) = names(nest.mod.list)
coef.nest.mods

imp.pred.df = function(guild){
  guild.mod.list = nest.mod.list[[guild]]
  imp.pred = data.frame(predicted = guild.mod.list[[11]]$fitted.values, imp = lc$imp)
  imp.pred = imp.pred[order(imp.pred$imp),]
  return(imp.pred)
}

can.pred.df = function(guild){
  guild.mod.list = nest.mod.list[[guild]]
  can.pred = data.frame(predicted = guild.mod.list[[7]]$fitted.values, can = lc$can)
  can.pred = can.pred[order(can.pred$can),]
  return(can.pred)
}

par(mar = c(5,5.2,2,1))
plot(predicted~imp, data = imp.pred.df('cavity.edificarian'), 
     ylim = c(0,.8), type = 'l', lwd = 2, bty = 'l', 
     xlab = '% Impervious', ylab = 'Relative abundance', cex.lab = 1.2)
  lines(predicted~imp, data = imp.pred.df('cavity.tree'), type = 'l', lwd = 2, col = 'blue')
  lines(predicted~imp, data = imp.pred.df('cup.ground'), type = 'l', lwd = 2, col = 'red')
  lines(predicted~imp, data = imp.pred.df('cup.shrub'), type = 'l', lwd = 2, col = 'orange')
  lines(predicted~imp, data = imp.pred.df('cup.tree'), type = 'l', lwd = 2, col = 'purple')
  legend('topleft', lwd = 2, legend = names(nest.mod.list), col = c(1, 'blue','red','orange','purple'), bty = 'n')


plot(predicted~can, data = can.pred.df('cavity.edificarian'), 
     ylim = c(0,.8), type = 'l', lwd = 2, bty = 'l', 
     xlab = '% Canopy', ylab = 'Relative abundance', cex.lab = 1.2)
lines(predicted~can, data = can.pred.df('cavity.tree'), type = 'l', lwd = 2, col = 'blue')
lines(predicted~can, data = can.pred.df('cup.ground'), type = 'l', lwd = 2, col = 'red')
lines(predicted~can, data = can.pred.df('cup.shrub'), type = 'l', lwd = 2, col = 'orange')
lines(predicted~can, data = can.pred.df('cup.tree'), type = 'l', lwd = 2, col = 'purple')
legend('topleft', lwd = 2, legend = names(nest.mod.list), col = c(1, 'blue','red','orange','purple'), bty = 'n')

#--------------------------------------------------------------------------------*
# Trophic relative abundance
#--------------------------------------------------------------------------------*

head(trophic.ra)
carnivore.mods = run.mods(trophic.ra,'carnivore')
frugivore.mods = run.mods(trophic.ra,'frugivore')
granivore.mods = run.mods(trophic.ra,'granivore')
insectivore.mods = run.mods(trophic.ra,'insectivore')
omnivore.mods = run.mods(trophic.ra,'omnivore')

trophic.mod.list = list(carnivore.mods, frugivore.mods, granivore.mods, insectivore.mods, omnivore.mods)
names(trophic.mod.list) = c('carnivore','frugivore','granivore', 'insectivore','omnivore')

aic.trophic.mods = list()
for (i in 1:length(trophic.mod.list)){
  aic.trophic.mods[[i]] = aictab(trophic.mod.list[[i]], modnames = names(trophic.mod.list[[i]]))
}

names(aic.trophic.mods) = names(trophic.mod.list)

aic.trophic.mods

coef.mods = function(guild.mod.list){
  modavg(guild.mod.list, 'imp', exclude = c('imp:can','I(imp^2)'))->modavg.imp
  mod.avg.imp = data.frame(modavg.imp[c(1, 3:7)])
  modavg(guild.mod.list, 'can', exclude = c('imp:can','I(can^2)'))->modavg.can
  mod.avg.can = data.frame(modavg.can[c(1, 3:7)])
  modavg(guild.mod.list, 'imp:can')->modavg.imp.can
  mod.avg.impcan = data.frame(modavg.imp.can[c(1, 3:7)])
  modavg(guild.mod.list, 'I(imp^2)')->modavg.imp2
  mod.avg.imp2 = data.frame(modavg.imp2[c(1, 3:7)])
  modavg(guild.mod.list, 'I(can^2)')->modavg.can2
  mod.avg.can2 = data.frame(modavg.can2[c(1, 3:7)])
  out.df = rbind(mod.avg.imp, mod.avg.can, mod.avg.impcan, mod.avg.imp2, mod.avg.can2)
  return(out.df)
}

coef.trophic.mods = list()
for (i in 1:length(trophic.mod.list)){
  coef.trophic.mods[[i]] = coef.mods(trophic.mod.list[[i]])
}

names(coef.trophic.mods) = names(trophic.mod.list)
coef.trophic.mods

imp.pred.df = function(guild){
  guild.mod.list = trophic.mod.list[[guild]]
  imp.pred = data.frame(predicted = guild.mod.list[[11]]$fitted.values, imp = lc$imp)
  imp.pred = imp.pred[order(imp.pred$imp),]
  return(imp.pred)
}

can.pred.df = function(guild){
  guild.mod.list = trophic.mod.list[[guild]]
  can.pred = data.frame(predicted = guild.mod.list[[7]]$fitted.values, can = lc$can)
  can.pred = can.pred[order(can.pred$can),]
  return(can.pred)
}

par(mar = c(5,5.2,2,1))
plot(predicted~imp, data = imp.pred.df('carnivore'), 
     ylim = c(0,.8), type = 'l', lwd = 2, bty = 'l', 
     xlab = '% Impervious', ylab = 'Relative abundance', cex.lab = 1.2)
lines(predicted~imp, data = imp.pred.df('frugivore'), type = 'l', lwd = 2, col = 'blue')
lines(predicted~imp, data = imp.pred.df('granivore'), type = 'l', lwd = 2, col = 'red')
lines(predicted~imp, data = imp.pred.df('insectivore'), type = 'l', lwd = 2, col = 'orange')
lines(predicted~imp, data = imp.pred.df('omnivore'), type = 'l', lwd = 2, col = 'purple')
legend('topleft', lwd = 2, cex = .8,legend = names(trophic.mod.list), col = c(1, 'blue','red','orange','purple'), bty = 'n')


plot(predicted~can, data = can.pred.df('carnivore'), 
     ylim = c(0,.8), type = 'l', lwd = 2, bty = 'l', 
     xlab = '% Canopy', ylab = 'Relative abundance', cex.lab = 1.2)
lines(predicted~can, data = can.pred.df('frugivore'), type = 'l', lwd = 2, col = 'blue')
lines(predicted~can, data = can.pred.df('granivore'), type = 'l', lwd = 2, col = 'red')
lines(predicted~can, data = can.pred.df('insectivore'), type = 'l', lwd = 2, col = 'orange')
lines(predicted~can, data = can.pred.df('omnivore'), type = 'l', lwd = 2, col = 'purple')
legend('topleft', lwd = 2, cex = .8, legend = names(trophic.mod.list), col = c(1, 'blue','red','orange','purple'), bty = 'n')

#--------------------------------------------------------------------------------*
# Foraging-Trophic relative abundance
#--------------------------------------------------------------------------------*

head(foraging_trophic.ra)
aerial_insectivore.mods = run.mods(foraging_trophic.ra,'aerial_insectivore')
foliage_insectivore.mods = run.mods(foraging_trophic.ra,'foliage_insectivore')
ground_insectivore.mods = run.mods(foraging_trophic.ra,'ground_insectivore')
bark_insectivore.mods = run.mods(foraging_trophic.ra,'bark_insectivore')

trophic.mod.list = list(aerial_insectivore.mods, foliage_insectivore.mods, ground_insectivore.mods, bark_insectivore.mods)
names(trophic.mod.list) = c('aerial_insectivore','foliage_insectivore','ground_insectivore', 'bark_insectivore')

aic.trophic.mods = list()
for (i in 1:length(trophic.mod.list)){
  aic.trophic.mods[[i]] = aictab(trophic.mod.list[[i]], modnames = names(trophic.mod.list[[i]]))
}

names(aic.trophic.mods) = names(trophic.mod.list)

aic.trophic.mods

coef.mods = function(guild.mod.list){
  modavg(guild.mod.list, 'imp', exclude = c('imp:can','I(imp^2)'))->modavg.imp
  mod.avg.imp = data.frame(modavg.imp[c(1, 3:7)])
  modavg(guild.mod.list, 'can', exclude = c('imp:can','I(can^2)'))->modavg.can
  mod.avg.can = data.frame(modavg.can[c(1, 3:7)])
  modavg(guild.mod.list, 'imp:can')->modavg.imp.can
  mod.avg.impcan = data.frame(modavg.imp.can[c(1, 3:7)])
  modavg(guild.mod.list, 'I(imp^2)')->modavg.imp2
  mod.avg.imp2 = data.frame(modavg.imp2[c(1, 3:7)])
  modavg(guild.mod.list, 'I(can^2)')->modavg.can2
  mod.avg.can2 = data.frame(modavg.can2[c(1, 3:7)])
  out.df = rbind(mod.avg.imp, mod.avg.can, mod.avg.impcan, mod.avg.imp2, mod.avg.can2)
  return(out.df)
}

coef.trophic.mods = list()
for (i in 1:length(trophic.mod.list)){
  coef.trophic.mods[[i]] = coef.mods(trophic.mod.list[[i]])
}

names(coef.trophic.mods) = names(trophic.mod.list)
coef.trophic.mods

imp.pred.df = function(guild){
  guild.mod.list = trophic.mod.list[[guild]]
  imp.pred = data.frame(predicted = guild.mod.list[[11]]$fitted.values, imp = lc$imp)
  imp.pred = imp.pred[order(imp.pred$imp),]
  return(imp.pred)
}

can.pred.df = function(guild){
  guild.mod.list = trophic.mod.list[[guild]]
  can.pred = data.frame(predicted = guild.mod.list[[7]]$fitted.values, can = lc$can)
  can.pred = can.pred[order(can.pred$can),]
  return(can.pred)
}

par(mar = c(5,5.2,2,1))
plot(predicted~imp, data = imp.pred.df('aerial_insectivore'), 
     ylim = c(0,.5), type = 'l', lwd = 2, bty = 'l', 
     xlab = '% Impervious', ylab = 'Relative abundance', cex.lab = 1.2)
lines(predicted~imp, data = imp.pred.df('foliage_insectivore'), type = 'l', lwd = 2, col = 'blue')
lines(predicted~imp, data = imp.pred.df('ground_insectivore'), type = 'l', lwd = 2, col = 'red')
lines(predicted~imp, data = imp.pred.df('bark_insectivore'), type = 'l', lwd = 2, col = 'orange')
legend('topleft', lwd = 2, cex = .8, legend = names(trophic.mod.list), col = c(1, 'blue','red','orange'), bty = 'n')


par(mar = c(5,5.2,2,1))
plot(predicted~can, data = can.pred.df('aerial_insectivore'), 
     ylim = c(0,.5), type = 'l', lwd = 2, bty = 'l', 
     xlab = '% Canopy', ylab = 'Relative abundance', cex.lab = 1.2)
lines(predicted~can, data = can.pred.df('foliage_insectivore'), type = 'l', lwd = 2, col = 'blue')
lines(predicted~can, data = can.pred.df('ground_insectivore'), type = 'l', lwd = 2, col = 'red')
lines(predicted~can, data = can.pred.df('bark_insectivore'), type = 'l', lwd = 2, col = 'orange')
legend('topleft', lwd = 2, cex = .8, legend = names(trophic.mod.list), col = c(1, 'blue','red','orange'), bty = 'n')

#####################################################################################################################################
# Adjusted D2 function (roughly equivalent to R2, see Guisan
# & Zimmermann 2000):

####################
# WORKS TO HERE
####################

d2adj = function(df, mod.list, i){
  mod = mod.list[[i]]
  d2 = (mod$null.deviance - mod$deviance)/mod$null.deviance
  k = length(mod$coefficients)
  n = length(df[,1])
  D2adj = 1 - ((n-1)/ (n-k))*(1-d2)
  if(k == 1) D2adj = d2
  D2adj = round(D2adj,3)
  return(D2adj)
  }

# Model out function:
# Note: input is the model list function above, output is a
# list of model outputs [[1]] where the length that of the list
# of models and a model selection table [[2]]

mod.outs = function(df, response){
  # Run models across list:
    mod.list = run.mods(df, response)
  # Set outputs to extract from model list:
    aic_out = numeric()
    formulas = character()
    k = numeric()
    deviance = numeric()
    D2adj = numeric()
  # Extract outputs from each model:
  for(i in 1:length(mod.list)){
    formulas[i] = as.character(mod.list[[i]]$formula[3])
    aic_out[i] = AIC(mod.list[[i]])
    k[i] = length(mod.list[[i]]$coefficients)
    deviance[i] = mod.list[[i]]$deviance
    D2adj[i] = d2adj(df, mod.list, response, i)
  }
  # Create model table:
  out.df = data.frame(formulas, k, aic_out, deviance, D2adj)
  # Sort by AIC
  out.df = out.df[order(out.df$aic_out),]
  # Add delta AIC and model weights:
  out.df$dAIC = out.df$aic_out-min(out.df$aic_out)
  out.df$w = round(exp(-0.5*out.df$dAIC)/sum(exp(-0.5*out.df$dAIC)),3)
  # Output table, with columns sorted:
  out.df = out.df[,c(1:3, 6:7, 4, 5)]
  # Output list:
  out.list = list(mod.list, out.df)
  return(out.list)
}

# Function extract the coefficient and standard error from a
# model:

get.coef.se = function(predictor, i, mod.list){
  # Summary table for a given model:
    sum.table = data.frame(summary(mod.list[[i]])[12])
  # Extract to predictor variable
    sum.table = sum.table[row.names(sum.table) == predictor,]
  # If the predictor variable isn't present, replace with NA
    if (length(sum.table[,1]) == 0) sum.table[1,] <- NA
  # Return estimate and standard error
    out.df = sum.table[,1:2]
    names(out.df) = c('estimate','se')
    return(out.df)
    }

# Function to calculate model-averaged coefficients (and standard
# errors) for a given predictor variable:

ma.coefs.se = function(predictor, mod.list, m.table){  
  # For loop to extract model averaged coefficients and se
  ma.beta = numeric()
  ma.se = numeric()
  w = numeric()
  for (i in 1:length(mod.list)){
    w[i] = m.table[row.names(m.table) == i,'w']
    ma.beta[i] = get.coef.se(predictor, i, mod.list)[,1]*w[i]
    ma.se[i] = get.coef.se(predictor, i, mod.list)[,2]*w[i]
  }
  # Remove NA's and returned sum coefficient and se values:
  ma.beta = na.omit(ma.beta)
  ma.se = na.omit(ma.se)
  beta.est = sum(ma.beta)
  beta.se = sum(ma.se)
  # Return data frame of estimate and se:
  out.df = data.frame(predictor, beta.est, beta.se)
  return(out.df)
}

# Function to create a model-averaged global model across predictor variables:

ma.pred = function(mod.list,m.table){
  predictors = c('(Intercept)','imp','can','imp:can','I(imp^2)', 'I(can^2)')
  out.list = list()
  for (i in 1:length(predictors)){
    out.list[[i]] = ma.coefs.se(predictors[i], mod.list, m.table)
  }
  # Output as data frame:
  do.call('rbind', out.list)
}

# Wrapper function that outputs the model table and beta estimates:

summary.outs = function(df, response){
  m1 = mod.outs(df, response)
  mod.list = m1[[1]]
  m.table = m1[[2]]
  mod.averages = ma.pred(mod.list, m.table)
  list.out = list(m.table, mod.averages)
  return(list.out)
}

# Table equivalent to Hurlbert-Liang (top-half ... summary stats
# of the best models):

HL.table = function(df, response){
  mt = mod.outs(df, response)[[2]]
  model.rank = 1:5
  D2adj = mt$D2adj[1:5]
  aic = mt$aic_out[1:5]
  dAIC = mt$dAIC[1:5]
  w = mt$w[1:5]
  t1 = data.frame(model.rank, D2adj, aic, dAIC,w)
  t1.df = data.frame(t(t1))
  t1.df
}

##########################################
# ---- WORKING ABOVE THIS POINT ----
##########################################

ma.coefs.se('imp', test[[1]], test[[2]])

ma.pred(test[[1]], m.table = test[[2]])

summary.outs(nest.ra,'cup.tree')

HL.table(nest.ab,'cup.gen.c')


# Graphical output (function):

scatterout(nest.ab,'cup.gen.c')

scatterout = function(df, response){
  df1 = data.frame(pc.abund$imp, pc.abund$can, df[,response])
  names(df1) = c('imp','can','z')
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
}

scatterout(nest.ab,'cavity.c')

scatterout.lc(nest.ab,'cavity.c', 'imp')


scatterout.lc = function(df, response, lc.in){
  # Make data frame:
  df = data.frame(df[,response], df$t, pc.abund[,lc.in])
  names(df) = c('response','t', 'predictor')
  # Get modelled data and construct grid to plot:
  model = glm(cbind(response, t-response)~predictor, data = df, family = binomial(link = 'logit'))
  grid <- with(df, expand.grid(predictor = seq(min(predictor), max(predictor), length = 50)))
  grid$response <- stats::predict(model, newdata=grid)
  err <- stats::predict(model, newdata=grid, se = TRUE)
  grid$ucl <- err$fit + 1.96 * err$se.fit
  grid$lcl <- err$fit - 1.96 * err$se.fit
  # Plot label for x-axis:
  xlabel = if (lc.in == 'imp') '% Impervious' else '% Canopy'
  # Plot data:
  ggplot(df, aes(predictor, response)) + geom_point(color = 'darkgray', size = .95) + 
    geom_smooth(aes(ymin = lcl, ymax = ucl), data = grid, stat='identity', col = 1) +
    ylab('Abundance') + 
    xlab(xlabel) + 
    # Add themes:
    theme(axis.text = element_text(size=10, color = 1),
          axis.title.x = element_text(vjust = -.5),
          axis.title.y = element_text(vjust = 1),
          axis.title = element_text(size=12, vjust = -1),
          axis.line = element_line(colour = "black"),
          legend.position="none",
          panel.background = element_blank())
}

# Write a function that combines the above into a single plot:

plot.outs = function(response){
  plot.list = list(scatterout(response), 
                   scatterout.lc(response,'imp'), 
                   scatterout.lc(response,'can'))
  plot.layout = matrix(c(1,1,3,1,1,2), nrow = 2, byrow = T)
  multiplot(plotlist = plot.list, layout = plot.layout)
}

#--------------------------------------------------------------------------------*
# ---- BY   trophic and nest guilds GUILD ----
#================================================================================*
# Lots of them!

# Wrap through to write output:

out.fun = function(df, i){
  file.name = paste('output/table_outs/abundance2/',names(df[i]),sep = '')
  write.csv(summary.outs(df,names(df)[i])[[1]],
            paste(file.name,'_aic_table.csv',sep =''),
            row.names = F)
  write.csv(summary.outs(df,names(df)[i])[[2]], 
            paste(file.name,'_beta_estimates_table.csv',sep =''),
            row.names = F)
  write.csv(HL.table(df,names(df)[i]),
            paste(file.name,'_HLtable.csv',sep =''),
            row.names = F)
}

# Wrap through save output:

for (i in 2:7){
  out.fun(nest.ab, i)
}

for (i in 2:14){
  out.fun(troph.ab, i)
}

# Plot function across guilds:

plot.fun = function(df, i){
  out_dir = 'output/plots/relative_abundance2/'
  out_name = paste('rel_abund_',names(df[i]),'.pdf', sep ='')
  out = paste(out_dir, out_name, sep = '')
  pdf(out, width = 7.68, height = 4.8)
  plot.outs(df[,i])
  dev.off()
}


# Plot output:

for (i in 5:length(names(df))){
  plot.fun(i)
}






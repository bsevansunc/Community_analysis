#--------------------------------------------------------------------------------*
# ---- SETUP ----
#================================================================================*

# Load packages:

library(ggplot2)

# Run source script to prepare data for analysis:

source('analysis-scripts/pc_prep2.R')

#--------------------------------------------------------------------------------*
# ---- FUNCTIONS FOR ANALYSES ----
#================================================================================*

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

run.mods = function(response){
  out.list = list()
  for(i in 1:length(models)){
    out.list[[i]] = glm(as.formula(models[[i]]), data = lc)
    }
  return(out.list)
  }

# Adjusted D2 function (roughly equivalent to R2, see Guisan
# & Zimmermann 2000):

d2adj = function(mod.list, response, i){
  mod = mod.list[[i]]
  d2 = (mod$null.deviance - mod$deviance)/mod$null.deviance
  k = length(mod$coefficients)
  n = length(response)
  D2adj = 1 - ((n-1)/ (n-k))*(1-d2)
  if(k == 1) D2adj = d2
  return(D2adj)
}

# Model out function:
# Note: input is the model list function above, output is a
# list of model outputs [[1]] where the length that of the list
# of models and a model selection table [[2]]

mod.outs = function(response){
  # Run models across list:
    mod.list = run.mods(response)
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
      D2adj[i] = d2adj(mod.list, response, i)
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

get.coef.se = function(predictor, response, i){
  # Summary table for a given model:
    sum.table = data.frame(summary(mod.list[[i]])[12])
  # Extract to predictor variable
    sum.table = sum.table[row.names(sum.table) == predictor,]
  # If the predictor variable isn't present, replace with NA
    if (dim(sum.table)[1] == 0) sum.table[1,] = NA
  # Return estimate and standard error
    out.df = sum.table[,1:2]
      names(out.df) = c('estimate','se')
    return(out.df)
  }

# Function to calculate model-averaged coefficients (and standard
# errors) for a given predictor variable:

ma.coefs.se = function(predictor, response){  
  # For loop to extract model averaged coefficients and se
    ma.beta = numeric()
    ma.se = numeric()
    for (i in 1:length(mod.list)){
      w[i] = m.table[row.names(m.table) == i,'w']
      ma.beta[i] = get.coef.se(predictor, response, i)[,1]*w[i]
      ma.se[i] = get.coef.se(predictor, response, i)[,2]*w[i]
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

ma.pred = function(response){
  predictors = c('(Intercept)','imp','can','imp:can','I(imp^2)', 'I(can^2)')
  out.list = list()
  for (i in 1:length(predictors)){
    out.list[[i]] = ma.coefs.se(predictors[i], response)
  }
  # Output as data frame:
    do.call('rbind', out.list)
}
  
# Wrapper function that outputs the model table and beta estimates:

summary.outs = function(response){
    m1 = mod.outs(response)
    mod.list = m1[[1]]
    m.table = m1[[2]]
    mod.averages = ma.pred(response)
    list.out = list(m.table, mod.averages)
    return(list.out)
}

# Table equivalent to Hurlbert-Liang (top-half ... summary stats
# of the best models):

HL.table = function(response){
  mt = mod.outs(response)[[2]]
  model.rank = 1:5
    model.rank = round(model.rank,0)
  D2adj = mt$D2adj[1:5]
  aic = mt$aic_out[1:5]
  dAIC = mt$dAIC[1:5]
  w = mt$w[1:5]
  t1 = data.frame(model.rank, D2adj, aic, dAIC,w)
  t1.df = data.frame(t(t1))
  t1.df
}

# Graphical output (function):

scatterout = function(response){
  df1 = data.frame(lc$imp, lc$can, response)
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

#--------------------------------------------------------------------------------*
# ---- Species richness ----
#================================================================================*

# Get species richness by site:

  sr = specnumber(pc)

# AIC table (full):

  summary.outs(sr) 

# HL AIC table:

  HL.table(sr)

# Plot output:

  scatterout(sr)

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
    scale_fill_continuous(name = 'Richness',low = "yellow", high = "red", limits = c(0,23))+
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

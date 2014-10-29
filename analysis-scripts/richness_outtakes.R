

# Create a dataframe of expected richness values by potentential combinations of
# impervious and canopy:

expected.frame = function(response,imp, can){
  # Create modeled values of imp and can:
  imp = imp
  can = can
  # Make grid:
  grid = mesh(imp,can)
  # Exclude values for which the percent combination exceeds 100
  g1 = grid$x + grid$y
  g1 = ifelse(g1>100, NA, 1)
  # Project values:
  z = with(grid, expected.values(sr,imp,can))
  # Limit to potential combinations of imp and can:
  z = z*g1
  # Create a dataframe for output
  out.df = melt(z)
  out.df = na.omit(out.df)
  names(out.df) = c('imp','can','richness')
  return(out.df)
}

expected.frame(sr, lc$imp,lc$can)

R2 <- cor(sr,e1)^2

plot(lc$imp, sr)
plot(lc$can, sr)

plot(lc$imp, e1)
plot(lc$can, e1)

# Function that takes the model average frame, imp and canopy and
# calculates the expected values (global model):

expected.values = function(response, imp, can){
  mod = summary.outs(response)[[2]][,2]
  mod[1]+mod[2]*imp+mod[3]*can+mod[4]*imp*can+mod[5]*(imp^2)+mod[6]*(can^2)
}

# Function to determine the adjusted R2 averaged across the model:

adj.r2 = function(response, imp, can){
  aic.table = summary.outs(sr)[[1]]
  observed = response
  expected = expected.values(sr, imp, can)
  r2 = cor(observed, expected)^2
  k.ave = sum(aic.table$k*aic.table$w)
  ((length(response)-1)*r2-k.ave)/(length(response)-1-k.ave)
}

# Calculate variable weights:

var.weights = function(response){
  mt = = mod.outs(response)[[2]]
  vars = c('imp','can','imp:can','I(imp^2)', 'I(can^2)')
  mod.vars = numeric()
  for(i in 1:length(vars)){
    mod.vars  = grep
  }
}

grep('imp',mt[,1])

# Plotting from yesterday

# Construct output table:

mod.tab.sr = mod.table(mod.list.sr)

# Look at results

mod.tab.sr

summary(mod.list[[2]])

summary(lm(sr~imp + can + imp:can + I(imp^2), data = lc))

imp = seq(0,100, 1)
can = seq(0,100, 1)

coef(mod.list.sr[[2]])['imp']*mod.tab.sr[2,'w']





ma.coefs('imp', mod.list.sr, mod.tab.sr)

mod.predictions = predict(mod.i2.inter, type = 'response')

df2 = data.frame(mod.predictions, lc$can, lc$imp)
names(df2) = c('mod','can','imp')

plot(mod.predictions~lc$imp)
plot(mod.predictions~lc$can)

t = expected.frame(sr, 100, 100)

p1 = qplot(x = imp, y = can, data = t, color = richness) 
p1 + scale_colour_gradient(c(0,5), low = 'yellow', high = 'red')+
  geom_point(shape = 1)+
  xlab('% Impervious') + ylab('% Canopy') + 
  # Add themes:
  theme(axis.text = element_text(size=14, color = 1),
        axis.title.x = element_text(vjust = -.5),
        axis.title.y = element_text(vjust = .5),
        axis.title = element_text(size=18, vjust = -1),
        axis.line = element_line(colour = "black"),
        panel.background = element_blank())

scale_colour_gradient(na.value = 'gray')

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

# Old guild analyses ... there may be a useful plotting function in here.

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

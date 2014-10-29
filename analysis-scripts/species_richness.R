#--------------------------------------------------------------------------------*
# ---- SETUP ----
#================================================================================*

# Load packages:

library(vegan)
library(reshape)

# Run source script to prepare data for analysis:

source('analysis-scripts/pc_prep2.R')


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

df2 = data.frame(mod.predictions, lc$can, lc$imp)
  names(df2) = c('mod','can','imp')

plot(mod.predictions~lc$imp)
plot(mod.predictions~lc$can)

p1 = qplot(x = imp, y = can, data = df2, color = mod) 
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

#================================================================================*
# ---- SET-UP ----
#================================================================================*

# Load packages:

library(vegan)
library(lme4)
library(reshape2)
library(AICcmodavg)

# Get data:

h <- read.csv('Jackson community mixed model/MLM_dataset.csv')

#--------------------------------------------------------------------------------*
# ---- Format data to run  (Convert from wide to long) ----
#================================================================================*
# This section was DEFINITELY NOT PARSIMONIOUS. I used Wickhams rehape2 package.

# Remove sites with no species (note: made the formula generalizable):

h <- h[rowSums(h[,26:39]) != 0,]

# Make sepaarate environment and sp frames:

env = h[,1:25]
sp = h[,c(1,26:dim(h)[2])]

# Put in long format and merge frames:

h.long = melt(sp, id = 'Plot', variable_name = 'sp', value.name = 'PRESENCE')
  names(h.long)[2] = 'SPP'

hh = merge(env, h.long)  

# Normalize predictor variables:

for(i in 5:25){hh[,i] <- (hh[,i] - mean(hh[,i]))/sd(hh[,i])} 

# Add quadratic forms of predictor variables:

hh$Elevation2 <- hh$Elevation^2
hh$LITU2 <- hh$LITU^2
hh$Ca2 <- hh$Ca^2
hh$P2 <- hh$P^2
hh$Herb2=hh$Herb^2

#--------------------------------------------------------------------------------*
# ---- Model selection ----
#================================================================================*
#AIC Model Selection (try all possible models)

# They didn't really show the trial of "all possible models". Against, what I 
# consider good model-building behavior, I will use the dredge function in the 
# MuMIn package to run through all potential models (ugh).

options(na.action = 'na.fail')
global.mod <- glmer(PRESENCE~ (1|SPP)+Elevation+Elevation2+(0+Elevation|SPP)+
                    Herb+Herb2+(0+Herb|SPP)+TreeBA+(0+TreeBA|SPP)+ACSA
                    +(0+ACSA|SPP)+LITU+(0+LITU|SPP)+Ca+Ca2+(0+Ca|SPP)+P+(0+P|SPP)
                    +P2, family=binomial, data=hh)

require(MuMIn)

dredge(global.mod) 

#AIC = 849.6 (this is the best model)

best <- glmer(PRESENCE~ (1|SPP)+Elevation+Elevation2+(0+Elevation|SPP)+Herb+Herb2+
                (0+Herb|SPP)+TreeBA+ACSA+LITU+(0+LITU|SPP)+Ca+Ca2+(0+Ca|SPP)+
                P+(0+P|SPP), family=binomial, data=hh)

summary(best)
ranef(best) 

nsite <- dim(h)[1]
nx <- 10 #number of fixed effects 
nspp <- 14 #number of species

################################################
# analyze data: lmer (MLM)

z=glmer(PRESENCE~(1|SPP)+Elevation+Elevation2+(0+Elevation|SPP)+Herb+Herb2+(0+Herb|SPP)+TreeBA+ACSA+LITU+(0+LITU|SPP)+Ca+Ca2+(0+Ca|SPP)+P+(0+P|SPP),family=binomial,data=hh)

# compute ranef pvalues
#Without elevation
z1 <- glmer(PRESENCE~(1|SPP)+Elevation+Elevation2+Herb+Herb2+(0+Herb|SPP)+TreeBA+ACSA+LITU+(0+LITU|SPP)+Ca+Ca2+(0+Ca|SPP)+P+(0+P|SPP),family=binomial,data=hh)

#Without Herb
z2 <- glmer(PRESENCE~(1|SPP)+Elevation+Elevation2+(0+Elevation|SPP)+Herb+Herb2+TreeBA+ACSA+LITU+(0+LITU|SPP)+Ca+Ca2+(0+Ca|SPP)+P+(0+P|SPP),family=binomial,data=hh)

#Without LITU
z3 <- glmer(PRESENCE~(1|SPP)+Elevation+Elevation2+(0+Elevation|SPP)+Herb+Herb2+(0+Herb|SPP)+TreeBA+ACSA+LITU+Ca+Ca2+(0+Ca|SPP)+P+(0+P|SPP),family=binomial,data=hh)

#Without Ca
z4 <- glmer(PRESENCE~(1|SPP)+Elevation+Elevation2+(0+Elevation|SPP)+Herb+Herb2+(0+Herb|SPP)+TreeBA+ACSA+LITU+(0+LITU|SPP)+Ca+Ca2+P+(0+P|SPP),family=binomial,data=hh)

#Without P
z5 <- glmer(PRESENCE~(1|SPP)+Elevation+Elevation2+(0+Elevation|SPP)+Herb+Herb2+(0+Herb|SPP)+TreeBA+ACSA+LITU+(0+LITU|SPP)+Ca+Ca2+(0+Ca|SPP)+P,family=binomial,data=hh)

# I added this function to extract the deviance (LL code, as written (commented out), does not work)

mod.list= list(z1,z2, z3, z4, z5)

dev.mod = function(mod) getME(mod, name = 'devcomp')$cmp[['dev']]

# LL <- c(attributes(z)$deviance[1], attributes(z1)$deviance[1], attributes(z2)$deviance[1], attributes(z3)$deviance[1], attributes(z4)$deviance[1], attributes(z5)$deviance[1])

fun1 = function(x) 1-pchisq(dev.mod(x) - dev.mod(z))

LL = c(dev.mod(z), dev.mod(z1), dev.mod(z2), dev.mod(z3), dev.mod(z4), dev.mod(z5))

mlm.pvals <- c(1-pchisq(LL[2] - LL[1],1), 1-pchisq(LL[3] - LL[1],1), 1-pchisq(LL[4] - LL[1],1), 1-pchisq(LL[5] - LL[1],1), 1-pchisq(LL[6] - LL[1],1))

mlm.pvals 

# compute residual effect of random effects environmental variables
z.r <- glmer(PRESENCE~(1|SPP)+Elevation+Elevation2+Herb+Herb2+TreeBA+ACSA+LITU+Ca+Ca2+P, data=hh, family=binomial)

MLM.fitted = array(predict(z) - predict(z.r), c(nsite, nspp))

rownames(MLM.fitted)=c(1:54)
colnames(MLM.fitted)=c("ACRA","ARTR","CATH","CHMA","DILA","GALA","GEMA","LISU","POBI","SACA","THDI","TRGR","UVGR","VIRO")

# standardize over spp
MLM.fitted.standard <- MLM.fitted
for(j in 1:nspp) MLM.fitted.standard[,j] <- (MLM.fitted[,j]-mean(MLM.fitted[,j]))/sd(MLM.fitted[,j])

ss <- cor(MLM.fitted.standard)
U <- svd(ss)
mlm.fit <- MLM.fitted.standard %*% U$v
mlm.fit <- mlm.fit[,1:2]

# environmental variables (only those with random effects)
envir.vars <- cbind(hh$Elevation, hh$Herb, hh$LITU, hh$Ca, hh$P)
mlm.envir <- NULL
for(j in 1:5)
	mlm.envir <- cbind(mlm.envir, envir.vars[,j]*mlm.fit[,1],envir.vars[,j]*mlm.fit[,2])

envir.points <- t(array(colMeans(mlm.envir),c(2,dim(mlm.envir)[2]/2)))

# plot mlm
par(mfcol=c(1,1))
plot(-mlm.fit,xlab="PC1",ylab="PC2",type="n")
text(-mlm.fit,label=1:nsite,cex=.5)

arrow.coordMLM <- cbind(array(0,dim(envir.points)),-envir.points)

arrows(arrow.coordMLM[,1],arrow.coordMLM[,2],arrow.coordMLM[,3],arrow.coordMLM[,4], col="black", length=0.1)

text(1.3*-envir.points,label=c("Elevation", "Herb", "LITU", "Ca", "P"),cex=.7)


################################################
# analyze data: CCA
comm.matrix=array(hh$PRESENCE,c(nsite,nspp))
colnames(comm.matrix)=c("ACRA","ARTR","CATH","CHMA","DILA","GALA","GEMA","LISU","POBI","SACA","THDI","TRGR","UVGR","VIRO")
rownames(comm.matrix)=c(1:54)

z=glmer(PRESENCE~(1|SPP)+Elevation+(0+Elevation|SPP)+Herb+(0+Herb|SPP)+TreeBA+ACSA+LITU+(0+LITU|SPP)+Ca+(0+Ca|SPP)+P+(0+P|SPP),family=binomial,data=hh)

nx.cca=7

envir.matrix <- cbind(hh$Elevation,hh$Herb,hh$TreeBA,hh$ACSA,hh$LITU,hh$Ca,hh$P)
envir.matrix <- envir.matrix[1:nsite,]
colnames(envir.matrix)=c("Elevation","Herb","TreeBA","ACSA","LITU","Ca","P")

cca.fit <- cca(comm.matrix,envir.matrix, scale=TRUE)

cca.partial.pvals <- NULL
for(j in 1:nx.cca){
	# partial cca
	Z <- envir.matrix[,-j]
	cca.fit.partial <- cca(comm.matrix,envir.matrix[,j],Z,scale=TRUE)
	anova.cca <- permutest(cca.fit.partial,permutations=10000)
	pval <- mean(anova.cca$F.0 < anova.cca$F.perm)
	cca.partial.pvals <- c(cca.partial.pvals, pval)
}
cca.partial.pvals

################################################
# analyze data: RDA
rda.fit <- rda(comm.matrix,envir.matrix, scale=TRUE)

rda.partial.pvals <- NULL
for(j in 1:nx.cca){
	# partial rda
	Z <- envir.matrix[,-j]
	rda.fit.partial <- rda(comm.matrix,envir.matrix[,j],Z,scale=TRUE)
	anova.rda <- permutest(rda.fit.partial,permutations=10000)
	pval <- mean(anova.rda$F.0 < anova.rda$F.perm)
	rda.partial.pvals <- c(rda.partial.pvals, pval)
}
rda.partial.pvals

################################################
# analyze data: NMDS
nmds.fit <- metaMDS(comm.matrix, distance="bray", k=2, trymin=200, trymax=200, zerodist="add", trace=0)
nmds.fit$stress 

nmds.envir <- envfit(nmds.fit, envir.matrix, permutations=10000)
nmds.pvals <- nmds.envir$vectors[4]$pvals
nmds.pvals

################################################
# analyze data: lm 
z.lm <- lm(PRESENCE~SPP*Elevation+SPP*Herb+SPP*TreeBA+SPP*ACSA+SPP*LITU+SPP*Ca+SPP*P, data=hh)

MLM.fitted.lm <- array(z.lm$fitted.values,c(nsite,nspp))
colnames(MLM.fitted.lm)=c("ACRA","ARTR","CATH","CHMA","DILA","GALA","GEMA","LISU","POBI","SACA","THDI","TRGR","UVGR","VIRO")
rownames(MLM.fitted.lm)=c(1:54)
mlm.fit.lm <- rda(MLM.fitted.lm,scale=TRUE)

################################################
# plot results
quartz()
dev.set(3)
par(mfrow=c(2,2))

# plot mlm
plot(-mlm.fit,xlab="PC1",ylab="PC2",type="n")
text(-mlm.fit,label=1:nsite,cex=.5)
arrows(arrow.coordMLM[,1],arrow.coordMLM[,2],arrow.coordMLM[,3],arrow.coordMLM[,4], col="black", length=0.05)

text(1.3*-envir.points,label=c("Elevation", "Herb","LITU", "Ca", "P"),cex=.7)

plot(rda.fit,type="n")
text(rda.fit,display="lc",cex=0.5,col="black")
text(rda.fit,display="bp",cex=0.7,select=c(1,4:7))

plot(cca.fit,type="n",xlim=c(-8,5))
text(cca.fit,display="lc",cex=0.5,col="black")
text(cca.fit,display="bp",cex=0.7,select=c(1,2,5,6,7))

plot(nmds.fit,type="n")
text(nmds.fit,display=c("sites"),col="black",cex=0.5)

arrow.list=c(1,2,4,6) #choose significant vectors
arrow.rvalues=sqrt(nmds.envir$vector$r[arrow.list])
#correlations with species determine length of the arrow
arrow.endpt=nmds.envir$vectors$arrows[arrow.list,]*arrow.rvalues
arrow.endpt=0.7*arrow.endpt/max(arrow.endpt[1,])
#adjust length of arrows to fit plot
#0.8 is a scaling coefficient; change if desired
arrow.coord=cbind(0,0,arrow.endpt)
arrow.labels=row.names(arrow.coord)
arrows(arrow.coord[,1],arrow.coord[,2],arrow.coord[,3],arrow.coord[,4],col="black",length=0.05)
text(1.1*arrow.coord[,3:4],label=arrow.labels,cex=0.7)

###############################################
dev.set(4)
par(mfrow=c(1,2))

plot(mlm.fit.lm, type="n")
text(mlm.fit.lm,display="sites",cex=0.5,col="black")
plot(mlm.envir.lm, col="black",cex=0.7)
plot(rda.fit,type="n")
text(rda.fit,cex=0.5,col="black")
text(rda.fit,display="bp",cex=0.7)

################################################
# method comparisons
################################################

# access site scores
mlm.scores <- mlm.fit
rda.scores <- scores(rda.fit, display="lc")
cca.scores <- scores(cca.fit, display="lc")
nmds.scores <- nmds.fit$points

# procrustes
procrustes(mlm.scores, rda.scores, scale = TRUE, symmetric = TRUE) #0.1857
procrustes(mlm.scores, cca.scores, scale = TRUE, symmetric = TRUE) #0.1872
procrustes(mlm.scores, nmds.scores, scale = TRUE, symmetric = TRUE) #0.63

procrustes(rda.scores, cca.scores, scale = TRUE, symmetric = TRUE) #0.1373
procrustes(rda.scores, nmds.scores, scale = TRUE, symmetric = TRUE) #0.7258
procrustes(cca.scores, nmds.scores, scale = TRUE, symmetric = TRUE) #0.7075


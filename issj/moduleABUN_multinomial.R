##
## Demo fitting a multinomial-mixture model
##
## Here we use bird point count data collected using a "removal" protocol.
## the observations y[i] = (y[i,1],y[i,2],y[i,3]) are the numbers of birds removed
## in each of three 5 minute intervals at point i. We assume y[i] ~ multinomial(N[i], ...)
## and that N[i] ~ Poisson(lam[i]) or some other suitable distribution
## We use the oven bird data set that comes with the unmarked package
##
## To load the data do this:
library("unmarked")
data(ovendata)
# note: every covariate is a "siteCovs" with a single primary sample, even if
# the covariate effects "p"

ovenFrame <- unmarkedFrameMPois(ovendata.list$data,
    siteCovs=as.data.frame(scale(ovendata.list$covariates[,-1])),
    type = "removal")


fm0 <- multinomPois(~1 ~1,ovenFrame)
fm1 <- multinomPois(~ 1 ~ ufp + trba, ovenFrame)
fm2 <- multinomPois(~ 1 ~ ufp + trba + ufp*trba, ovenFrame)
fm3 <- multinomPois(~ ufp ~ ufp + trba, ovenFrame)
fm4 <- multinomPois(~ ufp ~ ufp + trba +ufp*trba, ovenFrame)

ms<- fitList( 
"lam(.)p(.)"                                = fm0,
"lam(ufp+trba)p(.)"                         = fm1,
"lam(ufp+trba+ufp*trba)p(.)"                = fm2,
"lam(ufp+trba)p(ufp)"                       = fm3,
"lam(ufp+trba+ufp*trba)p(ufp)"              = fm4)

# Rank them by AIC
(ms1 <- modSel(ms))

# Table with everything you could possibly need
coef(ms1)
toExport <- as(ms1, "data.frame")



##
## use gmultmix to control the abundance model and fit OPEN models
##
## numPrimary=1
## extra formula, order: lambda phi p
## DIFFERENT FORMAT AND ORDER
## LIKELIHOODS ARE SCALED DIFFERENTLY TOO (DON'T COMPARE ACROSS METHODS)
##
##
data(ovendata)

ovenFrame <- unmarkedFrameGMM(ovendata.list$data,
    siteCovs=as.data.frame(scale(ovendata.list$covariates[,-1])),numPrimary=1,
    type = "removal")

fm0 <- gmultmix(~1 , ~1,  ~1 ,data=ovenFrame)
fm1 <- gmultmix(~ ufp+trba , ~1 , ~ 1, data=ovenFrame)
fm2 <- gmultmix(~ ufp+trba+ufp*trba , ~1 , ~ 1, data=ovenFrame)
# maybe p also depends on understory foliage?
fm3 <- gmultmix(~ ufp +trba, ~1 , ~ ufp,data= ovenFrame)
fm4 <- gmultmix(~ ufp + trba+ufp*trba,  ~1 , ~ ufp,data= ovenFrame)

fm0nb <- gmultmix(~1 , ~1,  ~1 ,mixture="NB",data=ovenFrame)
fm1nb <- gmultmix(~ ufp+trba , ~1 , ~ 1, mixture="NB",data=ovenFrame)
fm2nb <- gmultmix(~ ufp+trba+ufp*trba , ~1 , ~ 1, mixture="NB",data=ovenFrame)
# maybe p also depends on understory foliage?
fm3nb <- gmultmix(~ ufp +trba, ~1 , ~ ufp,,mixture="NB",data= ovenFrame)
fm4nb <- gmultmix(~ ufp + trba+ufp*trba,  ~1 , ~ ufp,mixture="NB",data= ovenFrame)

gms<- fitList(
"lam(.)p(.)"                                = fm0,
"lam(ufp+trba)p(.)"                         = fm1,
"lam(ufp+trba+ufp*trba)p(.)"                = fm2,
"lam(ufp+trba)p(ufp)"                       = fm3,
"lam(ufp+trba+ufp*trba)p(ufp)"              = fm4,
"NB,lam(.)p(.)"                                = fm0nb,
"NB,lam(ufp+trba)p(.)"                         = fm1nb,
"NB,lam(ufp+trba+ufp*trba)p(.)"                = fm2nb,
"NB,lam(ufp+trba)p(ufp)"                       = fm3nb,
"NB,lam(ufp+trba+ufp*trba)p(ufp)"              = fm4nb)

##
##  Note: Dispersion parameter of NB is usually near the boundary 1/tau = 0
##  indicating no over-dispersion which is supported by the model-selection
##  results


# Rank them by AIC
(gms1 <- modSel(gms))

# Table with everything you could possibly need
coef(gms1)
# generates an error -- bug needs fixed
toExport <- as(gms1, "data.frame")  ## This doesn't work














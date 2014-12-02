library(unmarked)
source("issj/data_100m.R")    # load ISSJ data file 


## Format data using unmarkedFrameGDS()
# Note that it's unmarkedFrameGDS, not unmarkedFrameDS
db <- c(0, 100, 200, 300)

# Fall data
umfFall <- unmarkedFrameGDS(y = Xall.100m, siteCovs=covs, dist.breaks=db,
                            unitsIn="m", survey="point", numPrimary=1)
summary(umfFall)    

# Spring data
umfSpr <- unmarkedFrameGDS(y = Xall.spring.100m, siteCovs=covs, dist.breaks=db,
                           unitsIn="m", survey="point", numPrimary=1)
summary(umfSpr)    

## Fit models
## Note that the order of formulas is lambda, phi, p
## Fall models

fall <- list()
fall$Null <- gdistsamp(lambdaformula=~1, phiformula=~1, pformula=~1, 
                       umfFall, output="abund", K=200)
fall$Null_NB <- gdistsamp(~1, ~1, ~1, umfFall, output="abund", mixture="NB", 
                          K=200)
fall$Chap._NB <- gdistsamp(~chap, ~1, ~1, umfFall, output="abund", mixture="NB",
                           K=200)
fall$Chap2._NB <- gdistsamp(~chap + I(chap^2), ~1, ~1, umfFall, 
                            output="abund", mixture="NB", K=200)
fall$Elev._NB <- gdistsamp(~elev, ~1, ~1, umfFall, output="abund", mixture="NB",
                           K=200)
fall$Elev2._NB <- gdistsamp(~elev + I(elev^2), ~1, ~1, umfFall, output="abund", 
                            mixture="NB", K=200)
fall$Forest._NB <- gdistsamp(~forest, ~1, ~1, umfFall, output="abund", 
                             mixture="NB", K=200)
fall$Forest2._NB <- gdistsamp(~forest + I(forest^2), ~1, ~1, umfFall, 
                              output="abund", mixture="NB", K=200)
fall$.Forest_NB <- gdistsamp(~1, ~1, ~forest, umfFall, output="abund", 
                             mixture="NB", starts=c(0, 4.5, 0, 0), K=200)
fall$.Chap_NB <- gdistsamp(~1, ~1, ~chap, umfFall, mixture="NB", 
                           starts=c(0, 4.5, 0, 0), output="abund", K=200) 
fall$C2E._NB <- gdistsamp(~chap + I(chap^2) + elev, ~1, ~1, 
                          umfFall, mixture="NB", starts=c(1, 0.5, 0, 0, 4.5, 0), output="abund", 
                          K=200)
fall$C2F2._NB <- gdistsamp(~chap + I(chap^2) + forest + I(forest^2), ~1, ~1, 
                           umfFall, mixture="NB", starts=c(1, 0.5, 0, 0, 0, 4.6, 0), 
                           output="abund", K=200)
fall$C2E.F_NB <- gdistsamp(~chap + I(chap^2) + elev, ~1, ~forest,
                           umfFall, mixture="NB", starts=c(1, 0.5, 0, 0, 4.6, 0, 0), 
                           output="abund", K=200)
fall$C2E.C_NB <- gdistsamp(~chap + I(chap^2) + elev, ~1, ~chap,
                           umfFall, mixture="NB", starts=c(1, 0.5, 0, 0, 4.6, 0, 0), 
                           output="abund", K=200)

## Assemble the various model fits into a "fitList" and do model selection
fitsFall <- fitList(fits=fall)
(msFall <- modSel(fitsFall))

## stuff needed for the table
exportFall <- cbind(coef(msFall), msFall@Full$AIC) # Coefs and AIC

## Spring models

spr <- list()
spr$Null <- gdistsamp(~1, ~1, ~1, umfSpr, output="abund", K=200)
spr$Null_NB <- gdistsamp(~1, ~1, ~1, umfSpr, output="abund", mixture="NB", 
                         K=200)
spr$Chap._NB <- gdistsamp(~chap, ~1, ~1, umfSpr, output="abund", mixture="NB", 
                          K=200)
spr$Chap2._NB <- gdistsamp(~chap + I(chap^2), ~1, ~1, umfSpr, 
                           output="abund", mixture="NB", K=200)
spr$Elev._NB <- gdistsamp(~elev, ~1, ~1, umfSpr, output="abund", mixture="NB", 
                          K=200)
spr$Elev2._NB <- gdistsamp(~elev + I(elev^2), ~1, ~1, umfSpr, output="abund", 
                           mixture="NB", K=200)
spr$Forest._NB <- gdistsamp(~forest, ~1, ~1, umfSpr, output="abund", 
                            mixture="NB", K=200)
spr$Forest2._NB <- gdistsamp(~forest + I(forest^2), ~1, ~1, umfSpr, 
                             output="abund", mixture="NB", K=200)
spr$.Forest_NB <- gdistsamp(~1, ~1, ~forest, umfSpr, output="abund", 
                            mixture="NB", starts=c(0, 4.5, 0, 0), K=200)
spr$.Chap_NB <- gdistsamp(~1, ~1, ~chap, umfSpr, mixture="NB", 
                          starts=c(0, 4.5, 0, 0), output="abund", K=200) 
spr$C2E._NB <- gdistsamp(~chap + I(chap^2) + elev, ~1, ~1, 
                         umfSpr, mixture="NB", starts=c(1, 0.5, 0, 0, 4.5, 0), output="abund", K=200)
spr$C2E2._NB <- gdistsamp(~chap + I(chap^2) + elev + I(elev^2), ~1, ~1, 
                          umfSpr, mixture="NB", starts=c(1, 0.5, 0, 0, 0, 4.6, 0), 
                          output="abund", K=200)
spr$C2E2.F_NB <- gdistsamp(~chap + I(chap^2) + elev + I(elev^2), ~1, ~forest,
                           umfSpr, mixture="NB", starts=c(1, 0.5, 0, 0, 0, 4.6, 0, 0), 
                           output="abund", K=200)

## assemble model fits into a fitList and do model selection
fitsSpr <- fitList(fits=spr)
(msSpr <- modSel(fitsSpr))

## stuff needed for the table
exportSpr <- cbind(coef(msSpr), msSpr@Full$AIC) # Coefs and AIC

## save all of this stuff
save.image("issj.RData")


## GoF analysis

## define a fit statistic
freeTuke <- function(fm) {
  observed <- getY(fm@data)
  expected <- fitted(fm)
  sum((sqrt(observed) - sqrt(expected))^2)
}

# should set nsim=100 or more, but nsim=2 here for illustration
pbFall <- parboot(fall$C2E.C_NB, freeTuke, nsim=2, report=1)
pbSpr <- parboot(spr$C2E2.F_NB, freeTuke, nsim=2, report=1)



# Pop size calculation for the sample points
getN <- function(fm) sum(predict(fm, type="lambda")[,1])

getN(fall$C2E.C_NB)
getN(spr$C2E2.F_NB)

## extract the coefficients for the best model

betaspr<-coef(spr$C2E2.F_NB)
betaspr<-betaspr[c("lambda(Int)","lambda(chap)","lambda(elev)","lambda(I(chap^2))","lambda(I(elev^2))")]
beta<-coef(fall$C2E.C_NB)
betafall<-beta[c("lambda(Int)","lambda(chap)","lambda(elev)","lambda(I(chap^2))")]

## use those to verify what predict() is doing. predict is doing exactly this:
Xspr<-cbind(rep(1,307),covs[,1],covs[,2],covs[,1]^2,covs[,2]^2)
pspr<-Xspr%*%(betaspr)
sum(exp(pspr))
## sample area, in ha:
point.area<- pi*300*300/10000

# now compute the population size for the entire island ("gelev" "gfor" "gchap" are grid values)
# note the multiplier 9/point.area rescales density to units of animals/9 ha, the size
# of a 300 m grid cell

## grab the grid covariates (these are standardized already)
source("grid_covariates.R")

## compute the expected abundance at each grid cell
## sum(lamnew) is the expected value of total population size
newdata<-data.frame(elev=gelev,forest=gfor,chap=gchap)
Xspr<- cbind(rep(1,nrow(newdata)),gchap,gelev,gchap*gchap,gelev*gelev)
lamnew<- (9/point.area)*exp(Xspr%*%(betaspr))
Xfall<-cbind(rep(1,nrow(newdata)),gchap,gelev,gchap*gchap)
lamnew<-(9/point.area)*exp(Xfall%*%(betafall))

### Now we create functions to compute total population size to use in a bootstrap
### simulation

gridN<- function(fm) {
  beta<-coef(fm)
  beta<-beta[c("lambda(Int)","lambda(chap)","lambda(elev)","lambda(I(chap^2))","lambda(I(elev^2))")]
  sum(   (9/point.area)*exp(Xspr%*%(beta) ) )  
}
gridN<- function(fm) {
  beta<-coef(fm)
  beta<-beta[c("lambda(Int)","lambda(chap)","lambda(elev)","lambda(I(chap^2))")]
  sum(   (9/point.area)*exp(Xfall%*%(beta) ) )  
}

### Do the simulations. In practice, use nsim=100 or 200
pbN.spr<- parboot(spr$C2E2.F_NB,gridN,nsim=2,report=1)
pbN.fall<- parboot(fall$C2E.C_NB,gridN,nsim=2,report=1)

###
###
### predictions based on 1985 data
###
###
###
### computing the predictions 

betaspr<-coef(spr$C2E2.F_NB)
betaspr<-betaspr[c("lambda(Int)","lambda(chap)","lambda(elev)","lambda(I(chap^2))","lambda(I(elev^2))")]
Xnew<- cbind(rep(1,nrow(newdata)),gchap85,gelev,gchap85*gchap85,gelev*gelev)
lamnew<- (9/point.area)*exp(Xnew%*%(betaspr))

beta<-coef(fall$C2E.C_NB)
beta<-beta[c("lambda(Int)","lambda(chap)","lambda(elev)","lambda(I(chap^2))")]
Xnew<-cbind(rep(1,nrow(newdata)),gchap85,gelev,gchap85*gchap85)
lamnew<-(9/point.area)*exp(Xnew%*%(beta))

### Getting a bootstrap SE for the 1985 predictions

predict85spr<-function(fm){
  betaspr<-coef(fm)
  betaspr<-betaspr[c("lambda(Int)","lambda(chap)","lambda(elev)","lambda(I(chap^2))","lambda(I(elev^2))")]
  Xnew<- cbind(rep(1,nrow(newdata)),gchap85,gelev,gchap85*gchap85,gelev*gelev)
  lamnew<- (9/point.area)*exp(Xnew%*%(betaspr))
  sum(lamnew)
}

predict85fall<-function(fm){
  beta<-coef(fm)
  beta<-beta[c("lambda(Int)","lambda(chap)","lambda(elev)","lambda(I(chap^2))")]
  Xnew<-cbind(rep(1,nrow(newdata)),gchap85,gelev,gchap85*gchap85)
  lamnew<-(9/point.area)*exp(Xnew%*%(beta))
  sum(lamnew)
}

## Do the bootstrapping. In practice, use nsim=100 or more.
##
pbN.spr85<- parboot(spr$C2E2.F_NB,predict85spr,nsim=2,report=1)
pbN.fall85<- parboot(fall$C2E.C_NB,predict85fall,nsim=2,report=1)
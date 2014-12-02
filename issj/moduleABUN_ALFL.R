

# -----------------------------------------------------------------------
# ------------------------ Poisson regression  -------------------------
# -----------------------------------------------------------------------


# ------------------------- Simulate data ------------------------------

# PART 1: Set-up the data (simulate it here)
#
# Create a covariate called vegHt
nSites <- 100
set.seed(443)                # so that we all get the same values of vegHt
vegHt <- runif(nSites, 1, 3) # uniform from 1 to 3


# Suppose that expected population size increases with vegHt
# The relationship is described by an intercept of -3 and
#    a slope parameter of 2 on the log scale

lambda <- exp(-3 + 2*vegHt)


# Now we go to 100 sites and observe the # of individuals (perfectly)
# 
N <- rpois(nSites, lambda)


## PART 2: Fit some models
##
##
# We can fit a model that relates abundance to vegHt using the glm() function
#  with "family=Poisson":

fm.glm1 <- glm(N ~ vegHt, family=poisson)

##
## PART 3: Do some analysis of the results
##
plot(vegHt, N, xlab="Vegetation height", ylab="Abundance (N)")
glm1.est <- coef(fm.glm1)
plot(function(x) exp(-3 + 2*x), 1, 3, add=TRUE, lwd=2)
plot(function(x) exp(glm1.est[1] + glm1.est[2]*x), 1, 3, add=TRUE,
     lwd=2, col="blue")
legend(1, 15.9, c("Truth", "Estimate"), col=c("black", "blue"), lty=1,
       lwd=2)



# -----------------------------------------------------------------------
# ------ Imperfect observation of abundance using a binomial sampling model
## ----- suppose y is a BINOMIAL sample based on population size N and parameter
## ----- p = "detection probability". i.e.,:
## -----   y ~ binomial(N, p)
## ----- the canonical "point counting" model
## ----- In practice, we think p < 1, say p = 0.6. This is INDIVIDUAL-LEVEL
## ------ detection -- i.e., each individual is detected with probability p
##
## ----- In this case, N is a LATENT VARIABLE (i.e., unobserved).
# -----------------------------------------------------------------------



nVisits <- 4
p <- 0.6
y <- matrix(NA, nSites, nVisits)
for(i in 1:nSites) {
    y[i,] <- rbinom(nVisits, N[i], p)
}

# Format for unmarked and summarize
library(unmarked)
umf <- unmarkedFramePCount(y=y, siteCovs=as.data.frame(vegHt))
summary(umf)

# Fit a model and extract estimates
# Detection covariates follow first tilde, then come abundance covariates

fm.nmix1 <- pcount(~1 ~vegHt, data=umf)

# Note the following warning message:
#> fm.nmix1 <- pcount(~1 ~vegHt, data=umf)
#Warning message:
#In pcount(~1 ~ vegHt, data = umf) : K was not specified and was set to 116.
# K is upper limit of summation for calculating the likelihood

# Consider other abundance models: NB = negative binomial, ZIP = zero-inflated Poisson
## currently no others available. 
##
fm.nmix2<-  pcount(~1 ~vegHt, data=umf,mixture="NB")
fm.nmix3<-  pcount(~1 ~vegHt, data=umf,mixture="ZIP")
beta1 <- coef(fm.nmix1) 


#
# Note, estimates of detection coefficients are on the logit-scale
# When there are no covariates, we can back-transform using:
exp(beta1[3]) / (1+exp(beta1[3]))   # or
plogis(beta1[3])                    # or
backTransform(fm.nmix1, type="det") # estimate with SE

# When covariates are present we can do something like
#
plot(function(x) exp(beta1[1] + beta1[2]*x), 1, 3,
     xlab="vegetation height", ylab="Expected Abundance")


# Or suppose you want predictions for new values of vegHt, say 1.2 and 3.1
newdat <- data.frame(vegHt=c(1.2, 3.1))
predict(fm.nmix1, type="state", newdata=newdat)


#
# -----------------------------------------------------------------------
# ----------- Now let's start from scratch with real data --------------
# -----------------------------------------------------------------------
#



# PART 1: Set-up the data for analysis
#
#
# -------------------------- Format data ---------------------------------
# This a subset of point-count data from Chandler et al. (Auk 2009)
# alfl is Alder Flycatcher (Empidonax alnorum)

# Import data and check structure
#alfl.data <- read.csv("alfl05.csv", row.names=1)
alfl.data <- read.csv("http://sites.google.com/site/unmarkedinfo/home/webinars/2012-january/data/alfl05.csv?attredirects=0&d=1", row.names=1)

str(alfl.data)


# Pull out count matrix 
# No need to covert to binary as we did for occupancy model

alfl.y <- alfl.data[,c("alfl1", "alfl2", "alfl3")]

# Standardize site-covariates
woody.mean <- mean(alfl.data$woody)
woody.sd <- sd(alfl.data$woody)
woody.z <- (alfl.data$woody-woody.mean)/woody.sd

struct.mean <- mean(alfl.data$struct)
struct.sd <- sd(alfl.data$struct)
struct.z <- (alfl.data$struct-struct.mean)/struct.sd


# Create unmarkedFrame
library(unmarked)
alfl.umf <- unmarkedFramePCount(y=alfl.y,
    siteCovs=data.frame(woody=woody.z, struct=struct.z),
    obsCovs=list(time=alfl.data[,c("time.1", "time.2", "time.3")],
                 date=alfl.data[,c("date.1", "date.2", "date.3")]))
summary(alfl.umf)


# Here's an easy way to standardize covariates after making the UMF
obsCovs(alfl.umf) <- scale(obsCovs(alfl.umf))
summary(alfl.umf)


#
#
# PART 2: Fit some models
#
# -------------------------- Model fitting  -----------------------------

(fm1 <-  pcount(~1 ~1, alfl.umf))
backTransform(fm1, type="state")
backTransform(fm1, type="det")

(fm2 <- pcount(~date+time ~1, alfl.umf))
(fm3 <- pcount(~date+time ~woody, alfl.umf))
(fm4 <- pcount(~date+time ~woody+struct, alfl.umf))
(fm5 <- pcount(~date+time ~1, alfl.umf,mixture="NB"))
(fm6 <- pcount(~date+time ~1, alfl.umf,mixture="ZIP"))
(fm7 <- pcount(~date+time ~woody,alfl.umf,mixture="ZIP"))
(fm8 <- pcount(~date+time ~struct,alfl.umf,mixture="ZIP"))
(fm9 <- pcount(~date+time ~woody+struct, alfl.umf,mixture="ZIP"))
(fm10<- pcount(~date+time ~woody+struct, alfl.umf,mixture="NB"))

# -------------------------- Model selection -----------------------------

# Put the fitted models in a "fitList"
fms <- fitList("lam(.)p(.)"                    = fm1,
               "lam(.)p(date+time)"            = fm2,
               "lam(woody)p(date+time)"        = fm3,
               "lam(woody+struct)p(date+time)" = fm4,
               "lam(.)p(date+time)NB"          = fm5,
               "lam(.)p(date+time)ZIP"         = fm6,
               "lam(woody)p(date+time)ZIP"     = fm7,
               "lam(struct)p(date+time)ZIP"    = fm8,
               "lam(woody+struct)p(date+time)ZIP"=fm9,
               "lam(woody+struct)p(date+time)NB" =fm10)

# Rank them by AIC
(ms <- modSel(fms))

# Table with everything you could possibly need
coef(ms)
toExport <- as(ms, "data.frame")


#
#
# PART 3: Do some analysis of the results
#
#
# ---------------------------- Prediction --------------------------------


# Expected detection probability as function of time of day
# We standardized "time", so we predict over range of values on that scale
# We must fix "date" at some arbitrary value (let's use the mean)
newData1 <- data.frame(time=seq(-2.08, 1.86, by=0.1), date=0)
E.p <- predict(fm4, type="det", newdata=newData1, appendData=TRUE)
head(E.p)

par(mfrow=c(2,1))
# Plot it
plot(Predicted ~ time, E.p, type="l", ylim=c(0,1),
     xlab="time of day (standardized)",
     ylab="Expected detection probability")
lines(lower ~ time, E.p, type="l", col=gray(0.5))
lines(upper ~ time, E.p, type="l", col=gray(0.5))



# Expected abundance over range of "woody"
newData2 <- data.frame(woody=seq(-1.6, 2.38,,50),struct=seq(-1.8,3.2,,50))
E.N <- predict(fm4, type="state", newdata=newData2, appendData=TRUE)
head(E.N)

# Plot predictions with 95% CI
plot(Predicted ~ woody, E.N, type="l", ylim=c(-.1,max(E.N$Predicted)),
     xlab="woody vegetation (standardized)",
     ylab="Expected abundance, E[N]")
lines(lower ~ woody, E.N, type="l", col=gray(0.5))
lines(upper ~ woody, E.N, type="l", col=gray(0.5))

par(mfrow=c(1,1))
# Plot it again, but this time convert the x-axis back to original scale
plot(Predicted ~ woody, E.N, type="l", ylim=c(-.1,max(E.N$Predicted)),
     xlab="Percent cover - woody vegetation",
     ylab="Expected abundance, E[N]",
     xaxt="n")
xticks <- -1:2
xlabs <- xticks*woody.sd + woody.mean
axis(1, at=xticks, labels=round(xlabs, 1))
lines(lower ~ woody, E.N, type="l", col=gray(0.5))
lines(upper ~ woody, E.N, type="l", col=gray(0.5))


## Goodness-of-Fit
###
###
### Here's an example of a bootstrap GoF analysis.
### Best model is in "fm4" object
###

# Function returning three fit-statistics.
fitstats <- function(fm) {
    observed <- getY(fm@data)
    expected <- fitted(fm)
    resids <- residuals(fm)
    sse <- sum(resids^2)
    chisq <- sum((observed - expected)^2 / expected)
    freeTuke <- sum((sqrt(observed) - sqrt(expected))^2)
    out <- c(SSE=sse, Chisq=chisq, freemanTukey=freeTuke)
    return(out)
    }

(pb <- parboot(fm4, fitstats, nsim=100, report=1))
### To look at bootstrap distributions do this
###plot(pb, main="")
print(pb)

## Now lets bootstrap a summary statistic 
## This is not too meaningful right now but we will do a similar thing later in 
## a more relevant context

# Total population size (derived parameter)
Nhat <- function(fm) {
    N <- sum(predict(fm, type="state")$Predicted, na.rm=TRUE)
    }
    
(pb.N <- parboot(fm4, Nhat, nsim=25, report=5))
plot(pb.N)



























# Here's an example of model-averaging predictions
# See pg 150, section 4.2.1, of Burnham and Anderson (2002)
# This might be worthwhile since fm3 and fm4 had similar support

newData3 <- data.frame(woody=seq(-1.6, 2.38,,50), struct=seq(-1.8,3.2,,50))

## averages over _all_ models in the fit list "fms":
E.N.bar <- predict(fms, type="state", newdata=newData3,
                     appendData=TRUE)
head(E.N.bar)

# Plot it
plot(Predicted ~ woody, E.N.bar, type="l", ylim=c(-0.1, max(E.N$Predicted)),
     xlab="Percent cover - woody vegetation",
     ylab="Expected abundance, E[N]",
     xaxt="n")
xticks <- -1:2
xlabs <- xticks*woody.sd + woody.mean
axis(1, at=xticks, labels=round(xlabs, 1))
lines(lower ~ woody, E.N.bar, type="l", col=gray(0.5))
lines(upper ~ woody, E.N.bar, type="l", col=gray(0.5))


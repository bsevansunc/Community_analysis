# Working with distance sampling data
# From the vignette on distance sampling written by RBC

# Data might exist as a matrix of distance observations:

library(unmarked)
source("http://sites.google.com/site/unmarkedinfo/home/webinars/2012-january/data/utils.R?attredirects=0&d=1")
source("utils.R")

#
# PART 1: Set-up the data for analysis.
# we're using the toy simulated data set that comes with unmarked:
#
dists <- read.csv(system.file("csv", "distdata.csv", package="unmarked"))
dists
levels(dists$transect)

#
# here we illustrate adding a site where no detections occurred. 
# The data set has transects a-f but there was also a "g" transect.
# So lets add that to the list of levels
# we MUST ACCOUNT FOR THE 0 OBSERVATIONS!
#
levels(dists$transect) <- c(levels(dists$transect), "g")
levels(dists$transect)

#
# lets put the data into "multinomial format" using the formatDistData function
#
yDat <- formatDistData(dists, distCol="distance",
         transectNameCol="transect", dist.breaks=c(0, 5, 10, 15, 20))

# Transect-specific covariates:

(covs <- data.frame(canopyHt = c(5, 8, 3, 2, 4, 7, 5),
         habitat = c('A','A','A','A','B','B','B'), row.names=letters[1:7]))

##
# Lets make an unmarkedFrameDS
# "By organizing the data this way, the user does not need to repetitively 
# specify these arguments during each call to distsamp, thereby reducing 
# the potential for errors and facilitating data summary and manipulation."

 umf <- unmarkedFrameDS(y=as.matrix(yDat), siteCovs=covs, survey="line",
 dist.breaks=c(0, 5, 10, 15, 20), tlength=rep(100, 7),
 unitsIn="m")

# Note: no obsCovs as there usually is not (this would be nsites x nclasses)

#
#It is important that both transect lengths and distance break
#points are provided in the same units specified by unitsIn.
#
summary(umf)
hist(umf, xlab="distance (m)", main="", cex.lab=0.8, cex.axis=0.8)

#
# PART 2: FIT SOME MODELS NOW
#
#
hn_Null <- distsamp(~1~1, umf)
hn_Null <- distsamp(~1~1, umf, keyfun="halfnorm", output="density", unitsOut="ha")
haz_Null <- distsamp(~1~1, umf, keyfun="hazard")
hn_Hab.Ht <- distsamp(~canopyHt ~habitat, umf)

#
# USE OF OFFSET
#
# if units of varying length do this: 
length<-c(1,1,1,2,2,1,1)
siteCovs(umf)$length<-length

hn_Hab.Ht.2<- distsamp(~canopyHt ~ habitat + offset(log(length)),umf)
# A summary method shows extra details including the scale on which parameters were estimated and
# convergence results.

haz_Null

dsfitlist<- fitList(
"hm(.)lam(.)"               = hn_Null,
"haz(.)lam(.)"              = haz_Null,
"hn(canopy)lam(habitat)"    = hn_Hab.Ht)
#
# Rank them by AIC
#
(ms1 <- modSel(dsfitlist))
#
# Table with everything you could possibly need
#
coef(ms1)
toExport <- as(ms1, "data.frame")

##
## PART 3: SUMMARY ANALYSIS
## density = per ha (log scale) so lets back-transform
##
summary(hn_Null)
backTransform(hn_Null, type="state")       # animals / ha
exp(coef(hn_Null, type="state", altNames=TRUE))     # same
backTransform(hn_Null, type="det")         # half-normal SD
hist(hn_Null, xlab="Distance (m)")	# Only works when there are no detection covars

##
## EXAMPLE 2:
## Line transect example from the help file
## Shows how to estimate the effective strip width and probability of encounter.
## Shows also how to plot the fitted detection function
##
##
data(linetran)
ltUMF <- unmarkedFrameDS(y = as.matrix(linetran[,c("dc1", "dc2","dc3", "dc4","dc5")]),
   siteCovs = linetran[,c("Length", "area", "habitat")],
   dist.breaks = c(0, 5, 10, 15, 20),
   tlength = linetran$Length * 1000, survey = "line", unitsIn = "m")

ltUMF
summary(ltUMF)
hist(ltUMF)

# Half-normal detection function. Density output (log scale). No covariates.
(fm1 <- distsamp(~ 1 ~ 1, ltUMF))

# Some methods to use on fitted model
summary(fm1)
backTransform(fm1, type="state")       # animals / ha
exp(coef(fm1, type="state", altNames=TRUE))     # same
backTransform(fm1, type="det")         # half-normal SD
hist(fm1, xlab="Distance (m)")	# Only works when there are no detection covars

# Effective strip half-width
(eshw <- integrate(gxhn, 0, 20, sigma=10.9)$value)

# Detection probability
eshw / 20 # 20 is strip-width


# Halfnormal. Covariates affecting both density and and detection.
(fm2 <- distsamp(~area + habitat ~ habitat, ltUMF))

# Hazard-rate detection function.
(fm3 <- distsamp(~ 1 ~ 1, ltUMF, keyfun="hazard"))

# Plot detection function.
fmhz.shape <- exp(coef(fm3, type="det"))
fmhz.scale <- exp(coef(fm3, type="scale"))
plot(function(x) gxhaz(x, shape=fmhz.shape, scale=fmhz.scale), 0, 25,
	xlab="Distance (m)", ylab="Detection probability")



##
##
## Point transect example from the helpfile
## Not much new here except you have to remember to adjust for the geometry of the circle
## for some of the calculations
##

data(pointtran)
ptUMF <- 
	unmarkedFrameDS(y = as.matrix(pointtran[,c("dc1","dc2","dc3","dc4","dc5")]),
	siteCovs = pointtran[,c("area","habitat")],
	dist.breaks = seq(0, 25, by=5), survey = "point", unitsIn = "m")
	

# Half-normal.
(fmp1 <- distsamp(~ 1 ~ 1, ptUMF))
hist(fmp1, ylim=c(0, 0.07), xlab="Distance (m)")

# effective radius
sig <- exp(coef(fmp1, type="det"))
ea <- 2*pi * integrate(grhn, 0, 25, sigma=sig)$value # effective area
sqrt(ea / pi) # effective radius

# detection probability
ea / (pi*25^2)


(fmp1 <- distsamp(~ 1 ~ 1, ptUMF))
fm2 <- distsamp(~1 ~area, ptUMF)
fms <- fitList(fmp1, fm2)
modSel(fms)





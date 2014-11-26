#--------------------------------------------------------------------------------*
# ---- SET-UP ----
#================================================================================*

# Load packages:


# Run source script to prepare data for analysis:

source('analysis-scripts/pc_prep2.R')

#--------------------------------------------------------------------------------*
# ---- ORDINATION BY SPECIES ----
#================================================================================*

t1 = cbind(lc, counts.sp)#[,-c(2:4)]

ord2 <- cca(counts.sp ~ imp+can+I(imp^2), data=lc)

# Ordination of imp and can only

ord2

plot(ord2, type = "n", xlim = c(-3,3), ylim =c(-3.5,3.5))
#points(ord2, display = "sites", cex = 0.5, pch=19, col="red")
text(ord2, display = "spec", cex=0.7, col="darkgreen")
ord2.fit <- envfit(ord2 ~ imp+can, data=lc, perm=1000)
plot(ord2.fit)

# Upper left quadrant, high can, low imp:

plot(ord2, type = "n", xlim = c(-1.5,0),ylim = c(0,1))
#points(ord2, display = "sites", cex = 0.5, pch=19, col="red")
text(ord2, display = "spec", cex=0.7, col="darkgreen")
#plot(ord2.fit)

# Lower left quadrant, low can, low imp:

plot(ord2, type = "n", xlim = c(-2,0),ylim = c(-2.5,0))
#points(ord2, display = "sites", cex = 0.5, pch=19, col="red")
text(ord2, display = "spec", cex=0.7, col="darkgreen")

# Upper right quadrant, "high" can, high imp:

plot(ord2, type = "n", xlim = c(0,2),ylim = c(0,1))
#points(ord2, display = "sites", cex = 0.5, pch=19, col="red")
text(ord2, display = "spec", cex=0.7, col="darkgreen")

# Lower right quadrant, low can, high imp:

plot(ord2, type = "n", xlim = c(0,5),ylim = c(-1,0))
#points(ord2, display = "sites", cex = 0.5, pch=19, col="red")
text(ord2, display = "spec", cex=0.7, col="darkgreen")

anova(ord2)

# Explained variance

ord2 # Ack! 0.083!

cca(counts.sp ~ imp, data=lc) # .045 ... can explains .038 that imp does not

cca(counts.sp ~ can, data=lc) # .044 ... imp explains .039 that can does not

#--------------------------------------------------------------------------------*
# ---- ORDINATION BY GUILD ----
#================================================================================*

# Explore explained variance

ord.nest = cca(nest.counts ~ imp + can + imp:can, data = lc)
ord.troph = cca(trophic.counts ~ imp + can + imp:can, data = lc)
ord.ftroph = cca(foraging_trophic.counts ~ imp + can + imp:can, data = lc)

ord.troph # .1034

cca(trophic ~ imp, data=lc) # .0861, .0173 explained by can not explained by imp

cca(trophic ~ can, data=lc) # .0652, .0382 explained by imp not explained by can

#

ord.nest # .2005

cca(nest ~ imp, data=lc) # .1826, .0179 explained by can not explained by imp

cca(nest ~ can, data=lc) # .1447, .0558 explained by imp not explained by can

#

ord.mig # .05710

cca(mig ~ imp, data=lc) # .0526, .0045 explained by can not explained by imp

cca(mig ~ can, data=lc) # .0196, .0375 explained by imp not explained by can

plot(ord.troph, type = "n",xlim = c(-.75,.75),ylim = c(-.75,.75))
ord.fit = envfit(ord.troph ~ imp + can + imp:can, data=lc, perm=1000)
plot(ord.fit, col =1)
text(ord.troph, display = "spec", cex=0.75, col="darkgreen")

plot(ord.nest, type = "n",xlim = c(-.75,.75),ylim = c(-.75,.75))
ord.fit = envfit(ord.nest ~ imp + can, data=lc, perm=1000)
plot(ord.fit, col =1)
text(ord.nest, display = "spec", cex=1, col="red")

plot(ord.mig, type = "n",xlim = c(-.15,.15),ylim = c(-.15,.15))
ord.fit = envfit(ord.mig ~ imp + can, data=lc, perm=1000)
plot(ord.fit, col =1)
text(ord.mig, display = "spec", cex=1.5, col="blue")


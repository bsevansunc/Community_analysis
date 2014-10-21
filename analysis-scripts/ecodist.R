
library(ecodist)

# edist

# jdist

jdist = distance(iris[,1:4], method = 'jaccard')

iris.bc2 <- distance(iris[,1:4], "bray-curtis")

# spec


# 

mantel(jdist~edist)

mgram(jdist, edist, breaks, nclass, stepsize, nperm = 1000,
      mrank = FALSE, nboot = 500, pboot = 0.9, cboot = 0.95,
      alternative = "two.sided", trace = FALSE)

# generate a simple surface
x <- matrix(1:10, nrow=10, ncol=10, byrow=FALSE)
y <- matrix(1:10, nrow=10, ncol=10, byrow=TRUE)
z <- x + 3*y
image(z)

# analyze the pattern of z across space

space <- cbind(as.vector(x), as.vector(y))
z <- as.vector(z)
space.d <- distance(space, "eucl")
z.d <- distance(z, "eucl")
z.mgram <- mgram(z.d, space.d, nperm=0)
plot(z.mgram)


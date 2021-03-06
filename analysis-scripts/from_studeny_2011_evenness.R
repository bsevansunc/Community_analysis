##  Supplementary material to "Goodness-of-fit measures of evenness: a new tool 
##  for exploring changes in community structure"
##  R code for simulations of Tokeshi's niche models

# The following R functions provide a discrete algorithm to generate a
# realisation of Tokeshi's niche models with a fixed number of species S
# and total abundance N. The latter also determines the total available
# niche space. The functions have to be given the model specifications
# (S, N and k for the power fraction model) as well as one initial
# value b (the first niche division point); b can also be generated randomly.

#### Power fraction model ###################

powerfrac <- function(b,S,N,k)
  ### sampling function for the power fraction model
{
  while(length(b)< S-1)
  {
    p <- c(1,b,N)
    L <- length(p)
    dis <- p[2:L]-p[1:(L-1)]
    if(k==0)
    {
      ind <- which(dis==1)
      if(length(ind)>0){
        dis[ind] <- 0
        dis[-ind] <- 1}else{
          dis <- rep(1, times=length(dis))}
      v <- sample(c(b,N), 1, prob=dis)
    }else{
      dis <- (dis-1)^k
      v <- sample(c(b,N),1,prob=dis)
      
      ## here a random point is chosen, the next breaking point
      ## will lie in between this point and the next lower breaking point
      
    }
    ## the probabiliy of choosing a point is dependant on its
    ## distance from its lower neighbour
    ## the -1 was included to exclude points which are direct
    ## neighbours and don't allow to choose another point in between
    
    index <- which(p==v)
    a1 <- p[index-1]+1
    a2 <- p[index]-1
    
    if(a1==a2){b.new <- a1}else{
      b.new <- sample(c(a1:a2),1)}
    ## here the new point is randomly chosen
    b <- sort(c(b, b.new))
    ## and included in the set of breaking points
  }
  res <- sort(c(b,N)-c(0,b), decreasing=T)
  return(res)
}


########## Dominance preemption ################################################

dompreem <- function(b,S,N){
  ## function which gives discrete equivalent for the dominence preemption model
  ## because of the algorithm it is less unevenness than would be the case of the
  ## original model
  
  if(b>N-S){
    return("initial b too big")}else{
      while (length(b)< S-1){
        if(N-max(b)>S-length(b)+1){
          k <- runif(1,0.5,1)
          b.new <- floor(k*(N-S+length(b)-max(b)))+max(b)
          b <- c(b,b.new)
          res <- sort(c(b,N)-c(0,b), decreasing=T)
        }else{
          b.rem <- rep(1, times=S-length(b))
          b <- c(b,b.rem)
          res <- b}
      }
      return(res)}
}

############### Dominance decay ################################################

domdecay <- function(b,S,N){
  res <- c(b,N)-c(0,b)
  while(length(b)<S-1){
    d <- c(0,b,N)
    ind <- which(res==max(res))
    if(length(ind)==1){
      b.new <- sample((d[ind]+1):(d[ind+1]-1),1)
      b <- sort(c(b,b.new), decreasing=F)
      res <- c(b,N)-c(0,b)
    }else{
      ind <- sample(ind,1)
      b.new <- sample((d[ind]+1):(d[ind+1]-1),1)
      b <- sort(c(b,b.new), decreasing=F)
      res <- c(b,N)-c(0,b)
    }
  }
  res <- sort(res, decreasing=T)
  return(res)
}


# The following R code generates 300 random realisations of Tokeshi's power
# fraction model for various parameter values k (analogously the dominance
# decay and the dominance preemption models can be simulated using the R
# functions above). A representation of each power fraction model is derived
# by averaging over the 300 random outcomes:

k <- c(0,0.4,0.7,1)

N <- 50000
S <- 100

set.seed(12323)

b.pow <- sample(2:(N-1),300*length(k), replace=T)
### initial b values for the 500 samples for each k

samples.pow <- array(0,dim=c(300,S,length(k)))

seeds.pow <- sample(1:100000, length(k), replace=T)

for(j in 1:length(k))
{
  set.seed(seeds.pow[j])
  for(i in 1:nrow(samples.pow))
  {
    samples.pow[i,,j] <- powerfrac(b.pow[i*j],S,N,k[j])
  }
}

models.pow <- matrix(0, nrow=length(k), ncol=S)
### matrix containing the k SADs generated by averaging over the samples

for(j in 1:length(k))
{
  models.pow[j,] <- apply(samples.pow[,,j],2,mean)
}
row.names(models.pow) <- paste("k=", k)

abd.pow <- models.pow
rel.abd.pow <- models.pow/N      ### SADs from the power-fraction models


### Simulations for the dominance preemption model #############################

samples.preemp <- matrix(0, 300, S)
set.seed(12345)

b.preemp <- ceiling(runif(300,0.5,1)*(N-S+1))

for(i in 1:300){
  samples.preemp[i,] <- dompreem(b.preemp[i],S,N)
}

abd.preemp <- apply(samples.preemp, 2, mean)

rel.abd.preemp <- abd.preemp/N

######## Simulations for the dominance decay model #############################

samples.dec <- matrix(0, 300, S)
set.seed(1356)

b.dec <- sample(1:49999,300,replace=T)

for(i in 1:300){
  samples.dec[i,] <- domdecay(b.dec[i],S,N)
}

abd.dec <- apply(samples.dec, 2, mean)

rel.abd.dec <- abd.dec/N

######## SADs from the simulated Tokeshi models ################################

abd <- rbind(abd.preemp, abd.pow, abd.dec)
row.names(abd)<- c("abd.preemp", "k=0", "k=0.4", "k=0.7", "k=1", "abd.dec")

rel.abd <- rbind(rel.abd.preemp, rel.abd.pow, rel.abd.dec)
row.names(rel.abd)<- c("abd.preemp", "k=0", "k=0.4", "k=0.7", "k=1", "abd.dec")


######### Goodness-of-fit measures ##############################################
######### (for scenario 1 as discussed in the paper) ############################
#################################################################################

gof <- function(lambda,S, rel.abd){
  res <- rep(0,times=length(lambda))
  
  for(i in 1:length(lambda)){
    res[i] <- 1/(lambda[i]*(lambda[i]+1))*sum(rel.abd*((rel.abd*S)^lambda[i]-1))
  }
  
  return(res)
}

### Set range of parameter in Gof from -5 to 5
lambda <- seq(-2.99, 3, by=0.02)

### calculate GoF-measure for -5 < lambda < 5 for the different SADs in abd

gof.Tok <- t(apply(rel.abd, 1, gof, lambda=lambda, S=S))

######### Evenness profile plots for Tokeshi's models ###########################
######### (plotted on a logarithmic scale #######################################
palette(rainbow(10))

plot(lambda,log10(gof.Tok[1,]), ylim=range(log10(gof.Tok)), xlim=c(-3,3), yaxt="n", type="l",
     col=1, ylab="Index", xlab=expression(lambda))
axis(2, at=c(-1,0,1,2,3,4), labels=c(0.1,1,10,100,1000,10000))
for(i in 2:nrow(abd)){
  lines(lambda,log10(gof.Tok[i,]), ylim=c(-1,3), col=i)
}

########### Standardized evenness index based on goodness-of-fit ################

maxi <- matrix(1/(lambda*(1+lambda))*(S^lambda -1), nrow=nrow(gof.Tok), ncol=length(lambda), byrow=T)

gof.stand.Tok <- 1- gof.Tok/maxi


plot(lambda[lambda>=-1],gof.stand.Tok[1,lambda>=-1], ylim=c(0,1), xlim=c(-1,1),
     type="l",col=1, ylab="Index", xlab=expression(lambda))
for(i in 2:nrow(abd)){
  lines(lambda[lambda>=-1], gof.stand.Tok[i,lambda>=-1], col=i)
}

################################################################################

################################################################################
######## Sampling properties of the gof measures ###############################
######## (scenario 2) ##########################################################

####### power fraction model with k=0, 0.4 and 1 ###############################

abd.pow.0 <- abd[2,]
rel.abd.pow.0 <- rel.abd[2,]

abd.pow.0.4   <- abd[3,]
rel.abd.pow.0.4 <- rel.abd[3,]

abd.pow.1   <- abd[5,]
rel.abd.pow.1 <- rel.abd[5,]

rel.abd.pow.samp <- rbind(rel.abd.pow.0, rel.abd.pow.0.4, rel.abd.pow.1)

###### gof measures for true SADs ##############################################

gof.Tok <- gof(rel.abd.pow.0.4, lambda=lambda, S=100)

##### Hill's numbers ###########################################################

hill <- function(beta,rel.abd){
  res <- rep(0,times=length(beta))
  
  for(i in 1:length(beta)){
    res[i] <- (sum(rel.abd^beta[i]))^(1/(1-beta[i]))
  }
  return(res)
}

beta <- seq(0,5,by=0.03)

Hill <- hill(rel.abd.pow.0.4, beta=beta)
Hill.2 <- hill(rel.abd.pow.samp,  beta=2)


#################################################################################
############### sample for power fraction model with k=0 ########################



set.seed(2) #for pow frac 0.4
#set.seed(1) #for pow frac 0
#set.seed(19) #for pow frac 1

obs <- t(rmultinom(1, 500, abd.pow.0.4))

S_obs <- sum(obs!=0)

obs.rank <- sort(obs, decreasing=T)

S1 <- 80
S2 <- 90
S3 <- 100
S4 <- 150
S5 <- 200

S <- c(S_obs, S1, S2, S3, S4, S5)  ### assumed number of species

#eps <- 0.02     ### for k=0
eps <- 0.08    ### for k=0.4
#eps <- 0.2     ### for k=1

Gof.Tok.samp <- vector("list",length(S))
Hill.samp  <- vector("list",length(S))
Hill.samp.2  <- vector("list",length(S))

for(l in 1:length(S)){
  #for(l in 1:length(eps)){
  
  
  if(S[l]<=100){
    x <- obs.rank[1:S[l]]
  }else{
    x <- c(obs.rank, rep(0, times=S[l]-100))
  }
  
  resamp <- x
  
  resamp[resamp==0] <- eps
  rsums <- sum(resamp)
  rel.obs <- resamp/rsums
  
  
  
  gof.obs.S <- gof(lambda, S[l], rel.obs)
  
  hill.obs.S<- hill(beta, rel.obs)
  
  Hill.samp.2[[l]] <- hill(2, rel.obs)
  Hill.samp[[l]] <- hill.obs.S
  Gof.Tok.samp[[l]] <- gof.obs.S
  
}

#################################################################################
######### plots for scenario 2 (for k=0.4) #################################################
palette(rainbow(12, alpha=.7))

### evenness profiles 

plot(lambda, log10(gof.Tok),ylim = c(-2,3),xlim=c(-3,3), type="l", lty=4 ,yaxt="n", xlab=expression(lambda), ylab="Index" ) 
axis(2, at=c(-1,0,1,2,3), labels=c(0.1,1,10,100,1000))

for(l in 1:length(S)){
  lines(lambda, log10(Gof.Tok.samp[[l]]), col=2*l)
}

legend("bottomright", c("true profile", paste("S=",S_obs, "[obs]"), "S=80","S=90", "S=100 [true]", "S=150", "S=200"), col=c("black",seq(2,12,2)), lty=c(4,rep(1,6)))



#### Hill's diversity numbers ##################################################

plot(beta, Hill, ylim=c(0,200), xlim=c(0,5),  type="l", lty=4,  xlab="a", ylab="Index")

for(l in 1:length(S)){
  lines(beta, Hill.samp[[l]], col=2*l)
}
legend("topright", c("true profile", paste("S=",S_obs, "[obs]"), "S=80","S=90", "S=100 [true]", "S=150", "S=200"), col=c("black",seq(2,12,2)), lty=c(4,rep(1,6)))

### Hill's scaled diversity numbers (evenness) #################################

### scaled by N(0)=S ###########################################################

plot(beta, (Hill/Hill[1])^2, ylim=c(0,1), xlim=c(0,5), xlab="a", ylab="Index", type="l", lty=4)

for(l in 1:length(S)){
  lines(beta, (Hill.samp[[l]]/Hill.samp[[l]][1])^2, col=2*l)
}

legend("topright", c("true profile", paste("S=",S_obs, "[obs]"), "S=80","S=90", "S=100 [true]", "S=150", "S=200"), col=c("black",seq(2,12,2)), lty=c(4,rep(1,6)))


### scaled by N(2) #############################################################

plot(beta, (Hill/Hill.2)^2, ylim=c(0,1), xlim=c(2,5), xlab="a", ylab="Index", type="l", lty=4)

for(l in 1:length(S)){
  lines(beta, (Hill.samp[[l]]/Hill.samp.2[[l]])^2, col=2*l)
}

legend("topright", c("true profile", paste("S=",S_obs, "[obs]"), "S=80","S=90", "S=100 [true]", "S=150", "S=200"), col=c("black",seq(2,12,2)), lty=c(4,rep(1,6)))


### Relative and relative logarithmic evenness #################################

plot(beta, (Hill-1)/(100-1), xlim=c(0,5), ylim=c(0,1), xlab="a", ylab="Index", type = "l", lty=4)        

for(l in 1:length(S)){
  lines(beta, (Hill.samp[[l]]-1)/(S[l]-1), col=2*l)
}
legend("topright", c("true profile", paste("S=",S_obs, "[obs]"), "S=80","S=90", "S=100 [true]", "S=150", "S=200"), col=c("black",seq(2,12,2)), lty=c(4,rep(1,6)))


plot(beta, log(Hill)/log(100), xlim=c(0,5), ylim=c(0,1),  xlab="a", ylab="Index",type = "l", lty=4)        
for(l in 1:length(S)){
  lines(beta, log(Hill.samp[[l]])/log(S[l]), col=2*l)
}
legend("bottomright", c("true profile", paste("S=",S_obs, "[obs]"), "S=80","S=90", "S=100 [true]", "S=150", "S=200"), col=c("black",seq(2,12,2)), lty=c(4,rep(1,6)))

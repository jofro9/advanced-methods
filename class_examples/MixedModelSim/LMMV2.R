####################################################
## Purpose: Generate longitudinal data determined ##
##          by linear mixed effects models        ##
## By: Kevin Josey                                ##
####################################################

library(lme4)
library(MASS)

###Consider the matrix framework for a single subject/cluster in a linear mixed effects model
### Y_i = X_iB+Z_ib_i+E_i, where b_i~MVN(0,G) and E_i~MVN(0,R)
###This code generates a simulation of Y_i, the vector of observations on a subject/cluster

##Step 1: Create the G matrix from a correlation matrix
# we will need gmat() to generate covariance matrix for random effects
# we may also want to use un() from the MultiMod.R script
# sig2 is the variance(s) of the random effects.  It can be a vector (see example below).
#rho is the correlation between the random effects.

gmat <- function(sig2, rho, p = length(sig2)){
  
  if(length(sig2) == 1)
    sig2 <- rep(sig2, p)
  else if (length(sig2) != p)
    stop("length(sig2) != p")
  
  if(length(rho) != 1)
    stop("length(rho) != 1")
  
  R <- matrix(rho, p, p)  ##This is not R in the G and R notation in the mixed effects model above.  It is a generic matrix used for making G.
  diag(R) <- 1
  D <- diag(sqrt(sig2), p, p)
  G <- D %*% R %*% D  ##This does the matrix math of the relationship between the covariance/variance and the correlation.  
  return(G)
  
}

# data dimensions
n1 <- 1000 # sample size for group 1
n2 <- 1000 # sample size for group 2
n <- n1 + n2 # total sample size
p <- 4 # number of repeated measurements

# design matrix
grp <- factor( rep(1:2, times = c(n1*p,n2*p)) ) # identify treatment assignment
time <- rep(0:(p-1), times = n) # time variable
X <- model.matrix(~ grp*time) # reference cell coding
id <- rep(1:n, each = p) # id indexing clusters

# fixed effect coefs
beta <- c(1, -1, 2, 0.5) 
# c(int, grp2, time, grp2*time)

# generate random correlated random effects
Sig0 <- gmat(sig2 = c(1,2), rho = 0.2, p = 2) #This is G in the notation above for mixed effects models.
b <- mvrnorm(n, mu = c(0,0), Sigma = Sig0)  #This is b in the notation above for mixed effects models.
b0 <- rep(d[,1], each = p)
b1 <- rep(d[,2], each = p)
e <- rnorm(n*p, 0, sqrt(4)) #These are your error terms.  You set 4, which is sigma^2_e
# note that the errors are independent
# between cluster heterogeneity is expressed using random effects

y_n <- c( X %*% beta + b0 + b1*time + e ) # generate outcomes

# fit model
rand.fit <- lmer(y_n ~ grp*time + (1 + time|id)) # using lme4 package
summary(rand.fit)

###########################################################
## Purpose: Generate binary longitudinal data determined ##
##          by generalized linear mixed effects models   ##
## By: Kevin Josey                                       ##
###########################################################

library(lme4)
library(MASS)

# we will need gmat() to generate covariance matrix for random effects
# we may also want to use un() from the MultiMod.R script
gmat <- function(sig2, rho, p = length(sig2)){
  
  if(length(sig2) == 1)
    sig2 <- rep(sig2, p)
  else if (length(sig2) != p)
    stop("length(sig2) != p")
  
  if(length(rho) != 1)
    stop("length(rho) != 1")
  
  R <- matrix(rho, p, p)
  diag(R) <- 1
  D <- diag(sqrt(sig2), p, p)
  V <- D %*% R %*% D
  return(V)
  
}

# data dimensions
n1 <- 100 # sample size for group 1
n2 <- 100 # sample size for group 2
n <- n1 + n2 # total sample size
p <- 4 # number of repeated measurements

# design matrix
grp <- factor( rep(1:2, times = c(n1*p,n2*p)) ) # identify treatment assignment
time <- rep(0:(p-1), times = n) # time variable
X <- model.matrix(~ grp*time) # reference cell coding
id <- rep(1:n, each = p) # id indexing clusters

# fixed effect coefs
beta <- c(1, 1, -2, -0.5) 
# c(int, grp2, time, grp2*time)

# generate correlated random effects for each ISU
Sig0 <- gmat(sig2 = c(1,2), rho = 0.2, p = 2)
d <- mvrnorm(n, mu = c(0,0), Sigma = Sig0)
d0 <- rep(d[,1], each = p)
d1 <- rep(d[,2], each = p)
# there is no error component in GzLMM

eta <- c( X %*% beta + d0 + d1*time ) # linear predictor
theta <- plogis(eta) # inverse logit function
y_b <- rbinom(n*p, 1, theta) # generate outcomes

# fit model
rand.logit.fit <- glmer(y_b ~ grp*time + (1 + time|id), family = binomial(link = "logit")) # using lme4 package
summary(rand.logit.fit)

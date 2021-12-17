############################################################
## Program Name: PrelimModels_Testing_Video7.R         	  ##
## Purpose: Investigate independent effects of PSA for    ##
##          age and computing DIC and posterior probabilities## 
## Created by: Nichole Carlson                             ##
############################################################


# Dependencies
library(readr)
library(rjags)
library(mcmcse)
library(tidyverse)

setwd("/Users/nichole/Repositories/BIOS6624/Bios6624ClassExamples/PSAExample/")
psaclean2 <- read_csv("DataProcessed/psaclean2.csv")

###########################################################################
#################LOG PSA AND AGE ANALYSIS##################################
###########################################################################

## This script will work for any linear model we wish to fit with JAGS
## In addition, we can use this script to generate models for the univariable interaction hypothesis

# y is the outcome in your linear regression model
y <- c(psaclean2$lpsa)

# Now we set up the design matrix in our linear regression
# For the first question in worksheet 3 we are looking at an intercept and age
X <- model.matrix(~ age, data = psaclean2) 

#These are variables needed in Rjags
#N is the number of observations in your dataset and p is the number of columns in the design matrix
N <- nrow(X)
p <- ncol(X)

## Hyperparameters for prior distributions. These change depending on your assumptions
## For this problem the model error variance is an inverse gamma and the coefficients
## are normal distribution with a mean and varinace, which can be written flexibly as a multivariate normal
## with a mean vector (m below) and a variance-covariance (R) below.  If we are going to model each coefficient
## independently, which would be typical.  We would have an independence matrix as is programmed below.
## Note jags works with precision rather than variance so we want to program the precision matrix rather than
## the var-covariance matrix.  So, really small values in this matrix are like being non-informative and having "large"
## variances in a prior.

a <- 0.001 # gamma shape for a non-informative prior
b <- 0.001 # gamma rate for a non-informative prior
m <- c(0.5, rep(0, (p-1))) # mvnorm mean  (mean in the prior on the regression coefficients)
R <- matrix(0, p, p) # mvnorm covariance (1/var-cov in the prior on the regression coefficients)
diag(R) <- 0.000001 # note that JAGS uses dispersion matrix (scalars) rather 

# create data list to pass to JAGS
jags_dat <- list(y = y, X = X, N = N, p = p,
                 a = a, b = b, m = m, R = R)

# Initialize the sampler (it's a Gibbs sampler so no need for several adaptations)
# In this code the linMod.jags is what programs up the Bayesian model formulation including the likelihood
# and the priors.  Go take a look at that now!

mod <- jags.model("Code/jags/linMod.jags", data = jags_dat, n.adapt = 1000, n.chains = 2)

# Sample observations from the posterior distributions
iter <- 25000 # number of draws for the MCMC
samples <- coda.samples(mod, variable.names = c("beta", "sigma2"), n.iter = iter) # generates matrix of draws from posterior dist.

## HPDI function which calculates credible intervals using highest posterior densities

hpd <- function(x, alpha = 0.05){
 
  n = length(x)
  m = round(n * alpha)
  x = sort(x)
  y = x[(n - m + 1):n] - x[1:m]
  z = min(y)
  k = which(y == z)[1]
  c(x[k], x[n - m + k])
  
}


## Generate output table when the final model has been selected to summarize findings

draws <- as.matrix(samples[[1]])
draws <- data.frame(draws[5001:iter,])
## Generate the boundaries of the posterior probabilities that I would like for beta1
## 5% boundaries: need log(1.05) since I am modeling on the log scale.  lower bound log(0.95).

nsample<- nrow(draws) 

inbounds05 <- draws %>% filter(beta.2.<log(1.05) & beta.2.>log(0.95))
postprob05<- 1 - nrow(inbounds05)/nsample
print("postprob05")
print(postprob05)

inbounds10 <- draws %>% filter(beta.2.<log(1.10) & beta.2.>log(0.90))
postprob10<- 1- nrow(inbounds10)/nsample
print("postprob10")
print(postprob10)


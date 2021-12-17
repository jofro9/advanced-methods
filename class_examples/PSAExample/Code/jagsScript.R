###########################################################
## Program Name: ModelBuilding.R      		               ##
## Purpose: Investigate independent effects of PSA.      ## 
##		        through multivariable regression 	         ##
## Created by: Nichole Carlson	(R code by: Kevin Josey) ##
###########################################################

# Dependencies
library(readr)
library(rjags)
library(mcmcse)

setwd("~/Github/Bios6624ClassExamples/PSAExample/Code/")
psaclean2 <- read_csv("~/Github/Bios6624ClassExamples/PSAExample/DataProcessed/psaclean2.csv")

iter <- 10000 # number of draws

## This script will work for any linear model we wish to fit with JAGS
## In addition, we can use this script to generate models for the univariable interaction hypothesis

y <- c(psaclean2$lpsa)
# change the formula to test the different models
X <- model.matrix(~ cavol + wt + age + bph + cappen + svi + grade6 + grade7, data = psaclean2) 
N <- nrow(X)
p <- ncol(X)

## Hyperparameters for prior distributions

a <- 2 # invgamma shape
b <- 1 # invgamma rate
m <- rep(0, p) # mvnorm mean
R <- matrix(0, p, p) # mvnorm covariance
diag(R) <- 0.001 # note that JAGS uses dispersion matrix (scalars)

# create data list to pass to JAGS
jags_dat <- list(y = y, X = X, N = N, p = p,
                 a = a, b = b, m = m, R = R)

# Initialize the sampler (it's a Gibbs sampler so no need for several adaptations)

mod <- jags.model("jags/linMod.jags", data = jags_dat, n.adapt = 1000, n.chains = 2)

# Sample observations from the posterior distributions

samples <- coda.samples(mod, variable.names = c("beta", "sigma2"), n.iter = iter) # generates matrix of draws from posterior dist.
samples_dic <- dic.samples(mod, n.iter = iter, type = "pD") # determine DIC of model fit (requires two chains)

## Diagnostics

# determine if the support of the parameter space is explored well (more adaptations)
plot(samples)

# determine amount of autocorrelation between draws (thinning parameter)
par(mfrow = c(ceiling(ncol(as.matrix(samples))/2), 2), mar = rep(1, 4))
apply(as.matrix(samples), 2, acf) # again the Gibb's sampler should have almost no autocorrelation

## HPDI function

hpd <- function(x, alpha = 0.05){
  
  n = length(x)
  m = round(n * alpha)
  x = sort(x)
  y = x[(n - m + 1):n] - x[1:m]
  z = min(y)
  k = which(y == z)[1]
  c(x[k], x[n - m + k])
  
}


## Generate output table when the final model has been selected

draws <- as.matrix(samples)

out_mat <- matrix("", nrow = ncol(draws), ncol = 5)
out_mat[,1] <- c(colnames(X), "sigma2") # names
out_mat[,2] <- apply(draws, 2, function(x) round(mcse(x)$est, 3)) # batch mean
out_mat[,3] <- apply(draws, 2, function(x) round(mcse(x)$se, 3)) # MCSE
out_mat[,4] <- round(apply(draws, 2, sd), 3) # Std. Dev.
out_mat[,5] <- apply(draws, 2, function(x)
  paste("(", round(hpd(x)[1], 3), ", ", round(hpd(x)[2], 3), ")", sep = "")) # HPDI

colnames(out_mat) <- c("Parameter", "Estimate", "MCSE", "Std. Dev.", "95% HPDI")

write_csv(as.data.frame(out_mat), "~/Github/Bios6624ClassExamples/PSAExample/Output/bayes-est.csv")

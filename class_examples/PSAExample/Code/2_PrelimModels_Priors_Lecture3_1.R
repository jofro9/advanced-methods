############################################################
## Program Name: PrelimModels_Priors_Lecture3_1.R      	  ##
## Purpose: Investigate independent effects of PSA for    ##
##           each variable separately and different priors## 
## Created by: Nichole Carlson and Kevin Josey            ##
############################################################

###INSTALL JAGS, WHICH IS A STAND ALONE PROGRAM
##For Macs go get homebrew at https://brew.sh/
### Then at your command line type "brew install jags"

##FOR ALL THIS PDF MIGHT BE OF USE
###https://faculty.washington.edu/jmiyamot/p548/installing.jags.pdf


# Dependencies
library(readr)
library(rjags)
library(mcmcse)
library(tidyverse)

setwd("/Users/nichole/Repositories/BIOS6624/Bios6624ClassExamples/PSAExample/")
psaclean2 <- read.csv("DataProcessed/psaclean2.csv")

###########################################################################
#################LOG PSA AND AGE ANALYSIS##################################
###########################################################################

####Run a linear model just to orient us using standard methods
agemodel <-lm(lpsa ~ age, data=psaclean2)

sink("Output/AgeLinearModelOutput.txt")
print(agemodel)
summary.lm(agemodel)
confint(agemodel)
sink()
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

draws <- as.matrix(samples)

out_mat <- matrix("", nrow = ncol(draws), ncol = 5)
out_mat[,1] <- c(colnames(X), "sigma2") # names
out_mat[,2] <- apply(draws, 2, function(x) round(mcse(x)$est, 3)) # batch mean
out_mat[,3] <- apply(draws, 2, function(x) round(mcse(x)$se, 3)) # MCSE
out_mat[,4] <- round(apply(draws, 2, sd), 3) # Std. Dev.
out_mat[,5] <- apply(draws, 2, function(x)
  paste("(", round(hpd(x)[1], 3), ", ", round(hpd(x)[2], 3), ")", sep = "")) # HPDI

colnames(out_mat) <- c("Parameter", "Estimate", "MCSE", "Std. Dev.", "95% HPDI")

write_csv(as.data.frame(out_mat), "Output/Priors_Video4/bayes-est-prior1-agemodel.csv")

# determine if the support of the parameter space is explored well (more adaptations)
plot(samples)

####additional code for checking MCMC model assumptions using ggplot
# acf plots
par(mfrow = c(ceiling(ncol(as.matrix(samples))/2), 2), mar = rep(1, 4))
apply(as.matrix(samples), 2, acf) # again the Gibb's sampler should have almost no autocorrelation


chain1 = samples[[1]]
chain2 = samples[[2]]

# make trace plots
colnames(chain1) = c(paste0("beta", 0:1), "sigma2")
colnames(chain2) = c(paste0("beta", 0:1), "sigma2")

chain1 = chain1 %>% as.matrix() %>% as_tibble() %>% mutate(iteration = row_number(), chain = 1)
chain2 = chain2 %>% as.matrix() %>% as_tibble() %>% mutate(iteration = row_number(), chain = 2)

samples = rbind(chain1, chain2) %>% mutate(chain = factor(chain)) 


# after the first 1000 iterations
samples %>% 
  filter(iteration >= 1000) %>%
  gather(parameter, value, beta0:sigma2) %>%
  ggplot(aes(iteration, value, group = chain, color = chain)) + geom_line(alpha = 0.4) +
  facet_wrap(~parameter, scales = "free", ncol = 5)


# first 100 iterations
# what happens to this plot if you have bad starting values?
samples %>% 
  filter(iteration < 50) %>%
  gather(parameter, value, beta0:sigma2) %>%
  ggplot(aes(iteration, value, group = chain, color = chain)) + 
  geom_line(alpha = 0.4) +
  facet_wrap(~parameter, scales = "free", ncol = 5)



#### density plots
samples %>% 
  filter(iteration >= 1000) %>%
  gather(parameter, value, beta0:sigma2) %>%
  ggplot(aes(x = value, group = chain, color = chain)) + 
  geom_density() +
  facet_wrap(~parameter, scales = "free", ncol = 5)


#### ACF plots
# sigma2
bacf <- acf(filter(samples, iteration > 1000, chain == 1)$sigma2, plot = FALSE)
bacfdf <- with(bacf, data.frame(lag, acf))

ggplot(data = bacfdf, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = lag, yend = 0))


# beta1
bacf <- acf(filter(samples, iteration > 1000, chain == 1)$beta1, plot = FALSE)
bacfdf <- with(bacf, data.frame(lag, acf))

ggplot(data = bacfdf, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = lag, yend = 0))




########################################################################################
########################################################################################
############################Vague priors################################################
a <- 0.001 # gamma shape for a non-informative prior
b <- 0.001 # gamma rate for a non-informative prior
m <- c(0.5, rep(0, (p-1))) # mvnorm mean  (mean in the prior on the regression coefficients)
R <- matrix(0, p, p) # mvnorm covariance (1/var-cov in the prior on the regression coefficients)
diag(R) <- 0.0001 # note that JAGS uses dispersion matrix (scalars) rather 

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


## Generate output table when the final model has been selected to summarize findings

draws <- as.matrix(samples)

out_mat <- matrix("", nrow = ncol(draws), ncol = 5)
out_mat[,1] <- c(colnames(X), "sigma2") # names
out_mat[,2] <- apply(draws, 2, function(x) round(mcse(x)$est, 3)) # batch mean
out_mat[,3] <- apply(draws, 2, function(x) round(mcse(x)$se, 3)) # MCSE
out_mat[,4] <- round(apply(draws, 2, sd), 3) # Std. Dev.
out_mat[,5] <- apply(draws, 2, function(x)
  paste("(", round(hpd(x)[1], 3), ", ", round(hpd(x)[2], 3), ")", sep = "")) # HPDI

colnames(out_mat) <- c("Parameter", "Estimate", "MCSE", "Std. Dev.", "95% HPDI")

write_csv(as.data.frame(out_mat), "Output/Priors_Video4/bayes-est-prior2-agemodel.csv")

# determine if the support of the parameter space is explored well (more adaptations)
plot(samples)

########################################################################################
########################################################################################
############################Informative/optimistic priors################################################
a <- 0.001 # gamma shape for a non-informative prior
b <- 0.001 # gamma rate for a non-informative prior
m <- c(0.5, rep(10, (p-1))) # mvnorm mean  (mean in the prior on the regression coefficients)
R <- matrix(0, p, p) # mvnorm covariance (1/var-cov in the prior on the regression coefficients)
diag(R) <- 0.01 # note that JAGS uses dispersion matrix (scalars) rather 

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


## Generate output table when the final model has been selected to summarize findings

draws <- as.matrix(samples)

out_mat <- matrix("", nrow = ncol(draws), ncol = 5)
out_mat[,1] <- c(colnames(X), "sigma2") # names
out_mat[,2] <- apply(draws, 2, function(x) round(mcse(x)$est, 3)) # batch mean
out_mat[,3] <- apply(draws, 2, function(x) round(mcse(x)$se, 3)) # MCSE
out_mat[,4] <- round(apply(draws, 2, sd), 3) # Std. Dev.
out_mat[,5] <- apply(draws, 2, function(x)
  paste("(", round(hpd(x)[1], 3), ", ", round(hpd(x)[2], 3), ")", sep = "")) # HPDI

colnames(out_mat) <- c("Parameter", "Estimate", "MCSE", "Std. Dev.", "95% HPDI")

write_csv(as.data.frame(out_mat), "Output/Priors_Video4/bayes-est-prior4-agemodel.csv")

# determine if the support of the parameter space is explored well (more adaptations)
plot(samples)

########################################################################################
########################################################################################
############################Extreme Priors ################################################
a <- 0.001 # gamma shape for a non-informative prior
b <- 0.001 # gamma rate for a non-informative prior
m <- c(10, rep(10, (p-1))) # mvnorm mean  (mean in the prior on the regression coefficients)
R <- matrix(0, p, p) # mvnorm covariance (1/var-cov in the prior on the regression coefficients)
diag(R) <- 1 # note that JAGS uses dispersion matrix (scalars) rather 

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


## Generate output table when the final model has been selected to summarize findings

draws <- as.matrix(samples)

out_mat <- matrix("", nrow = ncol(draws), ncol = 5)
out_mat[,1] <- c(colnames(X), "sigma2") # names
out_mat[,2] <- apply(draws, 2, function(x) round(mcse(x)$est, 3)) # batch mean
out_mat[,3] <- apply(draws, 2, function(x) round(mcse(x)$se, 3)) # MCSE
out_mat[,4] <- round(apply(draws, 2, sd), 3) # Std. Dev.
out_mat[,5] <- apply(draws, 2, function(x)
  paste("(", round(hpd(x)[1], 3), ", ", round(hpd(x)[2], 3), ")", sep = "")) # HPDI

colnames(out_mat) <- c("Parameter", "Estimate", "MCSE", "Std. Dev.", "95% HPDI")

write_csv(as.data.frame(out_mat), "Output/Priors_Video4/bayes-est-prior5-agemodel.csv")

# determine if the support of the parameter space is explored well (more adaptations)
plot(samples)
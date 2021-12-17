###########################################################
## Program Name: DICModels_LectureWeek4.R      		       ##
## Purpose: Investigate independent effects of PSA for   ##
##           each variable separately and compare models with DIC.## 
## Created by: Nichole Carlson and Kevin Josey           ##
###########################################################

# Dependencies
library(readr)
library(rjags)
library(mcmcse)

setwd("/Users/nichole/Repositories/BIOS6624/Bios6624ClassExamples/PSAExample/")
psaclean2 <- read_csv("DataProcessed/psaclean2.csv")

###########################################################################
#################LOG PSA AND AGE ANALYSIS##################################
###########################################################################

##Model with Age in the analysis
# y is the outcome in your linear regression model
y <- c(psaclean2$lpsa)

# Now we set up the design matrix in our linear regression
# For the first question in worksheet 3 we are looking at an intercept and age
X <- model.matrix(~ age, data = psaclean2) 

#These are variables needed in Rjags
#N is the number of observations in your dataset and p is the number of columns in the design matrix
N <- nrow(X)
p <- ncol(X)

a <- 0.001 # gamma shape for a non-informative prior
b <- 0.001 # gamma rate for a non-informative prior
m <- rep(0, p) # mvnorm mean  (mean in the prior on the regression coefficients)
R <- matrix(0, p, p) # mvnorm covariance (1/var-cov in the prior on the regression coefficients)
diag(R) <- 0.001 # note that JAGS uses dispersion matrix (scalars) rather 

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

##To be used later
samples_dic <- dic.samples(mod, n.iter = iter, type = "pD") # determine DIC of model fit (requires two chains)

###########################################################
#########Now fit an intercept only model for comparison####
###########################################################
# Now we set up the design matrix in our linear regression
#Intercept only model
Xint <- as.matrix(rep(1,nrow(psaclean2)) )

N <- nrow(Xint)
p <- ncol(Xint)

a <- 0.001 # gamma shape for a non-informative prior
b <- 0.001 # gamma rate for a non-informative prior
m <- rep(0, p) # mvnorm mean  (mean in the prior on the regression coefficients)
R <- matrix(0, p, p) # mvnorm covariance (1/var-cov in the prior on the regression coefficients)
diag(R) <- 0.001 # note that JAGS uses dispersion matrix (scalars) rather 

# create data list to pass to JAGS
jags_dat_int <- list(y = y, X = Xint, N = N, p = p,
                 a = a, b = b, m = m, R = R)

# Initialize the sampler (it's a Gibbs sampler so no need for several adaptations)
# In this code the linMod.jags is what programs up the Bayesian model formulation including the likelihood
# and the priors.  Go take a look at that now!

mod.int <- jags.model("Code/jags/linMod.jags", data = jags_dat_int, n.adapt = 1000, n.chains = 2)

# Sample observations from the posterior distributions
iter <- 25000 # number of draws for the MCMC
samples.int <- coda.samples(mod.int, variable.names = c("beta", "sigma2"), n.iter = iter) # generates matrix of draws from posterior dist.

##To be used later
samples_dic_int <- dic.samples(mod.int, n.iter = iter, type = "pD") # determine DIC of model fit (requires two chains)




###########################################################################
#################LOG PSA AND GLEASON SCORE ANALYSIS########################
###########################################################################
####Run a linear model just to orient us using standard methods
psaclean2$gleason<-as.factor(psaclean2$gleason)
gleasonmodel <-lm(lpsa ~ gleason, data=psaclean2)

sink("Output/GleasonLinearModelOutput.txt")
print(gleasonmodel)
summary.lm(gleasonmodel)
confint(gleasonmodel)
anova(gleasonmodel)
sink()

# y is the outcome in your linear regression model
y <- c(psaclean2$lpsa)

# Now we set up the design matrix in our linear regression
# For the first question in worksheet 3 we are looking at an intercept and age
X <- model.matrix(~ grade7 + grade8, data = psaclean2) 

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
## the var-covariance matrix.  So, really small values in this matrix are like being vague and having "large"
## variances in a prior.

a <- 0.001 # gamma shape for a non-informative prior
b <- 0.001 # gamma rate for a non-informative prior
m <- rep(0, p) # mvnorm mean  (mean in the prior on the regression coefficients)
R <- matrix(0, p, p) # mvnorm covariance (1/var-cov in the prior on the regression coefficients)
diag(R) <- 0.001 # note that JAGS uses dispersion matrix (scalars) rather 

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

##To be used later--don't run now
##samples_dic <- dic.samples(mod, n.iter = iter, type = "pD") # determine DIC of model fit (requires two chains)

## Diagnostics

# determine if the support of the parameter space is explored well (more adaptations)
plot(samples)

# determine amount of autocorrelation between draws (thinning parameter)
par(mfrow = c(ceiling(ncol(as.matrix(samples))/2), 2), mar = rep(1, 4))
apply(as.matrix(samples), 2, acf) # We will learn later how to assess if this graph says our MCMC worked correctly.

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

write_csv(as.data.frame(out_mat), "Output/bayes-est-worksheet3-gleasonmodel.csv")

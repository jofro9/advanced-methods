###########################################################
## Program Name: code.R                                  ##
## Purpose: Project 1, HAART study, investigating        ##
##          the effects of hard drug use                 ## 
## Created by: Nichole Carlson and Joe Froelicher        ##
###########################################################

# Packages
library(janitor)
library(ggplot2)
library(mcmcse)
library(naniar)
library(readr)
library(rjags)
library(tidyverse)

data = read.csv("../data_raw/hiv_6624_final.csv") %>% clean_names()

#################
# Data Cleaning #
#################

# take out variables not used for regression for simplicity
drop = c(
  'x', 'hashv', 'hashf', 'income', 'hbp',
  'diab', 'liv34', 'kid', 'frp', 'fp',
  'tchol', 'trig', 'ldl', 'dyslip', 'cesd',
  'dkgrp', 'heropiate', 'idu', 'hivpos', 'art',
  'ever_art'
)

data = data[,!(names(data) %in% drop)]

# set bmi values of 999 to NA and remove unrealistic values
for (i in 1:dim(data)[1]) {
  if (!is.na(data$bmi[i])) {
    if (data$bmi[i] == 999 || data$bmi[i] == -1 || data$bmi[i] > 100) {
      data$bmi[i] = NA
    }
  }
}


# collapse adherence
for (i in 1:dim(data)[1]) {
  if (!is.na(data$adh[i])) {
    if (data$adh[i] == 1 || data$adh[i] == 2) {
      data$adh[i] = 0
      
    } else if (data$adh[i] == 3 || data$adh[i] == 4) {
      data$adh[i] = 1    
    }
  }
}

# collapse race
for (i in 1:dim(data)[1]) {
  if (!is.na(data$race[i])) {
    if (data$race[i] == 4) {
      data$race[i] = 2
      
    } else if (data$race[i] > 3) {
      data$race[i] = 4
    }
  }
}

# collapse education
for (i in 1:dim(data)[1]) {
  if (!is.na(data$educbas[i])) {
    if (data$educbas[i] == 2) {
      data$educbas[i] = 1
      
    } else if (data$educbas[i] == 3) {
      data$educbas[i] = 2
    } else {
      data$educbas[i] = 3
    }
  }
}

# collapse smoking status
for (i in 1:dim(data)[1]) {
  if (!is.na(data$smoke[i])) {
    if (data$smoke[i] == 1 || data$smoke[i] == 2) {
      data$smoke[i] = 0
      
    } else if (data$smoke[i] == 3) {
      data$smoke[i] = 1    
    }
  }
}

# Make the data set with change score, and baseline values (except adh for analysis)
row_to_keep = vector("logical", length = dim(data)[1])
for (i in 1:dim(data)[1]) {
  if (data$years[i] == 0 | data$years[i] == 2) {
    row_to_keep[i] = TRUE
  }
}

data = data[row_to_keep,]

# make wide to create change scores 
data = pivot_wider(
  data,
  id_cols = newid,
  names_from = years,
  values_from = c(
    "agg_ment", "agg_phys", "bmi", "smoke", "leu3n", "vload", "adh",
    "race", "educbas", "age", "years", "hard_drugs"
  )
)

# get adherence at two years
data$adh_0 = data$adh_2

# get change scores
data$leu3n_change = data$leu3n_2 - data$leu3n_0
data$log_vload_change = log10(data$vload_2) - log10(data$vload_0)
data$agg_ment_change = data$agg_ment_2 - data$agg_ment_0
data$agg_phys_change = data$agg_phys_2 - data$agg_phys_0

# remove unecessary columns
data = data[, -c(7, 9, 14, 17, 19, 21, 22, 23, 25)]


# filter out NA values
data = data[!is.na(data$adh_2), ]
data = data[!is.na(data$log_vload_change), ]
data = data[!is.na(data$leu3n_change), ]
data = data[!is.na(data$agg_ment_change), ]
data = data[!is.na(data$agg_phys_change), ]
data = data[!is.na(data$bmi_0), ]

########################
# Frequentist Analysis #
########################
set.seed(8675309)

# univariate  frequentist models
leu3n_uni = lm(leu3n_change ~ leu3n_0 + hard_drugs_0, data=data)
vload_uni = lm(log_vload_change ~ hard_drugs_0, data=data)
aggment_uni = lm(agg_ment_change ~ hard_drugs_0, data=data)
aggphys_uni = lm(agg_phys_change ~ hard_drugs_0, data=data)


# multivariate frequentist models
leu3n_full = lm(leu3n_change ~ leu3n_0 + age_0 + bmi_0 + smoke_0 + educbas_0 + adh_2 + hard_drugs_0, data=data)
vload_full = lm(log_vload_change ~ vload_0 + age_0 + bmi_0 + smoke_0 + educbas_0 + adh_2 + hard_drugs_0, data=data)
aggment_full = lm(agg_ment_change ~ agg_phys_0 + age_0 + bmi_0 + smoke_0 + educbas_0 + adh_2 + hard_drugs_0, data=data)
aggphys_full = lm(agg_phys_change ~ agg_ment_0 + age_0 + bmi_0 + smoke_0 + educbas_0 + adh_2 + hard_drugs_0, data=data)

#####################
# Bayesian Analysis #
#####################

## HPDI function which calculates credible intervals using highest posterior densities
hpd = function(x, alpha = 0.05){
  n = length(x)
  m = round(n * alpha)
  x = sort(x)
  y = x[(n - m + 1):n] - x[1:m]
  z = min(y)
  k = which(y == z)[1]
  c(x[k], x[n - m + k])
}

# significance of detecting 10% difference function from @willippit
p_sig = function(ref, draws){
  p_left = mean(min(ref * c(.9, 1.1)) > draws + ref)
  p_right = mean(ref+draws > max(ref*c(.9, 1.1)))
  return(c(p_left, p_right, p_left + p_right))
}

### leu3n univariate ###
# y is the outcome in your linear regression model
y = c(data$leu3n_change)

# Now we set up the design matrix in our linear regression
# For the first question in worksheet 3 we are looking at an intercept and age
X = model.matrix(~ hard_drugs_0, data=data) 
N = nrow(X)
p = ncol(X)
a = 0.001 # gamma shape for a non-informative prior
b = 0.001 # gamma rate for a non-informative prior
m = c(0.5 , rep(0, (p-1))) # mvnorm mean  (mean in the prior on the regression coefficients)
R = matrix(0, p, p) # mvnorm covariance (1/var-cov in the prior on the regression coefficients)
diag(R) = 0.000001 # note that JAGS uses dispersion matrix (scalars) rather 

# create data list to pass to JAGS
jags_dat = list(y = y, X = X, N = N, p = p, a = a, b = b, m = m, R = R)
mod = jags.model("linMod.jags", data = jags_dat, n.adapt = 1000, n.chains = 2)

# Sample observations from the posterior distributions
set.seed(8675309)
iter = 25000 # number of draws for the MCMC
samples = coda.samples(mod, variable.names = c("beta", "sigma2"), n.iter = iter) # generates matrix of draws from posterior dist.

## To be used later
set.seed(8675309)
samples_dic = dic.samples(mod, n.iter = iter, type = "pD") # determine DIC of model fit (requires two chains)

## Diagnostics

# determine if the support of the parameter space is explored well (more adaptations)
plot(samples)

# determine amount of autocorrelation between draws (thinning parameter)
par(mfrow = c(ceiling(ncol(as.matrix(samples))/2), 2), mar = rep(1, 4))
apply(as.matrix(samples), 2, acf) # We will learn later how to assess if this graph says our MCMC worked correctly.

## Generate output table when the final model has been selected to summarize findings
draws = as.matrix(samples)

out_mat = matrix("", nrow = ncol(draws), ncol = 7)
out_mat[,1] = c(colnames(X), "sigma2") # names
out_mat[,2] = apply(draws, 2, function(x) round(mcse(x)$est, 3)) # batch mean
out_mat[,3] = apply(draws, 2, function(x) round(mcse(x)$se, 3)) # MCSE
out_mat[,4] = round(apply(draws, 2, sd), 3) # Std. Dev.
out_mat[,5] = apply(draws, 2, function(x)
  paste("(", round(hpd(x)[1], 3), ", ", round(hpd(x)[2], 3), ")", sep = "")) # HPDI
out_mat[,6] = sum(samples_dic$deviance)+sum(samples_dic$penalty)
out_mat[,7] = p_sig(mean(data[data$hard_drugs_0 == 0,]$leu3n_change), draws[, 2])[3]

colnames(out_mat) = c("Parameter", "Estimate", "MCSE", "Std. Dev.", "95% HPDI", "DIC", "P(10% Change)")

write_csv(as.data.frame(out_mat), "output/leu3n_unimodel.csv")

### leu3n full model ###
# y is the outcome in your linear regression model
y = c(data$leu3n_change)

# Now we set up the design matrix in our linear regression
# For the first question in worksheet 3 we are looking at an intercept and age
X = model.matrix(~ leu3n_0 + age_0 + bmi_0 + smoke_0 + educbas_0 + adh_2 + hard_drugs_0, data=data) 
N = nrow(X)
p = ncol(X)
a = 0.001 # gamma shape for a non-informative prior
b = 0.001 # gamma rate for a non-informative prior
m = c(0.5 , rep(0, (p-1))) # mvnorm mean  (mean in the prior on the regression coefficients)
R = matrix(0, p, p) # mvnorm covariance (1/var-cov in the prior on the regression coefficients)
diag(R) = 0.000001 # note that JAGS uses dispersion matrix (scalars) rather 

# create data list to pass to JAGS
jags_dat = list(y = y, X = X, N = N, p = p, a = a, b = b, m = m, R = R)
mod = jags.model("linMod.jags", data = jags_dat, n.adapt = 1000, n.chains = 2)

# Sample observations from the posterior distributions
set.seed(8675309)
iter = 25000 # number of draws for the MCMC
samples = coda.samples(mod, variable.names = c("beta", "sigma2"), n.iter = iter) # generates matrix of draws from posterior dist.

## To be used later
set.seed(8675309)
samples_dic = dic.samples(mod, n.iter = iter, type = "pD") # determine DIC of model fit (requires two chains)

## Diagnostics

# determine if the support of the parameter space is explored well (more adaptations)
plot(samples)

# determine amount of autocorrelation between draws (thinning parameter)
par(mfrow = c(ceiling(ncol(as.matrix(samples))/2), 2), mar = rep(1, 4))
apply(as.matrix(samples), 2, acf) # We will learn later how to assess if this graph says our MCMC worked correctly.

## Generate output table when the final model has been selected to summarize findings
draws = as.matrix(samples)

out_mat = matrix("", nrow = ncol(draws), ncol = 7)
out_mat[,1] = c(colnames(X), "sigma2") # names
out_mat[,2] = apply(draws, 2, function(x) round(mcse(x)$est, 3)) # batch mean
out_mat[,3] = apply(draws, 2, function(x) round(mcse(x)$se, 3)) # MCSE
out_mat[,4] = round(apply(draws, 2, sd), 3) # Std. Dev.
out_mat[,5] = apply(draws, 2, function(x)
  paste("(", round(hpd(x)[1], 3), ", ", round(hpd(x)[2], 3), ")", sep = "")) # HPDI
out_mat[,6] = sum(samples_dic$deviance)+sum(samples_dic$penalty)
out_mat[,7] = p_sig(mean(data[data$hard_drugs_0 == 0,]$leu3n_change), draws[, 8])[3]

colnames(out_mat) = c("Parameter", "Estimate", "MCSE", "Std. Dev.", "95% HPDI", "DIC", "P(10% Change)")

write_csv(as.data.frame(out_mat), "output/leu3n_fullmodel.csv")

### vload univariate ###
# y is the outcome in your linear regression model
y = c(data$log_vload_change)

# Now we set up the design matrix in our linear regression
# For the first question in worksheet 3 we are looking at an intercept and age
X = model.matrix(~ hard_drugs_0, data=data) 
N = nrow(X)
p = ncol(X)
a = 0.001 # gamma shape for a non-informative prior
b = 0.001 # gamma rate for a non-informative prior
m = c(0.5 , rep(0, (p-1))) # mvnorm mean  (mean in the prior on the regression coefficients)
R = matrix(0, p, p) # mvnorm covariance (1/var-cov in the prior on the regression coefficients)
diag(R) = 0.000001 # note that JAGS uses dispersion matrix (scalars) rather 

# create data list to pass to JAGS
jags_dat = list(y = y, X = X, N = N, p = p, a = a, b = b, m = m, R = R)
mod = jags.model("linMod.jags", data = jags_dat, n.adapt = 1000, n.chains = 2)

# Sample observations from the posterior distributions
set.seed(8675309)
iter = 25000 # number of draws for the MCMC
samples = coda.samples(mod, variable.names = c("beta", "sigma2"), n.iter = iter) # generates matrix of draws from posterior dist.

## To be used later
set.seed(8675309)
samples_dic = dic.samples(mod, n.iter = iter, type = "pD") # determine DIC of model fit (requires two chains)

## Diagnostics

# determine if the support of the parameter space is explored well (more adaptations)
plot(samples)

# determine amount of autocorrelation between draws (thinning parameter)
par(mfrow = c(ceiling(ncol(as.matrix(samples))/2), 2), mar = rep(1, 4))
apply(as.matrix(samples), 2, acf) # We will learn later how to assess if this graph says our MCMC worked correctly.

## Generate output table when the final model has been selected to summarize findings
draws = as.matrix(samples)

out_mat = matrix("", nrow = ncol(draws), ncol = 7)
out_mat[,1] = c(colnames(X), "sigma2") # names
out_mat[,2] = apply(draws, 2, function(x) round(mcse(x)$est, 3)) # batch mean
out_mat[,3] = apply(draws, 2, function(x) round(mcse(x)$se, 3)) # MCSE
out_mat[,4] = round(apply(draws, 2, sd), 3) # Std. Dev.
out_mat[,5] = apply(draws, 2, function(x)
  paste("(", round(hpd(x)[1], 3), ", ", round(hpd(x)[2], 3), ")", sep = "")) # HPDI
out_mat[,6] = sum(samples_dic$deviance)+sum(samples_dic$penalty)
out_mat[,7] = p_sig(mean(data[data$hard_drugs_0 == 0,]$log_vload_change), draws[, 2])[3]

colnames(out_mat) = c("Parameter", "Estimate", "MCSE", "Std. Dev.", "95% HPDI", "DIC", "P(10% Change)")

write_csv(as.data.frame(out_mat), "output/vload_unimodel.csv")

### vload full model ###
# y is the outcome in your linear regression model
y = c(data$log_vload_change)

# Now we set up the design matrix in our linear regression
# For the first question in worksheet 3 we are looking at an intercept and age
X = model.matrix(~ log10(vload_0) + age_0 + bmi_0 + smoke_0 + educbas_0 + adh_2 + hard_drugs_0, data=data) 
N = nrow(X)
p = ncol(X)
a = 0.001 # gamma shape for a non-informative prior
b = 0.001 # gamma rate for a non-informative prior
m = c(0.5 , rep(0, (p-1))) # mvnorm mean  (mean in the prior on the regression coefficients)
R = matrix(0, p, p) # mvnorm covariance (1/var-cov in the prior on the regression coefficients)
diag(R) = 0.000001 # note that JAGS uses dispersion matrix (scalars) rather 

# create data list to pass to JAGS
jags_dat = list(y = y, X = X, N = N, p = p, a = a, b = b, m = m, R = R)
mod = jags.model("linMod.jags", data = jags_dat, n.adapt = 1000, n.chains = 2)
# Sample observations from the posterior distributions
set.seed(8675309)
iter = 25000 # number of draws for the MCMC
samples = coda.samples(mod, variable.names = c("beta", "sigma2"), n.iter = iter) # generates matrix of draws from posterior dist.

## To be used later
set.seed(8675309)
samples_dic = dic.samples(mod, n.iter = iter, type = "pD") # determine DIC of model fit (requires two chains)

## Diagnostics

# determine if the support of the parameter space is explored well (more adaptations)
plot(samples)

# determine amount of autocorrelation between draws (thinning parameter)
par(mfrow = c(ceiling(ncol(as.matrix(samples))/2), 2), mar = rep(1, 4))
apply(as.matrix(samples), 2, acf) # We will learn later how to assess if this graph says our MCMC worked correctly.

## Generate output table when the final model has been selected to summarize findings
draws = as.matrix(samples)

out_mat = matrix("", nrow = ncol(draws), ncol = 7)
out_mat[,1] = c(colnames(X), "sigma2") # names
out_mat[,2] = apply(draws, 2, function(x) round(mcse(x)$est, 3)) # batch mean
out_mat[,3] = apply(draws, 2, function(x) round(mcse(x)$se, 3)) # MCSE
out_mat[,4] = round(apply(draws, 2, sd), 3) # Std. Dev.
out_mat[,5] = apply(draws, 2, function(x)
  paste("(", round(hpd(x)[1], 3), ", ", round(hpd(x)[2], 3), ")", sep = "")) # HPDI
out_mat[,6] = sum(samples_dic$deviance)+sum(samples_dic$penalty)
out_mat[,7] = p_sig(mean(data[data$hard_drugs_0 == 0,]$log_vload_change), draws[, 8])[3]

colnames(out_mat) = c("Parameter", "Estimate", "MCSE", "Std. Dev.", "95% HPDI", "DIC", "P(10% Change)")

write_csv(as.data.frame(out_mat), "output/vload_fullmodel.csv")

### agg_ment uni model ###
# y is the outcome in your linear regression model
y = c(data$agg_ment_change)

# Now we set up the design matrix in our linear regression
# For the first question in worksheet 3 we are looking at an intercept and age
X = model.matrix(~ hard_drugs_0, data=data) 
N = nrow(X)
p = ncol(X)
a = 0.001 # gamma shape for a non-informative prior
b = 0.001 # gamma rate for a non-informative prior
m = c(0.5 , rep(0, (p-1))) # mvnorm mean  (mean in the prior on the regression coefficients)
R = matrix(0, p, p) # mvnorm covariance (1/var-cov in the prior on the regression coefficients)
diag(R) = 0.000001 # note that JAGS uses dispersion matrix (scalars) rather 

# create data list to pass to JAGS
jags_dat = list(y = y, X = X, N = N, p = p, a = a, b = b, m = m, R = R)
mod = jags.model("linMod.jags", data = jags_dat, n.adapt = 1000, n.chains = 2)
# Sample observations from the posterior distributions
set.seed(8675309)
iter = 25000 # number of draws for the MCMC
samples = coda.samples(mod, variable.names = c("beta", "sigma2"), n.iter = iter) # generates matrix of draws from posterior dist.

## To be used later
set.seed(8675309)
samples_dic = dic.samples(mod, n.iter = iter, type = "pD") # determine DIC of model fit (requires two chains)

## Diagnostics

# determine if the support of the parameter space is explored well (more adaptations)
plot(samples)

# determine amount of autocorrelation between draws (thinning parameter)
par(mfrow = c(ceiling(ncol(as.matrix(samples))/2), 2), mar = rep(1, 4))
apply(as.matrix(samples), 2, acf) # We will learn later how to assess if this graph says our MCMC worked correctly.

## Generate output table when the final model has been selected to summarize findings
draws = as.matrix(samples)

out_mat = matrix("", nrow = ncol(draws), ncol = 7)
out_mat[,1] = c(colnames(X), "sigma2") # names
out_mat[,2] = apply(draws, 2, function(x) round(mcse(x)$est, 3)) # batch mean
out_mat[,3] = apply(draws, 2, function(x) round(mcse(x)$se, 3)) # MCSE
out_mat[,4] = round(apply(draws, 2, sd), 3) # Std. Dev.
out_mat[,5] = apply(draws, 2, function(x)
  paste("(", round(hpd(x)[1], 3), ", ", round(hpd(x)[2], 3), ")", sep = "")) # HPDI
out_mat[,6] = sum(samples_dic$deviance)+sum(samples_dic$penalty)
out_mat[,7] = p_sig(mean(data[data$hard_drugs_0 == 0,]$agg_ment_change), draws[, 2])[3]

colnames(out_mat) = c("Parameter", "Estimate", "MCSE", "Std. Dev.", "95% HPDI", "DIC", "P(10% Change)")

write_csv(as.data.frame(out_mat), "output/agg_ment_unimodel.csv")

### agg_ment multivariate ###
# y is the outcome in your linear regression model
y = c(data$agg_ment_change)

# Now we set up the design matrix in our linear regression
# For the first question in worksheet 3 we are looking at an intercept and age
X = model.matrix(~ agg_ment_0 + age_0 + bmi_0 + smoke_0 + educbas_0 + adh_2 + hard_drugs_0, data=data) 
N = nrow(X)
p = ncol(X)
a = 0.001 # gamma shape for a non-informative prior
b = 0.001 # gamma rate for a non-informative prior
m = c(0.5 , rep(0, (p-1))) # mvnorm mean  (mean in the prior on the regression coefficients)
R = matrix(0, p, p) # mvnorm covariance (1/var-cov in the prior on the regression coefficients)
diag(R) = 0.000001 # note that JAGS uses dispersion matrix (scalars) rather 

# create data list to pass to JAGS
jags_dat = list(y = y, X = X, N = N, p = p, a = a, b = b, m = m, R = R)
mod = jags.model("linMod.jags", data = jags_dat, n.adapt = 1000, n.chains = 2)
# Sample observations from the posterior distributions
set.seed(8675309)
iter = 25000 # number of draws for the MCMC
samples = coda.samples(mod, variable.names = c("beta", "sigma2"), n.iter = iter) # generates matrix of draws from posterior dist.

## To be used later
set.seed(8675309)
samples_dic = dic.samples(mod, n.iter = iter, type = "pD") # determine DIC of model fit (requires two chains)

## Diagnostics

# determine if the support of the parameter space is explored well (more adaptations)
plot(samples)

# determine amount of autocorrelation between draws (thinning parameter)
par(mfrow = c(ceiling(ncol(as.matrix(samples))/2), 2), mar = rep(1, 4))
apply(as.matrix(samples), 2, acf) # We will learn later how to assess if this graph says our MCMC worked correctly.

## Generate output table when the final model has been selected to summarize findings
draws = as.matrix(samples)

out_mat = matrix("", nrow = ncol(draws), ncol = 7)
out_mat[,1] = c(colnames(X), "sigma2") # names
out_mat[,2] = apply(draws, 2, function(x) round(mcse(x)$est, 3)) # batch mean
out_mat[,3] = apply(draws, 2, function(x) round(mcse(x)$se, 3)) # MCSE
out_mat[,4] = round(apply(draws, 2, sd), 3) # Std. Dev.
out_mat[,5] = apply(draws, 2, function(x)
  paste("(", round(hpd(x)[1], 3), ", ", round(hpd(x)[2], 3), ")", sep = "")) # HPDI
out_mat[,6] = sum(samples_dic$deviance)+sum(samples_dic$penalty)
out_mat[,7] = p_sig(mean(data[data$hard_drugs_0 == 0,]$agg_ment_change), draws[, 8])[3]

colnames(out_mat) = c("Parameter", "Estimate", "MCSE", "Std. Dev.", "95% HPDI", "DIC", "P(10% Change)")

write_csv(as.data.frame(out_mat), "output/agg_ment_fullmodel.csv")

### agg_phys uni ###
# y is the outcome in your linear regression model
y = c(data$agg_phys_change)

# Now we set up the design matrix in our linear regression
# For the first question in worksheet 3 we are looking at an intercept and age
X = model.matrix(~ hard_drugs_0, data=data) 
N = nrow(X)
p = ncol(X)
a = 0.001 # gamma shape for a non-informative prior
b = 0.001 # gamma rate for a non-informative prior
m = c(0.5 , rep(0, (p-1))) # mvnorm mean  (mean in the prior on the regression coefficients)
R = matrix(0, p, p) # mvnorm covariance (1/var-cov in the prior on the regression coefficients)
diag(R) = 0.000001 # note that JAGS uses dispersion matrix (scalars) rather 

# create data list to pass to JAGS
jags_dat = list(y = y, X = X, N = N, p = p, a = a, b = b, m = m, R = R)
mod = jags.model("linMod.jags", data = jags_dat, n.adapt = 1000, n.chains = 2)
# Sample observations from the posterior distributions
set.seed(8675309)
iter = 25000 # number of draws for the MCMC
samples = coda.samples(mod, variable.names = c("beta", "sigma2"), n.iter = iter) # generates matrix of draws from posterior dist.

## To be used later
set.seed(8675309)
samples_dic = dic.samples(mod, n.iter = iter, type = "pD") # determine DIC of model fit (requires two chains)

## Diagnostics

# determine if the support of the parameter space is explored well (more adaptations)
plot(samples)

# determine amount of autocorrelation between draws (thinning parameter)
par(mfrow = c(ceiling(ncol(as.matrix(samples))/2), 2), mar = rep(1, 4))
apply(as.matrix(samples), 2, acf) # We will learn later how to assess if this graph says our MCMC worked correctly.

## Generate output table when the final model has been selected to summarize findings
draws = as.matrix(samples)

out_mat = matrix("", nrow = ncol(draws), ncol = 7)
out_mat[,1] = c(colnames(X), "sigma2") # names
out_mat[,2] = apply(draws, 2, function(x) round(mcse(x)$est, 3)) # batch mean
out_mat[,3] = apply(draws, 2, function(x) round(mcse(x)$se, 3)) # MCSE
out_mat[,4] = round(apply(draws, 2, sd), 3) # Std. Dev.
out_mat[,5] = apply(draws, 2, function(x)
  paste("(", round(hpd(x)[1], 3), ", ", round(hpd(x)[2], 3), ")", sep = "")) # HPDI
out_mat[,6] = sum(samples_dic$deviance)+sum(samples_dic$penalty)
out_mat[,7] = p_sig(mean(data[data$hard_drugs_0 == 0,]$agg_phys_change), draws[, 2])[3]

colnames(out_mat) = c("Parameter", "Estimate", "MCSE", "Std. Dev.", "95% HPDI", "DIC", "P(10% Change)")

write_csv(as.data.frame(out_mat), "output/agg_phys_unimodel.csv")

### agg_phys full model ###
# y is the outcome in your linear regression model
y = c(data$agg_phys_change)

# Now we set up the design matrix in our linear regression
# For the first question in worksheet 3 we are looking at an intercept and age
X = model.matrix(~ agg_phys_0 + age_0 + bmi_0 + smoke_0 + educbas_0 + adh_2 + hard_drugs_0, data=data) 
N = nrow(X)
p = ncol(X)
a = 0.001 # gamma shape for a non-informative prior
b = 0.001 # gamma rate for a non-informative prior
m = c(0.5 , rep(0, (p-1))) # mvnorm mean  (mean in the prior on the regression coefficients)
R = matrix(0, p, p) # mvnorm covariance (1/var-cov in the prior on the regression coefficients)
diag(R) = 0.000001 # note that JAGS uses dispersion matrix (scalars) rather 

# create data list to pass to JAGS
jags_dat = list(y = y, X = X, N = N, p = p, a = a, b = b, m = m, R = R)
mod = jags.model("linMod.jags", data = jags_dat, n.adapt = 1000, n.chains = 2)
# Sample observations from the posterior distributions
set.seed(8675309)
iter = 25000 # number of draws for the MCMC
samples = coda.samples(mod, variable.names = c("beta", "sigma2"), n.iter = iter) # generates matrix of draws from posterior dist.

## To be used later
set.seed(8675309)
samples_dic = dic.samples(mod, n.iter = iter, type = "pD") # determine DIC of model fit (requires two chains)

## Diagnostics

# determine if the support of the parameter space is explored well (more adaptations)
plot(samples)

# determine amount of autocorrelation between draws (thinning parameter)
par(mfrow = c(ceiling(ncol(as.matrix(samples))/2), 2), mar = rep(1, 4))
apply(as.matrix(samples), 2, acf) # We will learn later how to assess if this graph says our MCMC worked correctly.

## Generate output table when the final model has been selected to summarize findings
draws = as.matrix(samples)

out_mat = matrix("", nrow = ncol(draws), ncol = 7)
out_mat[,1] = c(colnames(X), "sigma2") # names
out_mat[,2] = apply(draws, 2, function(x) round(mcse(x)$est, 3)) # batch mean
out_mat[,3] = apply(draws, 2, function(x) round(mcse(x)$se, 3)) # MCSE
out_mat[,4] = round(apply(draws, 2, sd), 3) # Std. Dev.
out_mat[,5] = apply(draws, 2, function(x)
  paste("(", round(hpd(x)[1], 3), ", ", round(hpd(x)[2], 3), ")", sep = "")) # HPDI
out_mat[,6] = sum(samples_dic$deviance)+sum(samples_dic$penalty)
out_mat[,7] = p_sig(mean(data[data$hard_drugs_0 == 0,]$agg_phys_change), draws[, 8])[3]

colnames(out_mat) = c("Parameter", "Estimate", "MCSE", "Std. Dev.", "95% HPDI", "DIC", "P(10% Change)")

write_csv(as.data.frame(out_mat), "output/agg_phys_fullmodel.csv")

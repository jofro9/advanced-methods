# Dependencies
library(readr)
library(rjags)
library(mcmcse)
library(tidyverse)

setwd("/Users/juliawrobel/Documents/teaching/2019/nichole_class/Bios6624_Julia/Bios6624_Julia/PSAExample")
psaclean2 <- read_csv("DataProcessed/psaclean2.csv")

iter <- 10000 # number of draws

# fit model with multiple parameters
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
# adaptations are changees in step size
mod <- jags.model("Code/jags/linMod.jags", data = jags_dat, n.adapt = 0, n.chains = 2)

samples <- coda.samples(mod, variable.names = c("beta", "sigma2"), n.iter = iter)

plot(samples)

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



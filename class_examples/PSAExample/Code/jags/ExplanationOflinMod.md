% This is what defines your model in jags.
% Everthing goes inside of the model statement.
% First we define our normally distributed likelihood.
% eta is our mean model and tau is the precison of the model error (or 1/sigma^2e) in the way we normally write up regressions
% eta is the linear regresison model and here is written in matrix format.  Here things are defined for each observation
% in our dataset so we need the row indicator i.  Write down your likelihood to see what we are setting this up this way.

%Next we define our priors.  We need the model error precision defined.  Here is is set up with a gamma prior.
% Then we need our regression coefficients.  In this framework all the values will be sampled together (more on the specifics next week).  You can write other set ups where you list each regression coefficient separately.

%sigma2 is simply doing a variable transformation so that we can print out summaries on a different parameter scale
%(the model error variance scale like we might like)


model {
	
	# likelihood

	for (i in 1:N) {

		y[i] ~ dnorm(eta[i], tau)
		eta[i] <- X[i,] %*% beta

	}

	# priors

	beta[1:p] ~ dmnorm(m, R)
	tau ~ dgamma(a, b)
	sigma2 <- 1/tau

}
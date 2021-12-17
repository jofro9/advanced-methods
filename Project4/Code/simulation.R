library(car)
library(dplyr)
library(glmnet)
library(hdrm)
library(olsrr)
library(tseries)

########################
# function definitions #
########################

pSelection = function(sim_data, betas, p, ci, names) {
  # run linear regression
  regression_model = lm(y ~ ., data=data.frame(y = sim_data$y, sim_data$X))

  # p value selection function
  p_selected_model = ols_step_backward_p(regression_model, prem = 0.15)
  idx = na.remove(match(colnames(data.frame(t(coef(p_selected_model$model)))), names))

  # get betas
  betas[idx] = coef(p_selected_model$model)[2:length(coef(p_selected_model$model))]

  # get p's
  p[idx] = summary(p_selected_model$model)$coefficients[2:length(coef(p_selected_model$model)), 4]
  
  # get ci's
  ci[, idx] = t(confint(p_selected_model$model))[, 2:length(coef(p_selected_model$model))]

  return(list(betas, p, ci))
}

aicSelection = function(sim_data, betas, p, ci, names) {
  # run linear regression
  regression_model = lm(y ~ ., data=data.frame(y = sim_data$y, sim_data$X))
  
  # aic value selection function
  aic_selected_model = step(regression_model, k=2, trace = 0, direction = "backward")
  idx = na.remove(match(colnames(data.frame(t(coef(aic_selected_model)))), names))
  
  # get betas
  betas[idx] = coef(aic_selected_model)[2:length(coef(aic_selected_model))]
  
  # get p's
  p[idx] = summary(aic_selected_model)$coefficients[2:length(coef(aic_selected_model)), 4]

  # get ci's
  ci[, idx] = t(confint(aic_selected_model))[, 2:length(coef(aic_selected_model))]

  return(list(betas, p, ci))
}

bicSelection = function(sim_data, betas, p, ci, names) {
  # run linear regression
  regression_model = lm(y ~ ., data=data.frame(y = sim_data$y, sim_data$X))
  
  # bic value selection function
  bic_selected_model = step(regression_model, k=log(n), trace = 0, direction = "backward")
  idx = na.remove(match(colnames(data.frame(t(coef(bic_selected_model)))), names))

  # get betas
  betas[idx] = coef(bic_selected_model)[2:length(coef(bic_selected_model))]

  # get p's
  p[idx] = summary(bic_selected_model)$coefficients[2:length(coef(bic_selected_model)), 4]

  # get ci's
  ci[, idx] = t(confint(bic_selected_model))[, 2:length(coef(bic_selected_model))]

  return(list(betas, p, ci))
}

lassoSelectionCV = function(grid, sim_data, betas, p, ci) {
  # lasso value selection function
  lasso = cv.glmnet(x = sim_data$X, y = sim_data$y, lambda = grid)
  
  # cross validated
  lasso_mod = glmnet(x = sim_data$X, y = sim_data$y, lambda = lasso$lambda.min)
  
  # get betas
  vals = lasso_mod$beta
  idx = which(vals != 0)
  betas[idx] = vals[idx]
  
  # get p's
  lin_reg = lm(sim_data$y ~ sim_data$X[, idx])
  p[idx] = summary(lin_reg)$coefficients[2:length(coef(lin_reg)), 4]
  
  # get ci's
  ci[, idx] = t(confint(lin_reg))[, 2:length(coef(lin_reg))]

  return(list(betas, p, ci))
}

lassoSelectionFixed = function(sim_data, betas, p, ci, lambda) {
  # lasso fixed value selection function
  lasso_mod = glmnet(x = sim_data$X, y = sim_data$y, lambda = lambda)
  
  # get betas
  vals = lasso_mod$beta
  idx = which(vals != 0)
  betas[idx] = vals[idx]
  
  # get p's
  lin_reg = lm(sim_data$y ~ sim_data$X[, idx])
  p[idx] = summary(lin_reg)$coefficients[2:length(coef(lin_reg)), 4]
  
  # get ci's
  ci[, idx] = t(confint(lin_reg))[, 2:length(coef(lin_reg))]

  return(list(betas, p, ci))
}

elasticNetSelectionCV = function(grid, sim_data, betas, p, ci, alpha) {
  # elastic net value selection function
  elastic_net = cv.glmnet(x = sim_data$X, y = sim_data$y, lambda = grid, alpha = alpha)
  
  # cross validated
  enet_mod = glmnet(x = sim_data$X, y = sim_data$y, lambda = elastic_net$lambda.min, alpha = alpha)
  
  # get betas
  vals = enet_mod$beta
  idx = which(vals != 0)
  betas[idx] = vals[idx]
  
  # get p's
  lin_reg = lm(sim_data$y ~ sim_data$X[, idx])
  p[idx] = summary(lin_reg)$coefficients[2:length(coef(lin_reg)), 4]
  
  # get ci's
  ci[, idx] = t(confint(lin_reg))[, 2:length(coef(lin_reg))]

  return(list(betas, p, ci))
}

elasticNetSelectionFixed = function(sim_data, betas, p, ci, lambda, alpha) {
  # elastic net fixed value selection function
  enet_mod = glmnet(x = sim_data$X, y = sim_data$y, lambda = lambda, alpha = alpha)
  
  # get betas
  vals = enet_mod$beta
  idx = which(vals != 0)
  betas[idx] = vals[idx]
  
  # get p's
  lin_reg = lm(sim_data$y ~ sim_data$X[, idx])
  p[idx] = summary(lin_reg)$coefficients[2:length(coef(lin_reg)), 4]
  
  # get ci's
  ci[, idx] = t(confint(lin_reg))[, 2:length(coef(lin_reg))]

  return(list(betas, p, ci))
}

modelSelection = function(sim_data, names, n) {
  # create vectors to store each set of betas
  p_betas = vector("double", length(colnames(sim_data$X)))
  aic_betas = vector("double", length(colnames(sim_data$X)))
  bic_betas = vector("double", length(colnames(sim_data$X)))
  cv_lasso_betas = vector("double", length(colnames(sim_data$X)))
  fixed_lasso_betas = vector("double", length(colnames(sim_data$X)))
  cv_enet_betas = vector("double", length(colnames(sim_data$X)))
  fixed_enet_betas = vector("double", length(colnames(sim_data$X)))
  
  # create vectors to store each set of ps
  p_ps = vector("double", length(colnames(sim_data$X)))
  aic_ps = vector("double", length(colnames(sim_data$X)))
  bic_ps = vector("double", length(colnames(sim_data$X)))
  cv_lasso_ps = vector("double", length(colnames(sim_data$X)))
  fixed_lasso_ps = vector("double", length(colnames(sim_data$X)))
  cv_enet_ps = vector("double", length(colnames(sim_data$X)))
  fixed_enet_ps = vector("double", length(colnames(sim_data$X)))
  
  # create vectors to store each set of ps
  p_ci = matrix(0, nrow = 2, ncol = length(colnames(sim_data$X)))
  aic_ci = matrix(0, nrow = 2, ncol = length(colnames(sim_data$X)))
  bic_ci = matrix(0, nrow = 2, ncol = length(colnames(sim_data$X)))
  cv_lasso_ci = matrix(0, nrow = 2, ncol = length(colnames(sim_data$X)))
  fixed_lasso_ci = matrix(0, nrow = 2, ncol = length(colnames(sim_data$X)))
  cv_enet_ci = matrix(0, nrow = 2, ncol = length(colnames(sim_data$X)))
  fixed_enet_ci = matrix(0, nrow = 2, ncol = length(colnames(sim_data$X)))

  ### backward selection by p-value ###
  p_list = pSelection(sim_data, p_betas, p_ps, p_ci, names)
  p_betas = p_list[[1]]
  p_ps = p_list[[2]]
  p_ci = p_list[[3]]
  
  ### backward selection by AIC ###
  aic_list = aicSelection(sim_data, aic_betas, aic_ps, aic_ci, names)
  aic_betas = aic_list[[1]]
  aic_ps = aic_list[[2]]
  aic_ci = aic_list[[3]]
  
  ### backward selection by BIC ###
  bic_list = bicSelection(sim_data, bic_betas, bic_ps, bic_ci, names)
  bic_betas = bic_list[[1]]
  bic_ps = bic_list[[2]]
  bic_ci = bic_list[[3]]
  
  ### selection by lasso ###
  grid = 10^seq(10, -2, length = 100)
  cv_lasso_list = lassoSelectionCV(grid, sim_data, cv_lasso_betas, cv_lasso_ps, cv_lasso_ci)
  cv_lasso_betas = cv_lasso_list[[1]]
  cv_lasso_ps = cv_lasso_list[[2]]
  cv_lasso_ci = cv_lasso_list[[3]]

  fixed_lasso_list = lassoSelectionFixed(sim_data, fixed_lasso_betas, fixed_lasso_ps, fixed_lasso_ci, 0.2)
  fixed_lasso_betas = fixed_lasso_list[[1]]
  fixed_lasso_ps = fixed_lasso_list[[2]]
  fixed_lasso_ci = fixed_lasso_list[[3]]

  ### selection by elastic net ###
  cv_enet_list = elasticNetSelectionCV(grid, sim_data, cv_enet_betas, cv_enet_ps, cv_enet_ci, 0.5)
  cv_enet_betas = cv_enet_list[[1]]
  cv_enet_ps = cv_enet_list[[2]]
  cv_enet_ci = cv_enet_list[[3]]

  fixed_enet_list = elasticNetSelectionFixed(sim_data, fixed_enet_betas, fixed_enet_ps, fixed_enet_ci, 0.2, 0.5)
  fixed_enet_betas = fixed_enet_list[[1]]
  fixed_enet_ps = fixed_enet_list[[2]]
  fixed_enet_ci = fixed_enet_list[[3]]

  return(
    list(
      list(p_betas, aic_betas, bic_betas, cv_lasso_betas, fixed_lasso_betas, cv_enet_betas, fixed_enet_betas),
      list(p_ps, aic_ps, bic_ps, cv_lasso_ps, fixed_lasso_ps, cv_enet_ps, fixed_enet_ps),
      list(p_ci, aic_ci, bic_ci, cv_lasso_ci, fixed_lasso_ci, cv_enet_ci, fixed_enet_ci)
    )
  )
}

#######################
#     simulation      #
#######################
### Scenario 1a ###
n = 250
corr = NULL
rho = 0
N = 1000
p = 20
beta = vector("double", p)
beta[1:5] = c(0.5, 1., 1.5, 2., 2.5) / 3

### storage for betas ###
p_betas = matrix(data = 0, nrow = N, ncol = p)
aic_betas = matrix(data = 0, nrow = N, ncol = p)
bic_betas = matrix(data = 0, nrow = N, ncol = p)
cv_lasso_betas = matrix(data = 0, nrow = N, ncol = p)
fixed_lasso_betas = matrix(data = 0, nrow = N, ncol = p)
cv_enet_betas = matrix(data = 0, nrow = N, ncol = p)
fixed_enet_betas = matrix(data = 0, nrow = N, ncol = p)

### storage for p ###
p_p = matrix(data = 0, nrow = N, ncol = p)
aic_p = matrix(data = 0, nrow = N, ncol = p)
bic_p = matrix(data = 0, nrow = N, ncol = p)
cv_lasso_p = matrix(data = 0, nrow = N, ncol = p)
fixed_lasso_p = matrix(data = 0, nrow = N, ncol = p)
cv_enet_p = matrix(data = 0, nrow = N, ncol = p)
fixed_enet_p = matrix(data = 0, nrow = N, ncol = p)

### storage for ci ###
p_ci = matrix(data = 0, nrow = N * 2, ncol = p)
aic_ci = matrix(data = 0, nrow = N * 2, ncol = p)
bic_ci = matrix(data = 0, nrow = N * 2, ncol = p)
cv_lasso_ci = matrix(data = 0, nrow = N * 2, ncol = p)
fixed_lasso_ci = matrix(data = 0, nrow = N * 2, ncol = p)
cv_enet_ci = matrix(data = 0, nrow = N * 2, ncol = p)
fixed_enet_ci = matrix(data = 0, nrow = N * 2, ncol = p)

### simulation ###
j = 0
for (i in 1:N) {
  print(i)
  # Generate data
  sim_data = genData(
    n = n,
    p = p,
    p1 = 5,
    beta = beta,
    family = "gaussian",
    SNR = 1,
    signal = "heterogeneous",
    corr = corr,
    rho = rho
  )
  
  # get betas
  parameters_list = modelSelection(sim_data, colnames(sim_data$X), n)
  
  # store betas
  p_betas[i,] = parameters_list[[1]][[1]]
  aic_betas[i,] = parameters_list[[1]][[2]]
  bic_betas[i,] = parameters_list[[1]][[3]]
  cv_lasso_betas[i,] = parameters_list[[1]][[4]]
  fixed_lasso_betas[i,] = parameters_list[[1]][[5]]
  cv_enet_betas[i,] = parameters_list[[1]][[6]]
  fixed_enet_betas[i,] = parameters_list[[1]][[7]]
  
  # store p
  p_p[i,] = parameters_list[[2]][[1]]
  aic_p[i,] = parameters_list[[2]][[2]]
  bic_p[i,] = parameters_list[[2]][[3]]
  cv_lasso_p[i,] = parameters_list[[2]][[4]]
  fixed_lasso_p[i,] = parameters_list[[2]][[5]]
  cv_enet_p[i,] = parameters_list[[2]][[6]]
  fixed_enet_p[i,] = parameters_list[[2]][[7]]
  
  # store ci
  p_ci[c((i + j), (2 * i)),] = parameters_list[[3]][[1]]
  aic_ci[c((i + j), (2 * i)),] = parameters_list[[3]][[2]]
  bic_ci[c((i + j), (2 * i)),] = parameters_list[[3]][[3]]
  cv_lasso_ci[c((i + j), (2 * i)),] = parameters_list[[3]][[4]]
  fixed_lasso_ci[c((i + j), (2 * i)),] = parameters_list[[3]][[5]]
  cv_enet_ci[c((i + j), (2 * i)),] = parameters_list[[3]][[6]]
  fixed_enet_ci[c((i + j), (2 * i)),] = parameters_list[[3]][[7]]
  
  j = j + 1
}

# 1a
write.table(p_betas, "../data_raw/p_betas_1a.txt", sep="\t")
write.table(aic_betas, "../data_raw/aic_betas_1a.txt", sep="\t")
write.table(bic_betas, "../data_raw/bic_betas_1a.txt", sep="\t")
write.table(cv_lasso_betas, "../data_raw/cv_lasso_betas_1a.txt", sep="\t")
write.table(fixed_lasso_betas, "../data_raw/fixed_lasso_betas_1a.txt", sep="\t")
write.table(cv_enet_betas, "../data_raw/cv_enet_betas_1a.txt", sep="\t")
write.table(fixed_enet_betas, "../data_raw/fixed_enet_betas_1a.txt", sep="\t")

write.table(p_p, "../data_raw/p_p_1a.txt", sep="\t")
write.table(aic_p, "../data_raw/aic_p_1a.txt", sep="\t")
write.table(bic_p, "../data_raw/bic_p_1a.txt", sep="\t")
write.table(cv_lasso_p, "../data_raw/cv_lasso_p_1a.txt", sep="\t")
write.table(fixed_lasso_p, "../data_raw/fixed_lasso_p_1a.txt", sep="\t")
write.table(cv_enet_p, "../data_raw/cv_enet_p_1a.txt", sep="\t")
write.table(fixed_enet_p, "../data_raw/fixed_enet_p_1a.txt", sep="\t")

write.table(p_ci, "../data_raw/p_ci_1a.txt", sep="\t")
write.table(aic_ci, "../data_raw/aic_ci_1a.txt", sep="\t")
write.table(bic_ci, "../data_raw/bic_ci_1a.txt", sep="\t")
write.table(cv_lasso_ci, "../data_raw/cv_lasso_ci_1a.txt", sep="\t")
write.table(fixed_lasso_ci, "../data_raw/fixed_lasso_ci_1a.txt", sep="\t")
write.table(cv_enet_ci, "../data_raw/cv_enet_ci_1a.txt", sep="\t")
write.table(fixed_enet_ci, "../data_raw/fixed_enet_ci_1a.txt", sep="\t")

### Scenario 1b ###
n = 250
corr = "exchangeable"
rho = 0.4
N = 1000
p = 20
beta = vector("double", p)
beta[1:5] = c(0.5, 1., 1.5, 2., 2.5) / 3

### storage for betas ###
p_betas = matrix(data = 0, nrow = N, ncol = p)
aic_betas = matrix(data = 0, nrow = N, ncol = p)
bic_betas = matrix(data = 0, nrow = N, ncol = p)
cv_lasso_betas = matrix(data = 0, nrow = N, ncol = p)
fixed_lasso_betas = matrix(data = 0, nrow = N, ncol = p)
cv_enet_betas = matrix(data = 0, nrow = N, ncol = p)
fixed_enet_betas = matrix(data = 0, nrow = N, ncol = p)

### storage for p ###
p_p = matrix(data = 0, nrow = N, ncol = p)
aic_p = matrix(data = 0, nrow = N, ncol = p)
bic_p = matrix(data = 0, nrow = N, ncol = p)
cv_lasso_p = matrix(data = 0, nrow = N, ncol = p)
fixed_lasso_p = matrix(data = 0, nrow = N, ncol = p)
cv_enet_p = matrix(data = 0, nrow = N, ncol = p)
fixed_enet_p = matrix(data = 0, nrow = N, ncol = p)

### storage for ci ###
p_ci = matrix(data = 0, nrow = N * 2, ncol = p)
aic_ci = matrix(data = 0, nrow = N * 2, ncol = p)
bic_ci = matrix(data = 0, nrow = N * 2, ncol = p)
cv_lasso_ci = matrix(data = 0, nrow = N * 2, ncol = p)
fixed_lasso_ci = matrix(data = 0, nrow = N * 2, ncol = p)
cv_enet_ci = matrix(data = 0, nrow = N * 2, ncol = p)
fixed_enet_ci = matrix(data = 0, nrow = N * 2, ncol = p)

### simulation ###
j = 0
for (i in 1:N) {
  print(i)
  # Generate data
  sim_data = genData(
    n = n,
    p = p,
    p1 = 5,
    beta = beta,
    family = "gaussian",
    SNR = 1,
    signal = "heterogeneous",
    corr = corr,
    rho = rho
  )
  
  # get betas
  parameters_list = modelSelection(sim_data, colnames(sim_data$X), n)
  
  # store betas
  p_betas[i,] = parameters_list[[1]][[1]]
  aic_betas[i,] = parameters_list[[1]][[2]]
  bic_betas[i,] = parameters_list[[1]][[3]]
  cv_lasso_betas[i,] = parameters_list[[1]][[4]]
  fixed_lasso_betas[i,] = parameters_list[[1]][[5]]
  cv_enet_betas[i,] = parameters_list[[1]][[6]]
  fixed_enet_betas[i,] = parameters_list[[1]][[7]]
  
  # store p
  p_p[i,] = parameters_list[[2]][[1]]
  aic_p[i,] = parameters_list[[2]][[2]]
  bic_p[i,] = parameters_list[[2]][[3]]
  cv_lasso_p[i,] = parameters_list[[2]][[4]]
  fixed_lasso_p[i,] = parameters_list[[2]][[5]]
  cv_enet_p[i,] = parameters_list[[2]][[6]]
  fixed_enet_p[i,] = parameters_list[[2]][[7]]
  
  # store ci
  p_ci[c((i + j), (2 * i)),] = parameters_list[[3]][[1]]
  aic_ci[c((i + j), (2 * i)),] = parameters_list[[3]][[2]]
  bic_ci[c((i + j), (2 * i)),] = parameters_list[[3]][[3]]
  cv_lasso_ci[c((i + j), (2 * i)),] = parameters_list[[3]][[4]]
  fixed_lasso_ci[c((i + j), (2 * i)),] = parameters_list[[3]][[5]]
  cv_enet_ci[c((i + j), (2 * i)),] = parameters_list[[3]][[6]]
  fixed_enet_ci[c((i + j), (2 * i)),] = parameters_list[[3]][[7]]
  
  j = j + 1
}

# 1b
write.table(p_betas, "../data_raw/p_betas_1b.txt", sep="\t")
write.table(aic_betas, "../data_raw/aic_betas_1b.txt", sep="\t")
write.table(bic_betas, "../data_raw/bic_betas_1b.txt", sep="\t")
write.table(cv_lasso_betas, "../data_raw/cv_lasso_betas_1b.txt", sep="\t")
write.table(fixed_lasso_betas, "../data_raw/fixed_lasso_betas_1b.txt", sep="\t")
write.table(cv_enet_betas, "../data_raw/cv_enet_betas_1b.txt", sep="\t")
write.table(fixed_enet_betas, "../data_raw/fixed_enet_betas_1b.txt", sep="\t")

write.table(p_p, "../data_raw/p_p_1b.txt", sep="\t")
write.table(aic_p, "../data_raw/aic_p_1b.txt", sep="\t")
write.table(bic_p, "../data_raw/bic_p_1b.txt", sep="\t")
write.table(cv_lasso_p, "../data_raw/cv_lasso_p_1b.txt", sep="\t")
write.table(fixed_lasso_p, "../data_raw/fixed_lasso_p_1b.txt", sep="\t")
write.table(cv_enet_p, "../data_raw/cv_enet_p_1b.txt", sep="\t")
write.table(fixed_enet_p, "../data_raw/fixed_enet_p_1b.txt", sep="\t")

write.table(p_ci, "../data_raw/p_ci_1b.txt", sep="\t")
write.table(aic_ci, "../data_raw/aic_ci_1b.txt", sep="\t")
write.table(bic_ci, "../data_raw/bic_ci_1b.txt", sep="\t")
write.table(cv_lasso_ci, "../data_raw/cv_lasso_ci_1b.txt", sep="\t")
write.table(fixed_lasso_ci, "../data_raw/fixed_lasso_ci_1b.txt", sep="\t")
write.table(cv_enet_ci, "../data_raw/cv_enet_ci_1b.txt", sep="\t")
write.table(fixed_enet_ci, "../data_raw/fixed_enet_ci_1b.txt", sep="\t")

### Scenario 2a ###
n = 500
corr = NULL
rho = 0
N = 1000
p = 20
beta = vector("double", p)
beta[1:5] = c(0.5, 1., 1.5, 2., 2.5) / 3

### storage for betas ###
p_betas = matrix(data = 0, nrow = N, ncol = p)
aic_betas = matrix(data = 0, nrow = N, ncol = p)
bic_betas = matrix(data = 0, nrow = N, ncol = p)
cv_lasso_betas = matrix(data = 0, nrow = N, ncol = p)
fixed_lasso_betas = matrix(data = 0, nrow = N, ncol = p)
cv_enet_betas = matrix(data = 0, nrow = N, ncol = p)
fixed_enet_betas = matrix(data = 0, nrow = N, ncol = p)

### storage for p ###
p_p = matrix(data = 0, nrow = N, ncol = p)
aic_p = matrix(data = 0, nrow = N, ncol = p)
bic_p = matrix(data = 0, nrow = N, ncol = p)
cv_lasso_p = matrix(data = 0, nrow = N, ncol = p)
fixed_lasso_p = matrix(data = 0, nrow = N, ncol = p)
cv_enet_p = matrix(data = 0, nrow = N, ncol = p)
fixed_enet_p = matrix(data = 0, nrow = N, ncol = p)

### storage for ci ###
p_ci = matrix(data = 0, nrow = N * 2, ncol = p)
aic_ci = matrix(data = 0, nrow = N * 2, ncol = p)
bic_ci = matrix(data = 0, nrow = N * 2, ncol = p)
cv_lasso_ci = matrix(data = 0, nrow = N * 2, ncol = p)
fixed_lasso_ci = matrix(data = 0, nrow = N * 2, ncol = p)
cv_enet_ci = matrix(data = 0, nrow = N * 2, ncol = p)
fixed_enet_ci = matrix(data = 0, nrow = N * 2, ncol = p)

### simulation ###
j = 0
for (i in 1:N) {
  print(i)
  # Generate data
  sim_data = genData(
    n = n,
    p = p,
    p1 = 5,
    beta = beta,
    family = "gaussian",
    SNR = 1,
    signal = "heterogeneous",
    corr = corr,
    rho = rho
  )
  
  # get betas
  parameters_list = modelSelection(sim_data, colnames(sim_data$X), n)
  
  # store betas
  p_betas[i,] = parameters_list[[1]][[1]]
  aic_betas[i,] = parameters_list[[1]][[2]]
  bic_betas[i,] = parameters_list[[1]][[3]]
  cv_lasso_betas[i,] = parameters_list[[1]][[4]]
  fixed_lasso_betas[i,] = parameters_list[[1]][[5]]
  cv_enet_betas[i,] = parameters_list[[1]][[6]]
  fixed_enet_betas[i,] = parameters_list[[1]][[7]]
  
  # store p
  p_p[i,] = parameters_list[[2]][[1]]
  aic_p[i,] = parameters_list[[2]][[2]]
  bic_p[i,] = parameters_list[[2]][[3]]
  cv_lasso_p[i,] = parameters_list[[2]][[4]]
  fixed_lasso_p[i,] = parameters_list[[2]][[5]]
  cv_enet_p[i,] = parameters_list[[2]][[6]]
  fixed_enet_p[i,] = parameters_list[[2]][[7]]
  
  # store ci
  p_ci[c((i + j), (2 * i)),] = parameters_list[[3]][[1]]
  aic_ci[c((i + j), (2 * i)),] = parameters_list[[3]][[2]]
  bic_ci[c((i + j), (2 * i)),] = parameters_list[[3]][[3]]
  cv_lasso_ci[c((i + j), (2 * i)),] = parameters_list[[3]][[4]]
  fixed_lasso_ci[c((i + j), (2 * i)),] = parameters_list[[3]][[5]]
  cv_enet_ci[c((i + j), (2 * i)),] = parameters_list[[3]][[6]]
  fixed_enet_ci[c((i + j), (2 * i)),] = parameters_list[[3]][[7]]
  
  j = j + 1
}

# 2a
write.table(p_betas, "../data_raw/p_betas_2a.txt", sep="\t")
write.table(aic_betas, "../data_raw/aic_betas_2a.txt", sep="\t")
write.table(bic_betas, "../data_raw/bic_betas_2a.txt", sep="\t")
write.table(cv_lasso_betas, "../data_raw/cv_lasso_betas_2a.txt", sep="\t")
write.table(fixed_lasso_betas, "../data_raw/fixed_lasso_betas_2a.txt", sep="\t")
write.table(cv_enet_betas, "../data_raw/cv_enet_betas_2a.txt", sep="\t")
write.table(fixed_enet_betas, "../data_raw/fixed_enet_betas_2a.txt", sep="\t")

write.table(p_p, "../data_raw/p_p_2a.txt", sep="\t")
write.table(aic_p, "../data_raw/aic_p_2a.txt", sep="\t")
write.table(bic_p, "../data_raw/bic_p_2a.txt", sep="\t")
write.table(cv_lasso_p, "../data_raw/cv_lasso_p_2a.txt", sep="\t")
write.table(fixed_lasso_p, "../data_raw/fixed_lasso_p_2a.txt", sep="\t")
write.table(cv_enet_p, "../data_raw/cv_enet_p_2a.txt", sep="\t")
write.table(fixed_enet_p, "../data_raw/fixed_enet_p_2a.txt", sep="\t")

write.table(p_ci, "../data_raw/p_ci_2a.txt", sep="\t")
write.table(aic_ci, "../data_raw/aic_ci_2a.txt", sep="\t")
write.table(bic_ci, "../data_raw/bic_ci_2a.txt", sep="\t")
write.table(cv_lasso_ci, "../data_raw/cv_lasso_ci_2a.txt", sep="\t")
write.table(fixed_lasso_ci, "../data_raw/fixed_lasso_ci_2a.txt", sep="\t")
write.table(cv_enet_ci, "../data_raw/cv_enet_ci_2a.txt", sep="\t")
write.table(fixed_enet_ci, "../data_raw/fixed_enet_ci_2a.txt", sep="\t")

### Scenario 2b ###
n = 500
corr = "exchangeable"
rho = 0.4
N = 1000
p = 20
beta = vector("double", p)
beta[1:5] = c(0.5, 1., 1.5, 2., 2.5) / 3

### storage for betas ###
p_betas = matrix(data = 0, nrow = N, ncol = p)
aic_betas = matrix(data = 0, nrow = N, ncol = p)
bic_betas = matrix(data = 0, nrow = N, ncol = p)
cv_lasso_betas = matrix(data = 0, nrow = N, ncol = p)
fixed_lasso_betas = matrix(data = 0, nrow = N, ncol = p)
cv_enet_betas = matrix(data = 0, nrow = N, ncol = p)
fixed_enet_betas = matrix(data = 0, nrow = N, ncol = p)

### storage for p ###
p_p = matrix(data = 0, nrow = N, ncol = p)
aic_p = matrix(data = 0, nrow = N, ncol = p)
bic_p = matrix(data = 0, nrow = N, ncol = p)
cv_lasso_p = matrix(data = 0, nrow = N, ncol = p)
fixed_lasso_p = matrix(data = 0, nrow = N, ncol = p)
cv_enet_p = matrix(data = 0, nrow = N, ncol = p)
fixed_enet_p = matrix(data = 0, nrow = N, ncol = p)

### storage for ci ###
p_ci = matrix(data = 0, nrow = N * 2, ncol = p)
aic_ci = matrix(data = 0, nrow = N * 2, ncol = p)
bic_ci = matrix(data = 0, nrow = N * 2, ncol = p)
cv_lasso_ci = matrix(data = 0, nrow = N * 2, ncol = p)
fixed_lasso_ci = matrix(data = 0, nrow = N * 2, ncol = p)
cv_enet_ci = matrix(data = 0, nrow = N * 2, ncol = p)
fixed_enet_ci = matrix(data = 0, nrow = N * 2, ncol = p)

### simulation ###
j = 0
for (i in 1:N) {
  print(i)
  # Generate data
  sim_data = genData(
    n = n,
    p = p,
    p1 = 5,
    beta = beta,
    family = "gaussian",
    SNR = 1,
    signal = "heterogeneous",
    corr = corr,
    rho = rho
  )
  
  # get betas
  parameters_list = modelSelection(sim_data, colnames(sim_data$X), n)
  
  # store betas
  p_betas[i,] = parameters_list[[1]][[1]]
  aic_betas[i,] = parameters_list[[1]][[2]]
  bic_betas[i,] = parameters_list[[1]][[3]]
  cv_lasso_betas[i,] = parameters_list[[1]][[4]]
  fixed_lasso_betas[i,] = parameters_list[[1]][[5]]
  cv_enet_betas[i,] = parameters_list[[1]][[6]]
  fixed_enet_betas[i,] = parameters_list[[1]][[7]]
  
  # store p
  p_p[i,] = parameters_list[[2]][[1]]
  aic_p[i,] = parameters_list[[2]][[2]]
  bic_p[i,] = parameters_list[[2]][[3]]
  cv_lasso_p[i,] = parameters_list[[2]][[4]]
  fixed_lasso_p[i,] = parameters_list[[2]][[5]]
  cv_enet_p[i,] = parameters_list[[2]][[6]]
  fixed_enet_p[i,] = parameters_list[[2]][[7]]
  
  # store ci
  p_ci[c((i + j), (2 * i)),] = parameters_list[[3]][[1]]
  aic_ci[c((i + j), (2 * i)),] = parameters_list[[3]][[2]]
  bic_ci[c((i + j), (2 * i)),] = parameters_list[[3]][[3]]
  cv_lasso_ci[c((i + j), (2 * i)),] = parameters_list[[3]][[4]]
  fixed_lasso_ci[c((i + j), (2 * i)),] = parameters_list[[3]][[5]]
  cv_enet_ci[c((i + j), (2 * i)),] = parameters_list[[3]][[6]]
  fixed_enet_ci[c((i + j), (2 * i)),] = parameters_list[[3]][[7]]
  
  j = j + 1
}

# 2b
write.table(p_betas, "../data_raw/p_betas_2b.txt", sep="\t")
write.table(aic_betas, "../data_raw/aic_betas_2b.txt", sep="\t")
write.table(bic_betas, "../data_raw/bic_betas_2b.txt", sep="\t")
write.table(cv_lasso_betas, "../data_raw/cv_lasso_betas_2b.txt", sep="\t")
write.table(fixed_lasso_betas, "../data_raw/fixed_lasso_betas_2b.txt", sep="\t")
write.table(cv_enet_betas, "../data_raw/cv_enet_betas_2b.txt", sep="\t")
write.table(fixed_enet_betas, "../data_raw/fixed_enet_betas_2b.txt", sep="\t")

write.table(p_p, "../data_raw/p_p_2b.txt", sep="\t")
write.table(aic_p, "../data_raw/aic_p_2b.txt", sep="\t")
write.table(bic_p, "../data_raw/bic_p_2b.txt", sep="\t")
write.table(cv_lasso_p, "../data_raw/cv_lasso_p_2b.txt", sep="\t")
write.table(fixed_lasso_p, "../data_raw/fixed_lasso_p_2b.txt", sep="\t")
write.table(cv_enet_p, "../data_raw/cv_enet_p_2b.txt", sep="\t")
write.table(fixed_enet_p, "../data_raw/fixed_enet_p_2b.txt", sep="\t")

write.table(p_ci, "../data_raw/p_ci_2b.txt", sep="\t")
write.table(aic_ci, "../data_raw/aic_ci_2b.txt", sep="\t")
write.table(bic_ci, "../data_raw/bic_ci_2b.txt", sep="\t")
write.table(cv_lasso_ci, "../data_raw/cv_lasso_ci_2b.txt", sep="\t")
write.table(fixed_lasso_ci, "../data_raw/fixed_lasso_ci_2b.txt", sep="\t")
write.table(cv_enet_ci, "../data_raw/cv_enet_ci_2b.txt", sep="\t")
write.table(fixed_enet_ci, "../data_raw/fixed_enet_ci_2b.txt", sep="\t")

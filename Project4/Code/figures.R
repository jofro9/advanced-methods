library(ggplot2)
library(ggpubr)
library(gt)
library(kableExtra)
library(tidyverse)

########################
# function definitions #
########################

# true positive rate function
truePositiveRate = function(betas) {
  sum = 0
  for (i in 1:length(betas)) {
    if (betas[i] != 0) {
      sum = sum + 1
    }
  }
  
  return(sum / length(betas))
}

# False positive rate function
falsePositiveRate = function(betas) {
  sum = 0
  for (i in 1:dim(betas)[1]) {
    for (j in 1:dim(betas)[2]) {
      if (betas[i, j] != 0) {
        sum = sum + 1
      }
    }
  }
  
  return(sum / (dim(betas)[1] * dim(betas)[2]))
}

# false discovery rate function
falseDiscoveryRate = function(full_betas, last_15_betas) {
  sum = 0
  for (i in 1:dim(full_betas)[1]) {
    full_idx = which(full_betas != 0)
    last_15_idx = which(last_15_betas != 0)
    sum = sum + (length(last_15_idx) / length(full_idx))
  }
  
  return(sum / dim(full_betas)[1])
}

# type one error rate function
typeOneErrorRate = function(p) {
  numerator = 0
  denomenator = 0
  for (i in 1:dim(p)[1]) {
    for (j in 1:dim(p)[2]) {
      if (p[i, j] < 0.05 & p[i, j] != 0) {
        numerator = numerator + 1
      }
      
      if (p[i, j] != 0) {
        denomenator = denomenator + 1
      }
    }
  }
  
  return(numerator / denomenator)
}

# type two error rate function
typeTwoErrorRate = function(p) {
  sum = 0
  for (i in 1:length(p)) {
    if (p[i] != 0 & p[i] > 0.05) {
      sum = sum + 1
    }
  }
  idx = which(p != 0)
  
  return(sum / length(idx))
}

# bias for sig vars function
parameterBiasNonZero = function(observed, expected) {
  sum = 0
  for (i in 1:length(observed)) {
    sum = sum + (expected - observed[i])
  }

  return(sum / length(observed))
}

# bias for non-sig vars function
parameterBiasZero = function(observed, expected) {
  return(sum(observed) / (nrow(observed) * ncol(observed)))
}

# coverage function
coverageCI = function(lower, upper, betas, coefficient) {
  sum = 0
  for (i in 1:length(lower)) {
    if (coefficient >= lower[i] & coefficient <= upper[i]) {
      sum = sum + 1
    }
  }
  idx = which(lower != 0)
  return(sum / length(lower))
}

########
#  1a  #
########
########################
#     read in data     #
########################

p_betas = read.delim("../data_raw/p_betas_1a.txt")
p_ci = read.delim("../data_raw/p_ci_1a.txt")
p_p = read.delim("../data_raw/p_p_1a.txt")

aic_betas = read.delim("../data_raw/aic_betas_1a.txt")
aic_ci = read.delim("../data_raw/aic_ci_1a.txt")
aic_p = read.delim("../data_raw/aic_p_1a.txt")

bic_betas = read.delim("../data_raw/bic_betas_1a.txt")
bic_ci = read.delim("../data_raw/bic_ci_1a.txt")
bic_p = read.delim("../data_raw/bic_p_1a.txt")

cv_lasso_betas = read.delim("../data_raw/cv_lasso_betas_1a.txt")
cv_lasso_ci = read.delim("../data_raw/cv_lasso_ci_1a.txt")
cv_lasso_p = read.delim("../data_raw/cv_lasso_p_1a.txt")

fixed_lasso_betas = read.delim("../data_raw/fixed_lasso_betas_1a.txt")
fixed_lasso_ci = read.delim("../data_raw/fixed_lasso_ci_1a.txt")
fixed_lasso_p = read.delim("../data_raw/fixed_lasso_p_1a.txt")

cv_enet_betas = read.delim("../data_raw/cv_enet_betas_1a.txt")
cv_enet_ci = read.delim("../data_raw/cv_enet_ci_1a.txt")
cv_enet_p = read.delim("../data_raw/cv_enet_p_1a.txt")

fixed_enet_betas = read.delim("../data_raw/fixed_enet_betas_1a.txt")
fixed_enet_ci = read.delim("../data_raw/fixed_enet_ci_1a.txt")
fixed_enet_p = read.delim("../data_raw/fixed_enet_p_1a.txt")

########################
# calculate parameters #
########################
beta = vector("double", 20)
beta[1:5] = c(0.5, 1., 1.5, 2., 2.5) / 3

# true positive rate
tpr = matrix(0, nrow = 7, ncol = 5)
for (i in 1:dim(p_betas[, 1:5])[2]) {
  tpr[1, i] = truePositiveRate(p_betas[, i])
}

for (i in 1:dim(aic_betas[, 1:5])[2]) {
  tpr[2, i] = truePositiveRate(aic_betas[, i])
}

for (i in 1:dim(bic_betas[, 1:5])[2]) {
  tpr[3, i] = truePositiveRate(bic_betas[, i])
}

for (i in 1:dim(cv_lasso_betas[, 1:5])[2]) {
  tpr[4, i] = truePositiveRate(cv_lasso_betas[, i])
}

for (i in 1:dim(fixed_lasso_betas[, 1:5])[2]) {
  tpr[5, i] = truePositiveRate(fixed_lasso_betas[, i])
}

for (i in 1:dim(cv_enet_betas[, 1:5])[2]) {
  tpr[6, i] = truePositiveRate(cv_enet_betas[, i])
}

for (i in 1:dim(fixed_enet_betas[, 1:5])[2]) {
  tpr[7, i] = truePositiveRate(fixed_enet_betas[, i])
}

# false positive rate
fpr = vector("double", 7)
fpr[1] = falsePositiveRate(p_betas[, 5:15])
fpr[2] = falsePositiveRate(aic_betas[, 5:15])
fpr[3] = falsePositiveRate(bic_betas[, 5:15])
fpr[4] = falsePositiveRate(cv_lasso_betas[, 5:15])
fpr[5] = falsePositiveRate(fixed_lasso_betas[, 5:15])
fpr[6] = falsePositiveRate(cv_enet_betas[, 5:15])
fpr[7] = falsePositiveRate(fixed_enet_betas[, 5:15])

# false discovery rate
fdr = vector("double", 7)
fdr[1] = falseDiscoveryRate(p_betas, p_betas[, 5:15])
fdr[2] = falseDiscoveryRate(aic_betas, aic_betas[, 5:15])
fdr[3] = falseDiscoveryRate(bic_betas, bic_betas[, 5:15])
fdr[4] = falseDiscoveryRate(cv_lasso_betas, cv_lasso_betas[, 5:15])
fdr[5] = falseDiscoveryRate(fixed_lasso_betas, fixed_lasso_betas[, 5:15])
fdr[6] = falseDiscoveryRate(cv_enet_betas, cv_enet_betas[, 5:15])
fdr[7] = falseDiscoveryRate(fixed_enet_betas, fixed_enet_betas[, 5:15])

# type I error
toe = vector("double", 7)
toe[1] = typeOneErrorRate(p_p[, 5:15])
toe[2] = typeOneErrorRate(aic_p[, 5:15])
toe[3] = typeOneErrorRate(bic_p[, 5:15])
toe[4] = typeOneErrorRate(cv_lasso_p[, 5:15])
toe[5] = typeOneErrorRate(fixed_lasso_p[, 5:15])
toe[6] = typeOneErrorRate(cv_enet_p[, 5:15])
toe[7] = typeOneErrorRate(fixed_enet_p[, 5:15])

# type II error rate
tte = matrix(0, nrow = 7, ncol = 5)
for (i in 1:dim(p_p[, 1:5])[2]) {
  tte[1, i] = typeTwoErrorRate(p_p[, i])
}

for (i in 1:dim(aic_p[, 1:5])[2]) {
  tte[2, i] = typeTwoErrorRate(aic_p[, i])
}

for (i in 1:dim(bic_p[, 1:5])[2]) {
  tte[3, i] = typeTwoErrorRate(bic_p[, i])
}

for (i in 1:dim(cv_lasso_p[, 1:5])[2]) {
  tte[4, i] = typeTwoErrorRate(cv_lasso_p[, i])
}

for (i in 1:dim(fixed_lasso_p[, 1:5])[2]) {
  tte[5, i] = typeTwoErrorRate(fixed_lasso_p[, i])
}

for (i in 1:dim(cv_enet_p[, 1:5])[2]) {
  tte[6, i] = typeTwoErrorRate(cv_enet_p[, i])
}

for (i in 1:dim(fixed_enet_p[, 1:5])[2]) {
  tte[7, i] = typeTwoErrorRate(fixed_enet_p[, i])
}

# parameters 1-5 bias
bias = matrix(0, nrow = 7, ncol = 5)
for (i in 1:dim(p_betas[, 1:5])[2]) {
  bias[1, i] = parameterBiasNonZero(p_betas[, i], beta[i])
}

for (i in 1:dim(aic_betas[, 1:5])[2]) {
  bias[2, i] = parameterBiasNonZero(aic_betas[, i], beta[i])
}

for (i in 1:dim(bic_betas[, 1:5])[2]) {
  bias[3, i] = parameterBiasNonZero(bic_betas[, i], beta[i])
}

for (i in 1:dim(cv_lasso_betas[, 1:5])[2]) {
  bias[4, i] = parameterBiasNonZero(cv_lasso_betas[, i], beta[i])
}

for (i in 1:dim(fixed_lasso_betas[, 1:5])[2]) {
  bias[5, i] = parameterBiasNonZero(fixed_lasso_betas[, i], beta[i])
}

for (i in 1:dim(cv_enet_betas[, 1:5])[2]) {
  bias[6, i] = parameterBiasNonZero(cv_enet_betas[, i], beta[i])
}

for (i in 1:dim(fixed_enet_betas[, 1:5])[2]) {
  bias[7, i] = parameterBiasNonZero(fixed_enet_betas[, i], beta[i])
}

# parameters 5-20 average bias
bias0 = vector("double", 7)
bias0[1] = parameterBiasZero(p_betas[, 5:15])
bias0[2] = parameterBiasZero(aic_betas[, 5:15])
bias0[3] = parameterBiasZero(bic_betas[, 5:15])
bias0[4] = parameterBiasZero(cv_lasso_betas[, 5:15])
bias0[5] = parameterBiasZero(fixed_lasso_betas[, 5:15])
bias0[6] = parameterBiasZero(cv_enet_betas[, 5:15])
bias0[7] = parameterBiasZero(fixed_enet_betas[, 5:15])

# parameters 1-5 CI coverage
coverage = matrix(0, nrow = 7, ncol = 5)
for (i in 1:dim(p_ci[, 1:5])[2]) {
  coverage[1, i] = coverageCI(p_ci[seq(1, ncol(p_ci) - 1, 2), i], p_ci[seq(2, ncol(p_ci), 2), i], p_betas[, i], beta[i])
}

for (i in 1:dim(aic_ci[, 1:5])[2]) {
  coverage[2, i] = coverageCI(aic_ci[seq(1, ncol(aic_ci) - 1, 2), i], aic_ci[seq(2, ncol(aic_ci), 2), i], aic_betas[, i], beta[i])
}

for (i in 1:dim(bic_ci[, 1:5])[2]) {
  coverage[3, i] = coverageCI(bic_ci[seq(1, ncol(bic_ci) - 1, 2), i], bic_ci[seq(2, ncol(bic_ci), 2), i], bic_betas[, i], beta[i])
}

for (i in 1:dim(cv_lasso_ci[, 1:5])[2]) {
  coverage[4, i] = coverageCI(cv_lasso_ci[seq(1, ncol(cv_lasso_ci) - 1, 2), i], cv_lasso_ci[seq(2, ncol(cv_lasso_ci), 2), i], cv_lasso_betas[, i], beta[i])
}

for (i in 1:dim(fixed_lasso_ci[, 1:5])[2]) {
  coverage[5, i] = coverageCI(fixed_lasso_ci[seq(1, ncol(fixed_lasso_ci) - 1, 2), i], fixed_lasso_ci[seq(2, ncol(fixed_lasso_ci), 2), i], fixed_lasso_betas[, i], beta[i])
}

for (i in 1:dim(cv_enet_ci[, 1:5])[2]) {
  coverage[6, i] = coverageCI(cv_enet_ci[seq(1, ncol(cv_enet_ci) - 1, 2), i], cv_enet_ci[seq(2, ncol(cv_enet_ci), 2), i], cv_enet_betas[, 1:5], beta[i])
}

for (i in 1:dim(fixed_enet_ci[, 1:5])[2]) {
  coverage[7, i] = coverageCI(fixed_enet_ci[seq(1, ncol(fixed_enet_ci) - 1, 2), i], fixed_enet_ci[seq(2, ncol(fixed_enet_ci), 2), i], fixed_enet_betas[, i], beta[i])
}

################
# build table #
################

table_1a = data.frame(
  "LRT" = c(
    fdr[1], fpr[1], tpr[1, 1], tpr[1, 2], tpr[1, 3],
    tpr[1, 4], tpr[1, 5], coverage[1, 1], coverage[1, 2],
    coverage[1, 3], coverage[1, 4], coverage[1, 5],
    tte[1, 1], tte[1, 2], tte[1, 3], tte[1, 4], tte[1, 5],
    toe[1], bias[1, 1], bias[1, 2], bias[1, 3], bias[1, 4],
    bias[1, 5], bias0[1]
  ),
  "AIC" = c(
    fdr[2], fpr[2], tpr[2, 1], tpr[2, 2], tpr[2, 3],
    tpr[2, 4], tpr[2, 5], coverage[2, 1], coverage[2, 2],
    coverage[2, 3], coverage[2, 4], coverage[2, 5],
    tte[2, 1], tte[2, 2], tte[2, 3], tte[2, 4], tte[2, 5],
    toe[1], bias[2, 1], bias[2, 2], bias[2, 3], bias[2, 4],
    bias[2, 5], bias0[2]
  ),
  "BIC" = c(
    fdr[3], fpr[3], tpr[3, 1], tpr[3, 2], tpr[3, 3],
    tpr[3, 4], tpr[3, 5], coverage[3, 1], coverage[3, 2],
    coverage[3, 3], coverage[3, 4], coverage[3, 5],
    tte[3, 1], tte[3, 2], tte[3, 3], tte[3, 4], tte[3, 5],
    toe[1], bias[3, 1], bias[3, 2], bias[3, 3], bias[3, 4],
    bias[3, 5], bias0[3]
  ),
  "cv-LASSO" = c(
    fdr[4], fpr[4], tpr[4, 1], tpr[4, 2], tpr[4, 3],
    tpr[4, 4], tpr[4, 5], coverage[4, 1], coverage[4, 2],
    coverage[4, 3], coverage[4, 4], coverage[4, 5],
    tte[4, 1], tte[4, 2], tte[4, 3], tte[4, 4], tte[4, 5],
    toe[1], bias[4, 1], bias[4, 2], bias[4, 3], bias[4, 4],
    bias[4, 5], bias0[4]
  ),
  "fixed-LASSO" = c(
    fdr[5], fpr[5], tpr[5, 1], tpr[5, 2], tpr[5, 3],
    tpr[5, 4], tpr[5, 5], coverage[5, 1], coverage[5, 2],
    coverage[5, 3], coverage[5, 4], coverage[5, 5],
    tte[5, 1], tte[5, 2], tte[5, 3], tte[5, 4], tte[5, 5],
    toe[1], bias[5, 1], bias[5, 2], bias[5, 3], bias[5, 4],
    bias[5, 5], bias0[5]
  ),
  "cv-ENET" = c(
    fdr[6], fpr[6], tpr[6, 1], tpr[6, 2], tpr[6, 3],
    tpr[6, 4], tpr[6, 5], coverage[6, 1], coverage[6, 2],
    coverage[6, 3], coverage[6, 4], coverage[6, 5],
    tte[6, 1], tte[6, 2], tte[6, 3], tte[6, 4], tte[6, 5],
    toe[1], bias[6, 1], bias[6, 2], bias[6, 3], bias[6, 4],
    bias[6, 5], bias0[6]
  ),
  "fixed-ENET" = c(
    fdr[7], fpr[7], tpr[7, 1], tpr[7, 2], tpr[7, 3],
    tpr[7, 4], tpr[7, 5], coverage[7, 1], coverage[7, 2],
    coverage[7, 3], coverage[7, 4], coverage[7, 5],
    tte[7, 1], tte[7, 2], tte[7, 3], tte[7, 4], tte[7, 5],
    toe[1], bias[7, 1], bias[7, 2], bias[7, 3], bias[7, 4],
    bias[7, 5], bias0[7]
  )
)

rownames(table_1a) = c(
  "FDR",
  "FPR",
  "TPR: X1",
  "TPR: X2",
  "TPR: X3",
  "TPR: X4",
  "TPR: X5",
  "CI Coverage: X1",
  "CI Coverage: X2",
  "CI Coverage: X3",
  "CI Coverage: X4",
  "CI Coverage: X5",
  "Type II Error: X1",
  "Type II Error: X2",
  "Type II Error: X3",
  "Type II Error: X4",
  "Type II Error: X5",
  "Type I Error: X6-X20",
  "Bias: X1",
  "Bias: X2",
  "Bias: X3",
  "Bias: X4",
  "Bias: X5",
  "Bias: X6-X20"
)

# fdr figure case 1a
fdr_fig1_dat = data.frame(fdr)
fdr_fig1_dat$method = c("LRT", "AIC", "BIC", "cv.LASSO", "fixed.LASSO", "cv.ENET", "fixed.ENET")

fdr_fig1 = ggplot(fdr_fig1_dat, aes(y = fdr, x = method, fill=method)) +
  geom_bar(stat="identity") +
  theme(legend.position = "none") +
  ylim(0, 0.5)

# fpr figure case 1a
fpr_fig1_dat = data.frame(fpr)
fpr_fig1_dat$method = c("LRT", "AIC", "BIC", "cv.LASSO", "fixed.LASSO", "cv.ENET", "fixed.ENET")

fpr_fig1 = ggplot(fpr_fig1_dat, aes(y = fpr, x = method, fill=method)) +
  geom_bar(stat="identity") +
  theme(legend.position = "none") +
  ylim(0, 0.5)

# toe figure case 1a
toe_fig1_dat = data.frame(toe)
toe_fig1_dat$method = c("LRT", "AIC", "BIC", "cv.LASSO", "fixed.LASSO", "cv.ENET", "fixed.ENET")

toe_fig1 = ggplot(toe_fig1_dat, aes(y = toe, x = method, fill=method)) +
  geom_bar(stat="identity") +
  theme(legend.position = "none") +
  ylim(0, 1)

# bias0 figure case 1a
bias0_fig1_dat = data.frame(bias0)
bias0_fig1_dat$method = c("LRT", "AIC", "BIC", "cv.LASSO", "fixed.LASSO", "cv.ENET", "fixed.ENET")

bias0_fig1 = ggplot(bias0_fig1_dat, aes(y = bias0, x = method, fill=method)) +
  geom_bar(stat="identity") +
  theme(legend.position = "none") +
  ylim(0, 0.1)

# tpr figure case 1a
tpr_fig1_dat = data.frame(tpr)
tpr_fig1_dat$id = 1:nrow(tpr_fig1_dat)
tpr_fig1_dat$method = c("LRT", "AIC", "BIC", "cv-LASSO", "fixed-LASSO", "cv-ENET", "fixed-ENET")
tpr_fig1_dat = tpr_fig1_dat %>%
  pivot_longer(
    cols = c(X1, X2, X3, X4, X5), 
    names_to = "method",
    values_to = "tpr",
    names_repair = "unique"
  )
colnames(tpr_fig1_dat) = c("id", "method", "variable", "tpr")

tpr_fig1 = ggplot(tpr_fig1_dat, aes(y = tpr, x = method, color=variable)) +
  geom_bar(aes(fill = variable), stat="identity", position = "dodge") +
  theme(legend.position = "none") +
  ylim(0, 1)

# coverage figure case 1a
coverage_fig1_dat = data.frame(coverage)
coverage_fig1_dat$id = 1:nrow(coverage_fig1_dat)
coverage_fig1_dat$method = c("LRT", "AIC", "BIC", "cv-LASSO", "fixed-LASSO", "cv-ENET", "fixed-ENET")
coverage_fig1_dat = coverage_fig1_dat %>%
  pivot_longer(
    cols = c(X1, X2, X3, X4, X5), 
    names_to = "method",
    values_to = "tpr",
    names_repair = "unique"
  )
colnames(coverage_fig1_dat) = c("id", "method", "variable", "coverage")

coverage_fig1 = ggplot(coverage_fig1_dat, aes(y = coverage, x = method, color=variable)) +
  geom_bar(aes(fill = variable), stat="identity", position = "dodge") +
  theme(legend.position = "none") +
  ylim(0, 1)

# type II error figure case 1a
tte_fig1_dat = data.frame(tte)
tte_fig1_dat$id = 1:nrow(tte_fig1_dat)
tte_fig1_dat$method = c("LRT", "AIC", "BIC", "cv-LASSO", "fixed-LASSO", "cv-ENET", "fixed-ENET")
tte_fig1_dat = tte_fig1_dat %>%
  pivot_longer(
    cols = c(X1, X2, X3, X4, X5), 
    names_to = "method",
    values_to = "tpr",
    names_repair = "unique"
  )
colnames(tte_fig1_dat) = c("id", "method", "variable", "tte")

tte_fig1 = ggplot(tte_fig1_dat, aes(y = tte, x = method, group = variable)) +
  geom_bar(aes(fill = variable), stat="identity", position = "dodge") +
  theme(legend.position = "none") +
  ylim(0, 0.6)

# bias figure case 1a
bias_fig1_dat = data.frame(bias)
bias_fig1_dat$id = 1:nrow(bias_fig1_dat)
bias_fig1_dat$method = c("LRT", "AIC", "BIC", "cv-LASSO", "fixed-LASSO", "cv-ENET", "fixed-ENET")
bias_fig1_dat = bias_fig1_dat %>%
  pivot_longer(
    cols = c(X1, X2, X3, X4, X5), 
    names_to = "method",
    values_to = "tpr",
    names_repair = "unique"
  )
colnames(bias_fig1_dat) = c("id", "method", "variable", "bias")

bias_fig1 = ggplot(bias_fig1_dat, aes(y = bias, x = method, group = variable)) +
  geom_bar(aes(fill = variable), stat="identity", position = "dodge") +
  theme(legend.position = "none") +
  ylim(0, 0.2)

########
#  1b  #
########
########################
#     read in data     #
########################

p_betas = read.delim("../data_raw/p_betas_1b.txt")
p_ci = read.delim("../data_raw/p_ci_1b.txt")
p_p = read.delim("../data_raw/p_p_1b.txt")

aic_betas = read.delim("../data_raw/aic_betas_1b.txt")
aic_ci = read.delim("../data_raw/aic_ci_1b.txt")
aic_p = read.delim("../data_raw/aic_p_1b.txt")

bic_betas = read.delim("../data_raw/bic_betas_1b.txt")
bic_ci = read.delim("../data_raw/bic_ci_1b.txt")
bic_p = read.delim("../data_raw/bic_p_1b.txt")

cv_lasso_betas = read.delim("../data_raw/cv_lasso_betas_1b.txt")
cv_lasso_ci = read.delim("../data_raw/cv_lasso_ci_1b.txt")
cv_lasso_p = read.delim("../data_raw/cv_lasso_p_1b.txt")

fixed_lasso_betas = read.delim("../data_raw/fixed_lasso_betas_1b.txt")
fixed_lasso_ci = read.delim("../data_raw/fixed_lasso_ci_1b.txt")
fixed_lasso_p = read.delim("../data_raw/fixed_lasso_p_1b.txt")

cv_enet_betas = read.delim("../data_raw/cv_enet_betas_1b.txt")
cv_enet_ci = read.delim("../data_raw/cv_enet_ci_1b.txt")
cv_enet_p = read.delim("../data_raw/cv_enet_p_1b.txt")

fixed_enet_betas = read.delim("../data_raw/fixed_enet_betas_1b.txt")
fixed_enet_ci = read.delim("../data_raw/fixed_enet_ci_1b.txt")
fixed_enet_p = read.delim("../data_raw/fixed_enet_p_1b.txt")

########################
# calculate parameters #
########################
beta = vector("double", 20)
beta[1:5] = c(0.5, 1., 1.5, 2., 2.5) / 3

# true positive rate
tpr = matrix(0, nrow = 7, ncol = 5)
for (i in 1:dim(p_betas[, 1:5])[2]) {
  tpr[1, i] = truePositiveRate(p_betas[, i])
}

for (i in 1:dim(aic_betas[, 1:5])[2]) {
  tpr[2, i] = truePositiveRate(aic_betas[, i])
}

for (i in 1:dim(bic_betas[, 1:5])[2]) {
  tpr[3, i] = truePositiveRate(bic_betas[, i])
}

for (i in 1:dim(cv_lasso_betas[, 1:5])[2]) {
  tpr[4, i] = truePositiveRate(cv_lasso_betas[, i])
}

for (i in 1:dim(fixed_lasso_betas[, 1:5])[2]) {
  tpr[5, i] = truePositiveRate(fixed_lasso_betas[, i])
}

for (i in 1:dim(cv_enet_betas[, 1:5])[2]) {
  tpr[6, i] = truePositiveRate(cv_enet_betas[, i])
}

for (i in 1:dim(fixed_enet_betas[, 1:5])[2]) {
  tpr[7, i] = truePositiveRate(fixed_enet_betas[, i])
}

# false positive rate
fpr = vector("double", 7)
fpr[1] = falsePositiveRate(p_betas[, 5:15])
fpr[2] = falsePositiveRate(aic_betas[, 5:15])
fpr[3] = falsePositiveRate(bic_betas[, 5:15])
fpr[4] = falsePositiveRate(cv_lasso_betas[, 5:15])
fpr[5] = falsePositiveRate(fixed_lasso_betas[, 5:15])
fpr[6] = falsePositiveRate(cv_enet_betas[, 5:15])
fpr[7] = falsePositiveRate(fixed_enet_betas[, 5:15])

# false discovery rate
fdr = vector("double", 7)
fdr[1] = falseDiscoveryRate(p_betas, p_betas[, 5:15])
fdr[2] = falseDiscoveryRate(aic_betas, aic_betas[, 5:15])
fdr[3] = falseDiscoveryRate(bic_betas, bic_betas[, 5:15])
fdr[4] = falseDiscoveryRate(cv_lasso_betas, cv_lasso_betas[, 5:15])
fdr[5] = falseDiscoveryRate(fixed_lasso_betas, fixed_lasso_betas[, 5:15])
fdr[6] = falseDiscoveryRate(cv_enet_betas, cv_enet_betas[, 5:15])
fdr[7] = falseDiscoveryRate(fixed_enet_betas, fixed_enet_betas[, 5:15])

# type I error
toe = vector("double", 7)
toe[1] = typeOneErrorRate(p_p[, 5:15])
toe[2] = typeOneErrorRate(aic_p[, 5:15])
toe[3] = typeOneErrorRate(bic_p[, 5:15])
toe[4] = typeOneErrorRate(cv_lasso_p[, 5:15])
toe[5] = typeOneErrorRate(fixed_lasso_p[, 5:15])
toe[6] = typeOneErrorRate(cv_enet_p[, 5:15])
toe[7] = typeOneErrorRate(fixed_enet_p[, 5:15])

# type II error rate
tte = matrix(0, nrow = 7, ncol = 5)
for (i in 1:dim(p_p[, 1:5])[2]) {
  tte[1, i] = typeTwoErrorRate(p_p[, i])
}

for (i in 1:dim(aic_p[, 1:5])[2]) {
  tte[2, i] = typeTwoErrorRate(aic_p[, i])
}

for (i in 1:dim(bic_p[, 1:5])[2]) {
  tte[3, i] = typeTwoErrorRate(bic_p[, i])
}

for (i in 1:dim(cv_lasso_p[, 1:5])[2]) {
  tte[4, i] = typeTwoErrorRate(cv_lasso_p[, i])
}

for (i in 1:dim(fixed_lasso_p[, 1:5])[2]) {
  tte[5, i] = typeTwoErrorRate(fixed_lasso_p[, i])
}

for (i in 1:dim(cv_enet_p[, 1:5])[2]) {
  tte[6, i] = typeTwoErrorRate(cv_enet_p[, i])
}

for (i in 1:dim(fixed_enet_p[, 1:5])[2]) {
  tte[7, i] = typeTwoErrorRate(fixed_enet_p[, i])
}

# parameters 1-5 bias
bias = matrix(0, nrow = 7, ncol = 5)
for (i in 1:dim(p_betas[, 1:5])[2]) {
  bias[1, i] = parameterBiasNonZero(p_betas[, i], beta[i])
}

for (i in 1:dim(aic_betas[, 1:5])[2]) {
  bias[2, i] = parameterBiasNonZero(aic_betas[, i], beta[i])
}

for (i in 1:dim(bic_betas[, 1:5])[2]) {
  bias[3, i] = parameterBiasNonZero(bic_betas[, i], beta[i])
}

for (i in 1:dim(cv_lasso_betas[, 1:5])[2]) {
  bias[4, i] = parameterBiasNonZero(cv_lasso_betas[, i], beta[i])
}

for (i in 1:dim(fixed_lasso_betas[, 1:5])[2]) {
  bias[5, i] = parameterBiasNonZero(fixed_lasso_betas[, i], beta[i])
}

for (i in 1:dim(cv_enet_betas[, 1:5])[2]) {
  bias[6, i] = parameterBiasNonZero(cv_enet_betas[, i], beta[i])
}

for (i in 1:dim(fixed_enet_betas[, 1:5])[2]) {
  bias[7, i] = parameterBiasNonZero(fixed_enet_betas[, i], beta[i])
}

# parameters 5-20 average bias
bias0 = vector("double", 7)
bias0[1] = parameterBiasZero(p_betas[, 5:15])
bias0[2] = parameterBiasZero(aic_betas[, 5:15])
bias0[3] = parameterBiasZero(bic_betas[, 5:15])
bias0[4] = parameterBiasZero(cv_lasso_betas[, 5:15])
bias0[5] = parameterBiasZero(fixed_lasso_betas[, 5:15])
bias0[6] = parameterBiasZero(cv_enet_betas[, 5:15])
bias0[7] = parameterBiasZero(fixed_enet_betas[, 5:15])

# parameters 1-5 CI coverage
coverage = matrix(0, nrow = 7, ncol = 5)
for (i in 1:dim(p_ci[, 1:5])[2]) {
  coverage[1, i] = coverageCI(p_ci[seq(1, ncol(p_ci) - 1, 2), i], p_ci[seq(2, ncol(p_ci), 2), i], p_betas[, i], beta[i])
}

for (i in 1:dim(aic_ci[, 1:5])[2]) {
  coverage[2, i] = coverageCI(aic_ci[seq(1, ncol(aic_ci) - 1, 2), i], aic_ci[seq(2, ncol(aic_ci), 2), i], aic_betas[, i], beta[i])
}

for (i in 1:dim(bic_ci[, 1:5])[2]) {
  coverage[3, i] = coverageCI(bic_ci[seq(1, ncol(bic_ci) - 1, 2), i], bic_ci[seq(2, ncol(bic_ci), 2), i], bic_betas[, i], beta[i])
}

for (i in 1:dim(cv_lasso_ci[, 1:5])[2]) {
  coverage[4, i] = coverageCI(cv_lasso_ci[seq(1, ncol(cv_lasso_ci) - 1, 2), i], cv_lasso_ci[seq(2, ncol(cv_lasso_ci), 2), i], cv_lasso_betas[, i], beta[i])
}

for (i in 1:dim(fixed_lasso_ci[, 1:5])[2]) {
  coverage[5, i] = coverageCI(fixed_lasso_ci[seq(1, ncol(fixed_lasso_ci) - 1, 2), i], fixed_lasso_ci[seq(2, ncol(fixed_lasso_ci), 2), i], fixed_lasso_betas[, i], beta[i])
}

for (i in 1:dim(cv_enet_ci[, 1:5])[2]) {
  coverage[6, i] = coverageCI(cv_enet_ci[seq(1, ncol(cv_enet_ci) - 1, 2), i], cv_enet_ci[seq(2, ncol(cv_enet_ci), 2), i], cv_enet_betas[, 1:5], beta[i])
}

for (i in 1:dim(fixed_enet_ci[, 1:5])[2]) {
  coverage[7, i] = coverageCI(fixed_enet_ci[seq(1, ncol(fixed_enet_ci) - 1, 2), i], fixed_enet_ci[seq(2, ncol(fixed_enet_ci), 2), i], fixed_enet_betas[, i], beta[i])
}

################
# build table #
################

table_1b = data.frame(
  "LRT" = c(
    fdr[1], fpr[1], tpr[1, 1], tpr[1, 2], tpr[1, 3],
    tpr[1, 4], tpr[1, 5], coverage[1, 1], coverage[1, 2],
    coverage[1, 3], coverage[1, 4], coverage[1, 5],
    tte[1, 1], tte[1, 2], tte[1, 3], tte[1, 4], tte[1, 5],
    toe[1], bias[1, 1], bias[1, 2], bias[1, 3], bias[1, 4],
    bias[1, 5], bias0[1]
  ),
  "AIC" = c(
    fdr[2], fpr[2], tpr[2, 1], tpr[2, 2], tpr[2, 3],
    tpr[2, 4], tpr[2, 5], coverage[2, 1], coverage[2, 2],
    coverage[2, 3], coverage[2, 4], coverage[2, 5],
    tte[2, 1], tte[2, 2], tte[2, 3], tte[2, 4], tte[2, 5],
    toe[1], bias[2, 1], bias[2, 2], bias[2, 3], bias[2, 4],
    bias[2, 5], bias0[2]
  ),
  "BIC" = c(
    fdr[3], fpr[3], tpr[3, 1], tpr[3, 2], tpr[3, 3],
    tpr[3, 4], tpr[3, 5], coverage[3, 1], coverage[3, 2],
    coverage[3, 3], coverage[3, 4], coverage[3, 5],
    tte[3, 1], tte[3, 2], tte[3, 3], tte[3, 4], tte[3, 5],
    toe[1], bias[3, 1], bias[3, 2], bias[3, 3], bias[3, 4],
    bias[3, 5], bias0[3]
  ),
  "cv-LASSO" = c(
    fdr[4], fpr[4], tpr[4, 1], tpr[4, 2], tpr[4, 3],
    tpr[4, 4], tpr[4, 5], coverage[4, 1], coverage[4, 2],
    coverage[4, 3], coverage[4, 4], coverage[4, 5],
    tte[4, 1], tte[4, 2], tte[4, 3], tte[4, 4], tte[4, 5],
    toe[1], bias[4, 1], bias[4, 2], bias[4, 3], bias[4, 4],
    bias[4, 5], bias0[4]
  ),
  "fixed-LASSO" = c(
    fdr[5], fpr[5], tpr[5, 1], tpr[5, 2], tpr[5, 3],
    tpr[5, 4], tpr[5, 5], coverage[5, 1], coverage[5, 2],
    coverage[5, 3], coverage[5, 4], coverage[5, 5],
    tte[5, 1], tte[5, 2], tte[5, 3], tte[5, 4], tte[5, 5],
    toe[1], bias[5, 1], bias[5, 2], bias[5, 3], bias[5, 4],
    bias[5, 5], bias0[5]
  ),
  "cv-ENET" = c(
    fdr[6], fpr[6], tpr[6, 1], tpr[6, 2], tpr[6, 3],
    tpr[6, 4], tpr[6, 5], coverage[6, 1], coverage[6, 2],
    coverage[6, 3], coverage[6, 4], coverage[6, 5],
    tte[6, 1], tte[6, 2], tte[6, 3], tte[6, 4], tte[6, 5],
    toe[1], bias[6, 1], bias[6, 2], bias[6, 3], bias[6, 4],
    bias[6, 5], bias0[6]
  ),
  "fixed-ENET" = c(
    fdr[7], fpr[7], tpr[7, 1], tpr[7, 2], tpr[7, 3],
    tpr[7, 4], tpr[7, 5], coverage[7, 1], coverage[7, 2],
    coverage[7, 3], coverage[7, 4], coverage[7, 5],
    tte[7, 1], tte[7, 2], tte[7, 3], tte[7, 4], tte[7, 5],
    toe[1], bias[7, 1], bias[7, 2], bias[7, 3], bias[7, 4],
    bias[7, 5], bias0[7]
  )
)

rownames(table_1b) = c(
  "FDR",
  "FPR",
  "TPR: X1",
  "TPR: X2",
  "TPR: X3",
  "TPR: X4",
  "TPR: X5",
  "CI Coverage: X1",
  "CI Coverage: X2",
  "CI Coverage: X3",
  "CI Coverage: X4",
  "CI Coverage: X5",
  "Type II Error: X1",
  "Type II Error: X2",
  "Type II Error: X3",
  "Type II Error: X4",
  "Type II Error: X5",
  "Type I Error: X6-X20",
  "Bias: X1",
  "Bias: X2",
  "Bias: X3",
  "Bias: X4",
  "Bias: X5",
  "Bias: X6-X20"
)

# fdr figure case 1b
fdr_fig2_dat = data.frame(fdr)
fdr_fig2_dat$method = c("LRT", "AIC", "BIC", "cv.LASSO", "fixed.LASSO", "cv.ENET", "fixed.ENET")

fdr_fig2 = ggplot(fdr_fig2_dat, aes(y = fdr, x = method, fill=method)) +
  geom_bar(stat="identity") +
  theme(legend.position = "none") +
  ylim(0, 0.5)

# fpr figure case 1b
fpr_fig2_dat = data.frame(fpr)
fpr_fig2_dat$method = c("LRT", "AIC", "BIC", "cv.LASSO", "fixed.LASSO", "cv.ENET", "fixed.ENET")

fpr_fig2 = ggplot(fpr_fig2_dat, aes(y = fpr, x = method, fill=method)) +
  geom_bar(stat="identity") +
  theme(legend.position = "none") +
  ylim(0, 0.5)

# toe figure case 1b
toe_fig2_dat = data.frame(toe)
toe_fig2_dat$method = c("LRT", "AIC", "BIC", "cv.LASSO", "fixed.LASSO", "cv.ENET", "fixed.ENET")

toe_fig2 = ggplot(toe_fig2_dat, aes(y = toe, x = method, fill=method)) +
  geom_bar(stat="identity") +
  theme(legend.position = "none") +
  ylim(0, 1)

# bias0 figure case 1b
bias0_fig2_dat = data.frame(bias0)
bias0_fig2_dat$method = c("LRT", "AIC", "BIC", "cv.LASSO", "fixed.LASSO", "cv.ENET", "fixed.ENET")

bias0_fig2 = ggplot(bias0_fig2_dat, aes(y = bias0, x = method, fill=method)) +
  geom_bar(stat="identity") +
  theme(legend.position = "none") +
  ylim(0, 0.1)

# tpr figure case 1b
tpr_fig2_dat = data.frame(tpr)
tpr_fig2_dat$id = 1:nrow(tpr_fig2_dat)
tpr_fig2_dat$method = c("LRT", "AIC", "BIC", "cv-LASSO", "fixed-LASSO", "cv-ENET", "fixed-ENET")
tpr_fig2_dat = tpr_fig2_dat %>%
  pivot_longer(
    cols = c(X1, X2, X3, X4, X5), 
    names_to = "method",
    values_to = "tpr",
    names_repair = "unique"
  )
colnames(tpr_fig2_dat) = c("id", "method", "variable", "tpr")

tpr_fig2 = ggplot(tpr_fig2_dat, aes(y = tpr, x = method, color=variable)) +
  geom_bar(aes(fill = variable), stat="identity", position = "dodge") +
  theme(legend.position = "none") +
  ylim(0, 1)

# coverage figure case 1b
coverage_fig2_dat = data.frame(coverage)
coverage_fig2_dat$id = 1:nrow(coverage_fig2_dat)
coverage_fig2_dat$method = c("LRT", "AIC", "BIC", "cv-LASSO", "fixed-LASSO", "cv-ENET", "fixed-ENET")
coverage_fig2_dat = coverage_fig2_dat %>%
  pivot_longer(
    cols = c(X1, X2, X3, X4, X5), 
    names_to = "method",
    values_to = "tpr",
    names_repair = "unique"
  )
colnames(coverage_fig2_dat) = c("id", "method", "variable", "coverage")

coverage_fig2 = ggplot(coverage_fig2_dat, aes(y = coverage, x = method, color=variable)) +
  geom_bar(aes(fill = variable), stat="identity", position = "dodge") +
  theme(legend.position = "none") +
  ylim(0, 1)

# type II error figure case 1b
tte_fig2_dat = data.frame(tte)
tte_fig2_dat$id = 1:nrow(tte_fig2_dat)
tte_fig2_dat$method = c("LRT", "AIC", "BIC", "cv-LASSO", "fixed-LASSO", "cv-ENET", "fixed-ENET")
tte_fig2_dat = tte_fig2_dat %>%
  pivot_longer(
    cols = c(X1, X2, X3, X4, X5), 
    names_to = "method",
    values_to = "tpr",
    names_repair = "unique"
  )
colnames(tte_fig2_dat) = c("id", "method", "variable", "tte")

tte_fig2 = ggplot(tte_fig2_dat, aes(y = tte, x = method, group = variable)) +
  geom_bar(aes(fill = variable), stat="identity", position = "dodge") +
  theme(legend.position = "none") +
  ylim(0, 0.6)

# bias figure case 1b
bias_fig2_dat = data.frame(bias)
bias_fig2_dat$id = 1:nrow(bias_fig2_dat)
bias_fig2_dat$method = c("LRT", "AIC", "BIC", "cv-LASSO", "fixed-LASSO", "cv-ENET", "fixed-ENET")
bias_fig2_dat = bias_fig2_dat %>%
  pivot_longer(
    cols = c(X1, X2, X3, X4, X5), 
    names_to = "method",
    values_to = "tpr",
    names_repair = "unique"
  )
colnames(bias_fig2_dat) = c("id", "method", "variable", "bias")

bias_fig2 = ggplot(bias_fig2_dat, aes(y = bias, x = method, group = variable)) +
  geom_bar(aes(fill = variable), stat="identity", position = "dodge") +
  theme(legend.position = "none") +
  ylim(0, 0.2)

########
#  2a  #
########
########################
#     read in data     #
########################

p_betas = read.delim("../data_raw/p_betas_2a.txt")
p_ci = read.delim("../data_raw/p_ci_2a.txt")
p_p = read.delim("../data_raw/p_p_2a.txt")

aic_betas = read.delim("../data_raw/aic_betas_2a.txt")
aic_ci = read.delim("../data_raw/aic_ci_2a.txt")
aic_p = read.delim("../data_raw/aic_p_2a.txt")

bic_betas = read.delim("../data_raw/bic_betas_2a.txt")
bic_ci = read.delim("../data_raw/bic_ci_2a.txt")
bic_p = read.delim("../data_raw/bic_p_2a.txt")

cv_lasso_betas = read.delim("../data_raw/cv_lasso_betas_2a.txt")
cv_lasso_ci = read.delim("../data_raw/cv_lasso_ci_2a.txt")
cv_lasso_p = read.delim("../data_raw/cv_lasso_p_2a.txt")

fixed_lasso_betas = read.delim("../data_raw/fixed_lasso_betas_2a.txt")
fixed_lasso_ci = read.delim("../data_raw/fixed_lasso_ci_2a.txt")
fixed_lasso_p = read.delim("../data_raw/fixed_lasso_p_2a.txt")

cv_enet_betas = read.delim("../data_raw/cv_enet_betas_2a.txt")
cv_enet_ci = read.delim("../data_raw/cv_enet_ci_2a.txt")
cv_enet_p = read.delim("../data_raw/cv_enet_p_2a.txt")

fixed_enet_betas = read.delim("../data_raw/fixed_enet_betas_2a.txt")
fixed_enet_ci = read.delim("../data_raw/fixed_enet_ci_2a.txt")
fixed_enet_p = read.delim("../data_raw/fixed_enet_p_2a.txt")

########################
# calculate parameters #
########################
beta = vector("double", 20)
beta[1:5] = c(0.5, 1., 1.5, 2., 2.5) / 3

# true positive rate
tpr = matrix(0, nrow = 7, ncol = 5)
for (i in 1:dim(p_betas[, 1:5])[2]) {
  tpr[1, i] = truePositiveRate(p_betas[, i])
}

for (i in 1:dim(aic_betas[, 1:5])[2]) {
  tpr[2, i] = truePositiveRate(aic_betas[, i])
}

for (i in 1:dim(bic_betas[, 1:5])[2]) {
  tpr[3, i] = truePositiveRate(bic_betas[, i])
}

for (i in 1:dim(cv_lasso_betas[, 1:5])[2]) {
  tpr[4, i] = truePositiveRate(cv_lasso_betas[, i])
}

for (i in 1:dim(fixed_lasso_betas[, 1:5])[2]) {
  tpr[5, i] = truePositiveRate(fixed_lasso_betas[, i])
}

for (i in 1:dim(cv_enet_betas[, 1:5])[2]) {
  tpr[6, i] = truePositiveRate(cv_enet_betas[, i])
}

for (i in 1:dim(fixed_enet_betas[, 1:5])[2]) {
  tpr[7, i] = truePositiveRate(fixed_enet_betas[, i])
}

# false positive rate
fpr = vector("double", 7)
fpr[1] = falsePositiveRate(p_betas[, 5:15])
fpr[2] = falsePositiveRate(aic_betas[, 5:15])
fpr[3] = falsePositiveRate(bic_betas[, 5:15])
fpr[4] = falsePositiveRate(cv_lasso_betas[, 5:15])
fpr[5] = falsePositiveRate(fixed_lasso_betas[, 5:15])
fpr[6] = falsePositiveRate(cv_enet_betas[, 5:15])
fpr[7] = falsePositiveRate(fixed_enet_betas[, 5:15])

# false discovery rate
fdr = vector("double", 7)
fdr[1] = falseDiscoveryRate(p_betas, p_betas[, 5:15])
fdr[2] = falseDiscoveryRate(aic_betas, aic_betas[, 5:15])
fdr[3] = falseDiscoveryRate(bic_betas, bic_betas[, 5:15])
fdr[4] = falseDiscoveryRate(cv_lasso_betas, cv_lasso_betas[, 5:15])
fdr[5] = falseDiscoveryRate(fixed_lasso_betas, fixed_lasso_betas[, 5:15])
fdr[6] = falseDiscoveryRate(cv_enet_betas, cv_enet_betas[, 5:15])
fdr[7] = falseDiscoveryRate(fixed_enet_betas, fixed_enet_betas[, 5:15])

# type I error
toe = vector("double", 7)
toe[1] = typeOneErrorRate(p_p[, 5:15])
toe[2] = typeOneErrorRate(aic_p[, 5:15])
toe[3] = typeOneErrorRate(bic_p[, 5:15])
toe[4] = typeOneErrorRate(cv_lasso_p[, 5:15])
toe[5] = typeOneErrorRate(fixed_lasso_p[, 5:15])
toe[6] = typeOneErrorRate(cv_enet_p[, 5:15])
toe[7] = typeOneErrorRate(fixed_enet_p[, 5:15])

# type II error rate
tte = matrix(0, nrow = 7, ncol = 5)
for (i in 1:dim(p_p[, 1:5])[2]) {
  tte[1, i] = typeTwoErrorRate(p_p[, i])
}

for (i in 1:dim(aic_p[, 1:5])[2]) {
  tte[2, i] = typeTwoErrorRate(aic_p[, i])
}

for (i in 1:dim(bic_p[, 1:5])[2]) {
  tte[3, i] = typeTwoErrorRate(bic_p[, i])
}

for (i in 1:dim(cv_lasso_p[, 1:5])[2]) {
  tte[4, i] = typeTwoErrorRate(cv_lasso_p[, i])
}

for (i in 1:dim(fixed_lasso_p[, 1:5])[2]) {
  tte[5, i] = typeTwoErrorRate(fixed_lasso_p[, i])
}

for (i in 1:dim(cv_enet_p[, 1:5])[2]) {
  tte[6, i] = typeTwoErrorRate(cv_enet_p[, i])
}

for (i in 1:dim(fixed_enet_p[, 1:5])[2]) {
  tte[7, i] = typeTwoErrorRate(fixed_enet_p[, i])
}

# parameters 1-5 bias
bias = matrix(0, nrow = 7, ncol = 5)
for (i in 1:dim(p_betas[, 1:5])[2]) {
  bias[1, i] = parameterBiasNonZero(p_betas[, i], beta[i])
}

for (i in 1:dim(aic_betas[, 1:5])[2]) {
  bias[2, i] = parameterBiasNonZero(aic_betas[, i], beta[i])
}

for (i in 1:dim(bic_betas[, 1:5])[2]) {
  bias[3, i] = parameterBiasNonZero(bic_betas[, i], beta[i])
}

for (i in 1:dim(cv_lasso_betas[, 1:5])[2]) {
  bias[4, i] = parameterBiasNonZero(cv_lasso_betas[, i], beta[i])
}

for (i in 1:dim(fixed_lasso_betas[, 1:5])[2]) {
  bias[5, i] = parameterBiasNonZero(fixed_lasso_betas[, i], beta[i])
}

for (i in 1:dim(cv_enet_betas[, 1:5])[2]) {
  bias[6, i] = parameterBiasNonZero(cv_enet_betas[, i], beta[i])
}

for (i in 1:dim(fixed_enet_betas[, 1:5])[2]) {
  bias[7, i] = parameterBiasNonZero(fixed_enet_betas[, i], beta[i])
}

# parameters 5-20 average bias
bias0 = vector("double", 7)
bias0[1] = parameterBiasZero(p_betas[, 5:15])
bias0[2] = parameterBiasZero(aic_betas[, 5:15])
bias0[3] = parameterBiasZero(bic_betas[, 5:15])
bias0[4] = parameterBiasZero(cv_lasso_betas[, 5:15])
bias0[5] = parameterBiasZero(fixed_lasso_betas[, 5:15])
bias0[6] = parameterBiasZero(cv_enet_betas[, 5:15])
bias0[7] = parameterBiasZero(fixed_enet_betas[, 5:15])

# parameters 1-5 CI coverage
coverage = matrix(0, nrow = 7, ncol = 5)
for (i in 1:dim(p_ci[, 1:5])[2]) {
  coverage[1, i] = coverageCI(p_ci[seq(1, ncol(p_ci) - 1, 2), i], p_ci[seq(2, ncol(p_ci), 2), i], p_betas[, i], beta[i])
}

for (i in 1:dim(aic_ci[, 1:5])[2]) {
  coverage[2, i] = coverageCI(aic_ci[seq(1, ncol(aic_ci) - 1, 2), i], aic_ci[seq(2, ncol(aic_ci), 2), i], aic_betas[, i], beta[i])
}

for (i in 1:dim(bic_ci[, 1:5])[2]) {
  coverage[3, i] = coverageCI(bic_ci[seq(1, ncol(bic_ci) - 1, 2), i], bic_ci[seq(2, ncol(bic_ci), 2), i], bic_betas[, i], beta[i])
}

for (i in 1:dim(cv_lasso_ci[, 1:5])[2]) {
  coverage[4, i] = coverageCI(cv_lasso_ci[seq(1, ncol(cv_lasso_ci) - 1, 2), i], cv_lasso_ci[seq(2, ncol(cv_lasso_ci), 2), i], cv_lasso_betas[, i], beta[i])
}

for (i in 1:dim(fixed_lasso_ci[, 1:5])[2]) {
  coverage[5, i] = coverageCI(fixed_lasso_ci[seq(1, ncol(fixed_lasso_ci) - 1, 2), i], fixed_lasso_ci[seq(2, ncol(fixed_lasso_ci), 2), i], fixed_lasso_betas[, i], beta[i])
}

for (i in 1:dim(cv_enet_ci[, 1:5])[2]) {
  coverage[6, i] = coverageCI(cv_enet_ci[seq(1, ncol(cv_enet_ci) - 1, 2), i], cv_enet_ci[seq(2, ncol(cv_enet_ci), 2), i], cv_enet_betas[, 1:5], beta[i])
}

for (i in 1:dim(fixed_enet_ci[, 1:5])[2]) {
  coverage[7, i] = coverageCI(fixed_enet_ci[seq(1, ncol(fixed_enet_ci) - 1, 2), i], fixed_enet_ci[seq(2, ncol(fixed_enet_ci), 2), i], fixed_enet_betas[, i], beta[i])
}

################
# build table #
################

table_2a = data.frame(
  "LRT" = c(
    fdr[1], fpr[1], tpr[1, 1], tpr[1, 2], tpr[1, 3],
    tpr[1, 4], tpr[1, 5], coverage[1, 1], coverage[1, 2],
    coverage[1, 3], coverage[1, 4], coverage[1, 5],
    tte[1, 1], tte[1, 2], tte[1, 3], tte[1, 4], tte[1, 5],
    toe[1], bias[1, 1], bias[1, 2], bias[1, 3], bias[1, 4],
    bias[1, 5], bias0[1]
  ),
  "AIC" = c(
    fdr[2], fpr[2], tpr[2, 1], tpr[2, 2], tpr[2, 3],
    tpr[2, 4], tpr[2, 5], coverage[2, 1], coverage[2, 2],
    coverage[2, 3], coverage[2, 4], coverage[2, 5],
    tte[2, 1], tte[2, 2], tte[2, 3], tte[2, 4], tte[2, 5],
    toe[1], bias[2, 1], bias[2, 2], bias[2, 3], bias[2, 4],
    bias[2, 5], bias0[2]
  ),
  "BIC" = c(
    fdr[3], fpr[3], tpr[3, 1], tpr[3, 2], tpr[3, 3],
    tpr[3, 4], tpr[3, 5], coverage[3, 1], coverage[3, 2],
    coverage[3, 3], coverage[3, 4], coverage[3, 5],
    tte[3, 1], tte[3, 2], tte[3, 3], tte[3, 4], tte[3, 5],
    toe[1], bias[3, 1], bias[3, 2], bias[3, 3], bias[3, 4],
    bias[3, 5], bias0[3]
  ),
  "cv-LASSO" = c(
    fdr[4], fpr[4], tpr[4, 1], tpr[4, 2], tpr[4, 3],
    tpr[4, 4], tpr[4, 5], coverage[4, 1], coverage[4, 2],
    coverage[4, 3], coverage[4, 4], coverage[4, 5],
    tte[4, 1], tte[4, 2], tte[4, 3], tte[4, 4], tte[4, 5],
    toe[1], bias[4, 1], bias[4, 2], bias[4, 3], bias[4, 4],
    bias[4, 5], bias0[4]
  ),
  "fixed-LASSO" = c(
    fdr[5], fpr[5], tpr[5, 1], tpr[5, 2], tpr[5, 3],
    tpr[5, 4], tpr[5, 5], coverage[5, 1], coverage[5, 2],
    coverage[5, 3], coverage[5, 4], coverage[5, 5],
    tte[5, 1], tte[5, 2], tte[5, 3], tte[5, 4], tte[5, 5],
    toe[1], bias[5, 1], bias[5, 2], bias[5, 3], bias[5, 4],
    bias[5, 5], bias0[5]
  ),
  "cv-ENET" = c(
    fdr[6], fpr[6], tpr[6, 1], tpr[6, 2], tpr[6, 3],
    tpr[6, 4], tpr[6, 5], coverage[6, 1], coverage[6, 2],
    coverage[6, 3], coverage[6, 4], coverage[6, 5],
    tte[6, 1], tte[6, 2], tte[6, 3], tte[6, 4], tte[6, 5],
    toe[1], bias[6, 1], bias[6, 2], bias[6, 3], bias[6, 4],
    bias[6, 5], bias0[6]
  ),
  "fixed-ENET" = c(
    fdr[7], fpr[7], tpr[7, 1], tpr[7, 2], tpr[7, 3],
    tpr[7, 4], tpr[7, 5], coverage[7, 1], coverage[7, 2],
    coverage[7, 3], coverage[7, 4], coverage[7, 5],
    tte[7, 1], tte[7, 2], tte[7, 3], tte[7, 4], tte[7, 5],
    toe[1], bias[7, 1], bias[7, 2], bias[7, 3], bias[7, 4],
    bias[7, 5], bias0[7]
  )
)

rownames(table_2a) = c(
  "FDR",
  "FPR",
  "TPR: X1",
  "TPR: X2",
  "TPR: X3",
  "TPR: X4",
  "TPR: X5",
  "CI Coverage: X1",
  "CI Coverage: X2",
  "CI Coverage: X3",
  "CI Coverage: X4",
  "CI Coverage: X5",
  "Type II Error: X1",
  "Type II Error: X2",
  "Type II Error: X3",
  "Type II Error: X4",
  "Type II Error: X5",
  "Type I Error: X6-X20",
  "Bias: X1",
  "Bias: X2",
  "Bias: X3",
  "Bias: X4",
  "Bias: X5",
  "Bias: X6-X20"
)

# fdr figure case 2a
fdr_fig3_dat = data.frame(fdr)
fdr_fig3_dat$method = c("LRT", "AIC", "BIC", "cv.LASSO", "fixed.LASSO", "cv.ENET", "fixed.ENET")

fdr_fig3 = ggplot(fdr_fig3_dat, aes(y = fdr, x = method, fill=method)) +
  geom_bar(stat="identity") +
  theme(legend.position = "none") +
  ylim(0, 0.5)

# fpr figure case 2a
fpr_fig3_dat = data.frame(fpr)
fpr_fig3_dat$method = c("LRT", "AIC", "BIC", "cv.LASSO", "fixed.LASSO", "cv.ENET", "fixed.ENET")

fpr_fig3 = ggplot(fpr_fig3_dat, aes(y = fpr, x = method, fill=method)) +
  geom_bar(stat="identity") +
  theme(legend.position = "none") +
  ylim(0, 0.5)

# toe figure case 2a
toe_fig3_dat = data.frame(toe)
toe_fig3_dat$method = c("LRT", "AIC", "BIC", "cv.LASSO", "fixed.LASSO", "cv.ENET", "fixed.ENET")

toe_fig3 = ggplot(toe_fig3_dat, aes(y = toe, x = method, fill=method)) +
  geom_bar(stat="identity") +
  theme(legend.position = "none") +
  ylim(0, 1)

# bias0 figure case 2a
bias0_fig3_dat = data.frame(bias0)
bias0_fig3_dat$method = c("LRT", "AIC", "BIC", "cv.LASSO", "fixed.LASSO", "cv.ENET", "fixed.ENET")

bias0_fig3 = ggplot(bias0_fig3_dat, aes(y = bias0, x = method, fill=method)) +
  geom_bar(stat="identity") +
  theme(legend.position = "none") +
  ylim(0, 0.1)

# tpr figure case 2a
tpr_fig3_dat = data.frame(tpr)
tpr_fig3_dat$id = 1:nrow(tpr_fig3_dat)
tpr_fig3_dat$method = c("LRT", "AIC", "BIC", "cv-LASSO", "fixed-LASSO", "cv-ENET", "fixed-ENET")
tpr_fig3_dat = tpr_fig3_dat %>%
  pivot_longer(
    cols = c(X1, X2, X3, X4, X5), 
    names_to = "method",
    values_to = "tpr",
    names_repair = "unique"
  )
colnames(tpr_fig3_dat) = c("id", "method", "variable", "tpr")

tpr_fig3 = ggplot(tpr_fig3_dat, aes(y = tpr, x = method, color=variable)) +
  geom_bar(aes(fill = variable), stat="identity", position = "dodge") +
  theme(legend.position = "none") +
  ylim(0, 1)

# coverage figure case 2a
coverage_fig3_dat = data.frame(coverage)
coverage_fig3_dat$id = 1:nrow(coverage_fig3_dat)
coverage_fig3_dat$method = c("LRT", "AIC", "BIC", "cv-LASSO", "fixed-LASSO", "cv-ENET", "fixed-ENET")
coverage_fig3_dat = coverage_fig3_dat %>%
  pivot_longer(
    cols = c(X1, X2, X3, X4, X5), 
    names_to = "method",
    values_to = "tpr",
    names_repair = "unique"
  )
colnames(coverage_fig3_dat) = c("id", "method", "variable", "coverage")

coverage_fig3 = ggplot(coverage_fig3_dat, aes(y = coverage, x = method, color=variable)) +
  geom_bar(aes(fill = variable), stat="identity", position = "dodge") +
  theme(legend.position = "none") +
  ylim(0, 1)

# type II error figure case 2a
tte_fig3_dat = data.frame(tte)
tte_fig3_dat$id = 1:nrow(tte_fig3_dat)
tte_fig3_dat$method = c("LRT", "AIC", "BIC", "cv-LASSO", "fixed-LASSO", "cv-ENET", "fixed-ENET")
tte_fig3_dat = tte_fig3_dat %>%
  pivot_longer(
    cols = c(X1, X2, X3, X4, X5), 
    names_to = "method",
    values_to = "tpr",
    names_repair = "unique"
  )
colnames(tte_fig3_dat) = c("id", "method", "variable", "tte")

tte_fig3 = ggplot(tte_fig3_dat, aes(y = tte, x = method, group = variable)) +
  geom_bar(aes(fill = variable), stat="identity", position = "dodge") +
  theme(legend.position = "none") +
  ylim(0, 0.6)

# bias figure case 2a
bias_fig3_dat = data.frame(bias)
bias_fig3_dat$id = 1:nrow(bias_fig3_dat)
bias_fig3_dat$method = c("LRT", "AIC", "BIC", "cv-LASSO", "fixed-LASSO", "cv-ENET", "fixed-ENET")
bias_fig3_dat = bias_fig3_dat %>%
  pivot_longer(
    cols = c(X1, X2, X3, X4, X5), 
    names_to = "method",
    values_to = "tpr",
    names_repair = "unique"
  )
colnames(bias_fig3_dat) = c("id", "method", "variable", "bias")

bias_fig3 = ggplot(bias_fig3_dat, aes(y = bias, x = method, group = variable)) +
  geom_bar(aes(fill = variable), stat="identity", position = "dodge") +
  theme(legend.position = "none") +
  ylim(0, 0.2)

########
#  2b  #
########
########################
#     read in data     #
########################

p_betas = read.delim("../data_raw/p_betas_2b.txt")
p_ci = read.delim("../data_raw/p_ci_2b.txt")
p_p = read.delim("../data_raw/p_p_2b.txt")

aic_betas = read.delim("../data_raw/aic_betas_2b.txt")
aic_ci = read.delim("../data_raw/aic_ci_2b.txt")
aic_p = read.delim("../data_raw/aic_p_2b.txt")

bic_betas = read.delim("../data_raw/bic_betas_2b.txt")
bic_ci = read.delim("../data_raw/bic_ci_2b.txt")
bic_p = read.delim("../data_raw/bic_p_2b.txt")

cv_lasso_betas = read.delim("../data_raw/cv_lasso_betas_2b.txt")
cv_lasso_ci = read.delim("../data_raw/cv_lasso_ci_2b.txt")
cv_lasso_p = read.delim("../data_raw/cv_lasso_p_2b.txt")

fixed_lasso_betas = read.delim("../data_raw/fixed_lasso_betas_2b.txt")
fixed_lasso_ci = read.delim("../data_raw/fixed_lasso_ci_2b.txt")
fixed_lasso_p = read.delim("../data_raw/fixed_lasso_p_2b.txt")

cv_enet_betas = read.delim("../data_raw/cv_enet_betas_2b.txt")
cv_enet_ci = read.delim("../data_raw/cv_enet_ci_2b.txt")
cv_enet_p = read.delim("../data_raw/cv_enet_p_2b.txt")

fixed_enet_betas = read.delim("../data_raw/fixed_enet_betas_2b.txt")
fixed_enet_ci = read.delim("../data_raw/fixed_enet_ci_2b.txt")
fixed_enet_p = read.delim("../data_raw/fixed_enet_p_2b.txt")

########################
# calculate parameters #
########################
beta = vector("double", 20)
beta[1:5] = c(0.5, 1., 1.5, 2., 2.5) / 3

# true positive rate
tpr = matrix(0, nrow = 7, ncol = 5)
for (i in 1:dim(p_betas[, 1:5])[2]) {
  tpr[1, i] = truePositiveRate(p_betas[, i])
}

for (i in 1:dim(aic_betas[, 1:5])[2]) {
  tpr[2, i] = truePositiveRate(aic_betas[, i])
}

for (i in 1:dim(bic_betas[, 1:5])[2]) {
  tpr[3, i] = truePositiveRate(bic_betas[, i])
}

for (i in 1:dim(cv_lasso_betas[, 1:5])[2]) {
  tpr[4, i] = truePositiveRate(cv_lasso_betas[, i])
}

for (i in 1:dim(fixed_lasso_betas[, 1:5])[2]) {
  tpr[5, i] = truePositiveRate(fixed_lasso_betas[, i])
}

for (i in 1:dim(cv_enet_betas[, 1:5])[2]) {
  tpr[6, i] = truePositiveRate(cv_enet_betas[, i])
}

for (i in 1:dim(fixed_enet_betas[, 1:5])[2]) {
  tpr[7, i] = truePositiveRate(fixed_enet_betas[, i])
}

# false positive rate
fpr = vector("double", 7)
fpr[1] = falsePositiveRate(p_betas[, 5:15])
fpr[2] = falsePositiveRate(aic_betas[, 5:15])
fpr[3] = falsePositiveRate(bic_betas[, 5:15])
fpr[4] = falsePositiveRate(cv_lasso_betas[, 5:15])
fpr[5] = falsePositiveRate(fixed_lasso_betas[, 5:15])
fpr[6] = falsePositiveRate(cv_enet_betas[, 5:15])
fpr[7] = falsePositiveRate(fixed_enet_betas[, 5:15])

# false discovery rate
fdr = vector("double", 7)
fdr[1] = falseDiscoveryRate(p_betas, p_betas[, 5:15])
fdr[2] = falseDiscoveryRate(aic_betas, aic_betas[, 5:15])
fdr[3] = falseDiscoveryRate(bic_betas, bic_betas[, 5:15])
fdr[4] = falseDiscoveryRate(cv_lasso_betas, cv_lasso_betas[, 5:15])
fdr[5] = falseDiscoveryRate(fixed_lasso_betas, fixed_lasso_betas[, 5:15])
fdr[6] = falseDiscoveryRate(cv_enet_betas, cv_enet_betas[, 5:15])
fdr[7] = falseDiscoveryRate(fixed_enet_betas, fixed_enet_betas[, 5:15])

# type I error
toe = vector("double", 7)
toe[1] = typeOneErrorRate(p_p[, 5:15])
toe[2] = typeOneErrorRate(aic_p[, 5:15])
toe[3] = typeOneErrorRate(bic_p[, 5:15])
toe[4] = typeOneErrorRate(cv_lasso_p[, 5:15])
toe[5] = typeOneErrorRate(fixed_lasso_p[, 5:15])
toe[6] = typeOneErrorRate(cv_enet_p[, 5:15])
toe[7] = typeOneErrorRate(fixed_enet_p[, 5:15])

# type II error rate
tte = matrix(0, nrow = 7, ncol = 5)
for (i in 1:dim(p_p[, 1:5])[2]) {
  tte[1, i] = typeTwoErrorRate(p_p[, i])
}

for (i in 1:dim(aic_p[, 1:5])[2]) {
  tte[2, i] = typeTwoErrorRate(aic_p[, i])
}

for (i in 1:dim(bic_p[, 1:5])[2]) {
  tte[3, i] = typeTwoErrorRate(bic_p[, i])
}

for (i in 1:dim(cv_lasso_p[, 1:5])[2]) {
  tte[4, i] = typeTwoErrorRate(cv_lasso_p[, i])
}

for (i in 1:dim(fixed_lasso_p[, 1:5])[2]) {
  tte[5, i] = typeTwoErrorRate(fixed_lasso_p[, i])
}

for (i in 1:dim(cv_enet_p[, 1:5])[2]) {
  tte[6, i] = typeTwoErrorRate(cv_enet_p[, i])
}

for (i in 1:dim(fixed_enet_p[, 1:5])[2]) {
  tte[7, i] = typeTwoErrorRate(fixed_enet_p[, i])
}

# parameters 1-5 bias
bias = matrix(0, nrow = 7, ncol = 5)
for (i in 1:dim(p_betas[, 1:5])[2]) {
  bias[1, i] = parameterBiasNonZero(p_betas[, i], beta[i])
}

for (i in 1:dim(aic_betas[, 1:5])[2]) {
  bias[2, i] = parameterBiasNonZero(aic_betas[, i], beta[i])
}

for (i in 1:dim(bic_betas[, 1:5])[2]) {
  bias[3, i] = parameterBiasNonZero(bic_betas[, i], beta[i])
}

for (i in 1:dim(cv_lasso_betas[, 1:5])[2]) {
  bias[4, i] = parameterBiasNonZero(cv_lasso_betas[, i], beta[i])
}

for (i in 1:dim(fixed_lasso_betas[, 1:5])[2]) {
  bias[5, i] = parameterBiasNonZero(fixed_lasso_betas[, i], beta[i])
}

for (i in 1:dim(cv_enet_betas[, 1:5])[2]) {
  bias[6, i] = parameterBiasNonZero(cv_enet_betas[, i], beta[i])
}

for (i in 1:dim(fixed_enet_betas[, 1:5])[2]) {
  bias[7, i] = parameterBiasNonZero(fixed_enet_betas[, i], beta[i])
}

# parameters 5-20 average bias
bias0 = vector("double", 7)
bias0[1] = parameterBiasZero(p_betas[, 5:15])
bias0[2] = parameterBiasZero(aic_betas[, 5:15])
bias0[3] = parameterBiasZero(bic_betas[, 5:15])
bias0[4] = parameterBiasZero(cv_lasso_betas[, 5:15])
bias0[5] = parameterBiasZero(fixed_lasso_betas[, 5:15])
bias0[6] = parameterBiasZero(cv_enet_betas[, 5:15])
bias0[7] = parameterBiasZero(fixed_enet_betas[, 5:15])

# parameters 1-5 CI coverage
coverage = matrix(0, nrow = 7, ncol = 5)
for (i in 1:dim(p_ci[, 1:5])[2]) {
  coverage[1, i] = coverageCI(p_ci[seq(1, ncol(p_ci) - 1, 2), i], p_ci[seq(2, ncol(p_ci), 2), i], p_betas[, i], beta[i])
}

for (i in 1:dim(aic_ci[, 1:5])[2]) {
  coverage[2, i] = coverageCI(aic_ci[seq(1, ncol(aic_ci) - 1, 2), i], aic_ci[seq(2, ncol(aic_ci), 2), i], aic_betas[, i], beta[i])
}

for (i in 1:dim(bic_ci[, 1:5])[2]) {
  coverage[3, i] = coverageCI(bic_ci[seq(1, ncol(bic_ci) - 1, 2), i], bic_ci[seq(2, ncol(bic_ci), 2), i], bic_betas[, i], beta[i])
}

for (i in 1:dim(cv_lasso_ci[, 1:5])[2]) {
  coverage[4, i] = coverageCI(cv_lasso_ci[seq(1, ncol(cv_lasso_ci) - 1, 2), i], cv_lasso_ci[seq(2, ncol(cv_lasso_ci), 2), i], cv_lasso_betas[, i], beta[i])
}

for (i in 1:dim(fixed_lasso_ci[, 1:5])[2]) {
  coverage[5, i] = coverageCI(fixed_lasso_ci[seq(1, ncol(fixed_lasso_ci) - 1, 2), i], fixed_lasso_ci[seq(2, ncol(fixed_lasso_ci), 2), i], fixed_lasso_betas[, i], beta[i])
}

for (i in 1:dim(cv_enet_ci[, 1:5])[2]) {
  coverage[6, i] = coverageCI(cv_enet_ci[seq(1, ncol(cv_enet_ci) - 1, 2), i], cv_enet_ci[seq(2, ncol(cv_enet_ci), 2), i], cv_enet_betas[, 1:5], beta[i])
}

for (i in 1:dim(fixed_enet_ci[, 1:5])[2]) {
  coverage[7, i] = coverageCI(fixed_enet_ci[seq(1, ncol(fixed_enet_ci) - 1, 2), i], fixed_enet_ci[seq(2, ncol(fixed_enet_ci), 2), i], fixed_enet_betas[, i], beta[i])
}

################
# build table #
################

table_2b = data.frame(
  "LRT" = c(
    fdr[1], fpr[1], tpr[1, 1], tpr[1, 2], tpr[1, 3],
    tpr[1, 4], tpr[1, 5], coverage[1, 1], coverage[1, 2],
    coverage[1, 3], coverage[1, 4], coverage[1, 5],
    tte[1, 1], tte[1, 2], tte[1, 3], tte[1, 4], tte[1, 5],
    toe[1], bias[1, 1], bias[1, 2], bias[1, 3], bias[1, 4],
    bias[1, 5], bias0[1]
  ),
  "AIC" = c(
    fdr[2], fpr[2], tpr[2, 1], tpr[2, 2], tpr[2, 3],
    tpr[2, 4], tpr[2, 5], coverage[2, 1], coverage[2, 2],
    coverage[2, 3], coverage[2, 4], coverage[2, 5],
    tte[2, 1], tte[2, 2], tte[2, 3], tte[2, 4], tte[2, 5],
    toe[1], bias[2, 1], bias[2, 2], bias[2, 3], bias[2, 4],
    bias[2, 5], bias0[2]
  ),
  "BIC" = c(
    fdr[3], fpr[3], tpr[3, 1], tpr[3, 2], tpr[3, 3],
    tpr[3, 4], tpr[3, 5], coverage[3, 1], coverage[3, 2],
    coverage[3, 3], coverage[3, 4], coverage[3, 5],
    tte[3, 1], tte[3, 2], tte[3, 3], tte[3, 4], tte[3, 5],
    toe[1], bias[3, 1], bias[3, 2], bias[3, 3], bias[3, 4],
    bias[3, 5], bias0[3]
  ),
  "cv-LASSO" = c(
    fdr[4], fpr[4], tpr[4, 1], tpr[4, 2], tpr[4, 3],
    tpr[4, 4], tpr[4, 5], coverage[4, 1], coverage[4, 2],
    coverage[4, 3], coverage[4, 4], coverage[4, 5],
    tte[4, 1], tte[4, 2], tte[4, 3], tte[4, 4], tte[4, 5],
    toe[1], bias[4, 1], bias[4, 2], bias[4, 3], bias[4, 4],
    bias[4, 5], bias0[4]
  ),
  "fixed-LASSO" = c(
    fdr[5], fpr[5], tpr[5, 1], tpr[5, 2], tpr[5, 3],
    tpr[5, 4], tpr[5, 5], coverage[5, 1], coverage[5, 2],
    coverage[5, 3], coverage[5, 4], coverage[5, 5],
    tte[5, 1], tte[5, 2], tte[5, 3], tte[5, 4], tte[5, 5],
    toe[1], bias[5, 1], bias[5, 2], bias[5, 3], bias[5, 4],
    bias[5, 5], bias0[5]
  ),
  "cv-ENET" = c(
    fdr[6], fpr[6], tpr[6, 1], tpr[6, 2], tpr[6, 3],
    tpr[6, 4], tpr[6, 5], coverage[6, 1], coverage[6, 2],
    coverage[6, 3], coverage[6, 4], coverage[6, 5],
    tte[6, 1], tte[6, 2], tte[6, 3], tte[6, 4], tte[6, 5],
    toe[1], bias[6, 1], bias[6, 2], bias[6, 3], bias[6, 4],
    bias[6, 5], bias0[6]
  ),
  "fixed-ENET" = c(
    fdr[7], fpr[7], tpr[7, 1], tpr[7, 2], tpr[7, 3],
    tpr[7, 4], tpr[7, 5], coverage[7, 1], coverage[7, 2],
    coverage[7, 3], coverage[7, 4], coverage[7, 5],
    tte[7, 1], tte[7, 2], tte[7, 3], tte[7, 4], tte[7, 5],
    toe[1], bias[7, 1], bias[7, 2], bias[7, 3], bias[7, 4],
    bias[7, 5], bias0[7]
  )
)

rownames(table_2b) = c(
  "FDR",
  "FPR",
  "TPR: X1",
  "TPR: X2",
  "TPR: X3",
  "TPR: X4",
  "TPR: X5",
  "CI Coverage: X1",
  "CI Coverage: X2",
  "CI Coverage: X3",
  "CI Coverage: X4",
  "CI Coverage: X5",
  "Type II Error: X1",
  "Type II Error: X2",
  "Type II Error: X3",
  "Type II Error: X4",
  "Type II Error: X5",
  "Type I Error: X6-X20",
  "Bias: X1",
  "Bias: X2",
  "Bias: X3",
  "Bias: X4",
  "Bias: X5",
  "Bias: X6-X20"
)

# fdr figure case 2b
fdr_fig4_dat = data.frame(fdr)
fdr_fig4_dat$method = c("LRT", "AIC", "BIC", "cv.LASSO", "fixed.LASSO", "cv.ENET", "fixed.ENET")

fdr_fig4 = ggplot(fdr_fig4_dat, aes(y = fdr, x = method, fill=method)) +
  geom_bar(stat="identity") +
  theme(legend.position = "none") +
  ylim(0, 0.5)

# fpr figure case 2b
fpr_fig4_dat = data.frame(fpr)
fpr_fig4_dat$method = c("LRT", "AIC", "BIC", "cv.LASSO", "fixed.LASSO", "cv.ENET", "fixed.ENET")

fpr_fig4 = ggplot(fpr_fig4_dat, aes(y = fpr, x = method, fill=method)) +
  geom_bar(stat="identity") +
  theme(legend.position = "none") +
  ylim(0, 0.5)

# toe figure case 2b
toe_fig4_dat = data.frame(toe)
toe_fig4_dat$method = c("LRT", "AIC", "BIC", "cv.LASSO", "fixed.LASSO", "cv.ENET", "fixed.ENET")

toe_fig4 = ggplot(toe_fig4_dat, aes(y = toe, x = method, fill=method)) +
  geom_bar(stat="identity") +
  theme(legend.position = "none") +
  ylim(0, 1)

# bias0 figure case 2b
bias0_fig4_dat = data.frame(bias0)
bias0_fig4_dat$method = c("LRT", "AIC", "BIC", "cv.LASSO", "fixed.LASSO", "cv.ENET", "fixed.ENET")

bias0_fig4 = ggplot(bias0_fig4_dat, aes(y = bias0, x = method, fill=method)) +
  geom_bar(stat="identity") +
  theme(legend.position = "none") +
  ylim(0, 0.1)

# tpr figure case 2b
tpr_fig4_dat = data.frame(tpr)
tpr_fig4_dat$id = 1:nrow(tpr_fig4_dat)
tpr_fig4_dat$method = c("LRT", "AIC", "BIC", "cv-LASSO", "fixed-LASSO", "cv-ENET", "fixed-ENET")
tpr_fig4_dat = tpr_fig4_dat %>%
  pivot_longer(
    cols = c(X1, X2, X3, X4, X5), 
    names_to = "method",
    values_to = "tpr",
    names_repair = "unique"
  )
colnames(tpr_fig4_dat) = c("id", "method", "variable", "tpr")

tpr_fig4 = ggplot(tpr_fig4_dat, aes(y = tpr, x = method, color=variable)) +
  geom_bar(aes(fill = variable), stat="identity", position = "dodge") +
  theme(legend.position = "none") +
  ylim(0, 1)

# coverage figure case 2b
coverage_fig4_dat = data.frame(coverage)
coverage_fig4_dat$id = 1:nrow(coverage_fig4_dat)
coverage_fig4_dat$method = c("LRT", "AIC", "BIC", "cv-LASSO", "fixed-LASSO", "cv-ENET", "fixed-ENET")
coverage_fig4_dat = coverage_fig4_dat %>%
  pivot_longer(
    cols = c(X1, X2, X3, X4, X5), 
    names_to = "method",
    values_to = "tpr",
    names_repair = "unique"
  )
colnames(coverage_fig4_dat) = c("id", "method", "variable", "coverage")

coverage_fig4 = ggplot(coverage_fig4_dat, aes(y = coverage, x = method, color=variable)) +
  geom_bar(aes(fill = variable), stat="identity", position = "dodge") +
  theme(legend.position = "none") +
  ylim(0, 1)

# type II error figure case 2b
tte_fig4_dat = data.frame(tte)
tte_fig4_dat$id = 1:nrow(tte_fig4_dat)
tte_fig4_dat$method = c("LRT", "AIC", "BIC", "cv-LASSO", "fixed-LASSO", "cv-ENET", "fixed-ENET")
tte_fig4_dat = tte_fig4_dat %>%
  pivot_longer(
    cols = c(X1, X2, X3, X4, X5), 
    names_to = "method",
    values_to = "tpr",
    names_repair = "unique"
  )
colnames(tte_fig4_dat) = c("id", "method", "variable", "tte")

tte_fig4 = ggplot(tte_fig4_dat, aes(y = tte, x = method, group = variable)) +
  geom_bar(aes(fill = variable), stat="identity", position = "dodge") +
  theme(legend.position = "none") +
  ylim(0, 0.6)

# bias figure case 2b
bias_fig4_dat = data.frame(bias)
bias_fig4_dat$id = 1:nrow(bias_fig4_dat)
bias_fig4_dat$method = c("LRT", "AIC", "BIC", "cv-LASSO", "fixed-LASSO", "cv-ENET", "fixed-ENET")
bias_fig4_dat = bias_fig4_dat %>%
  pivot_longer(
    cols = c(X1, X2, X3, X4, X5), 
    names_to = "method",
    values_to = "tpr",
    names_repair = "unique"
  )
colnames(bias_fig4_dat) = c("id", "method", "variable", "bias")

bias_fig4 = ggplot(bias_fig4_dat, aes(y = bias, x = method, group = variable)) +
  geom_bar(aes(fill = variable), stat="identity", position = "dodge") +
  theme(legend.position = "none") +
  ylim(0, 0.2)

# fdr plot
ggarrange(
  fdr_fig1,
  fdr_fig2,
  fdr_fig3,
  fdr_fig4,
  ncol = 2,
  nrow = 2,
  labels = c("1A", "1B", "2A", "2B"),
  label.x = 0.85,
  label.y = 0.975
)

# fpr plot
ggarrange(
  fpr_fig1,
  fpr_fig2,
  fpr_fig3,
  fpr_fig4,
  ncol = 2,
  nrow = 2,
  labels = c("1A", "1B", "2A", "2B"),
  label.x = 0.85,
  label.y = 0.975
)

# toe plot
ggarrange(
  toe_fig1,
  toe_fig2,
  toe_fig3,
  toe_fig4,
  ncol = 2,
  nrow = 2,
  labels = c("1A", "1B", "2A", "2B"),
  label.x = 0.85,
  label.y = 0.975
)

# bias0 plot
ggarrange(
  bias0_fig1,
  bias0_fig2,
  bias0_fig3,
  bias0_fig4,
  ncol = 2,
  nrow = 2,
  labels = c("1A", "1B", "2A", "2B"),
  label.x = 0.85,
  label.y = 0.975
)

# tpr plot
ggarrange(
  tpr_fig1,
  tpr_fig2,
  tpr_fig3,
  tpr_fig4,
  ncol = 2,
  nrow = 2,
  labels = c("1A", "1B", "2A", "2B"),
  common.legend = TRUE,
  legend = "bottom",
  label.x = 0.85,
  label.y = 0.975
)

# coverage plot
ggarrange(
  coverage_fig1,
  coverage_fig2,
  coverage_fig3,
  coverage_fig4,
  ncol = 2,
  nrow = 2,
  labels = c("1A", "1B", "2A", "2B"),
  common.legend = TRUE,
  legend = "bottom",
  label.x = 0.85,
  label.y = 0.975
)

# tte plot
ggarrange(
  tte_fig1,
  tte_fig2,
  tte_fig3,
  tte_fig4,
  ncol = 2,
  nrow = 2,
  labels = c("1A", "1B", "2A", "2B"),
  common.legend = TRUE,
  legend = "bottom",
  label.x = 0.85,
  label.y = 0.975
)

# bias plot
ggarrange(
  bias_fig1,
  bias_fig2,
  bias_fig3,
  bias_fig4,
  ncol = 2,
  nrow = 2,
  labels = c("1A", "1B", "2A", "2B"),
  common.legend = TRUE,
  legend = "bottom",
  label.x = 0.85,
  label.y = 0.975
)

# table 1a
table_1a %>% gt(rownames_to_stub = TRUE)

# table 1b
table_1b %>% gt(rownames_to_stub = TRUE)

# table 2a
table_2a %>% gt(rownames_to_stub = TRUE)

# table 2b
table_2b %>% gt(rownames_to_stub = TRUE)

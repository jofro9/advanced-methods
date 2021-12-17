#############################################################
## Program Name: Descriptives.R						                 ##
## Purpose: To read in the clean dataset and create graphs ##
##     		 and descriptive statistics				               ##
## Created by: Nichole Carlson	(R code by: Kevin Josey)	 ##
#############################################################

# Dependencies
library(readr)
library(dplyr)
library(psych)

###You will have to set your working directory####
setwd("/Users/nichole/Repositories/BIOS6624/Bios6624ClassExamples/PSAExample/")

### read in the first cleaned data set (psaclean.csv)

## read in the clean data from our data directory
psaclean <- read_csv("DataProcessed/psaclean.csv")
  
## create scatter plots with PSA as the outcome
png("Output/Descriptives-1.png",
    width = 1500, 
    height = 2000,
    res = 100, 
    units = "px")
par(mfrow = c(3,2))

plot(psa ~ cavol, data = psaclean, main = "psa vs. cavol")
y <- predict(lm(psa ~ cavol, data = psaclean))
x <- psaclean$cavol
lines(y = y, x = x, col = "navy")

plot(psa ~ wt, data = psaclean, main = "psa vs. wt")
y <- predict(lm(psa ~ wt, data = psaclean))
x <- psaclean$wt
lines(y = y, x = x, col = "navy")

plot(psa ~ age, data = psaclean, main = "psa vs. age")
y <- predict(lm(psa ~ age, data = psaclean))
x <- psaclean$age
lines(y = y, x = x, col = "navy")

plot(psa ~ bph, data = psaclean, main = "psa vs. bph")
y <- predict(lm(psa ~ bph, data = psaclean))
x <- psaclean$bph
lines(y = y, x = x, col = "navy")

plot(psa ~ cappen, data = psaclean, main = "psa vs. cappen")
y <- predict(lm(psa ~ cappen, data = psaclean))
x <- psaclean$cappen
lines(y = y, x = x, col = "navy")

plot.new()
legend("center", legend = c("Fitted Simple Regression Line"), col = c("navy"),
       lty = 1, lwd = 3, bty = "n", cex = 3)

dev.off()

## continuous variables
cont_desc <- psaclean %>%
  select(psa, lpsa, cavol, wt, age, bph, cappen) %>%
  sapply(describe) # depends on library(psych)

## categorical variables
cat_desc <- psaclean %>%
  select(svi, gleason) %>%
  lapply(table)

cat_desc_prop <- psaclean %>%
  select(svi, gleason) %>%
  lapply(table) %>%
  lapply(prop.table)

###Write to screen to see things #####
##Write to a file: Not pretty but it does the trick.  There are other Table 1 functions to use for real reports
sink("Output/DescriptiveTables1.txt")
print(cont_desc)
print(cat_desc)
print(cat_desc_prop)
sink()

### read in the final clean data from our data directory after addressing findings in initial descriptive statistics
psaclean2 <- read_csv("DataProcessed/psaclean2.csv")
  
## create scatter plots with PSA as the outcome
# png("C:/Repositories/Bios6623ClassExamples/PSA-Project/G1Analysis/Plots/Descriptives-2.png")
png("Output/Descriptives-2.png",
    width = 1500, 
    height = 2000,
    res = 100, 
    units = "px")
par(mfrow = c(3,2))

plot(psa ~ cavol, data = psaclean2, main = "psa vs. cavol")
y <- predict(lm(psa ~ cavol, data = psaclean2))
x <- psaclean2$cavol
lines(y = y, x = x, col = "navy")

plot(psa ~ wt, data = psaclean2, main = "psa vs. wt")
y <- predict(lm(psa ~ wt, data = psaclean2))
x <- psaclean2$wt
lines(y = y, x = x, col = "navy")

plot(psa ~ age, data = psaclean2, main = "psa vs. age")
y <- predict(lm(psa ~ age, data = psaclean2))
x <- psaclean2$age
lines(y = y, x = x, col = "navy")

plot(psa ~ bph, data = psaclean2, main = "psa vs. bph")
y <- predict(lm(psa ~ bph, data = psaclean2))
x <- psaclean2$bph
lines(y = y, x = x, col = "navy")

plot(psa ~ cappen, data = psaclean2, main = "psa vs. cappen")
y <- predict(lm(psa ~ cappen, data = psaclean2))
x <- psaclean2$cappen
lines(y = y, x = x, col = "navy")

plot.new()
legend("center", legend = c("Fitted Simple Regression Line"), col = c("navy"),
       lty = 1, lwd = 3, bty = "n", cex = 3)

dev.off()
# might be best to write a function for plotting these curves

## create scatter plots with log(PSA) as the outcome (log transforming the outcom => approximately normal)
# png("C:/Repositories/Bios6623ClassExamples/PSA-Project/G1Analysis/Plots/Descriptives-3.png")
png("Output/Descriptives-3.png",
    width = 1500, 
    height = 2000,
    res = 100, 
    units = "px")

par(mfrow = c(3,2))

plot(lpsa ~ cavol, data = psaclean2, main = "lpsa vs. cavol")
y <- predict(lm(lpsa ~ cavol, data = psaclean2))
x <- psaclean2$cavol
lines(y = y, x = x, col = "navy")

plot(lpsa ~ wt, data = psaclean2, main = "lpsa vs. wt")
y <- predict(lm(lpsa ~ wt, data = psaclean2))
x <- psaclean2$wt
lines(y = y, x = x, col = "navy")

plot(lpsa ~ age, data = psaclean2, main = "lpsa vs. age")
y <- predict(lm(lpsa ~ age, data = psaclean2))
x <- psaclean2$age
lines(y = y, x = x, col = "navy")

plot(lpsa ~ bph, data = psaclean2, main = "lpsa vs. bph")
y <- predict(lm(lpsa ~ bph, data = psaclean2))
x <- psaclean2$bph
lines(y = y, x = x, col = "navy")

plot(lpsa ~ cappen, data = psaclean2, main = "lpsa vs. cappen")
y <- predict(lm(lpsa ~ cappen, data = psaclean2))
x <- psaclean2$cappen
lines(y = y, x = x, col = "navy")

plot.new()
legend("center", legend = c("Fitted Simple Regression Line"), col = c("navy"),
       lty = 1, lwd = 3, bty = "n", cex = 3)

dev.off()

## Make descriptives for continuous variables
cont_desc <- psaclean2 %>%
  select(psa, lpsa, cavol, wt, age, bph, cappen) %>%
  sapply(describe) # depends on library(psych)

## Make descriptives for categorical variables
cat_desc <- psaclean2 %>%
  select(svi, gleason) %>%
  sapply(table)

# Make descriptives for categorical variable proportions
cat_desc_prop <- psaclean2 %>%
  select(svi, gleason) %>%
  lapply(table) %>%
  lapply(prop.table)

# Make pearson correlation coefficients
cor_mat <- psaclean2 %>%
  select(psa, lpsa, cavol, wt, age, bph, cappen) %>%
  corr.test(.) # depends on library(psych)

#Write to a file for future use: Not pretty but useful.
sink("Output/DescriptiveTablesClean.txt")
print(cont_desc)
print(cat_desc)
print(cat_desc_prop)
print(cor_mat)
sink()

## create scatter plots with log(PSA) as the outcome and an interaction for svi

svi <- psaclean2$svi # indicator for generating interaction plots

# png("C:/Repositories/Bios6623ClassExamples/PSA-Project/G1Analysis/Plots/Descriptives-4.png")
png("Output/Descriptives-4.png",
    width = 1500, 
    height = 2000,
    res = 100, 
    units = "px")
    
par(mfrow = c(3,2))
    
plot(lpsa ~ cavol, data = psaclean2, main = "lpsa vs. cavol",
     col = ifelse(svi == 0, "navy", "darkorange"))
y <- predict(lm(lpsa ~ cavol*svi, data = psaclean2))
x <- psaclean2$cavol
lines(y = y[svi == 0], x = x[svi == 0], col = "navy")
lines(y = y[svi == 1], x = x[svi == 1], col = "darkorange")

plot(lpsa ~ wt, data = psaclean2, main = "lpsa vs. wt",
     col = ifelse(svi == 0, "navy", "darkorange"))
y <- predict(lm(lpsa ~ wt*svi, data = psaclean2))
x <- psaclean2$wt
lines(y = y[svi == 0], x = x[svi == 0], col = "navy")
lines(y = y[svi == 1], x = x[svi == 1], col = "darkorange")

plot(lpsa ~ age, data = psaclean2, main = "lpsa vs. age",
     col = ifelse(svi == 0, "navy", "darkorange"))
y <- predict(lm(lpsa ~ age*svi, data = psaclean2))
x <- psaclean2$age
lines(y = y[svi == 0], x = x[svi == 0], col = "navy")
lines(y = y[svi == 1], x = x[svi == 1], col = "darkorange")

plot(lpsa ~ bph, data = psaclean2, main = "lpsa vs. bph",
     col = ifelse(svi == 0, "navy", "darkorange"))
y <- predict(lm(lpsa ~ bph*svi, data = psaclean2))
x <- psaclean2$bph
lines(y = y[svi == 0], x = x[svi == 0], col = "navy")
lines(y = y[svi == 1], x = x[svi == 1], col = "darkorange")

plot(lpsa ~ cappen, data = psaclean2, main = "lpsa vs. cappen",
     col = ifelse(svi == 0, "navy", "darkorange"))
y <- predict(lm(lpsa ~ cappen*svi, data = psaclean2))
x <- psaclean2$cappen
lines(y = y[svi == 0], x = x[svi == 0], col = "navy")
lines(y = y[svi == 1], x = x[svi == 1], col = "darkorange")

plot.new()
legend("center", legend = c("svi = 0", "svi = 1"), col = c("navy", "darkorange"),
       lty = 1, lwd = 3, bty = "n", cex = 3)
    
dev.off()

###Maybe make a table 1 by SVI for continuous variables.  We would need to put gleason as a facctor to add it.
psaclean2$gleason <-as.factor(psaclean2$gleason)
Table1bySVI <- CreateTableOne(vars = c("psa", "lpsa", "cavol", "wt", "age", "bph", "cappen","gleason"),strata=c("svi"),data=psaclean2)
sink("Output/DescriptiveTablesBySVI.txt")
print(Table1bySVI)
sink()
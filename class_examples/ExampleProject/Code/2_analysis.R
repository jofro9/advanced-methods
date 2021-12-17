################################################################################
## Project0: analysis.R
## Author: Student A
## Description: This file conducts the univariate and multivariable models
##              to assess whether groups differ on the two outcomes.
##              My analysis strategy is to use the year 1 data as the outcome,
##              because it is more powerful than the difference. I will adjust for 
##              baseline of the outcome for precision. 
################################################################################


##---------------------------------##
## Setup workspace
##---------------------------------##

## Load libraries
# install/load devtools for installing github packages
#   *** if using windows, you need to install Rtools from CRAN ***
#   *** https://cran.r-project.org/bin/windows/Rtools/ ***
if (!require("devtools")) { 
  install.packages("devtools") 
  library(devtools) 
}
install_github("dewittpe/qwraps")

## Note these will be custom to your approach in R
library(Hmisc)
library(xtable)
library(pwr)
library(ggplot2)
library(corrgram)
library(car)
library(qwraps)
library(multcomp)
library(lsmeans)

## Init functions


## Load data
load("/Users/Nichole/Repositories/Bios6623/DATA/analysis.Rdata")


##---------------------------------##
## Univariate comparisons
##---------------------------------##

#assess univariate relationships
lm.pbase   <- lm(pd1year ~ pdbase, data=df)
lm.psex    <- lm(pd1year ~ Sex, data=df)
lm.prace   <- lm(pd1year ~ Race, data=df)
lm.page    <- lm(pd1year ~ age, data=df)
lm.psmoker <- lm(pd1year ~ Smoker, data=df)
lm.ptrtgrp <- lm(pd1year ~ Trtgrp, data=df)
plot(lm.pbase)
plot(lm.psex)
plot(lm.prace)
plot(lm.page)
plot(lm.psmoker)
plot(lm.ptrtgrp)
summary(lm.pbase)
summary(lm.psex)
summary(lm.prace)
summary(lm.page)
summary(lm.psmoker)
summary(lm.ptrtgrp)


###Run univariates for the 2nd outcome: attachmentloss###
lm.abase   <- lm(attach1year ~ attachbase, data=df)
lm.asex    <- lm(attach1year ~ Sex, data=df)
lm.arace   <- lm(attach1year ~ Race, data=df)
lm.aage    <- lm(attach1year ~ age, data=df)
lm.asmoker <- lm(attach1year ~ Smoker, data=df)
lm.atrtgrp <- lm(attach1year ~ Trtgrp, data=df)
plot(lm.abase)
plot(lm.asex)
plot(lm.arace)
plot(lm.aage)
plot(lm.asmoker)
plot(lm.atrtgrp)
summary(lm.abase)
summary(lm.asex)
summary(lm.arace)
summary(lm.aage)
summary(lm.asmoker)
summary(lm.atrtgrp)




##---------------------------------##
## Analysis & Hypothesis Testing
##---------------------------------##

# pd1year model -- This is the base model, adjustment for other factor next
lm.p <- lm(pd1year ~ pdbase + Trtgrp, data=df, na.action=na.exclude)
# diagnostics
plot(lm.p, which=1:5)
summary(lm.p)
# contrasts
# usualcare vs low
K23 <- matrix(c(0, 0, -1, 1, 0, 0), 1)
t23 <- summary(glht(lm.p, linfct=K23))
# usualcare vs med
K24 <- matrix(c(0, 0, -1, 0, 1, 0), 1)
t24 <- summary(glht(lm.p, linfct=K24))
# usualcare vs high
K25 <- matrix(c(0, 0, -1, 0, 0, 1), 1)
t25 <- summary(glht(lm.p, linfct=K25)) 
t23;t24;t25

#means & CI's
lsmeans(lm.p, ~Trtgrp)
latex(lsmeans(lm.p, ~Trtgrp), digits=4, file=paste("/Users/Nichole/Repositories/Bios6623/BIOS6623Class/Project0/Reports/tab3part1", ".tex", sep=""))

#####Only Smoking is related to 1 year measures, so only adjust for that variable###
# pd1year model -- This is the base model, adjustment for other factor next
lm.p.adj <- lm(pd1year ~ pdbase + Trtgrp + Smoker, data=df, na.action=na.exclude)
# diagnostics
plot(lm.p.adj, which=1:5)
summary(lm.p.adj)
# contrasts
# usualcare vs low
t23.adj <- summary(glht(lm.p.adj, linfct=K23))
# usualcare vs med
t24.adj <- summary(glht(lm.p.adj, linfct=K24))
# usualcare vs high
t25.adj <- summary(glht(lm.p.adj, linfct=K25)) 
t23.adj;t24.adj;t25.adj

#means & CI's
lsmeans(lm.p.adj, ~Trtgrp)
latex(lsmeans(lm.p.adj, ~Trtgrp), digits=4, file=paste("/Users/Nichole/Repositories/Bios6623/BIOS6623Class/Project0/Reports/tab3part2", ".tex", sep=""))

##############
###Repeat the analysis for attachment loss####  
lm.a <- lm(attach1year ~ attachbase + Trtgrp, data=df, na.action=na.exclude)
plot(lm.a)
summary(lm.a)
# contrasts
# usualcare vs low
K23 <- matrix(c(0, 0, -1, 1, 0, 0), 1)
a23 <- summary(glht(lm.a, linfct=K23))
# usualcare vs med
K24 <- matrix(c(0, 0, -1, 0, 1, 0), 1)
a24 <- summary(glht(lm.a, linfct=K24))
# usualcare vs high
K25 <- matrix(c(0, 0, -1, 0, 0, 1), 1)
a25 <- summary(glht(lm.a, linfct=K25))
a23;a24;a25
lsmeans(lm.a, ~Trtgrp)
latex(lsmeans(lm.a, ~Trtgrp), digits=4, file=paste("/Users/Nichole/Repositories/Bios6623/BIOS6623Class/Project0/Reports/tab3part3", ".tex", sep=""))

###Smoking and gender will associated with attachment loss. Adjust final analysis
###for these additional variables
lm.a.adj <- lm(attach1year ~ attachbase + Trtgrp + Smoker + Sex, data=df, na.action=na.exclude)
plot(lm.a.adj)
summary(lm.a.adj)
# contrasts
# usualcare vs low
K23 <- matrix(c(0, 0, -1, 1, 0, 0), 1)
a23 <- summary(glht(lm.a, linfct=K23))
# usualcare vs med
K24 <- matrix(c(0, 0, -1, 0, 1, 0), 1)
a24 <- summary(glht(lm.a, linfct=K24))
# usualcare vs high
K25 <- matrix(c(0, 0, -1, 0, 0, 1), 1)
a25 <- summary(glht(lm.a, linfct=K25))
a23;a24;a25
lsmeans(lm.a.adj, ~Trtgrp)
latex(lsmeans(lm.a.adj, ~Trtgrp), digits=4, file=paste("/Users/Nichole/Repositories/Bios6623/BIOS6623Class/Project0/Reports/tab3part4", ".tex", sep=""))

# manually throwing together reg & contrast table
xtable(lm.p)
xtable(lm.p.adj)
xtable(lm.a)
xtable(lm.a.adj)



################################################################################
##  END OF FILE ##  END OF FILE ##  END OF FILE ## END OF FILE ## END OF FILE ## 
################################################################################

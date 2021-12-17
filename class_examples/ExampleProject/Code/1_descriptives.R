################################################################################
## Bios6623 - Project0: analysis.R
## Author: Student A
## Description: This analysis file performs descriptives statistics and graphics
##              for Project 0. We make a table by treatment group.
##              I investigated in analysis on difference in yr1 and baseline
##              So tables and graphics were made for each.
################################################################################
###SOME FUNCTIONS ARE NOT RUNNING YET!!!####

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
## Exploratory & Descriptives
##---------------------------------##

## View summary of entire dataset
describe(df)


## Check missing data by treatment group
with(subset(df, is.na(df$pd1year)==TRUE), table(trtgrp))
with(subset(df, is.na(df$pddiff)==TRUE), table(trtgrp))
with(subset(df, is.na(df$attach1year)==TRUE), table(trtgrp))
with(subset(df, is.na(df$attachdiff)==TRUE), table(trtgrp))
##These are the same value as expected since missing in diff is only due to year 1 data

## View joint sample size by treatment & sex
with(subset(df, is.na(df$pd1year)==TRUE), table(trtgrp, sex))
with(subset(df, is.na(df$attach1year)==TRUE), table(trtgrp, sex))
## results: table shows maximum of 10 obs missing in group 5.  this limits the
##          available degrees of freedom for the regression.  any additional 
##          variables may/will reduce this number to less than ~10 per category.

## missing data comparisons
### copy pasted from print(tab1)... latex function for another day
tab2 <- tableone(vars   = c("Sex", "White", "Smoker", "age"), 
                 by     = "missing", 
                 data   = df, 
                 tests  = c("chisq.test", "t.test"),
                 digits = getOption("qwraps", 3))
print(tab2)
latex(tab2$tab,
      file = "/Users/Nichole/Repositories/Bios6623/BIOS6623Class/Project0/Reports/tabletwo.tex",
      title = "",
      ctable = TRUE,
      cgroup = tab2$cgrp,
      n.cgroup = tab2$ncgrp,
      colhead  = NULL,
      rgroup   = tab2$rgrp,
      n.rgroup = tab2$nrgrp,
      rowname  = tab2$rwnm,
      caption  = "This is the caption",
      label    = "tab:tableone",
      col.just = rep("r",ncol(tab2$tab)))


## scatterplot matrix
scatterplot.matrix(~Sex+Race+Smoker+age+attachbase+attach1year|Trtgrp, data=df)

### copy pasted from print(tab1)... latex function for another day
tab1 <- tableone(vars   = c("Sex", "White", "Smoker", "age", "sites"), 
                 by     = "Trtgrp", 
                 data   = df, 
                 tests  = c("chisq.test", "aov"),
                 digits = 2)
print(tab1)
latex(tab1$tab,
      file = "./tex/tableone.tex",
      title = "",
      ctable = TRUE,
      cgroup = tab1$cgrp,
      n.cgroup = tab1$ncgrp,
      colhead  = NULL,
      rgroup   = tab1$rgrp,
      n.rgroup = tab1$nrgrp,
      rowname  = tab1$rwnm,
      caption  = "This is the caption",
      label    = "tab:tableone",
      col.just = rep("r",ncol(tab1$tab)))






################################################################################
##  END OF FILE ##  END OF FILE ##  END OF FILE ## END OF FILE ## END OF FILE ## 
################################################################################

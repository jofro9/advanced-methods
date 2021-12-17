################################################################################
## File: data_read.R
## Author: Student A
## Date: XXX
## Description: This file reads in the raw data received from the investigator.
##              I create labels and add log transformed and difference variables.
################################################################################

###This sets my working directory
#setwd("[Your directory here]")

## Read csv file
df.orig <- read.csv(file="/Users/Nichole/Repositories/Bios6623/DATA/RawData/Project0_dental_data.csv",header=TRUE, sep=",")

## Set up dataframe, add factor versions.
df <- data.frame("id"          = df.orig$id,
                 "trtgrp"      = df.orig$trtgroup,
                 "Trtgrp"      = factor(df.orig$trtgroup,
                                        levels=c(1,2,3,4,5),
                                        labels=c("Placebo",
                                                 "Usual Care",
                                                 "Low",
                                                 "Medium",
                                                 "High")),
                 "sex"         = df.orig$sex,
                 "Sex"         = factor(df.orig$sex, 
                                        levels=c(1,2), 
                                        labels=c("Male", "Female")),
                 "race"        = df.orig$race,
                 "Race"        = factor(df.orig$race, 
                                        levels=c(5,4,2,1), 
                                        labels=c("White",
                                                 "Asian",
                                                 "African American",
                                                 "Native American")),
                 "White"       = factor(ifelse(df.orig$race==5,1,0),
                                        levels=c(0,1),
                                        labels=c("No", "Yes")),                                        
                 "age"         = df.orig$age,
                 "smoker"      = df.orig$smoker,
                 "Smoker"      = factor(df.orig$smoker, 
                                        levels=c(0,1),
                                        labels=c("No", "Yes")),
                 "sites"       = df.orig$sites,
                 "attachbase"  = df.orig$attachbase,
                 "attach1year" = df.orig$attach1year,
                 "attachdiff"  = df.orig$attach1year-df.orig$attachbase,
                 "lattachbase" = log(df.orig$attachbase),
                 "lattach1year" = log(df.orig$attach1year),
                 "pdbase"      = df.orig$pdbase,
                 "pd1year"     = df.orig$pd1year,
                 "pddiff"      = df.orig$pd1year-df.orig$pdbase,
                 "lpdbase"     = log(df.orig$pdbase),
                 "lpd1year"    = log(df.orig$pd1year),
                 "missing"     = factor(ifelse(is.na(df.orig$pd1year),1,0),
                                        levels=c(0,1),
                                        labels=c("No","Yes"))
                                        
                 )


## Review data
## esp. NA's
factorvars <- sapply(df, is.factor)
sapply(df[, factorvars], table, useNA="always")

numvars <- sapply(df, is.numeric)
sapply(df[, numvars], summary, useNA="always")

str(df)

## Save analysis dataframes
save(df, file="/Users/Nichole/Repositories/Bios6623/DATA/Project0Data/analysis.Rdata")












################################################################################
##  END OF FILE ##  END OF FILE ##  END OF FILE ## END OF FILE ## END OF FILE ## 
################################################################################

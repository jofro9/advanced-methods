####################################################################
#Program name: 2_VariableSelection.R
#
#Purpose: To apply different variable selection approaches to the 
#         an example dataset.
#Create by: Nichole Carlson
#
######################################################################

# Dependencies
library(readr)
library(dplyr)
library(psych)
library(olsrr)  #This allows for p-value backward selection and best subset selection using AIC--but inefficient for simulation for AIC.
library(MASS)  #This allows for AIC best subset, stepAIC.  This is backward selection using AIC.  Not really the goal of the AIC approach for this project.
library(glmnet)  #This does LASSO and elastic net
library(bestglm)  #This fits all possible models and picks one with lowest AIC.

###You will have to set your working directory####
setwd("/Users/nichole/Repositories/BIOS6624/Bios6624ClassExamples/VariableSelectionExample/")

### read in the first cleaned data set (psaclean.csv)

## read in the clean data from our data directory.  This has an outlier removed in weight and log transform variables avaiable.
psaclean <- read_csv("DataProcessed/psaclean2.csv")

##Variable selection with backward selection using p-values (not AIC).  You could also do backward using AIC.

##Give our model with all predictors in the model (so it know the maximum number of variables)
# Treats a factor variable correctly as a package.
###Output from the models is in the output folder with the file name below. Watch how the variables are removed by largest p-value.

bwd.model<-lm(lpsa~cavol+wt+age+svi+cappen,data=psaclean2)
sink("Output/backward_output.txt")
bwd.steps<-ols_step_backward_p(bwd.model,prem=0.15,details = TRUE)
plot(bwd.steps)
sink()


####Variable selection using best subset with AIC
###FYI, there are many packages to do this.  

bwd.model<-lm(lpsa~cavol+wt+age+svi+cappen,data=psaclean2)
bwd.steps<-ols_step_all_possible(bwd.model,details = TRUE)
plot(bwd.steps)

#Find minimum AIC
min.aic <- filter(bwd.steps,rank(aic)==1)

#Fit model with the predictors in min.aic predictor column.
final.model<-lm(lpsa~cavol+wt+svi,data=psaclean2)
summary(final.model)

##This doesn't automate well.  


####Variable selection using best subset with AIC
###Using bestglm function

#Create dataset with only the variables that I want.  Will not deal with a categorical variable like gleason for now.

psa.allvar<- psaclean2 %>% dplyr::select(cavol,wt,age,bph,svi,cappen,lpsa)

psa.aic<-bestglm(psa.allvar,IC="AIC")


####Variable selection using best subset with AIC
###Using stepAIC function from the MASS package
psa.model<-lm(lpsa~cavol+wt+age+svi+cappen,data=psaclean2)
psa.step<-stepAIC(psa.model,direction="both",k=2)
summary(psa.step)  #prints final model

###If I use k=log(n), then BIC
psa.step.bic<-stepAIC(psa.model,direction="both",k=log(96))  #n=96
summary(psa.step.bic)


##########################################################
####Now run a LASSO.  Note I am leaving out gleason for now because it isn't relevant to your project.
###FYI - my first version had alpha=0 doing ridge regression.
x<- psaclean2 %>% dplyr::select(cavol,wt,age,bph,svi,cappen)
x<-as.matrix(x)
y<- psaclean2$lpsa
psa.lasso<-glmnet(x,y,family="gaussian",alpha=1,standardize = TRUE)
plot(psa.lasso,label=TRUE)  #look at the behavior of the variables over the L1 norm (different values of lambda). 
coeff<-coef.glmnet(psa.lasso)  # prints the coefficients.  You can pick a certain value of lambda to print the coefficients for by setting s=the value you want.
coeff09<-coef.glmnet(psa.lasso,s=0.09)

###To understand the effect of lambda, set a small value close to zero and a larger value
#Watch the size of the coefficients.
psa.lasso4<-glmnet(x,y,family="gaussian",alpha=1,standardize = TRUE,lambda=0.75)
coeff.lasso4<-coef.glmnet(psa.lasso4)

psa.lasso05<-glmnet(x,y,family="gaussian",alpha=1,standardize = TRUE,lambda=0.05)
coeff.lasso05<-coef.glmnet(psa.lasso05)

###How do I choose lambda in practice.  GCV (or AIC/BIC).
psa_lassocv<-cv.glmnet(x,y,alpha=1,family="gaussian",standardize = TRUE, nfolds = 5)  # I chose 5 since the dataset is a bit small.

plot(psa_lassocv)

###Select off the lambda that minimizes MSE
lambda_cv <- psa_lassocv$lambda.min

###Grab my coefficients from this model.
coeff.min<-coef.glmnet(psa.lasso,s=psa_lassocv$lambda.min)

###Note these coefficients are from standardize variables



##########################################################
####Now run an elastic net.  Note I am leaving out gleason for now because it isn't relevant to your project.

psa.enet<-glmnet(x,y,family="gaussian",alpha=0.5,standardize = TRUE)  
plot(psa.enet,label=TRUE) #Here we can see variables and movement among the variables staying and leaving the model 
coeff.enet<-coef.glmnet(psa.enet)



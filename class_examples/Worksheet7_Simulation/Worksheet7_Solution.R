################################################
###Worksheet7_Solution.R
###
###Created by: Nichole Carlson
###Purpose: To simulate some linear regression models and compute
###         properties of the estimation
################################################3


results<-NULL


for (i in 1:20) {
  x1<-c(rep(0,50),rep(1,50))
  age<-rnorm(100,39,5)
  age.c<-age-39.66
  x.age.c<-x1*age.c

  mean1<-1+1*x1+0.5*age.c+0.25*x.age.c
  error1<-rnorm(100,0,1)

  y<-mean1+error1

  fit1<-lm(y~x1+age.c+x.age.c)
  temp.ci<-confint.lm(fit1)
  temp.pvalue<-summary(fit1)$coefficients[,4]
  temp.results<-cbind(t(fit1$coefficients),t(temp.ci[,1]),t(temp.ci[,2]),t(temp.pvalue))
  results<-rbind(results,temp.results)
}

results.sum<-matrix(nrow=20,ncol=12)

##Bias calc
results.sum[,1]<-results[,1] - 1
results.sum[,2]<-results[,2] - 1
results.sum[,3]<-results[,3] - 0.5
results.sum[,4]<-results[,4] - 0.25

##CI Calc
results.sum[,5]<-(results[,5]<=1 & results[,9]<=1) 
results.sum[,6]<-(results[,6]<=1 & results[,10]<=1)
results.sum[,7]<-(results[,7]<=0.5 & results[,11]<=0.5)
results.sum[,8]<-(results[,8]<=0.25 & results[,12]<=0.25)

##Power calc
results.sum[,9]<-results[,13]<=0.05
results.sum[,10]<-results[,14]<=0.05
results.sum[,11]<-results[,15]<=0.05
results.sum[,12]<-results[,16]<=0.05



bias<-apply(results.sum[,1:4],2,mean)
coverage<-1-apply(results.sum[,5:8],2,sum)/20
power<-apply(results.sum[,9:12],2,sum)/20
  


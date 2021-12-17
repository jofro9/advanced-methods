********************************************************;
*Program Name: MultivariableAnalysis-G1.sas		;
*Purpose: Investigate independent effects of PSA.       ;
*		through multivariable regression        ;
*Create by: Nichole Carlson				;
********************************************************;

/*Read in the final clean data from our data directory after addressing findings in initial descriptive statistics*/
PROC IMPORT OUT= WORK.psaclean2 
            DATAFILE= "C:\Repositories\Bios6623ClassExamples\PSA-Project
\G2Analysis\Data\psaclean2.csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
RUN;

/***open ods (output delivery system) to write all my figures       ***/
/***and tables to a directory and file of my choosing for editing   ***/
/***later							***/
ods listing close;
ods html path='C:\Repositories\Bios6623ClassExamples\PSA-Project\G2Analysis\Reports' file = 'MultiAnalysisFullModel.html';

/***full multivariable regression models for all variables vs. PSA****/
PROC MCMC data=psaclean2 nbi=1000 nmc=10000 plots=all DIC;  
  parms betaint 0 betacavol 0 betawt 0 betaage 0 betabph 0 betacap 0 betasvi 0 betagd6v8 0 betagd7v8 0; 
  parms sigma2 1;
  prior betain betacavol betawt betaage betapbh betacap betasvi betagd6v8 betagd7v8~ normal(mean = 0, var = 1000);
  prior sigma2 ~ igamma(shape = 2.001,scale = 1.001);
  mu = betaint + betacavol*cavol + betawt*wt + betaage*age + betabph*bph + betacap*cappen + betasvi*svi + betagd6v8*grade6 + betagd7v8*grade7;
  model lpsa ~ normal(mu, var = sigma2);
  title "Model 1: Model of LOG PSA and All variables";
RUN;title;

PROC MCMC data=psaclean2 nbi=1000 nmc=10000 plots=all DIC;  
  parms betaint 0 betawt 0 betaage 0 betabph 0 betacap 0 betasvi 0 betagd6v8 0 betagd7v8 0; 
  parms sigma2 1;
  prior betain betawt betaage betabph betacap betasvi betagd6v8 betagd7v8~ normal(mean = 0, var = 1000);
  prior sigma2 ~ igamma(shape = 2.001,scale = 1.001);
  mu = betaint + betawt*wt + betaage*age + betabph*bph + betacap*cappen + betasvi*svi + betagd6v8*grade6 + betagd7v8*grade7;
  model lpsa ~ normal(mu, var = sigma2);
  title "Model 1-1: Model of LOG PSA and All variables with cavol removed";
RUN;title;

PROC MCMC data=psaclean2 nbi=1000 nmc=10000 plots=all DIC;  
  parms betaint 0 betacavol 0 betaage 0 betabph 0 betacap 0 betasvi 0 betagd6v8 0 betagd7v8 0; 
  parms sigma2 1;
  prior betain betacavol betaage betabph betacap betasvi betagd6v8 betagd7v8~ normal(mean = 0, var = 1000);
  prior sigma2 ~ igamma(shape = 2.001,scale = 1.001);
  mu = betaint + betacavol*cavol + betaage*age + betabph*bph + betacap*cappen + betasvi*svi + betagd6v8*grade6 + betagd7v8*grade7;
  model lpsa ~ normal(mu, var = sigma2);
  title "Model 1-2: Model of LOG PSA and All variables with weight removed";
RUN;title;

PROC MCMC data=psaclean2 nbi=1000 nmc=10000 plots=all DIC;  
parms betaint 0 betacavol 0 betawt 0 betabph 0 betacap 0 betasvi 0 betagd6v8 0 betagd7v8 0; 
  parms sigma2 1;
  prior betain betacavol betawt betabph betacap betasvi betagd6v8 betagd7v8~ normal(mean = 0, var = 1000);
  prior sigma2 ~ igamma(shape = 2.001,scale = 1.001);
  mu = betaint + betacavol*cavol + betawt*wt + betabph*bph + betacap*cappen + betasvi*svi + betagd6v8*grade6 + betagd7v8*grade7;
  model lpsa ~ normal(mu, var = sigma2);
  title "Model 1-3: Model of LOG PSA and All variables with age removed";
RUN;title;

PROC MCMC data=psaclean2 nbi=1000 nmc=10000 plots=all DIC;  
 parms betaint 0 betacavol 0 betawt 0 betaage 0 betacap 0 betasvi 0 betagd6v8 0 betagd7v8 0; 
  parms sigma2 1;
  prior betain betacavol betawt betaage betacap betasvi betagd6v8 betagd7v8~ normal(mean = 0, var = 1000);
  prior sigma2 ~ igamma(shape = 2.001,scale = 1.001);
  mu = betaint + betacavol*cavol + betawt*wt + betaage*age + betabph*bph + betacap*cappen + betasvi*svi + betagd6v8*grade6 + betagd7v8*grade7;
  model lpsa ~ normal(mu, var = sigma2);
 title "Model 1-4: Model of LOG PSA and All variables with BPH removed";
RUN;title;

PROC MCMC data=psaclean2 nbi=1000 nmc=10000 plots=all DIC;  
 parms betaint 0 betacavol 0 betawt 0 betaage 0 betabph 0 betasvi 0 betagd6v8 0 betagd7v8 0; 
  parms sigma2 1;
  prior betain betacavol betawt betaage betabph betasvi betagd6v8 betagd7v8~ normal(mean = 0, var = 1000);
  prior sigma2 ~ igamma(shape = 2.001,scale = 1.001);
  mu = betaint + betacavol*cavol + betawt*wt + betaage*age + betabph*bph + betasvi*svi + betagd6v8*grade6 + betagd7v8*grade7;
  model lpsa ~ normal(mu, var = sigma2);
title "Model 1-5: Model of LOG PSA and All variables with CAPPEN removed";
RUN;title;

PROC MCMC data=psaclean2 nbi=1000 nmc=10000 plots=all DIC;  
  parms betaint 0 betacavol 0 betawt 0 betaage 0 betabph 0 betacap 0 betagd6v8 0 betagd7v8 0; 
  parms sigma2 1;
  prior betain betacavol betawt betaage betapbh betacap betagd6v8 betagd7v8~ normal(mean = 0, var = 1000);
  prior sigma2 ~ igamma(shape = 2.001,scale = 1.001);
  mu = betaint + betacavol*cavol + betawt*wt + betaage*age + betabph*bph + betacap*cappen + betagd6v8*grade6 + betagd7v8*grade7;
  model lpsa ~ normal(mu, var = sigma2);
  title "Model 1: Model of LOG PSA and All variables";
title "Model 1-6: Model of LOG PSA and All variables with SVI removed";
RUN;title;

PROC MCMC data=psaclean2 nbi=1000 nmc=10000 plots=all DIC;  
  parms beta0 0 beta1 0 beta2 0 beta3 0 beta4 0 beta5 0 beta6 0 ; 
  parms sigma2 1;
  prior beta0 beta1 beta2 beta3 beta4 beta5 beta6~ normal(mean = 0, var = 1000);
  prior sigma2 ~ igamma(shape = 2.001,scale = 1.001);
  mu = beta0 + beta1*cavol + beta2*wt + beta3*age + beta4*bph + beta5*cappen + beta6*svi;
  model lpsa ~ normal(mu, var = sigma2);
  title "Model 1-7: Model of LOG PSA and All variables with Gleason removed";
RUN;title;


ods html close;
ods listing;

/***Investigate the joint effects of cappen and other factors to understand***/
/***why loses significance												****/
/***Note: I might also check the categorical variables since they are not***/
/*** in the correlation matrix.  Those are commented out below and      ***/
/*** did not influence the coefficient on cappen like cavol did         ***/
/*** so I ruled them out as influencers of the cappen change           ****/
ods listing close;
ods html path='C:\Repositories\Bios6623ClassExamples\PSA-Project\G2Analysis\Reports' file = 'MultiAnalysisCappenInvest.html';


PROC MCMC data=psaclean2 nbi=1000 nmc=10000 plots=all DIC;  
  parms beta0 0 beta1 0 beta2 0; 
  parms sigma2 1;
  prior beta0 beta1 beta2 ~ normal(mean = 0, var = 1000);
  prior sigma2 ~ igamma(shape = 2.001,scale = 1.001);
  mu = beta0 + beta1*cavol + beta2*cappen;
  model lpsa ~ normal(mu, var = sigma2);
  title "Model 2: Log PSA and Cancer Volume and Capsular Penetration";
RUN;title;

PROC MCMC data=psaclean2 nbi=1000 nmc=10000 plots=all DIC;  
  parms beta0 0 beta1 0; 
  parms sigma2 1;
  prior beta0 beta1 ~ normal(mean = 0, var = 1000);
  prior sigma2 ~ igamma(shape = 2.001,scale = 1.001);
  mu = beta0 + beta1*cappen;
  model lpsa ~ normal(mu, var = sigma2);
  title "Model 3: Log PSA and capsular penetration";
RUN;title;


PROC MCMC data=psaclean2 nbi=1000 nmc=10000 plots=all DIC;  
  parms beta0 0 beta1 0; 
  parms sigma2 1;
  prior beta0 beta1 ~ normal(mean = 0, var = 1000);
  prior sigma2 ~ igamma(shape = 2.001,scale = 1.001);
  mu = beta0 + beta1*cavol;
  model lpsa ~ normal(mu, var = sigma2);
  title "Model 4: Log PSA and Cancer Volume";
RUN;title;

/*
PROC MCMC data=psaclean2 nbi=1000 nmc=10000 plots=all DIC;  
  parms beta0 0 beta1 0; 
  parms sigma2 1;
  prior beta0 beta1 ~ normal(mean = 0, var = 1000);
  prior sigma2 ~ igamma(shape = 2.001,scale = 1.001);
  mu = beta0 + beta1*svi;
  model lpsa ~ normal(mu, var = sigma2);
  title "Log PSA and seminal vesical invasion (yes=1)";
RUN;title;

PROC MCMC data=psaclean2 nbi=1000 nmc=10000 plots=all DIC;  
  parms beta0 0 beta1 0 beta2 0; 
  parms sigma2 1;
  prior beta0 beta1 beta2~ normal(mean = 0, var = 1000);
  prior sigma2 ~ igamma(shape = 2.001,scale = 1.001);
  mu = beta0 + beta1*grade6 + beta2*grade7;
  model lpsa ~ normal(mu, var = sigma2);
  title "Log PSA and gleason score";
RUN;title;
*/

/****Investigate why prostate weight loses signficance in the multivarible model***/

PROC MCMC data=psaclean2 nbi=1000 nmc=10000 plots=all DIC;  
  parms beta0 0 beta1 0 beta2 0 beta3 0; 
  parms sigma2 1;
  prior beta0 beta1 beta2 beta3 ~ normal(mean = 0, var = 1000);
  prior sigma2 ~ igamma(shape = 2.001,scale = 1.001);
  mu = beta0 + beta1*wt + beta2*age + beta3*bph;
  model lpsa ~ normal(mu, var = sigma2);
  title "Model 5: Model of Log PSA and wt age bph";
RUN;title;


PROC MCMC data=psaclean2 nbi=1000 nmc=10000 plots=all DIC;  
  parms beta0 0 beta1 0 ; 
  parms sigma2 1;
  prior beta0 beta1 ~ normal(mean = 0, var = 1000);
  prior sigma2 ~ igamma(shape = 2.001,scale = 1.001);
  mu = beta0 + beta1*wt ;
  model lpsa ~ normal(mu, var = sigma2);
  title "Model 6: Model of Log PSA and wt ";
RUN;title;

PROC MCMC data=psaclean2 nbi=1000 nmc=10000 plots=all DIC;  
  parms beta0 0 beta2 0 beta3 0; 
  parms sigma2 1;
  prior beta0 beta2 beta3 ~ normal(mean = 0, var = 1000);
  prior sigma2 ~ igamma(shape = 2.001,scale = 1.001);
  mu = beta0 +  beta2*age + beta3*bph;
  model lpsa ~ normal(mu, var = sigma2);
  title "Model 7: Model of Log PSA and age bph";
RUN;title;


PROC MCMC data=psaclean2 nbi=1000 nmc=10000 plots=all DIC;  
  parms beta0 0 beta1 0 beta2 0; 
  parms sigma2 1;
  prior beta0 beta1 beta2 ~ normal(mean = 0, var = 1000);
  prior sigma2 ~ igamma(shape = 2.001,scale = 1.001);
  mu = beta0 + beta1*wt + beta2*svi;
  model lpsa ~ normal(mu, var = sigma2);
  title "Model 8: Model of Log PSA and wt and SVI";
RUN;title;
	

PROC MCMC data=psaclean2 nbi=1000 nmc=10000 plots=all DIC;  
  parms beta0 0 beta1 0 beta2 0 beta3 0; 
  parms sigma2 1;
  prior beta0 beta1 beta2 beta3 ~ normal(mean = 0, var = 1000);
  prior sigma2 ~ igamma(shape = 2.001,scale = 1.001);
  mu = beta0 + beta1*wt + beta2*grade6 + beta3*grade7;
  model lpsa ~ normal(mu, var = sigma2);
  title "Model 9: Model of Log PSA and wt age bph";
RUN;title;

ods html close;
ods listing;


/***for 9/18 start on interaction work--will be ahead of lectures***/

/***for 9/20 model building***/

/***Last week put together a comprehensive analysis story from all of the output
and work on previous lectures***/

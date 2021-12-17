Model Building ? Group 2
********************************************************;
*Program Name: ModelBuilding-G2.sas		;
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
ods html path='C:\Repositories\Bios6623ClassExamples\PSA-Project\G2Analysis\Reports' file = 'ModelBuilding-G2.html';

/***full multivariable regression models for all variables vs. PSA****/
PROC MCMC data=psaclean2 nbi=1000 nmc=10000 plots=all DIC;  
  parms betaint 0 betacavol 0 betawt 0 betaage 0 betabph 0 betacap 0 betasvi 0 betagd6v8 0 betagd7v8 0; 
  parms sigma2 1;
  prior betaint betacavol betawt betaage betapbh betacap betasvi betagd6v8 betagd7v8~ normal(mean = 0, var = 1000);
  prior sigma2 ~ igamma(shape = 2.001,scale = 1.001);
  mu = betaint + betacavol*cavol + betawt*wt + betaage*age + betabph*bph + betacap*cappen + betasvi*svi + betagd6v8*grade6 + betagd7v8*grade7;
  model lpsa ~ normal(mu, var = sigma2);
  title "Model 1: Model of LOG PSA and All variables";
RUN;title;

PROC MCMC data=psaclean2 nbi=1000 nmc=10000 plots=all DIC;  
  parms betaint 0 betawt 0 betaage 0 betabph 0 betacap 0 betasvi 0 betagd6v8 0 betagd7v8 0; 
  parms sigma2 1;
  prior betaint betawt betaage betabph betacap betasvi betagd6v8 betagd7v8~ normal(mean = 0, var = 1000);
  prior sigma2 ~ igamma(shape = 2.001,scale = 1.001);
  mu = betaint + betawt*wt + betaage*age + betabph*bph + betacap*cappen + betasvi*svi + betagd6v8*grade6 + betagd7v8*grade7;
  model lpsa ~ normal(mu, var = sigma2);
  title "Model 1-1: Model of LOG PSA and All variables with cavol removed";
RUN;title;

PROC MCMC data=psaclean2 nbi=1000 nmc=10000 plots=all DIC;  
  parms betaint 0 betacavol 0 betaage 0 betabph 0 betacap 0 betasvi 0 betagd6v8 0 betagd7v8 0; 
  parms sigma2 1;
  prior betaint betacavol betaage betabph betacap betasvi betagd6v8 betagd7v8~ normal(mean = 0, var = 1000);
  prior sigma2 ~ igamma(shape = 2.001,scale = 1.001);
  mu = betaint + betacavol*cavol + betaage*age + betabph*bph + betacap*cappen + betasvi*svi + betagd6v8*grade6 + betagd7v8*grade7;
  model lpsa ~ normal(mu, var = sigma2);
  title "Model 1-2: Model of LOG PSA and All variables with weight removed";
RUN;title;

PROC MCMC data=psaclean2 nbi=1000 nmc=10000 plots=all DIC;  
parms betaint 0 betacavol 0 betawt 0 betabph 0 betacap 0 betasvi 0 betagd6v8 0 betagd7v8 0; 
  parms sigma2 1;
  prior betaint betacavol betawt betabph betacap betasvi betagd6v8 betagd7v8~ normal(mean = 0, var = 1000);
  prior sigma2 ~ igamma(shape = 2.001,scale = 1.001);
  mu = betaint + betacavol*cavol + betawt*wt + betabph*bph + betacap*cappen + betasvi*svi + betagd6v8*grade6 + betagd7v8*grade7;
  model lpsa ~ normal(mu, var = sigma2);
  title "Model 1-3: Model of LOG PSA and All variables with age removed";
RUN;title;

PROC MCMC data=psaclean2 nbi=1000 nmc=10000 plots=all DIC;  
 parms betaint 0 betacavol 0 betawt 0 betaage 0 betacap 0 betasvi 0 betagd6v8 0 betagd7v8 0; 
  parms sigma2 1;
  prior betaint betacavol betawt betaage betacap betasvi betagd6v8 betagd7v8~ normal(mean = 0, var = 1000);
  prior sigma2 ~ igamma(shape = 2.001,scale = 1.001);
  mu = betaint + betacavol*cavol + betawt*wt + betaage*age + betabph*bph + betacap*cappen + betasvi*svi + betagd6v8*grade6 + betagd7v8*grade7;
  model lpsa ~ normal(mu, var = sigma2);
 title "Model 1-4: Model of LOG PSA and All variables with BPH removed";
RUN;title;

PROC MCMC data=psaclean2 nbi=1000 nmc=10000 plots=all DIC;  
 parms betaint 0 betacavol 0 betawt 0 betaage 0 betabph 0 betasvi 0 betagd6v8 0 betagd7v8 0; 
  parms sigma2 1;
  prior betaint betacavol betawt betaage betabph betasvi betagd6v8 betagd7v8~ normal(mean = 0, var = 1000);
  prior sigma2 ~ igamma(shape = 2.001,scale = 1.001);
  mu = betaint + betacavol*cavol + betawt*wt + betaage*age + betabph*bph + betasvi*svi + betagd6v8*grade6 + betagd7v8*grade7;
  model lpsa ~ normal(mu, var = sigma2);
title "Model 1-5: Model of LOG PSA and All variables with CAPPEN removed";
RUN;title;

PROC MCMC data=psaclean2 nbi=1000 nmc=10000 plots=all DIC;  
  parms betaint 0 betacavol 0 betawt 0 betaage 0 betabph 0 betacap 0 betagd6v8 0 betagd7v8 0; 
  parms sigma2 1;
  prior betaint betacavol betawt betaage betabph betacap betagd6v8 betagd7v8~ normal(mean = 0, var = 1000);
  prior sigma2 ~ igamma(shape = 2.001,scale = 1.001);
  mu = betaint + betacavol*cavol + betawt*wt + betaage*age + betabph*bph + betacap*cappen + betagd6v8*grade6 + betagd7v8*grade7;
  model lpsa ~ normal(mu, var = sigma2);
title "Model 1-6: Model of LOG PSA and All variables with SVI removed";
RUN;title;

PROC MCMC data=psaclean2 nbi=1000 nmc=10000 plots=all DIC;  
  parms betaint 0 betacavol 0 betawt 0 betaage 0 betabph 0 betacap 0 betasvi 0; 
  parms sigma2 1;
  prior betaint betacavol betawt betaage betabph betacap betasvi ~ normal(mean = 0, var = 1000);
  prior sigma2 ~ igamma(shape = 2.001,scale = 1.001);
  mu = betaint + betacavol*cavol + betawt*wt + betaage*age + betabph*bph + betacap*cappen + betasvi*svi;
  model lpsa ~ normal(mu, var = sigma2);
  title "Model 1-7: Model of LOG PSA and All variables with Gleason removed";
RUN;title;

/***Based on DIC it looks like neither Age nor CAPPEN are needed***/
/***Since these two variables are not highly correlated, to save***/
/***run time, I am going to test whether both can be removed. ***/

PROC MCMC data=psaclean2 nbi=1000 nmc=10000 plots=all DIC diagnostics=geweke;  
  parms betaint 0 betacavol 0 betawt 0 betabph 0 betasvi 0 betagd6v8 0 betagd7v8 0; 
  parms sigma2 1;
  prior betaint betacavol betawt betabph betasvi betagd6v8 betagd7v8~ normal(mean = 0, var = 1000);
  prior sigma2 ~ igamma(shape = 2.001,scale = 1.001);
  mu = betaint + betacavol*cavol + betawt*wt + betabph*bph + betasvi*svi + betagd6v8*grade6 + betagd7v8*grade7;
  model lpsa ~ normal(mu, var = sigma2);
  title "Model 2 - Model Building: Age, cappen removed";
RUN;title;

/**DIC was 228, so okay to remove them***/
/**Repeate the process as above ****/

PROC MCMC data=psaclean2 nbi=1000 nmc=10000 plots=all DIC diagnostics=geweke;  
  parms betaint 0 betawt 0 betabph 0 betasvi 0 betagd6v8 0 betagd7v8 0; 
  parms sigma2 1;
  prior betaint betawt betabph betasvi betagd6v8 betagd7v8~ normal(mean = 0, var = 1000);
  prior sigma2 ~ igamma(shape = 2.001,scale = 1.001);
  mu = betaint + betawt*wt + betabph*bph + betasvi*svi + betagd6v8*grade6 + betagd7v8*grade7;
  model lpsa ~ normal(mu, var = sigma2);
  title "Model 3a - Model Building: cavol removed";
RUN;title; /**244 DIC--strong predictor***/


PROC MCMC data=psaclean2 nbi=1000 nmc=10000 plots=all DIC diagnostics=geweke;  
  parms betaint 0 betacavol 0 betabph 0 betasvi 0 betagd6v8 0 betagd7v8 0; 
  parms sigma2 1;
  prior betaint betacavol betabph betasvi betagd6v8 betagd7v8~ normal(mean = 0, var = 1000);
  prior sigma2 ~ igamma(shape = 2.001,scale = 1.001);
  mu = betaint + betacavol*cavol + betabph*bph + betasvi*svi + betagd6v8*grade6 + betagd7v8*grade7;
  model lpsa ~ normal(mu, var = sigma2);
  title "Model 3b - Model Building: wt removed";
RUN;title; /**DIC: 229.6**/


PROC MCMC data=psaclean2 nbi=1000 nmc=10000 plots=all DIC diagnostics=geweke;  
  parms betaint 0 betacavol 0 betawt 0 betasvi 0 betagd6v8 0 betagd7v8 0; 
  parms sigma2 1;
  prior betaint betacavol betawt betasvi betagd6v8 betagd7v8~ normal(mean = 0, var = 1000);
  prior sigma2 ~ igamma(shape = 2.001,scale = 1.001);
  mu = betaint + betacavol*cavol + betawt*wt + betasvi*svi + betagd6v8*grade6 + betagd7v8*grade7;
  model lpsa ~ normal(mu, var = sigma2);
  title "Model 3c - Model Building: bph removed";
RUN;title; /**228.8***/


PROC MCMC data=psaclean2 nbi=1000 nmc=10000 plots=all DIC diagnostics=geweke;  
  parms betaint 0 betacavol 0 betawt 0 betabph 0 betagd6v8 0 betagd7v8 0; 
  parms sigma2 1;
  prior betaint betacavol betawt betabph betagd6v8 betagd7v8~ normal(mean = 0, var = 1000);
  prior sigma2 ~ igamma(shape = 2.001,scale = 1.001);
  mu = betaint + betacavol*cavol + betawt*wt + betabph*bph + betagd6v8*grade6 + betagd7v8*grade7;
  model lpsa ~ normal(mu, var = sigma2);
  title "Model 3d - Model Building: svi removed";
RUN;title; /**DIC 234.6***/


PROC MCMC data=psaclean2 nbi=1000 nmc=10000 plots=all DIC diagnostics=geweke;  
  parms betaint 0 betacavol 0 betawt 0 betabph 0 betasvi 0 ; 
  parms sigma2 1;
  prior betaint betacavol betawt betabph betasvi ~ normal(mean = 0, var = 1000);
  prior sigma2 ~ igamma(shape = 2.001,scale = 1.001);
  mu = betaint + betacavol*cavol + betawt*wt + betabph*bph + betasvi*svi;
  model lpsa ~ normal(mu, var = sigma2);
  title "Model 3e - Model Building: gleason removed";
RUN;title;  /**DIC 230***/

/**If rigid on DIC, I would stop here***/
/**But the BPH variable is negligibly better***/
/**Remove and assess other variables ***/


PROC MCMC data=psaclean2 nbi=1000 nmc=10000 plots=all DIC diagnostics=geweke;  
  parms betaint 0  betawt 0 betasvi 0 betagd6v8 0 betagd7v8 0; 
  parms sigma2 1;
  prior betaint betawt betasvi betagd6v8 betagd7v8~ normal(mean = 0, var = 1000);
  prior sigma2 ~ igamma(shape = 2.001,scale = 1.001);
  mu = betaint + betawt*wt + betasvi*svi + betagd6v8*grade6 + betagd7v8*grade7;
  model lpsa ~ normal(mu, var = sigma2);
  title "Model 4a - Model Building: bph (and cavol) removed";
RUN;title; /**243 DIC***/

PROC MCMC data=psaclean2 nbi=1000 nmc=10000 plots=all DIC diagnostics=geweke;  
  parms betaint 0 betacavol 0 betasvi 0 betagd6v8 0 betagd7v8 0; 
  parms sigma2 1;
  prior betaint betacavol betasvi betagd6v8 betagd7v8~ normal(mean = 0, var = 1000);
  prior sigma2 ~ igamma(shape = 2.001,scale = 1.001);
  mu = betaint + betacavol*cavol + betasvi*svi + betagd6v8*grade6 + betagd7v8*grade7;
  model lpsa ~ normal(mu, var = sigma2);
  title "Model 4b - Model Building: bph (and weight) removed";
RUN;title; /**DIC 239***/


PROC MCMC data=psaclean2 nbi=1000 nmc=10000 plots=all DIC diagnostics=geweke;  
  parms betaint 0 betacavol 0 betawt 0 betagd6v8 0 betagd7v8 0; 
  parms sigma2 1;
  prior betaint betacavol betawt betagd6v8 betagd7v8~ normal(mean = 0, var = 1000);
  prior sigma2 ~ igamma(shape = 2.001,scale = 1.001);
  mu = betaint + betacavol*cavol + betawt*wt + betagd6v8*grade6 + betagd7v8*grade7;
  model lpsa ~ normal(mu, var = sigma2);
  title "Model 4c - Model Building: bph and svi removed";
RUN;title;  /**DIC 234**/


PROC MCMC data=psaclean2 nbi=1000 nmc=10000 plots=all DIC diagnostics=geweke;  
  parms betaint 0 betacavol 0 betawt 0 betasvi 0 ; 
  parms sigma2 1;
  prior betaint betacavol betawt betasvi ~ normal(mean = 0, var = 1000);
  prior sigma2 ~ igamma(shape = 2.001,scale = 1.001);
  mu = betaint + betacavol*cavol + betawt*wt + betasvi*svi;
  model lpsa ~ normal(mu, var = sigma2);
  title "Model 4d - Model Building: bph and gleason removed";
RUN;title; /***DIC 232**/

/**All the variables with higher DIC when removed, so stop here and interpret model**/
/**Final model here again for interpreting***/
PROC MCMC data=psaclean2 nbi=1000 nmc=10000 plots=all DIC diagnostics=geweke;  
  parms betaint 0 betacavol 0 betawt 0 betasvi 0 betagd6v8 0 betagd7v8 0; 
  parms sigma2 1;
  prior betaint betacavol betawt betasvi betagd6v8 betagd7v8~ normal(mean = 0, var = 1000);
  prior sigma2 ~ igamma(shape = 2.001,scale = 1.001);
  mu = betaint + betacavol*cavol + betawt*wt + betasvi*svi + betagd6v8*grade6 + betagd7v8*grade7;
  model lpsa ~ normal(mu, var = sigma2);
  title "Model Final - Model Building";
RUN;title;

ods html close;
ods listing;


/***for 9/18 start on interaction work--will be ahead of lectures***/

/***for 9/20 model building***/

/***Last week put together a comprehensive analysis story from all of the output
and work on previous lectures***/

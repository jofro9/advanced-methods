********************************************************;
*Program Name: PrelimAnalysis_Worksheet3.sas						;
*Purpose: To read in the clean dataset and create crude regression models;
*Create by: Nichole Carlson								;
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
/***later															***/
ods listing close;
ods html path='C:\Repositories\Bios6623ClassExamples\PSA-Project\G2Analysis\Reports' file = 'PrelimAnalysisPSA.html';

/***create crude regression models for AGE and GLEASON SCORE vs. LOG(PSA)****/
PROC MCMC data=psaclean2 nbi=1000 nmc=10000 plots=all DIC;  
  parms beta0 0 beta1 0; 
  parms sigma2 1;
  prior beta0 beta1 ~ normal(mean = 0, var = 1000);
  prior sigma2 ~ igamma(shape = 2.001,scale = 1.001);
  mu = beta0 + beta1*age;
  model lpsa ~ normal(mu, var = sigma2);
  title "Log PSA and age";
RUN;title;

PROC MCMC data=psaclean2 nbi=1000 nmc=10000 plots=all DIC;  
  parms beta0 0 beta1 0 beta2; 
  parms sigma2 1;
  prior beta0 beta1 beta2~ normal(mean = 0, var = 1000);
  prior sigma2 ~ igamma(shape = 2.001,scale = 1.001);
  mu = beta0 + beta1*grade6 + beta2*grade7;
  model lpsa ~ normal(mu, var = sigma2);
  title "Log(PSA) and gleason score";
RUN;title;

**MODELS FOR THE OTHER VARIABLES****
PROC MCMC data=psaclean2 nbi=1000 nmc=10000 plots=all DIC;  
  parms beta0 0 beta1 0; 
  parms sigma2 1;
  prior beta0 beta1 ~ normal(mean = 0, var = 1000);
  prior sigma2 ~ igamma(shape = 2.001,scale = 1.001);
  mu = beta0 + beta1*cavol;
  model lpsa ~ normal(mu, var = sigma2);
  title "Log PSA and Cancer Volume";
RUN;title;

PROC MCMC data=psaclean2 nbi=1000 nmc=10000 plots=all DIC;  
  parms beta0 0 beta1 0; 
  parms sigma2 1;
  prior beta0 beta1 ~ normal(mean = 0, var = 1000);
  prior sigma2 ~ igamma(shape = 2.001,scale = 1.001);
  mu = beta0 + beta1*wt;
  model lpsa ~ normal(mu, var = sigma2);
  title "Log PSA and weight";
RUN;title;



PROC MCMC data=psaclean2 nbi=1000 nmc=10000 plots=all DIC;  
  parms beta0 0 beta1 0; 
  parms sigma2 1;
  prior beta0 beta1 ~ normal(mean = 0, var = 1000);
  prior sigma2 ~ igamma(shape = 2.001,scale = 1.001);
  mu = beta0 + beta1*bph;
  model lpsa ~ normal(mu, var = sigma2);
  title "Log PSA and hyperplasia";
RUN;title;

PROC MCMC data=psaclean2 nbi=1000 nmc=10000 plots=all DIC;  
  parms beta0 0 beta1 0; 
  parms sigma2 1;
  prior beta0 beta1 ~ normal(mean = 0, var = 1000);
  prior sigma2 ~ igamma(shape = 2.001,scale = 1.001);
  mu = beta0 + beta1*cappen;
  model lpsa ~ normal(mu, var = sigma2);
  title "Log PSA and capsular penetration";
RUN;title;

PROC MCMC data=psaclean2 nbi=1000 nmc=10000 plots=all DIC;  
  parms beta0 0 beta1 0; 
  parms sigma2 1;
  prior beta0 beta1 ~ normal(mean = 0, var = 1000);
  prior sigma2 ~ igamma(shape = 2.001,scale = 1.001);
  mu = beta0 + beta1*svi;
  model lpsa ~ normal(mu, var = sigma2);
  title "Log PSA and seminal vesical invasion (yes=1)";
RUN;title;

**A MODEL NEEDED FOR TESTING SIG OF VARIABLES
PROC MCMC data=psaclean2 nbi=1000 nmc=10000 plots=all DIC;  
  parms beta0 0; 
  parms sigma2 1;
  prior beta0 ~ normal(mean = 0, var = 1000);
  prior sigma2 ~ igamma(shape = 2.001,scale = 1.001);
  mu = beta0;
  model lpsa ~ normal(mu, var = sigma2);
  title "Log PSA only: obtain DIC with no covariates for reference";
RUN;title;

ods html close;
ods listing;



********************************************************;
*Program Name: InteractionAnalysis-G2.sas		;
*Purpose: To read in the clean dataset and create graphs;
*			and descriptive statistics	;
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
/***later															***/
ods listing close;
ods html path='C:\Repositories\Bios6623ClassExamples\PSA-Project\G2Analysis\Reports' file = 'òInteractionAnalysisPSA.html';

/***create crude interaction regression models for all variables vs. PSA****/
/***I also removed the interaction for each variable so I can assess DIC***/

/******Cancer Volume and SVI***/
PROC MCMC data=psaclean2 nbi=1000 nmc=10000 plots=all DIC diagnostics=geweke;  
  parms beta0 0 beta1 0 beta2 0 beta3 0; 
  parms sigma2 1;
  prior beta0 beta1 beta2 beta3 ~ normal(mean = 0, var = 1000);
  prior sigma2 ~ igamma(shape = 2.001,scale = 1.001);
  mu = beta0 + beta1*cavol + beta2*svi + beta3*cavolsviint;
  model lpsa ~ normal(mu, var = sigma2);
  title "Log PSA and Cancer Volume-SVI Int"ù;
RUN;title;


PROC MCMC data=psaclean2 nbi=1000 nmc=10000 plots=all DIC diagnostics=geweke;  
  parms beta0 0 beta1 0 beta2 0; 
  parms sigma2 1;
  prior beta0 beta1 beta2 ~ normal(mean = 0, var = 1000);
  prior sigma2 ~ igamma(shape = 2.001,scale = 1.001);
  mu = beta0 + beta1*cavol + beta2*svi;
  model lpsa ~ normal(mu, var = sigma2);
  title "Log PSA and Cancer Volume-SVI";
RUN;title;


/****Prostate Wt and SVI Models *****/
PROC MCMC data=psaclean2 nbc=1000 nmc=10000 plots=all DIC diagnostics=geweke;  
  parms beta0 0 beta1 0 beta2 0 beta3 0; 
  parms sigma2 1;
  prior beta0 beta1 beta2 beta3 ~ normal(mean = 0, var = 1000);
  prior sigma2 ~ igamma(shape = 2.001,scale = 1.001);
  mu = beta0 + beta1*wt + beta2*svi + beta3*wtsviint;
  model lpsa ~ normal(mu, var = sigma2);
  title "Log PSA and weight-SVI int";
RUN;title;


PROC MCMC data=psaclean2 nbc=1000 nmc=10000 plots=all DIC diagnostics=geweke;  
  parms beta0 0 beta1 0 beta2 0; 
  parms sigma2 1;
  prior beta0 beta1 beta2 ~ normal(mean = 0, var = 1000);
  prior sigma2 ~ igamma(shape = 2.001,scale = 1.001);
  mu = beta0 + beta1*wt + beta2*svi;
  model lpsa ~ normal(mu, var = sigma2);
  title "Log PSA and weight-SVI";
RUN;title;

/****Age and SVI models*****/
PROC MCMC data=psaclean2 nbi=1000 nmc=10000 plots=all DIC diagnostics=geweke;  
  parms beta0 0 beta1 0 beta2 0 beta3 0; 
  parms sigma2 1;
  prior beta0 beta1 beta2 beta3~ normal(mean = 0, var = 1000);
  prior sigma2 ~ igamma(shape = 2.001,scale = 1.001);
  mu = beta0 + beta1*age + beta2*svi + beta3*agesviint;
  model lpsa ~ normal(mu, var = sigma2);
  title "Log PSA and age-SVI int"ù;
RUN;title;


PROC MCMC data=psaclean2 nbi=1000 nmc=10000 plots=all DIC diagnostics=geweke;  
  parms beta0 0 beta1 0 beta2 0; 
  parms sigma2 1;
  prior beta0 beta1 beta2~ normal(mean = 0, var = 1000);
  prior sigma2 ~ igamma(shape = 2.001,scale = 1.001);
  mu = beta0 + beta1*age + beta2*svi;
  model lpsa ~ normal(mu, var = sigma2);
  title "Log PSA and age-SVI"ù;
RUN;title;

/****BPH and SVI models****/
PROC MCMC data=psaclean2 nbc=1000 nmc=10000 plots=all DIC diagnostics=geweke;  
  parms beta0 0 beta1 0 beta2 0 beta3 0; 
  parms sigma2 1;
  prior beta0 beta1 beta2 beta3~ normal(mean = 0, var = 1000);
  prior sigma2 ~ igamma(shape = 2.001,scale = 1.001);
  mu = beta0 + beta1*bph + beta2*svi + beta3*bphsviint;
  model lpsa ~ normal(mu, var = sigma2);
  title "Log PSA and hyperplasia-SVI Int";
RUN;title;


PROC MCMC data=psaclean2 nbi=1000 nmc=10000 plots=all DIC diagnostics=geweke;  
  parms beta0 0 beta1 0 beta2 0; 
  parms sigma2 1;
  prior beta0 beta1 beta2~ normal(mean = 0, var = 1000);
  prior sigma2 ~ igamma(shape = 2.001,scale = 1.001);
  mu = beta0 + beta1*bph + beta2*svi;
  model lpsa ~ normal(mu, var = sigma2);
  title "Log PSA and hyperplasia-SVI";
RUN;title;

/****Capsular Penetration and SVI models*****/
PROC MCMC data=psaclean2 nbi=1000 nmc=10000 plots=all DIC diagnostics=geweke;  
  parms beta0 0 beta1 0 beta2 0 beta3 0; 
  parms sigma2 1;
  prior beta0 beta1 beta2 beta3~ normal(mean = 0, var = 1000);
  prior sigma2 ~ igamma(shape = 2.001,scale = 1.001);
  mu = beta0 + beta1*cappen + beta2*svi + beta3*cappensviint;
  model lpsa ~ normal(mu, var = sigma2);
  title "Log PSA and capsular penetration-SVI Int";
RUN;title;


PROC MCMC data=psaclean2 nbi=1000 nmc=10000 plots=all DIC diagnostics=geweke;  
  parms beta0 0 beta1 0 beta2 0; 
  parms sigma2 1;
  prior beta0 beta1 beta2 ~ normal(mean = 0, var = 1000);
  prior sigma2 ~ igamma(shape = 2.001,scale = 1.001);
  mu = beta0 + beta1*cappen + beta2*svi;
  model lpsa ~ normal(mu, var = sigma2);
  title "Log PSA and capsular penetration-SVI";
RUN;title;


/****gleason and SVI models***/
PROC MCMC data=psaclean2 nbi=1000 nmc=10000 plots=all DIC diagnostics=geweke;  
  parms beta0 0 beta1 0 beta2 0 beta3 0 beta4 0 beta5 0; 
  parms sigma2 1;
  prior beta0 beta1 beta2 beta3 beta4 beta5 ~ normal(mean = 0, var = 1000);
  prior sigma2 ~ igamma(shape = 2.001,scale = 1.001);
  mu = beta0 + beta1*grade6 + beta2*grade7 + beta3*svi + beta4*gd6sviint + beta5*gd7sviint;
  model lpsa ~ normal(mu, var = sigma2);
  title "Log PSA and gleason score-SVI Int";
RUN;title;

PROC MCMC data=psaclean2 nbi=1000 nmc=10000 plots=all DIC diagnostics=geweke;  
  parms beta0 0 beta1 0 beta2 0 beta3 0; 
  parms sigma2 1;
  prior beta0 beta1 beta2 beta3~ normal(mean = 0, var = 1000);
  prior sigma2 ~ igamma(shape = 2.001,scale = 1.001);
  mu = beta0 + beta1*grade6 + beta2*grade7 + beta3*svi;
  model lpsa ~ normal(mu, var = sigma2);
  title "Log PSA and gleason score-SVI";
RUN;title;

ods html close;
ods listing;



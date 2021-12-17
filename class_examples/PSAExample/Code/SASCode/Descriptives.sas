********************************************************;
*Program Name: Desccriptives-G1.sas						;
*Purpose: To read in the clean dataset and create graphs;
*			and descriptive statistics					;
*Create by: Nichole Carlson								;
********************************************************;

/*Read in the clean data from our data directory*/
PROC IMPORT OUT= WORK.psaclean 
            DATAFILE= "C:\Repositories\Bios6623ClassExamples\PSA-Project
\G1Analysis\Data\psaclean.csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
RUN;

/***open ods (output delivery system) to write all my figures       ***/
/***and tables to a directory and file of my choosing for editing   ***/
/***later															***/
ods listing close;
ods html path='C:\Repositories\Bios6623ClassExamples\PSA-Project\G1Analysis\Reports' file = 'Descriptives.html';

/***create scatter plots with PSA and log(PSA) as the outcome****/
PROC GPLOT data=psaclean;
  plot psa*(cavol wt age bph cappen);
  symbol I=r1 value=dot color=black; 
RUN;

PROC MEANS data=psaclean;
	var psa lpsa cavol wt age bph cappen;
	title 'Descriptive Statistics for Continuous Variables in PSA-Project';
RUN;title;

PROC FREQ data=psaclean;
	tables svi gleason;
	title 'Descriptive Statistics for the Categorical Variables in PSA-Project';
RUN;title;

ods html close;
ods listing;

/*Read in the final clean data from our data directory after addressing findings in initial descriptive statistics*/
PROC IMPORT OUT= WORK.psaclean2 
            DATAFILE= "C:\Repositories\Bios6623ClassExamples\PSA-Project
\G1Analysis\Data\psaclean2.csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
RUN;

/***open ods (output delivery system) to write all my figures       ***/
/***and tables to a directory and file of my choosing for editing   ***/
/***later															***/
ods listing close;
ods html path='C:\Repositories\Bios6623ClassExamples\PSA-Project\G1Analysis\Reports' file = 'DescriptivesClean.html';

/***create scatter plots with PSA and log(PSA) as the outcome with data issues addressed****/
PROC GPLOT data=psaclean2;
  plot psa*(cavol wt age bph cappen);
  symbol I=r1 value=dot color=black; 
RUN;

PROC GPLOT data=psaclean2;
  plot lpsa*(cavol wt age bph cappen);
  symbol I=r1 value=dot color=black; 
RUN;

PROC MEANS data=psaclean2;
	var psa lpsa cavol wt age bph cappen;
	title 'Descriptive Statistics for Continuous Variables in PSA-Project';
	title2 'After data cleaning';
RUN;title;title2;

PROC FREQ data=psaclean2;
	tables svi gleason;
	title 'Descriptive Statistics for the Categorical Variables in PSA-Project';
	title2 'After data cleaning';
RUN;title;title2;


PROC CORR data=psaclean2;
	var psa lpsa cavol wt age bph cappen;
	title 'Correlations between all variables: After data cleaning';
RUN;title;

ods html close;
ods listing;

ods listing close;
ods html path='C:\Repositories\Bios6623ClassExamples\PSA-Project\G1Analysis\Reports' file = 'Descriptives.html';

/***create scatter plots with PSA and log(PSA) as the outcome****/
PROC GPLOT data=psaclean;
  plot psa*(cavol wt age bph cappen);
  symbol I=r1 value=dot color=black; 
RUN;

PROC MEANS data=psaclean;
	var psa lpsa cavol wt age bph cappen;
	title 'Descriptive Statistics for Continuous Variables in PSA-Project';
RUN;title;

PROC FREQ data=psaclean;
	tables svi gleason;
	title 'Descriptive Statistics for the Categorical Variables in PSA-Project';
RUN;title;

ods html close;
ods listing;

/*Read in the final clean data from our data directory after addressing findings in initial descriptive statistics*/
PROC IMPORT OUT= WORK.psaclean2 
            DATAFILE= "C:\Repositories\Bios6623ClassExamples\PSA-Project
\G1Analysis\Data\psaclean2.csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
RUN;

/***Descriptives for SVI interaction hypothesis***/
ods listing close;
ods html path='C:\Repositories\Bios6623ClassExamples\PSA-Project\G1Analysis\Reports' file = 'DescriptivesCleanInteraction.html';

PROC SORT data=psaclean2;
    by svi;
RUN;

/***create scatter plots with PSA and log(PSA) as the outcome with data issues addressed****/
PROC GPLOT data=psaclean2;
  plot lpsa*(cavol wt age bph cappen)=svi;
  symbol1 I=r1 value=dot color=black; 
  symbol2 I=r1 value=dot color=red;
RUN;

PROC MEANS data=psaclean2;
	by svi;
	var psa lpsa cavol wt age bph cappen;
	title 'Descriptive Statistics for Continuous Variables in PSA-Project';
	title2 '˜by SVI';
RUN;title;title2;

PROC FREQ data=psaclean2;
	tables svi*gleason;
	title 'Descriptive Statistics for the Categorical Variables in PSA-Project';
RUN;title;


PROC CORR data=psaclean2;
	by SVI;
	var psa lpsa cavol wt age bph cappen;
	title 'Correlations between all variables: By SVI';
RUN;title;

ods html close;
ods listing;


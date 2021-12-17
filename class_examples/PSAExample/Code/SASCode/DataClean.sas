********************************************************;
*Program Name: DataClean-G1.sas							;
*Purpose: To read in the dataset and label the variables;
*			Create any new variables.					;
*Output: A Clean .csv file is exported for reading in to;
*			all other analysis files.					;
*Create by: Nichole Carlson								;
********************************************************;

***READ IN THE RAW DATA AND CREATE VARIABLE TRANSFORMS****;
***Note: This is the place where I would come back to ****;
***      for any data cleaning, outlier removal more  ****;
***      transforms.  This keeps all data cleaning    ****;
***		 and variable decisions in one place.		  ****;
***
***;
**********************************************************;
PROC IMPORT OUT= WORK.PSA 
            DATAFILE= "C:\Repositories\Bios6623ClassExamples\PSA-Project
\G1Analysis\DataRaw\prostate.csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
RUN;

DATA psaclean;set psa;
	lpsa = log(psa);
	grade6 = 0;
	if gleason = 6 then grade6 = 1;
	grade7 = 0;
	if gleason = 7 then grade7 = 1;
	grade8 = 0;
	if gleason = 8 then grade8 = 1;

	gd6ageint = age*grade6;
	gd7ageint = age*grade7;
	gd8ageint = age*grade8;
RUN;

***Export the clean dataset for use in other analysis programs****;
PROC EXPORT DATA= WORK.PSACLEAN 
            OUTFILE= "C:\Repositories\Bios6623ClassExamples\PSA-Project\
G1Analysis\Data\psaclean.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;

***Further data cleaning based on reviewing descriptive statistics***;
DATA psaclean2;set psaclean;
	if wt > 400 then delete;
	
      cavolsviint = cavol*svi;
	wtsviint = wt*svi;
	agesviint = age*svi;
	bphsviint = bph*svi;
	cappensviint = cappen*svi;
	gd6sviint = grade6*svi;
	gd7sviint = grade7*svi;
	gd8sviint = grade8*svi;
	
RUN;

***Export the clean dataset for use in other analysis programs****;
PROC EXPORT DATA= WORK.PSACLEAN2 
            OUTFILE= "C:\Repositories\Bios6623ClassExamples\PSA-Project\
G1Analysis\Data\psaclean2.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;


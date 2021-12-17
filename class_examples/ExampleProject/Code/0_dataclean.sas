/*################################################################################
## File: data_read.SAS
## Author: Student A
## Date: XXX
## Description: This file reads in the raw data received from the investigator.
##              I create labels and add log transformed and difference variables.
################################################################################*/

***CREATE THE LIBRARY FOR WRITING THE ANALYSIS DATASET***;
LIBNAME datadir “E:\Users\Nichole\Repositories\Bios6623\DATA\Project0Data\”;

* Read csv file using import command. I copied the code here from a saved SAS file;

PROC IMPORT OUT= WORK.dental 
     DATAFILE= "E:\Users\Nichole\Repositories\Bios6623\DATA\RawDa
ta\Project0_dental_data.csv" TERMSTR=CR DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
RUN;

* Create new variables;
DATA DENTAL;SET DENTAL;
	White = Race;  *create a white indicator for investigating race/ethnicity;
	if Race = 5 then White = 1;
	else if Race < 5 then White = 0; **check this code**;

	attachdiff = attach1year - attachbase;
	lattachbase = log(attachbase);
	lattach1year = log(attach1year);
	pddiff = pd1year - pdbase;
	lpdbase = log(pdbase);
	lpd1year = log(pd1year);
	missing = 
	
**Create the formats for categorical variables for interpreting output;
PROC FORMAT;
     VALUE Trtgrp 1=‘Placebo’ 2=‘Usual Care’ 3=‘Low’ 4=‘Medium’ 5=‘High’;
     VALUE Sex 1=‘Male’ 2=‘Female’;
     VALUE Race 5=‘White’ 4=‘Asian’ 2=‘African American’ 1=‘Native American’;
     VALUE White 1=‘Yes’ 0=‘No’;
     VALUE Smoker 1=‘Yes’ 0=‘No’;
RUN;
 
**Set the formats and write out the permanent SAS dataset for analysis;
DATA datadir.DENTAL;SET DENTAL;
    FORMAT Trtgrp Trtgrp.;
	FORMAT Sex Sex.;
	FORMAT Race Race.;
	FORMAT White White.;
	FORMAT Smoker Smoker.;
RUN;


/*################################################################################
##  END OF FILE ##  END OF FILE ##  END OF FILE ## END OF FILE ## END OF FILE ## 
################################################################################*/

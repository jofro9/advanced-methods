/*################################################################################
## Bios6623 - Project0: descriptives.R
## Author: Student A
## Description: This analysis file performs descriptives statistics and graphics
##              for Project 0. We make a table by treatment group.
##              I investigated in analysis on difference in yr1 and baseline
##              So tables and graphics were made for each.
## NOTE: This code generates analysis, but it doesn’t produce pretty tables for using ##       in reports.
################################################################################*/

***CREATE THE LIBRARY FOR READING IN THE ANALYSIS DATASET***;
LIBNAME datadir “E:\Users\Nichole\Repositories\Bios6623\DATA\Project0Data\”;

DATA DENTAL;SET datadir.dental;RUN;

/*##CODE HAS NOT BEEN CHECKED YET!!!####*/


/*#————————————————##
## Exploratory & Descriptives
##————————————————##*/

*Check summary measures for all continuous variables;
PROC UNIVARIATE DATA=DENTAL PLOTS;
	VAR SITES AGE ATTACHBASE ATTACH1YEAR ATTACHDIFF LATTACHBASE LATTACH1YEAR PDBASE PD1YEAR PDDIFF LPDBASE LPD1YEAR;
TITLE ’UNIVARIATE STATS FOR EACH VARIABLE OVERALL’;
RUN;TITLE;

**Check summary measures by group. This will also tell us about missing values**;
PROC SORT DATA=DENTAL;BY TRTGRP;RUN;
PROC UNIVARIATE DATA=DENTAL PLOTS;
	BY TRTGRP;
	VAR SITES AGE ATTACHBASE ATTACH1YEAR ATTACHDIFF LATTACHBASE LATTACH1YEAR PDBASE PD1YEAR PDDIFF LPDBASE LPD1YEAR;
TITLE ’UNIVARIATE STATS FOR EACH VARIABLE BY TREATMENT GROUP’;
RUN;TITLE;

**Check counts for categorical variables overall;
PROC FREQ DATA=DENTAL MISSING;
	TABLE TRTGRP SEX RACE WHITE SMOKER;
	TITLE ‘COUNTS FOR CATEGORICAL VARIABLES’;
RUN;TITLE;

**Repeat tables for categorical variables by treatment group;
PROC FREQ DATA=DENTAL MISSING;
	TABLE TRTGRP * (SEX RACE WHITE SMOKER);
	TITLE ‘COUNTS FOR CATEGORICAL VARIABLES BY TREATMENT GROUP’;
RUN;TITLE;

/**ADD SUMMARY INFO NOTES AFTER RUNNING CODE**/
/*################################################################################
##  END OF FILE ##  END OF FILE ##  END OF FILE ## END OF FILE ## END OF FILE ##
################################################################################*/

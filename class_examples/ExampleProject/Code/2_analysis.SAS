/*################################################################################
## Project0: analysis.SAS
## Author: Student A
## Description: This file conducts the univariate and multivariable models
##              to assess whether groups differ on the two outcomes.
##              My analysis strategy is to use the year 1 data as the outcome,
##              because it is more powerful than the difference. I will adjust for
##              baseline of the outcome for precision.
################################################################################*/
***CODE NOT CHECKED YET!***;

***CREATE THE LIBRARY FOR READING IN THE ANALYSIS DATASET***;
LIBNAME datadir “E:\Users\Nichole\Repositories\Bios6623\DATA\Project0Data\”;

DATA DENTAL;SET datadir.dental;RUN;

/*##————————————————##
## Univariate comparisons for 1 year pocket depth as outcome
##————————————————##*/

PROC GLM DATA=DENTAL;
	MODEL PD1YEAR = PDBASE/S;
RUN;

PROC GPLOT DATA=DENTAL;
	PLOT PD1YEAR * PDBASE;
RUN;

PROC GLM DATA=DENTAL;
	CLASS SEX;
	MODEL PD1YEAR = SEX/S;
RUN;

**Check distributions of outcome by sex;
**sort first;

PROC SORT DATA=DENTAL;BY SEX;RUN;

PROC UNIVARIATE DATA=DENTAL PLOT;
	BY SEX;
	VAR PD1YEAR; RUN;

PROC GLM DATA=DENTAL;
	CLASS RACE;
	MODEL PD1YEAR = RACE/S;
RUN;

**Check distributions of outcome by RACE;
**sort first;

PROC SORT DATA=DENTAL;BY RACE;RUN;

PROC UNIVARIATE DATA=DENTAL PLOT;
	BY RACE;
	VAR PD1YEAR; RUN;


PROC GLM DATA=DENTAL;
	MODEL PD1YEAR = AGE/S;
RUN;

PROC GPLOT DATA=DENTAL;
	PLOT PD1YEAR*AGE;
RUN;

PROC GLM DATA=DENTAL;
	MODEL PD1YEAR = SMOKER/S;
RUN;

**Check distributions of outcome by SMOKER;
**sort first;

PROC SORT DATA=DENTAL;BY SMOKER;RUN;

PROC UNIVARIATE DATA=DENTAL PLOT;
	BY SMOKER;
	VAR PD1YEAR; RUN;


PROC GLM DATA=DENTAL;
	MODEL PD1YEAR = TRTGRP/S;
RUN;

**Check distributions of outcome by TRTGRP;
**sort first;

PROC SORT DATA=DENTAL;BY TRTGRP;RUN;

PROC UNIVARIATE DATA=DENTAL PLOT;
	BY TRTGRP;
	VAR PD1YEAR; RUN;

/**************************************************/
***REPEAT UNIVARIATE FOR ATTACHMENT LOSS AT 1 YEAR;
/**************************************************/
PROC GLM DATA=DENTAL;
	MODEL ATTACH1YEAR = PDBASE/S;
RUN;

PROC GPLOT DATA=DENTAL;
	PLOT ATTACH1YEAR * PDBASE;
RUN;

PROC GLM DATA=DENTAL;
	CLASS SEX;
	MODEL ATTACH1YEAR = SEX/S;
RUN;

**Check distributions of outcome by sex;
**sort first;

PROC SORT DATA=DENTAL;BY SEX;RUN;

PROC UNIVARIATE DATA=DENTAL PLOT;
	BY SEX;
	VAR ATTACH1YEAR; RUN;

PROC GLM DATA=DENTAL;
	CLASS RACE;
	MODEL ATTACH1YEAR = RACE/S;
RUN;

**Check distributions of outcome by RACE;
**sort first;

PROC SORT DATA=DENTAL;BY RACE;RUN;

PROC UNIVARIATE DATA=DENTAL PLOT;
	BY RACE;
	VAR ATTACH1YEAR; RUN;


PROC GLM DATA=DENTAL;
	MODEL ATTACH1YEAR = AGE/S;
RUN;

PROC GPLOT DATA=DENTAL;
	PLOT ATTACH1YEAR*AGE;
RUN;

PROC GLM DATA=DENTAL;
	MODEL ATTACH1YEAR = SMOKER/S;
RUN;

**Check distributions of outcome by SMOKER;
**sort first;

PROC SORT DATA=DENTAL;BY SMOKER;RUN;

PROC UNIVARIATE DATA=DENTAL PLOT;
	BY SMOKER;
	VAR ATTACH1YEAR; RUN;


PROC GLM DATA=DENTAL;
	MODEL ATTACH1YEAR = TRTGRP/S;
RUN;

**Check distributions of outcome by TRTGRP;
**sort first;

PROC SORT DATA=DENTAL;BY TRTGRP;RUN;

PROC UNIVARIATE DATA=DENTAL PLOT;
	BY TRTGRP;
	VAR ATTACH1YEAR; RUN;


***ADD SOME TESTING OF DIFFERENCES IN VARIABLES ACROSS GROUPS.  IN OBS TRIALS THIS WOULD HELP IDENTIFY CONFOUNDERS.  HERE WE LOOK FOR TREATMENT IMBALANCES TO PUT IN THE DISCUSSION SECTION OR PERHAPS PUT IN A MODEL TO ADJUST***;

****stopped here updating code with SAS code
##---------------------------------##
## Analysis & Hypothesis Testing
##---------------------------------##

# pd1year model -- This is the base model, adjustment for other factor next
lm.p <- lm(pd1year ~ pdbase + Trtgrp, data=df, na.action=na.exclude)
# diagnostics
plot(lm.p, which=1:5)
summary(lm.p)
# contrasts
# usualcare vs low
K23 <- matrix(c(0, 0, -1, 1, 0, 0), 1)
t23 <- summary(glht(lm.p, linfct=K23))
# usualcare vs med
K24 <- matrix(c(0, 0, -1, 0, 1, 0), 1)
t24 <- summary(glht(lm.p, linfct=K24))
# usualcare vs high
K25 <- matrix(c(0, 0, -1, 0, 0, 1), 1)
t25 <- summary(glht(lm.p, linfct=K25))
t23;t24;t25

#means & CI's
lsmeans(lm.p, ~Trtgrp)
latex(lsmeans(lm.p, ~Trtgrp), digits=4, file=paste("/Users/Nichole/Repositories/Bios6623/BIOS6623Class/Project0/Reports/tab3part1", ".tex", sep=""))

#####Only Smoking is related to 1 year measures, so only adjust for that variable###
# pd1year model -- This is the base model, adjustment for other factor next
lm.p.adj <- lm(pd1year ~ pdbase + Trtgrp + Smoker, data=df, na.action=na.exclude)
# diagnostics
plot(lm.p.adj, which=1:5)
summary(lm.p.adj)
# contrasts
# usualcare vs low
t23.adj <- summary(glht(lm.p.adj, linfct=K23))
# usualcare vs med
t24.adj <- summary(glht(lm.p.adj, linfct=K24))
# usualcare vs high
t25.adj <- summary(glht(lm.p.adj, linfct=K25))
t23.adj;t24.adj;t25.adj

#means & CI's
lsmeans(lm.p.adj, ~Trtgrp)
latex(lsmeans(lm.p.adj, ~Trtgrp), digits=4, file=paste("/Users/Nichole/Repositories/Bios6623/BIOS6623Class/Project0/Reports/tab3part2", ".tex", sep=""))

##############
###Repeat the analysis for attachment loss####
lm.a <- lm(attach1year ~ attachbase + Trtgrp, data=df, na.action=na.exclude)
plot(lm.a)
summary(lm.a)
# contrasts
# usualcare vs low
K23 <- matrix(c(0, 0, -1, 1, 0, 0), 1)
a23 <- summary(glht(lm.a, linfct=K23))
# usualcare vs med
K24 <- matrix(c(0, 0, -1, 0, 1, 0), 1)
a24 <- summary(glht(lm.a, linfct=K24))
# usualcare vs high
K25 <- matrix(c(0, 0, -1, 0, 0, 1), 1)
a25 <- summary(glht(lm.a, linfct=K25))
a23;a24;a25
lsmeans(lm.a, ~Trtgrp)
latex(lsmeans(lm.a, ~Trtgrp), digits=4, file=paste("/Users/Nichole/Repositories/Bios6623/BIOS6623Class/Project0/Reports/tab3part3", ".tex", sep=""))

###Smoking and gender will associated with attachment loss. Adjust final analysis
###for these additional variables
lm.a.adj <- lm(attach1year ~ attachbase + Trtgrp + Smoker + Sex, data=df, na.action=na.exclude)
plot(lm.a.adj)
summary(lm.a.adj)
# contrasts
# usualcare vs low
K23 <- matrix(c(0, 0, -1, 1, 0, 0), 1)
a23 <- summary(glht(lm.a, linfct=K23))
# usualcare vs med
K24 <- matrix(c(0, 0, -1, 0, 1, 0), 1)
a24 <- summary(glht(lm.a, linfct=K24))
# usualcare vs high
K25 <- matrix(c(0, 0, -1, 0, 0, 1), 1)
a25 <- summary(glht(lm.a, linfct=K25))
a23;a24;a25
lsmeans(lm.a.adj, ~Trtgrp)
latex(lsmeans(lm.a.adj, ~Trtgrp), digits=4, file=paste("/Users/Nichole/Repositories/Bios6623/BIOS6623Class/Project0/Reports/tab3part4", ".tex", sep=""))

# manually throwing together reg & contrast table
xtable(lm.p)
xtable(lm.p.adj)
xtable(lm.a)
xtable(lm.a.adj)



################################################################################
##  END OF FILE ##  END OF FILE ##  END OF FILE ## END OF FILE ## END OF FILE ##
################################################################################

# Advanced Statistcal Methods and Analysis

### PI: Nichole Carlson  
### Analyst: Joseph Froelicher

This is the repository for the project work for Project 4 Bios6624 for Joseph Froelicher.

The purpose of this analysis is to investigate several common variable selection techniques used primarily in linear regression. We are most interested in comparing the reliability of these variable selection methods when it comes to several measures of precision and accuracy, which will be detailed below. Most importantly, we want to know if variable selection techniques can pick out variables that are already known to be significant. And how that relates to the effect size of the variable coefficients. Additionally, we are also interested in how well variable selection techniques ignore variables that are not statistically significant.

To investigate variable selection techniques, the analysis team will perform a simulation of data. There are two primary simulations for this analysis. The first is to simulate data using a sample size of 250 and the second is to simulated data using a sample size of 500. For each sample size, we will simulate both correlated and uncorrelated data. The details of those simulations are outlined below.

It is hypothesized that all variable selection methods will select variables with larger effect sizes more accurately than those with small effect sizes. It is also apparent that more conservative selection methods will wrongly select insignificant variables less frequently that the more anti-conservative selection techniques. Lastly, it is hypothesized that methods that use penalized regression techniques will select insignificant variables at a much higher rate than other selection methods.

Details about the folder structure for a project:

File | Description
---|---------------------------------------------------------------------
Background	| contains the background information for the analysis
Code		| contains all R scripts for this project
DataRaw		| contain all raw data provided by investigators
DataProcessed	| contains the processed data used for analysis
Reports		| contains all output, rmarkdown files and report
Figures     | contains figures used to generate reports


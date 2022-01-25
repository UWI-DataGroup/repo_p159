** HEADER -----------------------------------------------------
**  DO-FILE METADATA
    //  algorithm name					MIsurvival.do
    //  project:				            BNR
    //  analysts:				       	Christina Howitt
    //  algorithm task			            Prepare data for MI survival analysis using 10 years of BNR data


    ** General algorithm set-up
    version 16
    clear all
    macro drop _all
    set more 1
    set linesize 120

    ** Set working directories: this is for DATASET and LOGFILE import and export
    ** DATASETS to encrypted SharePoint folder
    local datapath "X:\The University of the West Indies\DataGroup - repo_data\data_p159"
     ** LOGFILES to unencrypted OneDrive folder (.gitignore set to IGNORE log files on PUSH to GitHub)
    local logpath X:/OneDrive - The University of the West Indies/repo_datagroup/repo_p159
    ** GRAPHS to project output folder
    local outputpath "X:\The University of the West Indies\DataGroup - DG_Projects\PROJECT_p159"

    ** Close any open log file and open a new log file
    capture log close
    log using "`logpath'\MIsurvival_prep", replace

** HEADER -----------------------------------------------------

**-------------------------------------------------------------
** PART ONE: STANDARD POPULATION FOR INCIDENCE 
**-------------------------------------------------------------

** Standard World population (WHO 2002 standard)
** REFERENCE
** Age-standardization of rates: A new WHO-standard, 2000. Ahmad OB, Boschi-Pinto C, Lopez AD, Murray CJL, Lozano R, Inoue M.
** GPE Discussion paper Series No: 31 EIP/GPE/EBD. World Health Organization

		drop _all
		input age10 pop
            1	17550
            2	17070
            3	16150
            4	14760
            5	12630
            6	9920
            7	6680
            8	3730
            9	1545
		end

** TEN age groups in TEN-year bands. 

label define age10_lab  1 "0-9" 2 "10-19"  3 "20-29"	///
		            4 "30-39"  5 "40-49"  6 "50-59"	///
				7 "60-69"  8 "70-79"  9 "80 & over" , modify
label values age10 age10_lab

save "`datapath'\version01\1-input\who2000_10-1.dta", replace

*/


**-------------------------------------------------------------
** PART TWO: PREPARE BNR DATA FOR ANALYSIS
**-------------------------------------------------------------

* ----------------------
** STEP 1. Load Data
* ----------------------
/*
use "`datapath'\version01\1-input\heart_2009-2019_v8_anonymisedFUdata_clean.dta", clear

** -------------------------------------------------------------------------------------------------------------------- 
** DATA checking     
** -------------------------------------------------------------------------------------------------------------------- 

*Reduce dataset to what is needed as there are 1290 variables
*keep anon_pid age sex etype dom t_ami tom toh doh hage hfu1date disdt dodtod dodyear deathdate dobabs sudd vstatus f1vstatus age5 age10 age_10 death dcat f2vstatus

count   // there are 3831 observations

** Are they all confirmed cases? 
count if dom==. // 131 are missing dates of MI. What to do with them? DROP!

** Which variable is correct for vital status at 28-day FU? Different results for death vs f1vstatus. 
tab death, miss

/*

 death |      Freq.     Percent        Cum.
-----------------+-----------------------------------
 dead by 28 days |      2,291       59.80       59.80
alive at 28 days |      1,394       36.39       96.19
               . |        146        3.81      100.00
-----------------+-----------------------------------
           Total |      3,831      100.0

*/

tab f1vstatus, miss // use for 28 days
/*
      Vital |
  Status at |
 28-Day F/U |      Freq.     Percent        Cum.
------------+-----------------------------------
      Alive |      1,284       33.52       33.52
       Dead |      2,154       56.23       89.74
    Unknown |         85        2.22       91.96
          . |        308        8.04      100.00
------------+-----------------------------------
      Total |      3,831      100.00
*/

tab vstatus, miss // again this does not tie in with previous results. useless
/*
     Vital |
     Status |      Freq.     Percent        Cum.
------------+-----------------------------------
      Alive |      1,524       39.78       39.78
   Deceased |        735       19.19       58.97
          . |      1,572       41.03      100.00
------------+-----------------------------------
      Total |      3,831      100.00
*/

tab f2vstatus, miss


/*
      Vital |
status at 1 |
       year |      Freq.     Percent        Cum.
------------+-----------------------------------
      alive |        962       25.11       25.11
       dead |        600       15.66       40.77
         99 |         87        2.27       43.04   what does this mean? Can't get hold of them?
          . |      2,182       56.96      100.00
------------+-----------------------------------
      Total |      3,831      100.00
*/

** SO SHOULD I DERIVE SURVIVAL TIME USING TIME OF MI AND TIME OF DEATH OR LOSS TO FOLLOW UP?
** Lots of variables called  fu2... labelled as 28 day follow up. Do they mean 1 yr follow up?  

** Decision made on 9th July 2021 with Natasha to recreate dataset as problems were caused by merge

use "X:\The University of the West Indies\DataGroup - repo_data\data_p159\heart_2009-2019_v8_anonymisedFUdata_Stata_v16_clean(13-Jul-2021)", clear 

** Keep only bare essentials for survival analysis
tempfile `heart28' `2015heart1yr' `2016heart1yr' `2017heart1yr' `2018heart1yr' `2019heart1yr'
keep anon_pid age sex dom deathdate f1vstatus  
save `heart28', replace

keep anon_pid age sex dom deathdate fu1doa f1vstatus fu2doa f2vstatus

order fu1doa, before(f1vstatus)
*/


******************************************************************************************************
*** NEW DATASET CREATED ON 13-JULY-2021 AS A RESULT OF DATA CLEANING
******************************************************************************************************
** Load clean dataset
use "X:\The University of the West Indies\DataGroup - repo_data\data_p159\heart_2009-2019_v8_anonymisedFUdata_Stata_v16_clean(13-Jul-2021)", clear


**-------------------------------------------------
** PREPARE DATA FOR ANALYSIS
**-------------------------------------------------

keep abstracted anon_pid age sex dom deathdate fu1doa f1vstatus fu2doa f2vstatus hfu1date

** check variables for consistency
codebook anon_pid 
codebook age 
codebook sex // 1 coded as 99 and 1 as missing. Change for consistency
      recode sex 99=.
codebook dom // 1 is missing, so have to drop
      drop if dom==.
codebook deathdate // missing .:  1,567/4,248

** Create variable for year of MI
gen year = year(dom)
label variable year "year of MI event"
order year, after(dom)

** we are only using full years, so 2009 to be dropped
drop if year == 2009
sort dom 
codebook dom

** Create a variable for survival time
gen time = .
replace time = deathdate - dom 
replace time = fu2doa - dom if deathdate==.
label variable time "survival time (days)"
** data checks
list anon_pid if f2vstatus==2 & deathdate==. // there are 42 people with missing death dates but vital status is dead. Ashley to provide death dates
      replace time = .z if f2vstatus==2 & deathdate==.

count if dom !=. & deathdate==. & f1vstatus==. & f2vstatus==. & fu1doa==. & fu2doa==. // 124 people in this category
      replace time = .z if dom !=. & deathdate==. & f1vstatus==. & f2vstatus==. & fu1doa==. & fu2doa==.

codebook time
sort time   // lots of missing values of time because no dates of abstraction for follow ups, but vital status states alive. Need some guidance here.


**SAVE DATA FOR ANALYSES
save "`datapath'\version01\2-working\MI_survival", replace

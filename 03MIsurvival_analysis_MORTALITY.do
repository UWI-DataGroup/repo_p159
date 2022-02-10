** HEADER -----------------------------------------------------
**  DO-FILE METADATA
    //  algorithm name					MIsurvival.do
    //  project:				            BNR
    //  analysts:				       	Christina Howitt
    //  algorithm task			            MI mortality analysis using 10 years of BNR data


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
    log using "`logpath'\MIsurvival_mort", replace


** Load clean dataset
use "`datapath'\version02\1-input\heart_2009-2019_v8_anonymisedFUdata_Stata_v16_clean(16-Jul-2021)", clear


**-------------------------------------------------
** PREPARE DATA FOR ANALYSIS
**-------------------------------------------------

keep abstracted death anon_pid age sex dom deathdate fu1doa f1vstatus fu2doa f2vstatus hfu1date death

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

** Age group preparation - 10 year bands
gen age10 = recode(age,9,19,29,39,49,59,69,79,200)
recode age10 9=1 19=2 29=3 39=4 49=5 59=6 69=7 79=8 200=9
label define age10      1 "0-9"    2 "10-19"  3 "20-29"	///
						4 "30-39"  5 "40-49"  6 "50-59"	///
						7 "60-69"  8 "70-79"  9 "80 & over" , modify
label values age10 age10
label variable age10 "age in 10yr bands"

**check vital status for consistency
codebook f1vstatus // 94 cases coded as 99 "unknown" and 342 as "."
      recode f1vstatus 99=.
codebook f2vstatus // 87 cases coded as 99 "unknown" and 2416 as "."
      recode f2vstatus 99=.

** data checks
count if dom!=. & deathdate==. & f1vstatus==. & f2vstatus==. & fu1doa==. & fu2doa==. // 188 people in this category
      * the variable "death" has 28-day vital status
      codebook death
      codebook f1vstatus // death and f1vstatus are coded differently. 
      *create new variable for death that has similar coding to f1vstatus in order to check for differences
      gen f1vstatus2=.
      replace f1vstatus2=1 if death==2
      replace f1vstatus2=2 if death==1
      tab f1vstatus2, miss
      gen f1_diff=0
      replace f1_diff=1 if f1vstatus != f1vstatus2
      label values f1vstatus2 f1vstatus_lab
      tab f1_diff // 457 are different

      list anon_pid f1vstatus f1vstatus2 deathdate if f1_diff==1, abbrev(12) divider

      * where f1vstatus is missing, but there is a value for death(f1vstatus2), we will use the value for death
      replace f1vstatus = f1vstatus2 if f1vstatus==. & f1vstatus2!=.
      count if dom!=. & deathdate==. & f1vstatus==. & f2vstatus==. & fu1doa==. & fu2doa==. // now only 13 people

      * remaining differences:
      gen f1_diff2=0
      replace f1_diff2=1 if f1vstatus != f1vstatus2
      list anon_pid f1vstatus f1vstatus2 deathdate if f1_diff2==1, abbrev(12) divider
      * note that the remaining inconsistencies are taken care of in the generation of 28-day mortality variable. 
           
** Create a variable for survival time
gen time = .
replace time = deathdate - dom 
replace time = 365 if f2vstatus == 1 // participants who are alive at 1-yr follow up have a survival time of 365 days
replace time = . if f1vstatus==2 & deathdate==. // 5 cases 
replace time = . if f2vstatus==2 & deathdate==. // 0 cases 
* censor at 365 days
replace time=365 if time>365 & time <.


**Create survival status 
gen case = 1
gen died = 0
replace died = 1 if deathdate != .
replace died = 1 if f1vstatus==2
replace died = 1 if f2vstatus==2
label variable died "Survival status (1=dead)"


**checking to see how many have vital status recorded as dead but are missing death date
gen died2=0
replace died2=1 if f1vstatus==2
replace died2=1 if f2vstatus==2
gen died_diff=0
replace died_diff=1 if died!=died2
tab died_diff  // 32 cases are recorded as dead but are missing death date



**-------------------------------------------------------------
** PART ONE: SURVIVAL TIME BASIC SUMMARIES 
**-------------------------------------------------------------
** How many died on day of event?
gen mort0=.
replace mort0=1 if time==0
replace mort0=0 if time>0
label define mort0 0 "Alive on day of event" 1 "Died on day of event"
label values mort0 mort0
label variable mort0 "vital status on day of event"
tab mort0, miss
tab mort0 sex, col chi
tab mort0 year, col
** 28 day mortality
gen mort28=.
replace mort28=1 if f1vstatus==2 
replace mort28=1 if time<28
replace mort28=0 if time>=28 & time<.
label define mort28 0 "Alive at 28 days" 1 "Died by 28 days"
label values mort28 mort28
label variable mort28 "vital status at 28 days"
tab mort28
tab mort28 sex, col chi
** 1yr mortality
tab died 
tab died sex, col chi 

**trends in mortality over time
tab year mort0, row
tab year mort28, row
tab year died, row

tab mort28 sex, col chi
tab died sex, col chi 

** SAVE DATASET FOR FURTHER ANALYSES
save "`datapath'\version02\2-working\heart2010-2019_mort", replace

**-------------------------------------------------------------
** PART TWO: CRUDE AND STANDARDISED MORTALITY 
**-------------------------------------------------------------

** ANALYSIS NOTES: for directly standardised rates (incidence and mortality), the distrate command will be used. It requires
** a row for every age group, so after collapsing the dataset, sort by age group and then determine 
** whether any groups are missing. If so, fill in the blanks.

** Prepare data for SEX-STRATIFIED mortality calculations - must be aggregated by the same 10-year age groups as in the standardized population dataset
* We use the Barbados population in 2015 according to WPP to calculate mortality
merge m:1 age10 sex using "`datapath'\version02\1-input\bb2015_sex"
tab _merge
drop if _merge!=3 
drop _merge

numlabel, add mask ("#",)

      preserve
            collapse (sum) mort28 died (mean) bb_pop, by(age10 sex)

            distrate mort28 bb_pop using "`datapath'\version02\1-input\who2000_10-1", by(sex) stand(age10) popstand(pop) mult(100000) format(%8.2f)
/*  +----------------------------------------------------------------------------------------------------+
  |       sex   mort28        N    crude   rateadj   lb_gam   ub_gam   se_gam    srr   lb_srr   ub_srr |
  |----------------------------------------------------------------------------------------------------|
  | 1 ,Female     1204   147779   814.73    372.00   349.62   395.73    11.61   1.00        .        . |
  |   2 ,Male     1252   137548   910.23    577.62   545.07   611.81    16.86   1.55     1.43     1.69 |
  +----------------------------------------------------------------------------------------------------+
*/

            distrate died bb_pop using "`datapath'\version02\1-input\who2000_10-1", by(sex) stand(age10) popstand(pop) mult(100000) format(%8.2f)

/*+--------------------------------------------------------------------------------------------------+
  |       sex   died        N    crude   rateadj   lb_gam   ub_gam   se_gam    srr   lb_srr   ub_srr |
  |--------------------------------------------------------------------------------------------------|
  | 1 ,Female   1293   147779   874.96    400.57   377.38   425.11    12.02   1.00        .        . |
  |   2 ,Male   1323   137548   961.85    609.13   575.73   644.15    17.29   1.52     1.40     1.65 |
  +--------------------------------------------------------------------------------------------------+ */

      restore

      
***************************************************************************************************************************************************
**    NEXT WE CALCULATE STANDARDIZED INCIDENCE BY YEAR IN WOMEN THEN MEN
***************************************************************************************************************************************************
** WOMEN
**2010 
preserve
      keep if year==2010
      keep if sex==1
      collapse (sum) died mort28 (mean) bb_pop, by(age10)
      sort age10 
      list
      expand 4 in 1
      list

      replace age10=1 in 7
      replace died=0 in 7
      replace mort28=0 in 7
      replace bb_pop=16280 in 7

      replace age10=2 in 8
      replace died=0 in 8
      replace mort28=0 in 8
      replace bb_pop=18550 in 8

      replace age10=3 in 9
      replace died=0 in 9
      replace mort28=0 in 9
      replace bb_pop=18614 in 9

      sort age10
      list 

            distrate died bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)
            distrate mort28 bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)
restore

**2011
preserve
      keep if year==2011
      keep if sex==1
      collapse (sum) died mort28 (mean) bb_pop, by(age10)
      sort age10 
      list 
      expand 3 in 1
      list

      replace age10=2 in 8
      replace died=0 in 8
      replace mort28=0 in 8
      replace bb_pop=18550 in 8

      replace age10=3 in 9
      replace died=0 in 9
      replace mort28=0 in 9
      replace bb_pop=18614 in 9

      sort age10
      list 

      distrate died bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)
      distrate mort28 bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)
restore

**2012
preserve
      keep if year==2012
      keep if sex==1
      collapse (sum) died mort28 (mean) bb_pop, by(age10)
      sort age10 
      list 
      expand 3 in 1
      list

      replace age10=1 in 8
      replace died=0 in 8
      replace mort28=0 in 8
      replace bb_pop=16280 in 8

      replace age10=3 in 9
      replace died=0 in 9
      replace mort28=0 in 9
      replace bb_pop=18614 in 9

      sort age10
      list 

            distrate died bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)
            distrate mort28 bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)
restore

**2013
preserve
      keep if year==2013
      keep if sex==1
      collapse (sum) died mort28 (mean) bb_pop, by(age10)
      sort age10 
      list 

      expand 4 in 1
      list

      replace age10=1 in 7
      replace died=0 in 7
      replace mort28=0 in 7
      replace bb_pop=16280 in 7

      replace age10=2 in 8
      replace died=0 in 8
      replace mort28=0 in 8
      replace bb_pop=18550 in 8

      replace age10=3 in 9
      replace died=0 in 9
      replace mort28=0 in 9
      replace bb_pop=18614 in 9

      sort age10
      list 

      distrate died bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)
      distrate mort28 bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)
restore

**2014
preserve
      keep if year==2014
      keep if sex==1
      collapse (sum) died mort28 (mean) bb_pop, by(age10)
      sort age10 
      list 
      expand 4 in 1
      list

      replace age10=1 in 7
      replace died=0 in 7
      replace mort28=0 in 7
      replace bb_pop=16280 in 7

      replace age10=2 in 8
      replace died=0 in 8
      replace mort28=0 in 8
      replace bb_pop=18550 in 8

      replace age10=3 in 9
      replace died=0 in 9
      replace mort28=0 in 9
      replace bb_pop=18614 in 9

      sort age10
      list 

            distrate died bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)
            distrate mort28 bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)
restore


**2015
preserve
      keep if year==2015
      keep if sex==1
      collapse (sum) died mort28 (mean) bb_pop, by(age10)
      sort age10 
      list
      expand 4 in 1    
      
      replace age10=1 in 7
      replace died=0 in 7
      replace mort28=0 in 7
      replace bb_pop=16280 in 7

      replace age10=2 in 8
      replace died=0 in 8
      replace mort28=0 in 8
      replace bb_pop=18550 in 8

      replace age10=3 in 9
      replace died=0 in 9
      replace mort28=0 in 9
      replace bb_pop=18614 in 9

      sort age10
      list 

            distrate died bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)
            distrate mort28 bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)
restore

**2016
preserve
      keep if year==2016
      keep if sex==1
      collapse (sum) died mort28 (mean) bb_pop, by(age10)
      sort age10 
      list
      expand 4 in 1    
      
      replace age10=1 in 7
      replace died=0 in 7
      replace mort28=0 in 7
      replace bb_pop=16280 in 7

      replace age10=2 in 8
      replace died=0 in 8
      replace mort28=0 in 8
      replace bb_pop=18550 in 8

      replace age10=3 in 9
      replace died=0 in 9
      replace mort28=0 in 9
      replace bb_pop=18614 in 9

      sort age10
      list 

      distrate died bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)
      distrate mort28 bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)
restore

**2017
preserve
      keep if year==2017
      keep if sex==1
      collapse (sum) died mort28 (mean) bb_pop, by(age10)
      sort age10 
      list
      expand 3 in 1    
      
      replace age10=1 in 8
      replace died=0 in 8
      replace mort28=0 in 8
      replace bb_pop=16280 in 8

      replace age10=2 in 9
      replace died=0 in 9
      replace mort28=0 in 9
      replace bb_pop=18550 in 9

      sort age10
      list 

      distrate died bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)
      distrate mort28 bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)
restore

**2018
preserve
      keep if year==2018
      keep if sex==1
      collapse (sum) died mort28 (mean) bb_pop, by(age10)
      sort age10 
      list
      expand 4 in 1    
      
      replace age10=1 in 7
      replace died=0 in 7
      replace mort28=0 in 7
      replace bb_pop=16280 in 7

      replace age10=2 in 8
      replace died=0 in 8
      replace mort28=0 in 8
      replace bb_pop=18550 in 8

      replace age10=3 in 9
      replace died=0 in 9
      replace mort28=0 in 9
      replace bb_pop=18614 in 9

      sort age10
      list 
restore

**2019
preserve
      keep if year==2019
      keep if sex==1
      collapse (sum) died mort28 (mean) bb_pop, by(age10)
      sort age10 
      list
      expand 3 in 1

      replace age10=1 in 8
      replace died=0 in 8
      replace mort28=0 in 8
      replace bb_pop=16280 in 8

      replace age10=2 in 9
      replace died=0 in 9
      replace mort28=0 in 9
      replace bb_pop=18550 in 9

      sort age10
      list 

      distrate died bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)
      distrate mort28 bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)
restore

**MEN
**2010
preserve
      keep if year==2010
      keep if sex==2
      collapse (sum) died mort28 (mean) bb_pop, by(age10)
      sort age10 
      list      
      expand 2 in 1

      replace age10=2 in 9
      replace died=0 in 9
      replace mort28=0 in 9
      replace bb_pop=19413 in 9
      sort age10
      list 

            distrate died bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)
            distrate mort28 bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)
restore

**2011
preserve
      keep if year==2011
      keep if sex==2
      collapse (sum) died mort28 (mean) bb_pop, by(age10)
      sort age10 
      list
      expand 3 in 1 

      replace age10=2 in 8
      replace died=0 in 8
      replace mort28=0 in 8
      replace bb_pop=19413 in 8

      replace age10=3 in 9
      replace died=0 in 9
      replace mort28=0 in 9
      replace bb_pop=18578 in 9

      sort age10
      list 

            distrate died bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)
            distrate mort28 bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)
restore

**2012
preserve
      keep if year==2012
      keep if sex==2
      collapse (sum) died mort28 (mean) bb_pop, by(age10)
      sort age10 
      list
      expand 2 in 1

      replace age10=1 in 9
      replace died=0 in 9
      replace mort28=0 in 9
      replace bb_pop=16834 in 9

      sort age10
      list 

            distrate died bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)
            distrate mort28 bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)
restore

**2013
preserve
      keep if year==2013
      keep if sex==2
      collapse (sum) died mort28 (mean) bb_pop, by(age10)
      sort age10 
      list
      expand 3 in 1

      replace age10=1 in 8
      replace died=0 in 8
      replace mort28=0 in 8
      replace bb_pop=16834 in 8

      replace age10=2 in 9
      replace died=0 in 9
      replace mort28=0 in 9
      replace bb_pop=19413 in 9
      
      sort age10
      list 

            distrate died bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)
            distrate mort28 bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)
restore

**2014
preserve
      keep if year==2014
      keep if sex==2
      collapse (sum) died mort28 (mean) bb_pop, by(age10)
      sort age10 
      list
      expand 4 in 1    
      
      replace age10=1 in 7
      replace died=0 in 7
      replace mort28=0 in 7
      replace bb_pop=16834 in 7

      replace age10=2 in 8
      replace died=0 in 8
      replace mort28=0 in 8
      replace bb_pop=19413 in 8

      replace age10=3 in 9
      replace died=0 in 9
      replace mort28=0 in 9
      replace bb_pop=18578 in 9

      sort age10
      list 

            distrate died bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)
            distrate mort28 bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)
restore

*2015
preserve
      keep if year==2015
      keep if sex==2
      collapse (sum) died mort28 (mean) bb_pop, by(age10)
      sort age10 
      list
      expand 3 in 1

      replace age10=1 in 8
      replace died=0 in 8
      replace mort28=0 in 8
      replace bb_pop=16834 in 8

      replace age10=2 in 9
      replace died=0 in 9
      replace mort28=0 in 9
      replace bb_pop=19413 in 9
      
      sort age10
      list 

            distrate died bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)
            distrate mort28 bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)
restore

*2016
preserve
      keep if year==2016
      keep if sex==2
      collapse (sum) died mort28 (mean) bb_pop, by(age10)
      sort age10 
      list
      expand 4 in 1    
      
      replace age10=1 in 7
      replace died=0 in 7
      replace mort28=0 in 7
      replace bb_pop=16834 in 7

      replace age10=2 in 8
      replace died=0 in 8
      replace mort28=0 in 8
      replace bb_pop=19413 in 8

      replace age10=3 in 9
      replace died=0 in 9
      replace mort28=0 in 9
      replace bb_pop=18578 in 9

      sort age10
      list 

            distrate died bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)
            distrate mort28 bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)
restore

*2017
preserve
      keep if year==2017
      keep if sex==2
      collapse (sum) died mort28 (mean) bb_pop, by(age10)
      sort age10 
      list
      expand 4 in 1    
      
      replace age10=1 in 7
      replace died=0 in 7
      replace mort28=0 in 7
      replace bb_pop=16834 in 7

      replace age10=2 in 8
      replace died=0 in 8
      replace mort28=0 in 8
      replace bb_pop=19413 in 8

      replace age10=3 in 9
      replace died=0 in 9
      replace mort28=0 in 9
      replace bb_pop=18578 in 9

      sort age10
      list 

            distrate died bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)
            distrate mort28 bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)
restore

*2018
preserve
      keep if year==2018
      keep if sex==2
      collapse (sum) died mort28 (mean) bb_pop, by(age10)
      sort age10 
      list
      expand 2 in 1

      replace age10=1 in 9
      replace died=0 in 9
      replace mort28=0 in 9
      replace bb_pop=16834 in 9
      sort age10
      list 

            distrate died bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)
            distrate mort28 bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)
restore

*2019
preserve
      keep if year==2019
      keep if sex==2
      collapse (sum) died mort28 (mean) bb_pop, by(age10)
      sort age10 
      list
      expand 2 in 1 

      replace age10=2 in 9
      replace died=0 in 9
      replace mort28=0 in 9
      replace bb_pop=19413 in 9 
      sort age10
      list 

      distrate died bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)
      distrate mort28 bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)

restore

*The results from the above analyses were copied and pasted into Excel and saved in X:\The University of the West Indies\DataGroup - repo_data\data_p159\version02\2-working\Mortality.xlsx
*This file will be used to produce graphics

**-------------------------------------------------------------
**    GRAPHIC: TRENDS IN AGE-STANDARDIZED MORTALITY OVER TIME
**------------------------------------------------------------- 

import excel "`datapath'\version02\2-working\Mortality2.xlsx", sheet("1yr") firstrow clear

tempfile women_mort 

** Lowess Smoothing
preserve 
      drop if sex==2
      lowess rateadj year, gen(incLW) bwidth(0.5) nograph
      lowess lowCI year, gen(lLow) bwidth(0.5) nograph
      lowess upCI year, gen(upLow) bwidth(0.5) nograph
      save `women_mort', replace
restore

 
drop if sex==1
lowess rateadj year, gen(incLW) bwidth(0.5) nograph
lowess lowCI year, gen(lLow) bwidth(0.5) nograph
lowess upCI year, gen(upLow) bwidth(0.5) nograph

append using `women_mort'

#delimit ; 
                  graph twoway 
                              (rarea lLow upLow year if sex==1, col("254 224 210%40") lw(none))
                              (scatter rateadj year if sex==1, lp("l") mc("222 45 38") lc("222 45 38")) 
                              (line incLW year if sex==1, lc("222 45 38") lw(0.4) lp("-"))

                              (rarea lLow upLow year if sex==2, col("222 235 247%40") lw(none))
                              (scatter rateadj year if sex==2, lp("l") mc("49 130 189") lc("49 130 189")) 
                              (line incLW year if sex==2, lc("49 130 189") lw(0.4) lp("-"))
                              ,
                              /// Format x axis
                              xlab(2010 "2010" 2011 "2011" 2012 "2012" 2013 "2013" 2014 "2014" 2015 "2015" 2016 "2016" 2017 "2017" 2018 "2018" 2019 "2019", angle(45))
                              /// change title 
                              xtitle (Year)


                              /// Format y axis
                              ylab(0(10)100, labs(small) nogrid angle(0))
                              ymtick(0(5)100)
                              ytitle("Mortality rate per 100,000 population")            

                              /// format legend 
                              legend (off)

                              /// Graph region and plot region. 
                              graphregion (c(gs16))
                              ysize(3)
                              
                              name(mortality)
                              ;
            #delimit cr 
            graph export "`outputpath'/graphs.png", replace
   
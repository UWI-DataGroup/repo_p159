** HEADER -----------------------------------------------------
**  DO-FILE METADATA
    //  algorithm name					MIsurvival.do
    //  project:				            BNR
    //  analysts:				       	Christina Howitt
    //  algorithm task			            MI survival analysis using 10 years of BNR data


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
    log using "`logpath'\MIsurvival", replace


** Load clean dataset
use "`datapath'\version02\1-input\heart_2009-2019_v8_anonymisedFUdata_Stata_v16_clean(16-Jul-2021)", clear


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
sort time   

**Create survival status and censor at 365 days
gen case = 1
gen died = 0
replace died = 1 if deathdate != .
label variable died "Survival status (1=dead)"
replace time=365 if time>365

** Age group preparation - 10 year bands
gen age10 = recode(age,9,19,29,39,49,59,69,79,200)
recode age10 9=1 19=2 29=3 39=4 49=5 59=6 69=7 79=8 200=9
label define age10      1 "0-9"    2 "10-19"  3 "20-29"	///
						4 "30-39"  5 "40-49"  6 "50-59"	///
						7 "60-69"  8 "70-79"  9 "80 & over" , modify
label values age10 age10
label variable age10 "age in 10yr bands"


**-------------------------------------------------------------
** PART ONE: SURVIVAL TIME BASIC SUMMARIES 
**-------------------------------------------------------------



**-------------------------------------------------------------
** SURVIVAL ANALYSIS
**-------------------------------------------------------------

** Declare data to be survival-time data
gen time1 = time + 1
stset time1, failure(died) id (anon_pid)

** KAPLAN MEIER curve 1: all cases
stdescribe
sts graph

#delimit ;
		sts graph, by(sex)   
			plotregion(c(gs16) ic(gs16) ilw(thin) lw(thin)) 
			graphregion(color(gs16) ic(gs16) ilw(thin) lw(thin)) 

				ylab(0(0.2)1, labs(small) nogrid angle(0)) 
				ytick (0(0.1)1) 
				ytitle("Proportion surviving") 

				xtitle (Time from Acute MI (days))
				xtick (0(50)400) xmtick (0(25)400)

			legend (position(6) cols(1)
			lab(1 "Female")
			lab(2 "Male")

            )
 ;
#delimit cr 

*/

** KAPLAN MEIER curve 2: only abstracted cases
keep if abstracted == 1
stdescribe
sts graph

#delimit ;
		sts graph, by(sex)   
			plotregion(c(gs16) ic(gs16) ilw(thin) lw(thin)) 
			graphregion(color(gs16) ic(gs16) ilw(thin) lw(thin)) 

				ylab(0.6(0.1)1, labs(small) nogrid angle(0)) 
				ytick (0.6(0.05)1) 
				ytitle("Proportion surviving") 

				xtitle (Time from Acute MI (days))
				xtick (0(50)400) xmtick (0(25)400)

			legend (position(6) cols(1)
			lab(1 "Female")
			lab(2 "Male")

            )
 ;
#delimit cr 

** MEDIAN SURVIVAL TIME
stci 
stci, by(sex)

stcox i.sex
stcox i.sex age
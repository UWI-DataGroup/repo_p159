
cls
** HEADER -----------------------------------------------------
**  DO-FILE METADATA
    //  algorithm name          20_analysis cancer.do
    //  project:                BNR
    //  analysts:               Jacqueline CAMPBELL / Kern ROCKE
    //  date first created      02-DEC-2019
    // 	date last modified      20-JAN-2021
    //  algorithm task          Analyzing combined cancer dataset: (1) Numbers (2) ASIRs (3) Survival
    //  status                  Completed
    //  objective               To have one dataset with cleaned and grouped 2013, 2014 data for inclusion in 2015 cancer report.
    //  methods                 See 30_report cancer.do for detailed methods of each statistic

    ** General algorithm set-up
    version 16.0
    clear all
    macro drop _all
    set more off

    ** Initialising the STATA log and allow automatic page scrolling
    capture {
            program drop _all
    	drop _all
    	log close
    	}

    ** Set working directories: this is for DATASET and LOGFILE import and export
    ** DATASETS to encrypted SharePoint folder
    local datapath "X:/The University of the West Indies/DataGroup - repo_data/data_p117"
    ** LOGFILES to unencrypted OneDrive folder (.gitignore set to IGNORE log files on PUSH to GitHub)
    local logpath X:/OneDrive - The University of the West Indies/repo_datagroup/repo_p117

    ** Close any open log file and open a new log file
    capture log close
    log using "`logpath'/20_analysis cancer.smcl", replace
** HEADER -----------------------------------------------------



***************************************************************************
* SECTION 1: NUMBERS 
*        (1.1) total number & number of multiple events
*        (1.2) DCOs
*    	 (1.3) tumours by age-group: 
*				NOTE: missing/unknown age (code 999) are 
*				to be included in the age group that has a median total if 
*			  	total number of unk age is small, i.e. 5 cases with unk age; 
*			  	if larger then they would be distributed amongst more than
*			  	one age groups with median totals (NS update on 14-Oct-2020)
****************************************************************************
 
** LOAD cancer incidence dataset INCLUDING DCOs
use "`datapath'\version02\3-output\2013_2014_2015_cancer_nonsurvival" ,clear

** CASE variable
*drop case
gen case=1
label var case "cancer patient (tumour)"
 
*************************************************
** (1.1) Total number of events & multiple events
*************************************************
count //2744
tab dxyr ,m
/*
DiagnosisYe |
         ar |      Freq.     Percent        Cum.
------------+-----------------------------------
       2013 |        852       31.05       31.05
       2014 |        857       31.23       62.28
       2015 |      1,035       37.72      100.00
------------+-----------------------------------
      Total |      2,744      100.00
*/
tab patient dxyr ,m //2691 patients & 53 MPs; 2015: 1011 patients & 24 MPs (Checked this)
/*
                |          DiagnosisYear
cancer patient |      2013       2014       2015 |     Total
---------------+---------------------------------+----------
       patient |       840        840      1,011 |     2,691 
separate event |        12         17         24 |        53 
---------------+---------------------------------+----------
         Total |       852        857      1,035 |     2,744
*/

** JC updated AR's 2008 code for identifying MPs
tab ptrectot ,m
tab ptrectot patient ,m
tab ptrectot dxyr ,m

tab eidmp dxyr,m

duplicates list pid, nolabel sepby(pid) 
duplicates tag pid, gen(mppid_analysis)
sort pid cr5id
count if mppid_analysis>0 //86
//list pid topography morph ptrectot eidmp cr5id icd10 dxyr if mppid_analysis>0 ,sepby(pid)
 
** Of 2691 patients, 53 had >1 tumour

** note: remember to check in situ vs malignant from behaviour (beh)
tab beh ,m // 3908 malignant; 134 in-situ; 18 uncertain/benign
/*
  Behaviour |      Freq.     Percent        Cum.
------------+-----------------------------------
  Malignant |      2,744      100.00      100.00
------------+-----------------------------------
      Total |      2,744      100.00
*/

*************************************************
** (1.2) DCOs - patients identified only at death
*************************************************
tab basis beh ,m
/*
                      | Behaviour
     BasisOfDiagnosis | Malignant |     Total
----------------------+-----------+----------
                  DCO |       216 |       216 
        Clinical only |        99 |        99 
Clinical Invest./Ult  |       115 |       115 
Exploratory surg./aut |        20 |        20 
Lab test (biochem/imm |         9 |         9 
        Cytology/Haem |       102 |       102 
           Hx of mets |        44 |        44 
        Hx of primary |     1,997 |     1,997 
        Autopsy w/ Hx |        19 |        19 
              Unknown |       123 |       123 
----------------------+-----------+----------
                Total |     2,744 |     2,744  
*/

tab basis dxyr ,m
/*

                      |          DiagnosisYear
     BasisOfDiagnosis |      2013       2014       2015 |     Total
----------------------+---------------------------------+----------
                  DCO |        43         39        134 |       216 
        Clinical only |        19         38         42 |        99 
Clinical Invest./Ult  |        50         29         36 |       115 
Exploratory surg./aut |        10          5          5 |        20 
Lab test (biochem/imm |         3          3          3 |         9 
        Cytology/Haem |        31         43         28 |       102 
           Hx of mets |        13         13         18 |        44 
        Hx of primary |       634        623        740 |     1,997 
        Autopsy w/ Hx |         6          9          4 |        19 
              Unknown |        43         55         25 |       123 
----------------------+---------------------------------+----------
                Total |       852        857      1,035 |     2,744 
*/
/* JC 03mar20 checked to see if any duplicated observations occurred but no, seems like a legitimate new prostate case
preserve
drop if dxyr!=2015 & siteiarc!=39
sort pid cr5id
quietly by pid cr5id :  gen duppidcr5id = cond(_N==1,0,_n)
sort pid cr5id
count if duppidcr5id>0 //0
list pid cr5id deathid eidmp ptrectot primarysite duppidcr5id if duppidcr5id>0
restore
*/
tab basis dxyr if patient==1
/*
                      |          DiagnosisYear
     BasisOfDiagnosis |      2013       2014       2015 |     Total
----------------------+---------------------------------+----------
                  DCO |        43         36        132 |       211 
        Clinical only |        19         36         41 |        96 
Clinical Invest./Ult  |        49         28         35 |       112 
Exploratory surg./aut |        10          5          5 |        20 
Lab test (biochem/imm |         3          3          3 |         9 
        Cytology/Haem |        30         43         28 |       101 
           Hx of mets |        13         13         18 |        44 
        Hx of primary |       624        614        720 |     1,958 
        Autopsy w/ Hx |         6          9          4 |        19 
              Unknown |        43         53         25 |       121 
----------------------+---------------------------------+----------
                Total |       840        840      1,011 |     2,691 
*/

//This section assesses DCO % in relation to tumour, patient and behaviour totals
**********
** 2015 **
**********
** As a percentage of all events: 12.95%
cii proportions 1035 134

** As a percentage of all events with known basis: 13.27%
cii proportions 1010 134

** As a percentage of all patients: 13.06%
cii proportions 1011 132

tab basis beh if dxyr==2015 ,m
/*
                      | Behaviour
     BasisOfDiagnosis | Malignant |     Total
----------------------+-----------+----------
                  DCO |       134 |       134 
        Clinical only |        42 |        42 
Clinical Invest./Ult  |        36 |        36 
Exploratory surg./aut |         5 |         5 
Lab test (biochem/imm |         3 |         3 
        Cytology/Haem |        28 |        28 
           Hx of mets |        18 |        18 
        Hx of primary |       740 |       740 
        Autopsy w/ Hx |         4 |         4 
              Unknown |        25 |        25 
----------------------+-----------+----------
                Total |     1,035 |     1,035
*/
** Below no longer applicable as non-malignant dx were removed from ds (23-Oct-2020)
** As a percentage for all those which were non-malignant: 0%
//cii proportions 18 0
 
** As a percentage of all malignant tumours: 12.95%
//cii proportions 1035 134

**********
** 2014 **
**********
** As a percentage of all events: 4.46%
cii proportions 857 39

** As a percentage of all events with known basis: 4.86%
cii proportions 802 39
 
** As a percentage of all patients: 4.29%
cii proportions 840 36

tab basis beh if dxyr==2014 ,m
/*
                      | Behaviour
     BasisOfDiagnosis | Malignant |     Total
----------------------+-----------+----------
                  DCO |        39 |        39 
        Clinical only |        38 |        38 
Clinical Invest./Ult  |        29 |        29 
Exploratory surg./aut |         5 |         5 
Lab test (biochem/imm |         3 |         3 
        Cytology/Haem |        43 |        43 
           Hx of mets |        13 |        13 
        Hx of primary |       623 |       623 
        Autopsy w/ Hx |         9 |         9 
              Unknown |        55 |        55 
----------------------+-----------+----------
                Total |       857 |       857
*/
** Below no longer applicable as non-malignant dx were removed from ds (23-Oct-2020)
** As a percentage for all those which were non-malignant: 0%
//cii proportions 23 0
 
** As a percentage of all malignant tumours: 4.58%
//cii proportions 874 40

**********
** 2013 **
**********
** As a percentage of all events: 5.05%
cii proportions 852 43

** As a percentage of all events with known basis: 5.32%
cii proportions 809 43
 
** As a percentage of all patients: 5.12%
cii proportions 840 43

tab basis beh if dxyr==2013 ,m
/*
                      | Behaviour
     BasisOfDiagnosis | Malignant |     Total
----------------------+-----------+----------
                  DCO |        43 |        43 
        Clinical only |        19 |        19 
Clinical Invest./Ult  |        50 |        50 
Exploratory surg./aut |        10 |        10 
Lab test (biochem/imm |         3 |         3 
        Cytology/Haem |        31 |        31 
           Hx of mets |        13 |        13 
        Hx of primary |       634 |       634 
        Autopsy w/ Hx |         6 |         6 
              Unknown |        43 |        43 
----------------------+-----------+----------
                Total |       852 |       852
*/
** Below no longer applicable as non-malignant dx were removed from ds (23-Oct-2020)
** As a percentage for all those which were non-malignant: 0%
//cii proportions 9 0
 
** As a percentage of all malignant tumours: 4.92%
//cii proportions 874 43


*************************
** Number of cases by sex
*************************
tab sex ,m

tab sex patient,m

** Mean age by sex overall (where sex: male=1, female=2)... BY TUMOUR
ameans age
ameans age if sex==1
ameans age if sex==2

 
** Mean age by sex overall (where sex: male=1, female=2)... BY PATIENT
preserve
keep if patient==1 //15 obs deleted
ameans age
ameans age if sex==1
ameans age if sex==2
restore
 
***********************************
** 1.4 Number of cases by age-group
***********************************
** Age labelling
gen age5 = recode(age,4,9,14,19,24,29,34,39,44,49,54,59,64,69,74,79,84,200)

recode age5 4=1 9=2 14=3 19=4 24=5 29=6 34=7 39=8 44=9 49=10 54=11 59=12 64=13 /// 
                        69=14 74=15 79=16 84=17 200=18

label define age5_lab 1 "0-4" 	 2 "5-9"    3 "10-14" ///
					  4 "15-19"  5 "20-24"  6 "25-29" ///
					  7 "30-34"  8 "35-39"  9 "40-44" ///
					 10 "45-49" 11 "50-54" 12 "55-59" ///
					 13 "60-64" 14 "65-69" 15 "70-74" ///
					 16 "75-79" 17 "80-84" 18 "85 & over", modify
label values age5 age5_lab
gen age_10 = recode(age5,3,5,7,9,11,13,15,17,200)
recode age_10 3=1 5=2 7=3 9=4 11=5 13=6 15=7 17=8 200=9

label define age_10_lab 1 "0-14"   2 "15-24"  3 "25-34" ///
                        4 "35-44"  5 "45-54"  6 "55-64" ///
                        7 "65-74"  8 "75-84"  9 "85 & over" , modify

label values age_10 age_10_lab

sort sex age_10

tab age_10 ,m
*/
** Save this new dataset without population data
label data "2013-2015 BNR-Cancer analysed data - Numbers"
note: TS This dataset does NOT include population data 
save "`datapath'\version02\2-working\2013_2014_2015_cancer_numbers", replace


* *********************************************
* ANALYSIS: SECTION 3 - cancer sites
* Covering:
*  3.1  Classification of cancer by site
*  3.2 	ASIRs by site; overall, men and women
* *********************************************
** NOTE: bb popn and WHO popn data prepared by IH are ALL coded M=2 and F=1
** Above note by AR from 2008 dofile

** Load the dataset
use "`datapath'\version02\2-working\2013_2014_2015_cancer_numbers", clear

****************************************************************************** 2015 ****************************************************************************************
drop if dxyr!=2015 //1709 deleted

count //1035

**********************************************************************************
** ASIR and 95% CI for Table 1 using AR's site groupings - using WHO World popn **
**					WORLD POPULATION PROSPECTS (WPP): 2015						**
*drop pfu
gen pfu=1 // for % year if not whole year collected; not done for cancer        **
**********************************************************************************
** NSobers confirmed use of WPP populations
labelbook sex_lab

********************************************************************
* (2.4c) IR age-standardised to WHO world popn - ALL TUMOURS: 2015
********************************************************************
** Using WHO World Standard Population
//tab siteiarc ,m

*drop _merge
merge m:m sex age_10 using "`datapath'\version02\2-working\pop_wpp_2015-10"
/*
    Result                           # of obs.
    -----------------------------------------
    not matched                             0
    matched                             1,035  (_merge==3)
    -----------------------------------------
*/
** No unmatched records

** SF requested by email on 16-Oct-2020 age and sex specific rates for top 10 cancers
/*
What is age-specific incidence rate? 
Age-specific rates provide information on the incidence of a particular event in an age group relative to the total number of people at risk of that event in the same age group.

What is age-standardised incidence rate?
The age-standardized incidence rate is the summary rate that would have been observed, given the schedule of age-specific rates, in a population with the age composition of some reference population, often called the standard population.
*/
preserve
collapse (sum) case (mean) pop_wpp, by(pfu age_10 sex siteiarc)
gen incirate=case/pop_wpp*100000
drop if siteiarc!=39 & siteiarc!=29 & siteiarc!=13 & siteiarc!=33 ///
		& siteiarc!=14 & siteiarc!=21 & siteiarc!=53 & siteiarc!=11 ///
		& siteiarc!=18 & siteiarc!=55
//by sex,sort: tab age_10 incirate ,m
sort siteiarc age_10 sex
//list incirate age_10 sex
//list incirate age_10 sex if siteiarc==13

format incirate %04.2f
gen year=2015
rename siteiarc cancer_site
rename incirate age_specific_rate
drop pfu case pop_wpp
order year cancer_site sex age_10 age_specific_rate
save "`datapath'\version02\2-working\2015_top10_age+sex_rates" ,replace
restore

** Check for missing age as these would need to be added to the median group for that site when assessing ASIRs to prevent creating an outlier
count if age==.|age==999 //0

** Below saved in pathway: 
//X:\The University of the West Indies\DataGroup - repo_data\data_p117\version02\2-working\WPP_population by sex_2013.txt
tab pop_wpp age_10  if sex==1 //female
tab pop_wpp age_10  if sex==2 //male

** Next, IRs for invasive tumours only
*preserve
	drop if age_10==.
	drop if beh!=3 //0 deleted
	
	collapse (sum) case (mean) pop_wpp, by(pfu age_10 sex)
	sort age sex
	** No missing age groups
		
	** -distrate is a user written command.
	** type -search distrate,net- at the Stata prompt to find and install this command

sort age_10
total pop_wpp

distrate case pop_wpp using "`datapath'\version02\2-working\who2000_10-2", 	///	
		         stand(age_10) popstand(pop) mult(100000) format(%8.2f)
** THIS IS FOR ALL SITES (INVASIVE TUMOURS ONLY) - STD TO WHO WORLD POPN

/*
  +-------------------------------------------------------------+
  | case        N    crude   rateadj   lb_gam   ub_gam   se_gam |
  |-------------------------------------------------------------|
  | 1035   285327   362.74    231.95   217.36   247.33     7.57 |
  +-------------------------------------------------------------+
*/
** JC update: Save these results as a dataset for reporting
gen cancer_site=1
gen year=1
matrix list r(adj)
matrix number = r(NDeath)
matrix asir = r(adj)
matrix ci_lower = r(lb_G)
matrix ci_upper = r(ub_G)
svmat number
svmat asir
svmat ci_lower
svmat ci_upper

collapse cancer_site year number asir ci_lower ci_upper
rename number1 number
rename asir1 asir 
rename ci_lower1 ci_lower
rename ci_upper1 ci_upper
replace asir=round(asir,0.01)
replace ci_lower=round(ci_lower,0.01)
replace ci_upper=round(ci_upper,0.01)
gen percent=number/1035*100
replace percent=round(percent,0.01)

label define cancer_site_lab 1 "all" 2 "prostate" 3 "breast" 4 "colon" 5 "rectum" 6 "corpus uteri" 7 "stomach" ///
							 8 "lung" 9 "multiple myeloma" 10 "non-hodgkin lymphoma" 11 "pancreas" ,modify
label values cancer_site cancer_site_lab
label define year_lab 1 "2015" 2 "2014" 3 "2013" ,modify
label values year year_lab
order cancer_site number percent asir ci_lower ci_upper
sort cancer_site number
save "`datapath'\version02\2-working\ASIRs" ,replace
restore


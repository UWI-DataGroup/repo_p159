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


/* GRAPHIC: TRENDS IN MORTALITY OVER TIME
			preserve

				drop _all  
				input year mort measure
						2010	63.95	1
						2011	59.86	1
						2012	54.95	1
						2013	59.66	1
						2014	49.76	1
						2015	40.74	1
						2016	53.53	1
						2017	45.18	1
						2018	51.97	1
						2019	50.09	1
						2010	69.48	2
						2011	64.97	2
						2012	60.89	2
						2013	67.33	2
						2014	60.24	2
						2015	48.46	2
						2016	59.00   2
						2017	53.32	2
						2018	60.66	2
						2019	57.04	2
						2010	70.06	3
						2011	65.65	3
						2012	65.84	3
						2013	70.17	3
						2014	62.68	3
						2015	50.00	3
						2016	60.59	3
						2017	56.53	3
						2018	65.01	3
						2019	60.33	3
				end

				#delimit ; 

				graph twoway 
						(connected mort year if measure==1, lp("l") mc("252 174 145") lc("252 174 145")) /// Mortality (Day 0)
						(connected mort year if measure==2, lp("-") mc("222 45 38") lc("222 45 38")) /// Mortality (Day 28)
						(connected mort year if measure==3, lp("-") mc("165 15 21") lc("165 15 21")) /// Mortality (1 Year)
						,
						/// Format x axis
						xlab(2010 "2010" 2011 "2011" 2012 "2012" 2013 "2013" 2014 "2014" 2015 "2015" 2016 "2016" 2017 "2017" 2018 "2018" 2019 "2019", angle(45))
						/// change title 
						xtitle (Year)


						/// Format y axis
						/// change label orientation
						ylab(30(10)80, labs(small) nogrid angle(0))
						/// change title. Can't work out how to move the title so it's not squished up next to the label (tried to make it better by making label text small)
						ytitle("Mortality (%)")            

						/// format legend 
						/// change position and stack labels vertically. The colours of my lines are different - not sure why but as it doesn't make a difference, I will leave it
						legend (position(6) cols(1)
						lab(1 "Mortality (Day 0)")
						lab(2 "Mortality (Day 28)")
						lab(3 "Mortality (1 Year)") 
						)

						/// Graph region and plot region. 
						graphregion (c(gs16))
						ysize(3)
						


						name(Mortality)
						;
			#delimit cr 

			restore	

*/


**-------------------------------------------------------------
** CRUDE AND STANDARDIZED INCIDENCE 
**-------------------------------------------------------------
** Age group preparation - 10 year bands
gen age_10 = recode(age,9,19,29,39,49,59,69,79,200)
recode age_10 9=1 19=2 29=3 39=4 49=5 59=6 69=7 79=8 200=9
label define age_10      1 "0-9"    2 "10-19"  3 "20-29"	///
						4 "30-39"  5 "40-49"  6 "50-59"	///
						7 "60-69"  8 "70-79"  9 "80 & over" , modify
label values age_10 age_10
label variable age_10 "age in 10yr bands"

* We use the Barbados population in 2015 according to WPP to calculate incidence
merge m:1 age_10 using "`datapath'\version02\1-input\bb2015"
drop _merge
drop if age_10==.

** Prepare data for incidence calculations - must be aggregated by the same 10-year age groups as in the standardized population dataset
collapse (sum) case bb_pop, by(age_10)

      **-------------------------------------------------------------
      **    Standardized incidence using WHO 2000 - 2025 population
      **-------------------------------------------------------------     
      sort age_10

      distrate case bb_pop using "`datapath'\version02\1-input\who2000_10-2", stand(age_10) popstand(pop) mult(100000) format(%8.2f)

      **-------------------------------------------------------------
      **    Age-specific incidence
      **------------------------------------------------------------- 
            *incidence in each 10-year age band per 100,000 population
                  generate inc = (case/bb_pop) * (10^5) 
                  label variable inc "Crude Incidence"
                  * Standard Error
                        gen se = ( (case^(1/2)) / bb_pop) * (10^5)
                  * Lower 95% CI
                        gen lower = ( (0.5 * invchi2(2*case, (0.05/2))) / bb_pop ) * (10^5)
                  * Upper 95% CI
                        gen upper = ( (0.5 * invchi2(2*(case+1), (1-(0.05/2)))) / bb_pop) * (10^5)

                  * Display the results           
                        label var bb_pop "BB population"
                        label var case "Cases"
                        label var se "SE"
                        label var lower "95% lo"
                        label var upper "95% hi"
                  
                        foreach var in inc se lower upper {
                                    format `var' %8.2f
                                    }

                        list age_10 case bb_pop inc se lower upper , noobs table

            /*Incidence standardized to WHO 2000 to 2025 population
            *MANUAL CALCULATION FIRST
            *merge with WHO standardized dataset (with proportion for each age stratum, rather than population size)
            merge 1:1 age10 using "X:\The University of the West Indies\DataGroup - repo_data\data_p159\version02\1-input\who2000_proportions.dta"
            list
            drop _merge
            list
            *rename who_pop population 
            * Weight stratum-specific rates according to standard population distribution and then sum them to create overall standardized rate
            generate product = inc*who_pop
            egen adj_rate = sum(product)
            list


/*
generate crMort = died/bb_pop 
label variable crMort "Crude Mortality"
distrate deaths population using "X:\OneDrive - The University of the West Indies\Christina_work\BNR_survivalanalyses\1962.dta", stand(age_category) popstand(pop) mult(100000) format(%8.2f)


 *dstdize case pop age10, by(nation) using (X:\The University of the West Indies\DataGroup - repo_data\data_p159\version02\1-input\who2000.dta)

*ATTEMPT TO USE DSTDIZE
            use `ami_incidence', clear
            generate nation=1
            save `ami_incidence', replace

            import excel "`datapath'\version02\1-input\dummy_nation.xlsx", sheet("Sheet1") firstrow clear
            destring age10, replace
            append using "`ami_incidence'"
            rename bb_pop population 

*/




*/*-------------------------------------------------------------
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

**Abstracted only
keep if abstracted==1
tab mort28 sex, col chi
tab died sex, col chi 




/* Prepare data for standardized mortality calculations - must be aggregated by the same 10-year age groups as in the standardized population dataset
*preserve
      collapse (sum) mort28 died bb_pop, by(age10)

      **-------------------------------------------------------------
      **    Standardized incidence using WHO 2000 - 2025 population
      **-------------------------------------------------------------     
      sort age10

      distrate mort28 bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)
      distrate mort28 bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)



*/
*/

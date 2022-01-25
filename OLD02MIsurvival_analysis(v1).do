** HEADER -----------------------------------------------------
**  DO-FILE METADATA
    //  algorithm name					MIsurvival.do
    //  project:				        BNR
    //  analysts:				       	Christina Howitt
    //  algorithm task			        MI survival analysis using 10 years of BNR data


    ** General algorithm set-up
    version 16
    clear all
    macro drop _all
    set more 1
    set linesize 120

    ** Set working directories: this is for DATASET and LOGFILE import and export
    ** DATASETS to encrypted SharePoint folder
    local datapath "X:\The University of the West Indies\DataGroup - repo_data\data_p116"
     ** LOGFILES to unencrypted OneDrive folder (.gitignore set to IGNORE log files on PUSH to GitHub)
    local logpath X:/OneDrive - The University of the West Indies/repo_datagroup/repo_p116
    ** GRAPHS to project output folder
    local outputpath "X:\The University of the West Indies\DataGroup - DG_Projects\PROJECT_p116\MI_survival\"

    ** Close any open log file and open a new log file
    capture log close
    log using "`logpath'\MIsurvival_analysis", replace


**Load dataset
use "`datapath'\version01\2-working\MI_survival", clear

**Create surviival status and censor at 365 days
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
gen age10 = recode(age,9,19,29,39,49,59,69,79,200)
recode age10 9=1 19=2 29=3 39=4 49=5 59=6 69=7 79=8 200=9
label define age10      1 "0-9"    2 "10-19"  3 "20-29"	///
						4 "30-39"  5 "40-49"  6 "50-59"	///
						7 "60-69"  8 "70-79"  9 "80 & over" , modify
label values age10 age10

merge m:m age10 using "`datapath'\version01\1-input\who2000_10-1"  // one participant has a missing age, so is not assigned an age group
drop _merge
drop if age==.

rename pop pop_who 


** CRUDE INCIDENCE 
preserve
    drop if age10==.
	collapse (sum) case (mean) pop_who, by(age10 sex)
	collapse (sum) case pop_who, by(age10)

    gen asir = (case / pop_who) * (10^5)
	label var asir "Age-stratified Incidence Rate"

    * Standard Error
	gen se = ( (case^(1/2)) / pop_who) * (10^5)

	* Lower 95% CI
	gen lower = ( (0.5 * invchi2(2*case, (0.05/2))) / pop_who ) * (10^5)
	* Upper 95% CI
	gen upper = ( (0.5 * invchi2(2*(case+1), (1-(0.05/2)))) / pop_who) * (10^5)

	* Display the results
	label var pop_who "P-Y"
	label var case "Cases"
	label var se "SE"
	label var lower "95% lo"
	label var upper "95% hi"
	foreach var in asir se lower upper {
			format `var' %8.2f
			}
	list age10 case pop_who asir se lower upper , noobs table

restore

/** INCIDENCE AGE-STANDARDIZED TO WHO WORLD POPULATION
preserve
	drop if age10==.
	collapse (sum) case (mean) pop_who, by(age10)
	sort age10
	distrate case pop_who using "`datapath'\version01\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)

restore 
*/


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


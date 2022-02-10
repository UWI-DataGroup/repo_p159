** HEADER -----------------------------------------------------
**  DO-FILE METADATA
    //  algorithm name					MIsurvival.do
    //  project:				            BNR
    //  analysts:				       	Christina Howitt
    //  algorithm task			            MI Incidence using 10 years of BNR data


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
    log using "`logpath'\MIsurvival_inc", replace


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


/**-------------------------------------------------------------
** PART ONE: CRUDE AND STANDARDIZED INCIDENCE 
**-------------------------------------------------------------
** Age group preparation - 10 year bands
gen age10 = recode(age,9,19,29,39,49,59,69,79,200)
recode age10 9=1 19=2 29=3 39=4 49=5 59=6 69=7 79=8 200=9
label define age10      1 "0-9"    2 "10-19"  3 "20-29"	///
						4 "30-39"  5 "40-49"  6 "50-59"	///
						7 "60-69"  8 "70-79"  9 "80 & over" , modify
label values age10 age10
label variable age10 "age in 10yr bands"



** Prepare data for incidence calculations - must be aggregated by the same 10-year age groups as in the standardized population dataset
preserve
      * We use the Barbados population in 2015 according to WPP to calculate incidence
      merge m:1 age10 using "`datapath'\version02\1-input\bb2015"
      drop _merge
      drop if age10==.

      collapse (sum) case (mean) bb_pop, by(age10)

      **-------------------------------------------------------------
      **    Standardized incidence using WHO 2000 - 2025 population
      **-------------------------------------------------------------     
      sort age10

      distrate case bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)

         /*               
                  +--------------------------------------------------------------+
                  | case        N     crude   rateadj   lb_gam   ub_gam   se_gam |
                  |--------------------------------------------------------------|
                  | 4063   285327   1423.98    828.26   801.75   855.54    13.64 |
                  +--------------------------------------------------------------+     */


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

                        list age10 case bb_pop inc se lower upper , noobs table

restore
            
** ANALYSIS NOTES: for directly standardised rates (incidence and mortality), the distrate command will be used. It requires
** a row for every age group, so after collapsing the dataset, sort by age group and then determine 
** whether any groups are missing. If so, fill in the blanks.            
            
            
            ***********************************************************************
             ** Standardized incidence by year (using WHO 2000 - 2025 population)
                  **2010
                  preserve
                        keep if year==2010
                        collapse (sum) case (mean) bb_pop, by(age10)
                        sort age10 
                        list 
                        expand 2 in 1
                        list

                        replace age10=2 in 9
                        replace case=0 in 9
                        replace bb_pop=37963 in 9

                        sort age10 
                        list 
                        
                        distrate case bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)

                     /*    +-------------------------------------------------------------+
                        | case        N    crude   rateadj   lb_gam   ub_gam   se_gam |
                        |-------------------------------------------------------------|
                        |  344   285327   120.56     69.97    62.38    78.35     4.00 |
                        +-------------------------------------------------------------+ */

                  restore

                  **2011
                  preserve
                        keep if year==2011 
                        collapse (sum) case (mean) bb_pop, by(age10)
                        sort age10
                        list 
                        expand 3 in 1 
                        list
                        replace age10=2 in 8
                        replace case=0 in 8
                        replace bb_pop=37963 in 8

                        replace age10=3 in 9
                        replace case=0 in 9
                        replace bb_pop=37192 in 9

                        sort age10
                        distrate case bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)

                          /*   +-------------------------------------------------------------+
                              | case        N    crude   rateadj   lb_gam   ub_gam   se_gam |
                              |-------------------------------------------------------------|
                              |  294   285327   103.04     59.92    52.90    67.73     3.71 |
                              +-------------------------------------------------------------+ */
  
                  restore

                  **2012
                  preserve
                  keep if year==2012
                        collapse (sum) case (mean) bb_pop, by(age10)
                        sort age10
                        list 
                        expand 2 in 1 
                        list

                        replace age10=1 in 9
                        replace case=0 in 9
                        replace bb_pop=33114 in 9

                        sort age10
                        list
                        
                        distrate case bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)

                        /*      +-------------------------------------------------------------+
                              | case        N    crude   rateadj   lb_gam   ub_gam   se_gam |
                              |-------------------------------------------------------------|
                              |  404   285327   141.59     82.32    74.09    91.34     4.32 |
                              +-------------------------------------------------------------+ */    
                  restore

                  **2013
                  preserve
                        keep if year==2013
                        collapse (sum) case (mean) bb_pop, by(age10)
                        sort age10
                        list 
                        expand 3 in 1
                        list

                        replace age10=1 in 8
                        replace case=0 in 8
                        replace bb_pop=33114 in 8
                        list

                        replace age10=2 in 9
                        replace case=0 in 9
                        replace bb_pop=37963 in 9

                        sort age10
                        list
                        distrate case bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)

                          /*    +-------------------------------------------------------------+
                              | case        N    crude   rateadj   lb_gam   ub_gam   se_gam |
                              |-------------------------------------------------------------|
                              |  352   285327   123.37     71.67    63.99    80.13     4.04 |
                              +-------------------------------------------------------------+ */
                       
                  restore

                  **2014
                  preserve 
                        keep if year==2014
                        collapse (sum) case (mean) bb_pop, by(age10)
                        sort age10
                        list 
                        expand 4 in 1
                        list
                        replace age10=1 in 7
                        replace case=0 in 7
                        replace bb_pop=33114 in 7

                        replace age10=2 in 8
                        replace case=0 in 8
                        replace bb_pop=37963 in 8

                        replace age10=3 in 9 
                        replace case=0 in 9
                        replace bb_pop=37192 in 9

                        sort age10 
                        list 
                        distrate case bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)
                  
                          /*    +-------------------------------------------------------------+
                              | case        N    crude   rateadj   lb_gam   ub_gam   se_gam |
                              |-------------------------------------------------------------|
                              |  410   285327   143.69     81.28    73.22    90.11     4.24 |
                              +-------------------------------------------------------------+   */
                  restore

                  **2015
                   preserve
                        keep if year==2015
                        collapse (sum) case (mean) bb_pop, by(age10)
                        sort age10
                        list
                        expand 3 in 1
                        list

                        replace age10=1 in 8
                        replace case=0 in 8
                        replace bb_pop=33114 in 8
                        list 
                              
                        replace age10=2 in 9
                        replace case=0 in 9
                        replace bb_pop=37963 in 9     
                        list

                        sort age10 
                        list 
                        distrate case bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)

                        /*      +-------------------------------------------------------------+
                              | case        N    crude   rateadj   lb_gam   ub_gam   se_gam |
                              |-------------------------------------------------------------|
                              |  324   285327   113.55     68.03    60.51    76.35     3.97 |
                              +-------------------------------------------------------------+   */

                  restore

                  **2016
                  preserve
                        keep if year==2016
                        collapse (sum) case (mean) bb_pop, by(age10)
                        sort age10
                        list

                        expand 4 in 1
                        list
                        replace age10=1 in 7
                        replace case=0 in 7
                        replace bb_pop=33114 in 7

                        replace age10=2 in 8
                        replace case=0 in 8
                        replace bb_pop=37963 in 8

                        replace age10=3 in 9 
                        replace case=0 in 9
                        replace bb_pop=37192 in 9

                        sort age10 
                        list 
                        distrate case bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)

                        /*      +-------------------------------------------------------------+
                              | case        N    crude   rateadj   lb_gam   ub_gam   se_gam |
                              |-------------------------------------------------------------|
                              |  438   285327   153.51     90.59    81.97    99.99     4.52 |
                              +-------------------------------------------------------------+  */

                        
                  restore

                  **2017
                  preserve
                        keep if year==2017
                        collapse (sum) case (mean) bb_pop, by(age10)
                        sort age10
                        list
                        expand 3 in 1
                        list

                        replace age10=1 in 8
                        replace case=0 in 8
                        replace bb_pop=33114 in 8
                        list 
                              
                        replace age10=2 in 9
                        replace case=0 in 9
                        replace bb_pop=37963 in 9     
                        list

                        sort age10 
                        list 
                        distrate case bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)

                        /*      +-------------------------------------------------------------+
                              | case        N    crude   rateadj   lb_gam   ub_gam   se_gam |
                              |-------------------------------------------------------------|
                              |  467   285327   163.67     95.90    87.03   105.55     4.65 |
                              +-------------------------------------------------------------+   */
                  
                  restore

                  **2018
                  preserve
                        keep if year==2018
                              collapse (sum) case (mean) bb_pop, by(age10)
                              sort age10
                              list
                              expand 2 in 1

                              replace age10=1 in 9
                              replace case=0 in 9
                              replace bb_pop=33114 in 9
                              list 

                              sort age10 
                              list 
                              distrate case bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)
                                /*    +-------------------------------------------------------------+
                                    | case        N    crude   rateadj   lb_gam   ub_gam   se_gam |
                                    |-------------------------------------------------------------|
                                    |  483   285327   169.28     99.21    90.17   109.02     4.73 |
                                    +-------------------------------------------------------------+  */                              
                  restore

                  **2019
                  preserve
                        keep if year==2019
                        collapse (sum) case (mean) bb_pop, by(age10)
                        sort age10 
                        list 
                        expand 2 in 1

                        replace age10=2 in 9
                        replace case=0 in 9
                        replace bb_pop=37963 in 9 

                        list 
                        sort age10 
                        list 
                        distrate case bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)

                          /*    +-------------------------------------------------------------+
                              | case        N    crude   rateadj   lb_gam   ub_gam   se_gam |
                              |-------------------------------------------------------------|
                              |  547   285327   191.71    109.37    99.96   119.55     4.92 |
                              +-------------------------------------------------------------      */

                             
                  restore                

                        **-------------------------------------------------------------
                        **    GRAPHIC: TRENDS IN AGE-STANDARDIZED INCIDENCE OVER TIME
                        **------------------------------------------------------------- 
                       preserve

                        drop _all
                        input year	case	crude	rateadj lowCI upCI SE
                        2010	344	120.56 69.97 62.38 78.35 4
                        2011	294	103.04 59.92 52.9	67.73	3.71
                        2012	404	141.59 82.32 74.09 91.34 4.32
                        2013	352	123.37 71.67 63.99 80.13 4.04
                        2014	410	143.69 81.28 73.22 90.11 4.24
                        2015	324	113.55 68.03 60.51 76.35 3.97
                        2016	438	153.51 90.59 81.97 99.99 4.52
                        2017	467	163.67 95.9	87.03	105.55 4.65
                        2018	483	169.28 99.21 90.17 109.02 4.73
                        2019	547	191.71 109.37 99.96 119.55 4.92

                        end

                        ** Lowess Smoothing
                        lowess rateadj year, gen(incLW) bwidth(0.6) nograph
                        lowess lowCI year, gen(lLow) bwidth(0.6) nograph
                        lowess upCI year, gen(upLow) bwidth(0.6) nograph

                                    #delimit ; 
                                          graph twoway 
                                                      (rarea lLow upLow year, col("181 215 244%40") lw(none))
                                                      (scatter rateadj year, lp("l") mc("71 129 179") lc("71 129 179")) 
                                                      (line incLW year, lc("71 129 179") lw(0.4) lp("-"))

                                                      ,
                                                      /// Format x axis
                                                      xlab(2010 "2010" 2011 "2011" 2012 "2012" 2013 "2013" 2014 "2014" 2015 "2015" 2016 "2016" 2017 "2017" 2018 "2018" 2019 "2019", angle(45))
                                                      /// change title 
                                                      xtitle (Year)


                                                      /// Format y axis
                                                      /// change label orientation
                                                      ylab(40(10)120, labs(small) nogrid angle(0))
                                                      ymtick(40(5)120)
                                                      ytitle("Incidence rate per 100,000 population")            

                                                      /// format legend 
                                                      legend (off)

                                                      /// Graph region and plot region. 
                                                      graphregion (c(gs16))
                                                      ysize(3)
                                                      
                                                      name(incidence)
                                                      ;
                                    #delimit cr 
                                    graph export "`outputpath'/graphs.png", replace
                        restore   




** Prepare data for SEX-STRATIFIED incidence calculations - must be aggregated by the same 10-year age groups as in the standardized population dataset
* We use the Barbados population in 2015 according to WPP to calculate incidence
merge m:1 age10 sex using "`datapath'\version02\1-input\bb2015_sex"
tab _merge
list anon_pid age sex if _merge==1  // 1 participant is missing age and another sex. These need to be dropped.
drop if age==.
drop if sex==.
drop _merge

preserve
      collapse (sum) case (mean) bb_pop, by(age10 sex)

      **-------------------------------------------------------------
      **    Standardized incidence using WHO 2000 - 2025 population
      **-------------------------------------------------------------     
      sort age10

      distrate case bb_pop using "`datapath'\version02\1-input\who2000_10-1", by(sex) stand(age10) popstand(pop) mult(100000) format(%8.2f)

    /*    +------------------------------------------------------------------------------------------------+
  |    sex   case        N     crude   rateadj   lb_gam   ub_gam   se_gam    srr   lb_srr   ub_srr |
  |------------------------------------------------------------------------------------------------|
  | Female   1903   137548   1383.52    822.77   785.06   862.06    19.48   1.00        .        . |
  |   Male   2159   147779   1460.97    860.10   822.18   899.52    19.58   1.05     0.98     1.12 |
  +------------------------------------------------------------------------------------------------+
  */

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

                        list age10 sex case bb_pop inc se lower upper , noobs table

restore

***********************************************************************
** Standardized incidence by year (using WHO 2000 - 2025 population) IN WOMEN 
***********************************************************************
numlabel, add mask ("#",)
**2010 
preserve
      keep if year==2010
      keep if sex==1
      collapse (sum) case (mean) bb_pop, by(age10)
      sort age10 
      list 
      expand 4 in 1
      list

      replace age10=1 in 7
      replace case=0 in 7
      replace bb_pop=16280 in 7

      replace age10=2 in 8
      replace case=0 in 8
      replace bb_pop=18550 in 8  

      replace age10=3 in 9
      replace case=0 in 9
      replace bb_pop=18614 in 9                  

      sort age10 
      list 

      distrate case bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)

      /*+-------------------------------------------------------------+
      | case        N    crude   rateadj   lb_gam   ub_gam   se_gam |
      |-------------------------------------------------------------|
      |  179   147779   121.13     60.00    50.90    70.57     4.88 |
      +-------------------------------------------------------------+ */


restore

**2011 
preserve
      keep if year==2011
      keep if sex==1
      collapse (sum) case (mean) bb_pop, by(age10)
      sort age10 
      list 
      expand 3 in 1
      list

      replace age10=2 in 8
      replace case=0 in 8
      replace bb_pop=18550 in 8

      replace age10=3 in 9
      replace case=0 in 9
      replace bb_pop=18614 in 9

      sort age10
      list
                              
      distrate case bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)

      /*    +------------------------------------------------------------+
      | case        N   crude   rateadj   lb_gam   ub_gam   se_gam |
      |------------------------------------------------------------|
      |  126   147779   85.26     41.55    34.02    50.57     4.09 |
      +------------------------------------------------------------+   */

restore

**2012 
preserve
      keep if year==2012
      keep if sex==1
      collapse (sum) case (mean) bb_pop, by(age10)
      sort age10 
      list 
      expand 3 in 1
      list 

      replace age10=1 in 8
      replace case=0 in 8
      replace bb_pop=16280 in 8

      replace age10=3 in 9
      replace case=0 in 9
      replace bb_pop=18614 in 9 

      sort age10
      list 
      distrate case bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)

      /*  +-------------------------------------------------------------+
            | case        N    crude   rateadj   lb_gam   ub_gam   se_gam |
            |-------------------------------------------------------------|
            |  201   147779   136.01     68.73    58.92    80.00     5.24 |
            +-------------------------------------------------------------+ */


restore 


      **2013
preserve
      keep if year==2013
      keep if sex==1
      collapse (sum) case (mean) bb_pop, by(age10)
      sort age10 
      list 
      expand 4 in 1
      list 

      replace age10=1 in 7
      replace case=0 in 7
      replace bb_pop=16280 in 7

      replace age10=2 in 8
      replace case=0 in 8
      replace bb_pop=18550 in 8

      replace age10=3 in 9
      replace case=0 in 9
      replace bb_pop=18614 in 9 

      sort age10
      list 
      distrate case bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)

      /*  +-------------------------------------------------------------+
      | case        N    crude   rateadj   lb_gam   ub_gam   se_gam |
      |-------------------------------------------------------------|
      |  161   147779   108.95     55.35    46.53    65.66     4.74 |
      +-------------------------------------------------------------+   */

      restore 


*2014
preserve
      keep if year==2014
      keep if sex==1
      collapse (sum) case (mean) bb_pop, by(age10)
      sort age10 
      list 
      expand 4 in 1
      list 

      replace age10=1 in 7
      replace case=0 in 7
      replace bb_pop=16280 in 7

      replace age10=2 in 8
      replace case=0 in 8
      replace bb_pop=18550 in 8

      replace age10=3 in 9
      replace case=0 in 9
      replace bb_pop=18614 in 9 

      sort age10
      list 
      distrate case bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)
      
            /* +-------------------------------------------------------------+
            | case        N    crude   rateadj   lb_gam   ub_gam   se_gam |
            |-------------------------------------------------------------|
            |  190   147779   128.57     62.70    53.47    73.38     4.94 |
            +-------------------------------------------------------------+   */
restore


      *2015
preserve
      keep if year==2015
      keep if sex==1
      collapse (sum) case (mean) bb_pop, by(age10)
      sort age10 
      list 
      expand 4 in 1
      list 

      replace age10=1 in 7
      replace case=0 in 7
      replace bb_pop=16280 in 7

      replace age10=2 in 8
      replace case=0 in 8
      replace bb_pop=18550 in 8

      replace age10=3 in 9
      replace case=0 in 9
      replace bb_pop=18614 in 9 

      sort age10
      list 
      distrate case bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)

            /* +-------------------------------------------------------------+
            | case        N    crude   rateadj   lb_gam   ub_gam   se_gam |
            |-------------------------------------------------------------|
            |  155   147779   104.89     53.53    44.95    63.60     4.62 |
            +-------------------------------------------------------------+ */
      
restore


*2016
preserve
      keep if year==2016
      keep if sex==1
      collapse (sum) case (mean) bb_pop, by(age10)
      sort age10 
      list 
      expand 4 in 1
      list 

      replace age10=1 in 7
      replace case=0 in 7
      replace bb_pop=16280 in 7

      replace age10=2 in 8
      replace case=0 in 8
      replace bb_pop=18550 in 8

      replace age10=3 in 9
      replace case=0 in 9
      replace bb_pop=18614 in 9 

      sort age10
      list 
      distrate case bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)

            /* +-------------------------------------------------------------+
            | case        N    crude   rateadj   lb_gam   ub_gam   se_gam |
            |-------------------------------------------------------------|
            |  188   147779   127.22     62.86    53.61    73.57     4.96 |
            +-------------------------------------------------------------+   */



restore

      *2017
preserve
      keep if year==2017
      keep if sex==1
      collapse (sum) case (mean) bb_pop, by(age10)
      sort age10 
      list 
      expand 3 in 1
      list

      replace age10=1 in 8
      replace case=0 in 8
      replace bb_pop=16280 in 8

      replace age10=2 in 9
      replace case=0 in 9
      replace bb_pop=18550 in 9

      sort age10
      list
                              
      distrate case bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)

      /*  +-------------------------------------------------------------+
            | case        N    crude   rateadj   lb_gam   ub_gam   se_gam |
            |-------------------------------------------------------------|
            |  214   147779   144.81     74.45    64.14    86.25     5.50 |
            +-------------------------------------------------------------+ */


restore 

      *2018
preserve
      keep if year==2018
      keep if sex==1
      collapse (sum) case (mean) bb_pop, by(age10)
      sort age10 
      list 
      expand 4 in 1
      list 

      replace age10=1 in 7
      replace case=0 in 7
      replace bb_pop=16280 in 7

      replace age10=2 in 8
      replace case=0 in 8
      replace bb_pop=18550 in 8

      replace age10=3 in 9
      replace case=0 in 9
      replace bb_pop=18614 in 9 

      sort age10
      list 
      distrate case bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)

      /*  +-------------------------------------------------------------+
      | case        N    crude   rateadj   lb_gam   ub_gam   se_gam |
      |-------------------------------------------------------------|
      |  213   147779   144.13     72.22    62.25    83.66     5.32 |
      +-------------------------------------------------------------+   */

restore 


*2019
preserve
      keep if year==2019
      keep if sex==1
      collapse (sum) case (mean) bb_pop, by(age10)
      sort age10 
      list 
      expand 3 in 1
      list

      replace age10=1 in 8
      replace case=0 in 8
      replace bb_pop=16280 in 8

      replace age10=2 in 9
      replace case=0 in 9
      replace bb_pop=18550 in 9

      sort age10
      list
                              
      distrate case bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)

            /* +-------------------------------------------------------------+
            | case        N    crude   rateadj   lb_gam   ub_gam   se_gam |
            |-------------------------------------------------------------|
            |  276   147779   186.77     88.90    78.03   101.19     5.77 |
            +-------------------------------------------------------------+   */
restore

                        
***********************************************************************
** Standardized incidence by year (using WHO 2000 - 2025 population) IN MEN 
***********************************************************************

*2010
preserve
      keep if year==2010
      keep if sex==2
      collapse (sum) case (mean) bb_pop, by(age10)
      sort age10 
      list 
      expand 2 in 1
      list

      replace age10=2 in 9
      replace case=0 in 9
      replace bb_pop=19413 in 9

      sort age10
      list
      distrate case bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)

  /*      +-------------------------------------------------------------+
            | case        N    crude   rateadj   lb_gam   ub_gam   se_gam |
            |-------------------------------------------------------------|
            |  165   137548   119.96     80.59    68.39    94.52     6.52 |
            +-------------------------------------------------------------+   */
restore

*2011
preserve
      keep if year==2011
      keep if sex==2
      collapse (sum) case (mean) bb_pop, by(age10)
      sort age10 
      list 
      expand 3 in 1
      list

      replace age10=2 in 8
      replace case=0 in 8
      replace bb_pop=19413 in 8

      replace age10=3 in 9
      replace case=0 in 9
      replace bb_pop=18578 in 9

      sort age10
      list
      
      distrate case bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)

       /*   +-------------------------------------------------------------+
            | case        N    crude   rateadj   lb_gam   ub_gam   se_gam |
            |-------------------------------------------------------------|
            |  168   137548   122.14     81.34    69.17    95.24     6.50 |
            +-------------------------------------------------------------+   */

restore


*2012
preserve
      keep if year==2012
      keep if sex==2
      collapse (sum) case (mean) bb_pop, by(age10)
      sort age10 
      list 
      expand 2 in 1
      list

      replace age10=1 in 9
      replace case=0 in 9
      replace bb_pop=16834 in 9

      sort age10
      list
     
      distrate case bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)

  /*      +-------------------------------------------------------------+
  | case        N    crude   rateadj   lb_gam   ub_gam   se_gam |
  |-------------------------------------------------------------|
  |  203   137548   147.58     97.45    84.12   112.49     7.09 |
  +-------------------------------------------------------------+ */


restore

*2013
preserve
      keep if year==2013
      keep if sex==2
      collapse (sum) case (mean) bb_pop, by(age10)
      sort age10 
      list 
      expand 3 in 1
      list

      replace age10=1 in 8
      replace case=0 in 8
      replace bb_pop=16834 in 8

      replace age10=2 in 9
      replace case=0 in 9
      replace bb_pop=19413 in 9

      sort age10
      list 

      distrate case bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)

 /*  +-------------------------------------------------------------+
  | case        N    crude   rateadj   lb_gam   ub_gam   se_gam |
  |-------------------------------------------------------------|
  |  191   137548   138.86     90.99    78.16   105.52     6.83 |
  +-------------------------------------------------------------+
*/

restore


*2014
preserve
      keep if year==2014
      keep if sex==2
      collapse (sum) case (mean) bb_pop, by(age10)
      sort age10 
      list 
      expand 4 in 1
      list

      replace age10=1 in 7
      replace case=0 in 7
      replace bb_pop=16834 in 7

      replace age10=2 in 8
      replace case=0 in 8
      replace bb_pop=19413 in 8

      replace age10=3 in 9
      replace case=0 in 9
      replace bb_pop=18578 in 9

      sort age10
      list 

      distrate case bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)
/*  +-------------------------------------------------------------+
  | case        N    crude   rateadj   lb_gam   ub_gam   se_gam |
  |-------------------------------------------------------------|
  |  220   137548   159.94    102.98    89.46   118.19     7.18 |
  +-------------------------------------------------------------+
*/
restore


*2015
preserve
      keep if year==2015
      keep if sex==2
      collapse (sum) case (mean) bb_pop, by(age10)
      sort age10 
      list 
      expand 3 in 1
      list

      replace age10=1 in 8
      replace case=0 in 8
      replace bb_pop=16834 in 8

      replace age10=2 in 9
      replace case=0 in 9
      replace bb_pop=19413 in 9

      sort age10
      list 

      distrate case bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)
/*
  +-------------------------------------------------------------+
  | case        N    crude   rateadj   lb_gam   ub_gam   se_gam |
  |-------------------------------------------------------------|
  |  169   137548   122.87     84.12    71.62    98.36     6.67 |
  +-------------------------------------------------------------+  */


restore

*2016
preserve
      keep if year==2016
      keep if sex==2
      collapse (sum) case (mean) bb_pop, by(age10)
      sort age10 
      list 
      expand 4 in 1
      list

      replace age10=1 in 7
      replace case=0 in 7
      replace bb_pop=16834 in 7

      replace age10=2 in 8
      replace case=0 in 8
      replace bb_pop=19413 in 8

      replace age10=3 in 9
      replace case=0 in 9
      replace bb_pop=18578 in 9

      sort age10
      list 

      distrate case bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)

/* +-------------------------------------------------------------+
  | case        N    crude   rateadj   lb_gam   ub_gam   se_gam |
  |-------------------------------------------------------------|
  |  250   137548   181.75    122.49   107.49   139.19     7.93 |
  +-------------------------------------------------------------+   */

restore


*2017
preserve
      keep if year==2017
      keep if sex==2
      collapse (sum) case (mean) bb_pop, by(age10)
      sort age10 
      list 
      expand 4 in 1

      replace age10=1 in 7
      replace case=0 in 7
      replace bb_pop=16834 in 7

      replace age10=2 in 8
      replace case=0 in 8
      replace bb_pop=19413 in 8

      replace age10=3 in 9
      replace case=0 in 9
      replace bb_pop=18578 in 9

      sort age10
      list 

      distrate case bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)

restore

*2018
preserve
      keep if year==2018
      keep if sex==2
      collapse (sum) case (mean) bb_pop, by(age10)
      sort age10 
      list 
      expand 2 in 1

      replace age10=1 in 9
      replace case=0 in 9
      replace bb_pop=16834 in 9

      sort age10
      list 

      distrate case bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)
/*  +-------------------------------------------------------------+
  | case        N    crude   rateadj   lb_gam   ub_gam   se_gam |
  |-------------------------------------------------------------|
  |  269   137548   195.57    129.98   114.54   147.12     8.16 |
  +-------------------------------------------------------------+
*/

restore

*2019
preserve
      keep if year==2019
      keep if sex==2
      collapse (sum) case (mean) bb_pop, by(age10)
      sort age10 
      list 
      expand 2 in 1

      replace age10=2 in 9
      replace case=0 in 9
      replace bb_pop=19413 in 9

      sort age10
      list 

      distrate case bb_pop using "`datapath'\version02\1-input\who2000_10-1", stand(age10) popstand(pop) mult(100000) format(%8.2f)

  /* +-------------------------------------------------------------+
  | case        N    crude   rateadj   lb_gam   ub_gam   se_gam |
  |-------------------------------------------------------------|
  |  271   137548   197.02    131.27   115.71   148.52     8.22 |
  +-------------------------------------------------------------+
*/

restore

*The results from the above analyses were copied and pasted into Excel and saved in X:\The University of the West Indies\DataGroup - repo_data\data_p159\version02\2-working\MIincidence.xlsx
*This file will be used to produce graphics

**-------------------------------------------------------------
**    GRAPHIC: TRENDS IN AGE-STANDARDIZED INCIDENCE OVER TIME
**------------------------------------------------------------- 

import excel "`datapath'\version02\2-working\MIincidence.xlsx", sheet("Sheet1") firstrow clear

tempfile women_inc 

** Lowess Smoothing
preserve 
      drop if sex==2
      lowess rateadj year, gen(incLW) bwidth(0.5) nograph
      lowess lowCI year, gen(lLow) bwidth(0.5) nograph
      lowess upCI year, gen(upLow) bwidth(0.5) nograph
      save `women_inc', replace
restore

 
drop if sex==1
lowess rateadj year, gen(incLW) bwidth(0.5) nograph
lowess lowCI year, gen(lLow) bwidth(0.5) nograph
lowess upCI year, gen(upLow) bwidth(0.5) nograph

append using `women_inc'



            #delimit ; 
                  graph twoway 
                              (rarea lLow upLow year if sex==1, col("254 224 210%60") lw(none))
                              (scatter rateadj year if sex==1, lp("l") mc("222 45 38") lc("222 45 38")) 
                              (line incLW year if sex==1, lc("222 45 38") lw(0.4) lp("-"))

                              (rarea lLow upLow year if sex==2, col("222 235 247%60") lw(none))
                              (scatter rateadj year if sex==2, lp("l") mc("49 130 189") lc("49 130 189")) 
                              (line incLW year if sex==2, lc("49 130 189") lw(0.4) lp("-"))
                              ,
                              /// Format x axis
                              xlab(2010 "2010" 2011 "2011" 2012 "2012" 2013 "2013" 2014 "2014" 2015 "2015" 2016 "2016" 2017 "2017" 2018 "2018" 2019 "2019", angle(45))
                              /// change title 
                              xtitle (Year)


                              /// Format y axis
                              /// change label orientation
                              ylab(20(20)160, labs(small) nogrid angle(0))
                              ymtick(40(10)160)
                              ytitle("Incidence rate per 100,000 population")            

                              /// format legend 
                              legend (off)

                              /// Graph region and plot region. 
                              graphregion (c(gs16))
                              ysize(3)
                              
                              name(incidence)
                              ;
            #delimit cr 
            graph export "`outputpath'/graphs.png", replace
*restore   

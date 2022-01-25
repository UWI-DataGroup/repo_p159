** HEADER -----------------------------------------------------
**  DO-FILE METADATA
    //  algorithm name					MIsurvival.do
    //  project:				        BNR
    //  analysts:				       	Christina Howitt
    //  algorithm task			        Prepare data required for MI survival analysis using 10 years of BNR data


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
** PART ONE: BARBADOS POPULATION IN 10-YR AGE GROUPS 
**-------------------------------------------------------------
** -----------------------------------------------------------
** In order to calculate crude incidence and mortality, we need the population size in Barbados in 10 year age bands. We used the WPP 2015 population estimates, downloaded from 
** https://population.un.org/wpp/Download/Standard/Population/
** The full population data file was downloaded and saved as"WPP2019_POP_F07_1_POPULATION_BY_AGE_BOTH_SEXES.xlsx" in "`datapath'\version02\1-input". The Barbados 2015 row was copied and 
** pasted into "WPP_BB_2015.xlsx" and the data was converted from 5 year bands into the 10 year bands entered below
** -----------------------------------------------------------
drop _all
input age10 bb_pop
1	33114
2	37963
3	37192
4	38155
5	40612
6	40457
7	28782
8	17047
9	12005
end

			label define age10_lab  1 "0-9"	   2 "10-19"  3 "20-29"	///
									4 "30-39"  5 "40-49"  6 "50-59"	///
									7 "60-69"  8 "70-79"  9 "80 & over" , modify
			label values age10 age10_lab

save "`datapath'\version02\1-input\bb2015", replace


drop _all 
input age10 sex bb_pop
1	1	16280
2	1	18550
3	1	18614
4	1	19485
5	1	21079
6	1	21543
7	1	15349
8	1	9616
9	1	7263
1	2	16834
2	2	19413
3	2	18578
4	2	18670
5	2	19533
6	2	18914
7	2	13433
8	2	7431
9	2	4742
end 

label values age10 age10_lab
label define sex 1 "female" 2 "male"
label values sex sex

save "`datapath'\version02\1-input\bb2015_sex", replace


**-------------------------------------------------------------
** PART TWO: STANDARD POPULATION FOR INCIDENCE 
**-------------------------------------------------------------
** -----------------------------------------------------------
** We use the World (WHO 2000-2025) Standard population (https://www.who.int/healthinfo/paper31.pdf) in order to standardize incidence and mortality estimates (saved in "`datapath'\
** version02\1-input" as "WHO2000.xlsx"). 
** -----------------------------------------------------------

drop _all
		input age5 pop
		1	8860
		2	8690
		3	8600
		4	8470
		5	8220
		6	7930
		7	7610
		8	7150
		9	6590
		10	6040
		11	5370
		12	4550
		13	3720
		14	2960
		15	2210
		16	1520
		17	910
		18	635
		end

** Age labelling
label define whoage5_lab  	1 "0-4"    2 "5-9"	  3 "10-14"	///
		            4 "15-19"  5 "20-24"  6 "25-29"	///
					7 "30-34"  8 "35-39"  9 "40-44"	///
		           10 "45-49" 11 "50-54" 12 "55-59"	///
				   13 "60-64" 14 "65-69" 15 "70-74"	///
		           16 "75-79" 17 "80-84" 18 "85 & over", modify
label values age5 whoage5_lab
label var age5 "WHO standard 5-year age-grouping (18 groups)"


** TEN age groups in TEN-year bands. 
** This is the standard for all standard population distributions
gen age10 = recode(age5,2,4,6,8,10,12,14,16,17)

recode age10 2=1 4=2 6=3 8=4 10=5 12=6 14=7 16=8 17=9
label values age10 age10_lab


collapse (sum) pop , by(age10)
	label data "WHO world standard million: 10-year age bands1"
	sort age10
save "`datapath'\version02\1-input\who2000_10-1", replace

/*drop _all
input age_10 who_pop
1	0.1755
2	0.1707
3	0.1615
4	0.1476
5	0.1263
6	0.0992
7	0.0668
8	0.0373
9	0.0135
end 

save "`datapath'\version02\1-input\who2000", replace
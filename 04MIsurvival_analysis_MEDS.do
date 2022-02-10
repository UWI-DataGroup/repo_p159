** HEADER -----------------------------------------------------
**  DO-FILE METADATA
    //  algorithm name					MIsurvival.do
    //  project:				         BNR
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
    log using "`logpath'\MIsurvival_meds", replace


** Load clean dataset
use "`datapath'\version02\1-input\heart_2009-2019_v8_anonymisedFUdata_Stata_v16_clean(16-Jul-2021)", clear


**-------------------------------------------------
** PREPARE DATA FOR ANALYSIS
**-------------------------------------------------

** Create variable for year of MI
drop year // this is an existing variable in the datset but there are a couple of differences between it and the derived version below, which we will consider more accurate
gen year = year(dom)
label variable year "year of MI event"
order year, after(dom)

** we are only using full years, so 2009 to be dropped
drop if year == 2009
sort dom 
codebook dom

/** Limit to alive at discharge and abstracted only
keep if vstatus==1 & abstracted==1

**********************************************************************************************
** THERAPEUTIC #1: EITHER ASPIRIN OR CLOPIDOGREL
**********************************************************************************************
* There are two variables that define aspirin use: aspdis (aspirin at discharge) & asp_use (use of aspirin)
/* Natasha's explanation is as follows:
        When a person is admitted we ask whether they are on aspirin chronically and whether they have had 
        aspirin in the acute aftermath of the event. So there we will say- acute use, chronic use and no 
        record of use.
        The discharge variable is separate, here the abstractor will look to see if the person has been 
        given aspirin in their list of take home medications.
        So a person may be on aspirin chronically (before they come in) and then tell the doc they have 
        enough at home so then this is not in their discharge list.

        For the purpose of analysis, in the past I have assumed that a person will receive secondary 
        prevention form aspirin if they are on aspirin chronically (before admission) OR if they are 
        discharged on aspirin. It's not a perfect system but knowing the system I believe that this 
        method prevents against undercount and given our poor documentation system undercount is a 
        far bigger problem.
*/

* create aspirin/clopidogrel secondary prevention variable
gen AC_sec = 0
* aspirin at discharge 
codebook aspdis
replace AC_sec = 1 if aspdis == 1
* chronic aspirin use before admission
codebook asp_use
replace AC_sec = 1 if asp_use == 2 | asp_use == 4
* clopidogrel at discharge
codebook pladis
replace AC_sec = 1 if pladis == 1
* chronic clopidogrel use before admission
codebook pla_use
replace AC_sec = 1 if pla_use == 2 | pla_use == 4
* aggrenox at discharge
codebook aggdis // no recorded use of aggrenox at discharge
* chronic aggrenox use
codebook agg_use // no recorded chronic use of aggrenox
* aspirin + dypridamol at discharge
codebook agg_dis
replace AC_sec = 1 if agg_dis == 1
* dipyridamole at discharge
codebook dipyrdis
replace AC_sec = 1 if dipyrdis == 1
* chronic use dipyridamole 
codebook dipyr_use // no recorded chronic use of dipyridamole

label variable AC_sec "aspirin/clopidogrel secondary prevention"
label define noyes 0 "No" 1 "Yes"
label values AC_sec noyes
tab AC_sec


**********************************************************************************************
** THERAPEUTIC #2: BETA BLOCKERS
**********************************************************************************************
* create BB secondary prevention variable
gen bb_sec = 0
* BB at discharge
codebook betadis
replace bb_sec = 1 if betadis == 1
** BB chronic use
codebook betachr
replace bb_sec = 1 if betachr == 1

label variable bb_sec "beta blocker use at discharge"
label values bb_sec noyes
tab bb_sec

**********************************************************************************************
** THERAPEUTIC #3: STATINS
**********************************************************************************************
* create statin secondary prevention variable
gen stat_sec = 0
* statins at discharge
codebook statdis
replace stat_sec = 1 if statdis == 1
** statin chronic use
codebook stat_use 
replace stat_sec = 1 if stat_use == 2 | stat_use == 4

label variable stat_sec "statin use at discharge"
label values stat_sec noyes
tab stat_sec

**********************************************************************************************
** OVERALL SECONDARY PREVENTION SCORE
**********************************************************************************************
egen score_sec = rowtotal (AC_sec bb_sec stat_sec)
tab score_sec 


** SAVE DATASET FOR FURTHER ANALYSES
keep anon_pid AC_sec bb_sec stat_sec score_sec
save "`datapath'\version02\2-working\heart2010-2019_meds", replace
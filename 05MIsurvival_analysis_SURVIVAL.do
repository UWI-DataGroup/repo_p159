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


** Load mortality dataset
use "`datapath'\version02\2-working\heart2010-2019_mort", clear
** merge with medication dataset
merge 1:1 anon_pid using "`datapath'\version02\2-working\heart2010-2019_meds"
tab _merge
drop _merge


**-------------------------------------------------------------
** SURVIVAL ANALYSIS
**-------------------------------------------------------------
gen time1 = time + 1
drop time
rename time1 time

** Declare data to be survival-time data
stset time, failure(died) id (anon_pid)


** KAPLAN MEIER curve: only abstracted cases
keep if abstracted == 1
stdescribe
sts graph

#delimit ;
		sts graph, by(sex)   
			plotregion(c(gs16) ic(gs16) ilw(thin) lw(thin)) 
			graphregion(color(gs16) ic(gs16) ilw(thin) lw(thin)) 

				ylab(0.5(0.1)1, labs(small) nogrid angle(0)) 
				ytick (0.5(0.05)1) 
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

stcox i.score_sec
stcox score_sec
stcox i.sex age i.score_sec
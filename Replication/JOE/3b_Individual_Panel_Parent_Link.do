*******************************************************************************
* Codes to replicate the sample in ABBL (2021)
* This version: January 2021
* Questions: jdlight@uchicago.edu
*******************************************************************************

*******************************************************************************
* This file links individuals to their father
*******************************************************************************

clear all

cd "$output"

* Load the individual file
u "$data\2019_FAMILY_ID.dta", clear

keep PID2 PID3 PID23 PID24 

* Generate uniqe identifier
gen persid = (PID2*1000) + PID3
gen persid_father = (PID23*1000) + PID24

keep persid persid_father

save "$temp\parent_id_links.dta", replace

cd "$programs"

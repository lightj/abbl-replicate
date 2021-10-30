*******************************************************************************
* Codes to replicate the sample in ABBL (2021)
* This version: January 2021
* Questions: jdlight@uchicago.edu
*******************************************************************************

*******************************************************************************
* This file runs the first stage regressions to get residuals for income, 
* consumption and assets
*******************************************************************************

cd "$output"

*******************************************************************************
*** 1. Get residuals
*******************************************************************************

u "$temp\Full_Panel_3b.dta",clear
set more off
local wave 2 /* Determines minimum # of consecutive waves */
local first_stage 1 /* 1 yes, 0 no ==> run first stage-regs? */ 
local equiv 0 /*Use equivalence scales? */
local dual_earners 1 /* 1 dual earners, 0 no restriction */

*** Merge parental wealth, income and consumption
drop persid 
ren persid_father persid
merge m:1 persid using "$temp\Full_Panel_3b_parents.dta"
gen merge_linked = 0
replace merge_linked =1 if _merge==3
drop _merge

*** Variables to see if correlated with xi (ex-post)
bysort person: egen max_rent = max(rent)
gen renter = 0
replace renter = 1 if max_rent>0
bysort person: egen home_in_wealth = mean(home_equity/total_wealth)
bysort person: egen medical_proportion = mean(any_medical /(non_durable + services + durable))
gen negative_equity = 0
replace negative_equity = 1 if home_equity < 0

*** If using equivalence scales
if `equiv'==1{
replace log_labor_y=log_labor_y-log(equivalence)
replace log_net_y=log_net_y-log(equivalence)
replace log_ass=log_ass-log(equivalence)
replace log_cons_1 =log_cons_1-log(equivalence)
replace log_cons_2 =log_cons_2-log(equivalence)
replace log_cons_ABB=log_cons_ABB-log(equivalence)
}

*** Basic sample selection
keep if year >= 2005 & year <=2017
sort person year
keep if age >=25 & age <= 60
egen max_select = max(sample_select), by(person)
drop if max_select==0

*** Drop households with missing covariate info
foreach j in person yb year age educ age kids state he whe race wrace ybw weduc {
forvalues yr = 1999(2)2017{
drop if `j' ==. & year == `yr'
}
}

*** Drop households with missing info on consumption, income and assets
foreach j in log_net_y log_labor_y log_cons_ABB log_cons_1 log_cons_2 log_ass {
forvalues yr = 1999(2)2017{
drop if `j' ==. & year == `yr'
}
}

*** Apply sample selection for dual earners
if `dual_earners'==1{
drop if ly==0
drop if wly==0
}

*** Keep unbalanced panel but drop re-entrants
tsset person year, delta(2)
by person, sort: replace numwav=_N
keep if numwav >= `wave'
gen temp=1
by person, sort: gen n=sum(temp)
forvalues j = 1(1)10{
drop if n!=1 & year-2!=L.year
}
drop temp 
drop n
by person, sort: replace numwav=_N
keep if numwav >= `wave'

*** Run first-stage regressions
if `first_stage'==0{

* First-stage regressions
local measure log_labor_y log_net_y log_cons_1 log_cons_2 log_cons_ABB log_ass 
foreach j in `measure'{
gen u_`j' = `j' 
}

}
gen log_share = log(food_share)
else if `first_stage'==1{

gen age2 = age^2
gen age3 = age^3 
gen age4 = age^4

/*
xi i.state i.race i.he*i.coh i.whe*i.coh
xtreg log_net_y age age2 age3, fe 

gen u_log_net_y =log_net_y - _b[_cons] - _b[age]*age -_b[age2]*age2 - _b[age3]*age3
*/


xi	i.he*i.yb i.whe*i.ybw i.state i.fsize i.kids i.race i.wrace i.year i.age i.age*i.year

* First-stage regressions
local measure log_labor_y log_net_y log_cons_1 log_cons_2 log_cons_ABB log_ass log_food food_share log_share
foreach j in `measure'{
reg `j' kidsout extra  _I* /* t1 t12 t13 t14 */
predict u_`j' if e(sample), res
predict hat_`j' if e(sample)

drop if u_`j' ==.
}

}

*** Keep unbalanced panel but drop re-entrants (need to re-run after first-stage regs)
tsset person year, delta(2)
by person, sort: replace numwav=_N
keep if numwav >= `wave'
gen temp=1
by person, sort: gen n=sum(temp)
forvalues j = 1(1)10{
drop if n!=1 & year-2!=L.year
}
drop temp 
drop n
by person, sort: replace numwav=_N
keep if numwav >= `wave'

*** Get descriptive summary tables
label variable level_food "Food"
label variable level_nd_excl_food "Non-durables (excl. food)"
label variable level_cons_1 "Total Non-durables"
label variable level_house "Home equity"
label variable level_wealth_excl_house "Wealth (excl. home)"
label variable level_ass "Total wealth"
label variable level_labor_y "Labor income"
label variable level_net_y "Net income"
label variable negative_equity "Negative Equity Dummy"
label variable age "Age"
label variable educ "Education"
label variable kids "Kids"
forvalues year = 2005(2)2017{
eststo es`year': quietly estpost summarize level_food level_nd_excl_food level_cons_1 level_house negative_equity level_wealth_excl_house level_ass level_labor_y level_net_y if year ==`year'
}
esttab es2005 es2007 es2009 es2011 es2013 es2015 es2017 using all_assets_sample.tex, collabels(none) label nostar replace main(mean %9.2fc) aux(sd %9.2fc) mtitle("2005" "2007" "2009" "2011" "2013" "2015" "2017")

*** Get balanced data summary tables
/*
cap drop tempwav
by person, sort: gen tempwav=_n
* Output summary tables
forvalues year = 1(1)7{
eststo rs`year': quietly estpost summarize age educ kids level_food level_nd_excl_food level_cons_1 level_house negative_equity level_wealth_excl_house level_ass level_labor_y level_net_y if numwav ==`year' & tempwav==1
}

esttab rs1 rs2 rs3 rs4 rs5 rs6 rs7 using main_bal_sample.tex, collabels(none) label nostar replace main(mean %9.2fc) aux(sd %9.2fc) mtitle("Waves 1" "Waves 2" "Waves 3" "Waves 4" "Waves 5" "Waves 6" "Waves 7")

forvalues bl = 2(1)7{
eststo bal`bl': estpost summarize numwav if numwav ==`bl' & tempwav==1}
esttab bal2 bal3 bal4 bal5 bal6 bal7 using bal_sample.tex, collabels(none) label nostar replace main(N %6.2f) 
*/

*** Get full panel
tsfill, full

*** Get HE and RACE dummies
gen hemax = 0
replace hemax = 1 if max(educ,weduc)==6
egen hemax2 = max(hemax), by(person)
replace hemax=hemax2
gen racemax = 0
replace racemax =1 if race ==1

*** Prep + Save data for MATLAB
sort person year
gen D=0
replace D=1 if sample_select !=. 
order person year age D u_log_net_y hemax yb u_log_labor_y u_log_cons_1 u_log_cons_2 u_log_cons_ABB u_log_ass hat_log_net_y hat_log_labor_y hat_log_cons_1 hat_log_cons_2 hat_log_cons_ABB hat_log_ass renter home_in_wealth medical_proportion negative_equity u_log_food u_food_share parent_asset parent_income parent_cons u_parent_asset u_parent_income u_parent_cons
keep person year age D u_log_net_y hemax yb u_log_labor_y u_log_cons_1 u_log_cons_2 u_log_cons_ABB u_log_ass hat_log_net_y hat_log_labor_y hat_log_cons_1 hat_log_cons_2 hat_log_cons_ABB hat_log_ass renter home_in_wealth medical_proportion negative_equity u_log_food u_food_share parent_asset parent_income parent_cons u_parent_asset u_parent_income u_parent_cons
foreach var in person year age D u_log_net_y hemax yb u_log_labor_y u_log_cons_1 u_log_cons_2 u_log_cons_ABB u_log_ass hat_log_net_y hat_log_labor_y hat_log_cons_1 hat_log_cons_2 hat_log_cons_ABB hat_log_ass renter home_in_wealth medical_proportion negative_equity u_log_food u_food_share parent_asset parent_income parent_cons u_parent_asset u_parent_income u_parent_cons{
replace `var' = 0 if `var' ==.
}
* outsheet person year age D u_log_net_y hemax yb u_log_labor_y u_log_cons_1 u_log_cons_2 u_log_cons_ABB u_log_ass hat_log_net_y hat_log_labor_y hat_log_cons_1 hat_log_cons_2 hat_log_cons_ABB hat_log_ass renter home_in_wealth medical_proportion negative_equity u_log_food u_food_share parent_asset parent_income parent_cons u_parent_asset u_parent_income u_parent_cons using "$output/data4est_vFINAL_WITH_DUAL_RESTRICTION.out", nonames replace

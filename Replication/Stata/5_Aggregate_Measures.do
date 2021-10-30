********************************************************************************************
* Update: Jan 2019
* Queries: jdlight@uchicago.edu
********************************************************************************************

********************************************************************************************
* This file does the basic aggregation of main income, consumption and asset measures
********************************************************************************************
*
*\\\ Structure of File: \\\*
*
* 1. Loads data from TAXSIM and a CPI price index 
* 2. Winsorizing as specified by user 
* 3. Consumption
* 3. Assets
* 5. Income 
* 
********************************************************************************************


clear all
cd "$temp"

/*  1 = Assigning missing for observations with wage lower then 1/2 state ****/
global min_wage     = 1
                                                                                                                                                            
/* User parameter to control dropping outliers */
global pec_drop     = 0.25  /* Drop observation in the bottom "pec_drop" percent of "jumps" */


********************************************************************************************
* 1. Load TAXSIM and Price Indicies
********************************************************************************************

u Full_Panel_2b.dta,clear

drop if sample_select==0

keep if year >= 1999 & year <= 2017
sort person year

* Drop if not in age range
keep if age >=25 & age <= 60

* Merge with spouse links
merge m:1 id year using "$temp\spouse_links.dta"
ren _merge spouse_link


* Merge with parental links
gen persid_store = persid /* Store for later use */
merge m:1 persid using "$temp\parent_id_links.dta"
ren _merge merge_parent_link
ren persid_father persid_father_of_head
drop persid
rename spouse_persid persid
merge m:1 persid using "$temp\parent_id_links.dta"
ren _merge merge_parent_link2
ren persid_father persid_father_of_spouse
drop persid

replace persid_father_of_head =. if persid_father_of_head ==0
replace persid_father_of_spouse =. if persid_father_of_spouse ==0

gen persid_father = persid_father_of_head
replace persid_father =persid_father_of_spouse if persid_father==.

gen persid = persid_store


* Merge with price index
merge m:1 year using "$data\CPI.dta"
keep if _merge == 3 // even years from the using need to be dropped 
drop _merge

* Merge with TAXSIM estimates
merge m:1 taxsimid using "$TAXSIM\TAXSIM_main.dta"
keep if _merge == 3 // even years from the using need to be dropped 
drop _merge
gen taxes_all = fiitax + siitax + fica/2
drop fiitax siitax fica frate srate ficar

* Merge with minimum wages
merge m:1 state year using "$data\minwage.dta"
keep if _merge == 3 // even years from the using need to be dropped 
drop _merge
drop min_fed_mw min_mw max_fed_mw max_mw


gen log_ass = log(total_wealth) - log(price) 

********************************************************************************************
* 2. Dropping some outliers before constructing aggregates
********************************************************************************************

* if hourly wage very small, measurement error
gen test = 0
 if $min_wage==1 {
  gen oly = ly
  gen owly = wly
  drop if (ly/hours)<0.5*mean_mw | (wly/hourw)<0.5*mean_mw

 * replace ly=.  if (ly/hours)<0.5*mean_mw
 * replace wly=. if (wly/hourw)<0.5*mean_mw
}


if $pec_drop>0 {
 
  gen log_ly  = log(ly) - log(lag_price)
  gen log_wly = log(wly) - log(lag_price)

  egen non_durable = rsum(food_home food_out food_deliv gasoline clothing), missing
  egen services 	= 	rsum(utilities telecom auto_ins parking bus_train taxi other_trans education ///
  childcare inst_medical doctor prescription health_ins trips other_rec), missing
  gen totcons = non_durable+services

  replace totcons = totcons/price 
  gen log_totcons = log(totcons)	  /* log consumption */ 
 
  *** Use raw moments (not residual moments) that are associated with measurement error
  tsset person year
  local npec = 100/${pec_drop}
  di `npec'

  /* Generate the first difference for logs and the interaction between first difference and lagged first difference (which would be large in absolute value for large values of transitory shocks or for measurement error */               
  foreach var of varlist log_* {
    gen d_`var' = `var' - l2.`var'
    gen d_`var'_lag = d_`var'*l2.d_`var'
  }

  /* Generate percentiles of the interacted difference by year */
  foreach var of varlist d_*_lag {
    egen pec_`var'=  xtile(`var'), by(year) n(`npec')
	}

  /* Assign missing values for the variable with the potential measurement error  */
  foreach var of varlist ly wly totcons {
    replace `var'=. if f2.pec_d_log_`var'_lag==1  /* Note that we are only assinging missing values to the year with the jump */
  }
  drop d_* pec_* log_* totcons non_durable services
}


********************************************************************************************
* 3. Consumption
********************************************************************************************
gen equivalence = 1 + (fsize-kids-1)*0.7 + 0.5*kids

replace fstmp =0 if fstmp ==.
gen m_non_durable=0
foreach var in food_home food_out food_deliv gasoline clothing{
replace m_non_durable = 1 if `var' ==.
}
egen non_durable = rsum(clothing food_home food_out food_deliv gasoline), missing
egen all_food = rsum(food_home food_out food_deliv), missing
egen non_durable_excl_food = rsum(clothing gasoline), missing

gen m_housing=0
foreach var in mortgage rent prop_tax homeinsure{
replace m_housing = 1 if `var' ==.
}
egen housing = rsum(mortgage rent prop_tax homeinsure), missing


gen m_durable=0
foreach var in /*home_repair home_furnish*/ auto_repair vehicle_loan vehicle_down vehicle_lease other_vehicle{
replace m_durable = 1 if `var' ==.
}
egen durable = rsum(/*home_repair home_furnish*/ auto_repair vehicle_loan vehicle_down vehicle_lease other_vehicle), missing


gen m_services=0
foreach var in utilities telecom auto_ins parking bus_train taxi other_trans education ///
childcare inst_medical doctor prescription   trips other_rec{
replace m_services = 1 if `var' ==.
}
egen services 	= 	rsum(utilities telecom auto_ins parking bus_train taxi other_trans education ///
childcare inst_medical doctor prescription health_ins trips other_rec), missing

gen log_nd_excl_food = log(non_durable_excl_food + services) - log(price)
gen level_nd_excl_food = (non_durable_excl_food + services)/(price)

* Generate alternative imputed housing measure
gen m_imputed_housing=0
gen imputed_housing = .
replace imputed_housing = rent if rent>0
replace imputed_housing = 0.06*house if rent==0 & house>0
replace m_imputed_housing = 1 if imputed_housing==.


* Various aggregate consumption measures 
* For consistency with survey timing footstamps need weighting by lag price

gen log_cons_1 = log(((non_durable + services + fstmp)/price) + (fstmp/lag_price))
gen level_cons_1 = (((non_durable + services + fstmp)/price) + (fstmp/lag_price))

gen log_cons_2 = log(((non_durable + services + durable)/price) + (fstmp/lag_price))
gen log_food = log((all_food/price)+(fstmp/lag_price))
gen level_food = ((all_food/price)+(fstmp/lag_price))

gen food_share = ((all_food/price)+(fstmp/lag_price))/(((non_durable + services + durable)/price) + (fstmp/lag_price))

egen ABB = rsum(food_home food_out food_deliv gasoline utilities auto_ins parking bus_train taxi other_trans education ///
childcare inst_medical doctor prescription health_ins), missing

gen log_cons_ABB = log((ABB/price)+ (fstmp/lag_price))

egen any_medical = rsum(inst_medical doctor prescription health_ins), missing

********************************************************************************************
* 3. Wealth
********************************************************************************************

label var total_wealth "Total wealth net of debt"

gen log_house = log(home_equity) - log(price)
gen level_house = (home_equity/price)

gen log_wealth_excl_house =log(total_wealth-home_equity) - log(price)
gen level_wealth_excl_house =(total_wealth-home_equity/price)

gen total_wealth2=total_wealth/price

gen bus_debt_temp = - bus_debt
gen other_debt_temp = -other_debt
gen real_estate_debt_temp = -real_estate_debt
egen wealth_ABB = rsum(cash other_ass stocks bus_ass bus_debt_temp ira other_debt_temp home_equity vehicles real_estate real_estate_debt_temp)
gen wealth_ABB2 =wealth_ABB/price

/* ^^^^ Constructed wealth variable, including equity. This imputed variable is constructed as
the sum of values of seven asset types (ER71429 (business assets), ER71435 (cash), ER71439 (other real estate), ER71445 (stocks), ER71447 (car),
ER71451 (other ass), ER71455 (annuity)) net of debt value (ER71431, ER71441, ER71459, ER71463, ER71467, ER71471,
ER71475, ER71479) plus value of home equity (ER71481). All missing data were assigned.. */

gen log_ass = log(total_wealth) - log(price)
gen level_ass = total_wealth/price
gen log_ass2 = log(wealth_ABB2)

********************************************************************************************
* 4. Income
********************************************************************************************

gen log_w=log(ly/weeks_emp)-log(price)

gen log_net_y = log((ly + wly + trhw + soc_sec + soc_secw + fstmp - taxes_all)) - log(lag_price)
gen level_net_y = ((ly + wly + trhw + soc_sec + soc_secw + fstmp - taxes_all))/(lag_price)

gen log_labor_y = log((ly + wly)) - log(lag_price)
gen level_labor_y = ((ly + wly)) /(lag_price)

********************************************************************************************
* 5. Other variables needed for regressions
********************************************************************************************

gen kidsout	= outkid==1       // dummy for kids out
gen extra=(tyoth)>0           // dummy for income recipient other than h/w



save "Full_Panel_3b.dta", replace
 
cd "$programs"


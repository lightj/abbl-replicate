********************************************************************************************
* Update: Jan 2019
* Queries: jdlight@uchicago.edu
********************************************************************************************

********************************************************************************************
* This file basically replicates Aggregate Measures 
* We will use the first version for the "young" households
* We can then merge in the parents from this version to the young
********************************************************************************************


clear all
cd "$temp"

/*  1 = Assigning missing for observations with wage lower then 1/2 state ****/
global min_wage     = 1
                                                                                                                                                            
/* User parameter to control dropping outlyers */
global pec_drop     = 0.25  /* Drop observation in the bottom "pec_drop" percent of "jumps" */


********************************************************************************************
* 1. Load TAXSIM and Price Indicies
********************************************************************************************

u Full_Panel_2b.dta,clear

drop if sample_select==0

keep if year >= 1999 & year <= 2017
sort person year

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

xi	i.yb i.year

* First-stage regressions
local measure log_ass log_net_y log_cons_1
foreach j in `measure'{
reg `j' _I*
predict u_`j' if e(sample), res
predict hat_`j' if e(sample)
}

ren log_ass parent_asset
ren log_net_y parent_income
ren log_cons_1 parent_cons

ren u_log_ass u_parent_asset
ren u_log_net_y u_parent_income
ren u_log_cons_1 u_parent_cons

keep persid* person year parent_asset parent_income parent_cons u_parent_asset u_parent_income u_parent_cons

collapse (mean) parent_asset=parent_asset parent_income=parent_income parent_cons=parent_cons u_parent_asset=u_parent_asset u_parent_income=u_parent_income u_parent_cons=u_parent_cons, by(persid)

save "Full_Panel_3b_parents.dta", replace
 
cd "$programs"

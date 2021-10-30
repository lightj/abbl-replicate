*******************************************************************************
* Codes to replicate the sample in ABBL (2021)
* This version: January 2021
* Questions: jdlight@uchicago.edu
*******************************************************************************

********************************************************************************************
* This file constructs the panel of households from the PSID, focusing on variables
* related to income, consumption and wealth
*
* Structure of File:
*
* 0. Load and prepare 2017
* 1. Load and prepare 2015
* 2. Load and prepare 2013
* 3. Load and prepare 2011
* 4. Load and prepare 2009
* 5. Load and prepare 2007 - append wealth file
* 6. Load and prepare 2005 - append wealth file
* 7. Load and prepare 2003 - append wealth file
* 8. Load and prepare 2001 - append wealth file
* 9. Load and prepare 1999 - append wealth file
********************************************************************************************

********************************************************************************************
* 2017
********************************************************************************************

clear all

use "$data\2017_H", clear

* Rename used variables
ren ER66002 id
ren ER66016 fsize
ren ER66017 age
ren ER66018 sex
ren ER66019 agew
ren ER66020 sexw
ren ER66021 kids
ren ER66022 newborn
ren ER66024 marit
ren ER70882 race1
ren ER70883 race2
ren ER70884 race3
ren ER70885 race4
ren ER70744 wrace1
ren ER70745 wrace2
ren ER70746 wrace3
ren ER70747 wrace4
ren ER67759 outkid
ren ER66031 house
ren ER71533 smsa
ren ER66007 fchg
ren ER66003 state
ren ER71570 weight 
ren ER66164 empst1
ren ER66165 empst2
ren ER66166 empst3
ren ER66439 wempst1
ren ER66440 wempst2
ren ER66441 wempst3
ren ER71538 educ
ren ER71539 weduc
ren ER66766 foodstamps
ren ER71227 weeks_emp
ren ER71248 weeks_empw
ren ER66231 exp_job_years
ren ER66232 exp_job_months
ren ER66233 exp_job_weeks
ren ER66369 unempl 
ren ER66362 layoff

ren ER71331 tr_tanf
ren ER71333 tr_ssi
ren ER71335 tr_welf
ren ER71337 tr_vap
ren ER71339 tr_pen
ren ER71341 tr_ann
ren ER71343 tr_ira
ren ER71345 tr_ret
ren ER71347 tr_unem
ren ER71349 tr_comp
ren ER71351 tr_child
ren ER71353 tr_alim
ren ER71355 tr_help
ren ER71357 tr_oth
ren ER71359 tr_misc

ren ER71361 tr_tanfw
ren ER71363 tr_ssiw
ren ER71365 tr_welfw
ren ER71367 tr_vapw
ren ER71369 tr_penw
ren ER71371 tr_annw
ren ER71373 tr_iraw
ren ER71375 tr_retw
ren ER71377 tr_unemw
ren ER71379 tr_compw
ren ER71381 tr_childw
ren ER71383 tr_alimw
ren ER71385 tr_helpw
ren ER71387 tr_othw
ren ER71389 tr_miscw

replace weeks_emp =. if weeks_emp == 98 | weeks_emp ==99
replace exp_job_years =. if exp_job_years == 98 | exp_job_years ==99
replace exp_job_months =. if exp_job_months == 98 | exp_job_months ==99
replace exp_job_weeks =. if exp_job_weeks == 98 | exp_job_weeks ==99

replace unempl =. if unempl == 8 | unempl == 9
replace unempl =0 if unempl == 5

replace layoff =. if layoff == 8 | layoff == 9
replace layoff =0 if layoff == 5

replace educ =. if educ == 99
replace weduc =. if weduc == 99

/// ASSETS ///

* Note that these are taken from the imputed values, _acc variables denote imputed entries

ren ER71429 bus_ass
ren ER71430 bus_ass_acc

ren ER71431 bus_debt
ren ER71432 bus_debt_acc

ren ER71435 cash
ren ER71436 cash_acc  

ren ER71439 real_estate
ren ER71440 real_estate_acc

ren ER71441 real_estate_debt
ren ER71442 real_estate_debt_acc

ren ER71445 stocks
ren ER71446 stocks_acc

ren ER71447 vehicles
ren ER71448 vehicles_acc

ren ER71451 other_ass
ren ER71452 other_ass_acc

ren ER71455 ira
ren ER71456 ira_acc

ren ER71459 credit_card 
ren ER71460 credit_card_acc

ren ER71463 stu_debt
ren ER71464 stu_debt_acc

ren ER71467 med_debt
ren ER71468 med_debt_acc

ren ER71471 legal_debt
ren ER71472 legal_debt_acc

ren ER71475 fam_debt
ren ER71476 fam_debt_acc

ren ER71479 other_debt 
ren ER71480 other_debt_acc

/* ^^^^ Note this includes student, legal, medical and fam debt */

ren ER71481 home_equity
ren ER71482 home_equity_acc

ren ER71485 total_wealth
ren ER71486 total_wealth_acc


/* ^^^^ Constructed wealth variable, including equity. This imputed variable is constructed as
the sum of values of seven asset types (ER71429 (business assets), ER71435 (cash), ER71439 (other real estate), ER71445 (stocks), ER71447 (car),
ER71451 (other ass), ER71455 (annuity)) net of debt value (ER71431, ER71441, ER71459, ER71463, ER71467, ER71471,
ER71475, ER71479) plus value of home equity (ER71481). All missing data were assigned.. */

ren ER68010 pension_plan
replace pension_plan=. if pension_plan>999999997
ren ER68227 pension_sp
replace pension_sp=. if pension_sp>999999997

/// EXPENDITURE ///

gen       fstmp=0
cap program drop doit
program def doit
        local i=66770
        while `i' < ER66782 {
            replace ER`i'=. if ER`i'>1
        	local i=`i'+1
       }
end
qui doit
egen nm=rsum(ER66770-ER66781)
replace fstmp=ER66767              if ER66768==6 
replace fstmp=ER66767*nm           if ER66768==5 
replace fstmp=ER66767*((26/12)*nm) if ER66768==4  
replace fstmp=ER66767*((52/12)*nm) if ER66768==3
replace fstmp=.					 if ER66767>999997  | ER66768>6 
drop nm



* Note that these are taken from the imputed values, _acc variables denote imputed entries
ren ER71488 food_home
ren ER71489 food_out
ren ER71490 food_deliv

ren ER71492 mortgage  
ren ER66053 mortgage_raw1
ren ER66074 mortgage_raw2
gen mortgage_acc = 0
replace mortgage_acc = 1 if mortgage_raw1 >= 99998 | mortgage_raw2 >= 99998

ren ER71494 rent
ren ER66090 rent_raw
gen rent_acc = 0
replace rent_acc = 1 if rent_raw >= 99998
 
ren ER71495 prop_tax 
ren ER66045 prop_tax_raw
gen prop_tax_acc = 0
replace prop_tax_acc = 1 if prop_tax_raw >= 99998 

ren ER71496 homeinsure
ren ER66047 homeinsure_raw
gen homeinsure_acc = 0
replace homeinsure_acc = 1 if homeinsure_raw>=9998

ren ER71497 utilities
ren ER71498 gas_raw
ren ER71499 electric_raw
ren ER71500 water_raw
ren ER71501 otheru_raw
gen utilities_acc = 0
replace utilities_acc = 1 if gas_raw >= 9999999 | electric_raw >= 9999999 | water_raw >= 9999999 | otheru_raw >= 9999999 

ren ER71502 telecom
ren ER66119 telecom_raw
gen telecom_acc = 0
replace telecom_acc = 1 if telecom_raw>=9998

ren ER71504 vehicle_loan
ren ER66869 vehicle_loan_raw1 
ren ER66893 vehicle_loan_raw2
ren ER66917 vehicle_loan_raw3
gen vehicle_loan_acc = 0
replace vehicle_loan_acc = 1 if vehicle_loan_raw1 >=999998 | vehicle_loan_raw2 >=999998 | vehicle_loan_raw3 >=999998 

ren ER71505 vehicle_down
ren ER66866 vehicle_down_raw1
ren ER66890 vehicle_down_raw2
ren ER66914 vehicle_down_raw3
gen vehicle_down_acc = 0
replace vehicle_down_acc = 1 if vehicle_down_raw1 >= 999998 | vehicle_down_raw2 >= 999998 | vehicle_down_raw3 >= 999998 

ren ER71506 vehicle_lease
ren ER66873 vehicle_lease_raw1
ren ER66897 vehicle_lease_raw2
ren ER66921 vehicle_lease_raw3
gen vehicle_lease_acc = 0
replace vehicle_lease_acc = 1 if vehicle_lease_raw1 >=999998 | vehicle_lease_raw2 >=999998 |vehicle_lease_raw3 >=999998 

ren ER71507 auto_ins
ren ER66926 auto_ins_raw
gen auto_ins_acc = 0
replace auto_ins_acc = 1 if auto_ins_raw >= 999998

ren ER71508 other_vehicle
ren ER66930 other_vehicle_raw
gen other_vehicle_acc = 0
replace other_vehicle_acc  = 1 if other_vehicle_raw >=999998

ren ER71510 gasoline
ren ER66931 gasoline_raw
gen gasoline_acc = 0
replace gasoline_acc = 1 if gasoline_raw >= 99998

ren ER71509 auto_repair
ren ER66929 auto_repair_raw
gen auto_repair_acc = 0
replace auto_repair_acc = 1 if auto_repair_raw >=99998

ren ER71511 parking
ren ER66932 parking_raw
gen parking_acc = 0
replace parking_acc = 1 if parking_raw >= 99998 

ren ER71512 bus_train
ren ER66933 bus_train_raw
gen bus_train_acc = 0
replace bus_train_acc = 1 if bus_train_raw >= 99998

ren ER71513 taxi
ren ER66934 taxi_raw
gen taxi_acc = 0
replace taxi_acc = 1 if taxi_raw >= 99998

ren ER71514 other_trans
ren ER66935 other_trans_raw
gen other_trans_acc = 0
replace other_trans_acc = 1 if other_trans_raw >=99998

ren ER71515 education
ren ER66937 education_raw1
ren ER66939 education_raw2
gen education_acc = 0
replace education_acc = 1 if education_raw1 >= 999998 | education_raw2 >= 999998

ren ER71516 childcare
ren ER66744 childcare_raw
gen childcare_acc = 0
replace childcare_raw = 1 if childcare_raw >= 999998

ren ER71518 inst_medical
ren ER70689 inst_medical_raw
gen inst_medical_acc = 0
replace inst_medical_acc = 1 if inst_medical_raw >= 999998

ren ER71519 doctor
ren ER70694 doctor_raw
gen doctor_acc = 0
replace doctor_acc = 1 if doctor_raw >= 999998

ren ER71520 prescription
ren ER70698 prescription_raw
gen prescription_acc = 0
replace prescription_acc = 1 if prescription_raw >= 999998

ren ER71521 health_ins
ren ER70683 health_ins_raw
gen health_ins_acc = 0
replace health_ins_acc = 1 if health_ins_raw >= 99998

ren ER71523 home_repair
ren ER66944 home_repair_raw
gen home_repair_acc = 0
replace home_repair_acc = 1 if home_repair_raw >= 999998

ren ER71524 home_furnish
ren ER66949 home_furnish_raw
gen home_furnish_acc = 0
replace home_furnish_acc = 1 if home_furnish_raw >= 999998

ren ER71525 clothing
ren ER66954 clothing_raw
gen clothing_acc = 0
replace clothing_acc = 1 if clothing_raw >= 999998

ren ER71526 trips
ren ER66959 trips_raw
gen trips_acc = 0
replace trips_acc = 1 if trips_raw >= 999998

ren ER71527 other_rec
ren ER66964 other_rec_raw
gen other_rec_acc = 0
replace other_rec_acc = 1 if other_rec >= 999998


/// INCOME ///

ren ER71426 y
gen truncy=y<1
replace y=1 if y<1


ren ER71330 tyhw
ren ER71391 trhw
ren ER71398 tyoth
ren ER71419 troth
ren ER71420 soc_sec
ren ER71422 soc_secw
ren ER71424 soc_secyoth



* Asset Income *

ren ER71294 renty
ren ER71295 renty_acc
ren ER71322 rentyw
ren ER71323 rentyw_acc

ren ER71296 divi
ren ER71297 divi_acc
ren ER71324 diviw
ren ER71325 diviw_acc

ren ER71298 interest
ren ER71299 interest_acc
ren ER71326 interestw
ren ER71327 interestw_acc

ren ER71300 trustinc
ren ER71301 trustinc_acc
ren ER71328 trustincw
ren ER71329 trustincw_acc

ren ER71275 businc
ren ER71303 busincw

* Generate labor and asset part of farm income *
gen farmlabor=0
replace farmlabor=0.5*ER71272 			if ER71272>0
replace farmlabor=. 					if ER71272>9999997

gen farmasset=0
replace farmasset=0.5*ER71272 			if ER71272>0
replace farmasset=ER71272 				if ER71272<0
replace farmasset=. 					if ER71272>9999997

* Hours Worked *
ren ER71233 hours
ren ER71254 hourw

* Income *
ren ER71293 wages
ren ER71321 wagesw
ren ER71274 y_business
ren ER71302 yw_business

egen ly = rsum(wages farmlabor y_business) 
egen wly= rsum(wagesw yw_business)


/// Prep variables to save ///

#delimit;

keep 

id fsize age sex agew sexw kids marit race* wrace*
outkid house smsa fchg state weight *educ newborn tr_*

bus_ass bus_debt cash real_estate real_estate_debt stocks vehicles other_ass
ira other_debt home_equity total_wealth

fstmp food_home food_out food_deliv mortgage rent prop_tax homeinsure utilities telecom
vehicle_loan vehicle_down vehicle_lease auto_ins other_vehicle gasoline auto_repair
parking bus_train taxi other_trans education childcare inst_medical doctor prescription
health_ins home_repair home_furnish clothing trips other_rec

y truncy tyhw trhw tyoth troth soc_sec soc_secw soc_secyoth
farmlabor farmasset hours* wages* y_business yw_business ly wly weeks_emp* hourw
exp_job_years exp_job_months exp_job_weeks empst* wempst* unempl layoff
rent* divi* interest* trust* bus*

; 

#delimit cr

gen year=2017

save "$temp\HH_2017.dta", replace


********************************************************************************************
* 2015 
********************************************************************************************
clear all

use "$data\2015_H", clear

* Rename used variables
ren ER60002 id
ren ER60016 fsize
ren ER60017 age
ren ER60018 sex
ren ER60019 agew
ren ER60020 sexw
ren ER60021 kids
ren ER60022 newborn
ren ER60024 marit
ren ER64810 race1
ren ER64811 race2
ren ER64812 race3
ren ER64813 race4
ren ER64671 wrace1
ren ER64672 wrace2
ren ER64673 wrace3
ren ER64674 wrace4
ren ER61706 outkid
ren ER60031 house
ren ER65454 smsa
ren ER60007 fchg
ren ER60003 state
ren ER65492 weight 
ren ER60163 empst1
ren ER60164 empst2
ren ER60165 empst3
ren ER60426 wempst1
ren ER60427 wempst2
ren ER60428 wempst3
ren ER65459 educ
ren ER65460 weduc
ren ER60719 foodstamps
ren ER65150 weeks_emp
ren ER65171 weeks_empw
ren ER60228 exp_job_years
ren ER60229 exp_job_months
ren ER60230 exp_job_weeks
ren ER60366 unempl 
ren ER60359 layoff

replace weeks_emp =. if weeks_emp == 98 | weeks_emp ==99
replace exp_job_years =. if exp_job_years == 98 | exp_job_years ==99
replace exp_job_months =. if exp_job_months == 98 | exp_job_months ==99
replace exp_job_weeks =. if exp_job_weeks == 98 | exp_job_weeks ==99

replace unempl =. if unempl == 8 | unempl == 9
replace unempl =0 if unempl == 5

replace layoff =. if layoff == 8 | layoff == 9
replace layoff =0 if layoff == 5

replace educ =. if educ == 99
replace weduc =. if weduc == 99

/// ASSETS ///

* Note that these are taken from the imputed values, _acc variables denote imputed entries

ren ER65352 bus_ass
ren ER65353 bus_ass_acc

ren ER65354 bus_debt
ren ER65355 bus_debt_acc

ren ER65358 cash
ren ER65359 cash_acc  

ren ER65362 real_estate
ren ER65363 real_estate_acc

ren ER65364 real_estate_debt
ren ER65365 real_estate_debt_acc

ren ER65368 stocks
ren ER65369 stocks_acc

ren ER65370 vehicles
ren ER65371 vehicles_acc

ren ER65374 other_ass
ren ER65375 other_ass_acc

ren ER65378 ira
ren ER65379 ira_acc

ren ER65382 credit_card 
ren ER65383 credit_card_acc

ren ER65386 stu_debt
ren ER65387 stu_debt_acc

ren ER65390 med_debt
ren ER65391 med_debt_acc

ren ER65394 legal_debt
ren ER65395 legal_debt_acc

ren ER65398 fam_debt
ren ER65399 fam_debt_acc

ren ER65402 other_debt 
ren ER65403 other_debt_acc

/* ^^^^ Note this includes student, legal, medical and fam debt */

ren ER65404 home_equity
ren ER65405 home_equity_acc

ren ER65408 total_wealth
ren ER65409 total_wealth_acc

/* ^^^^ Constructed wealth variable, including equity. This imputed variable is constructed as
the sum of values of seven asset types (ER65352, ER65358, ER65362, ER65368, ER65370,
ER65374, ER65378) net of debt value (ER65354, ER65364, ER65382, ER65386, ER65390, ER65394,
ER65398, ER65402) plus value of home equity (ER65404). All missing data were assigned. */

ren ER61956 pension_plan
replace pension_plan=. if pension_plan>999999997
ren ER62173 pension_sp
replace pension_sp=. if pension_sp>999999997

/// EXPENDITURE ///

gen       fstmp=0
cap program drop doit
program def doit
        local i=60723
        while `i' < 60735 {
            replace ER`i'=. if ER`i'>1
        	local i=`i'+1
       }
end
qui doit
egen nm=rsum(ER60723-ER60734)
replace fstmp=ER60720              if ER60721==6 
replace fstmp=ER60720*nm           if ER60721==5 
replace fstmp=ER60720*((26/12)*nm) if ER60721==4  
replace fstmp=ER60720*((52/12)*nm) if ER60721==3
replace fstmp=.					 if ER60720>999997  | ER60721>6 
drop nm



* Note that these are taken from the imputed values, _acc variables denote imputed entries
ren ER65411 food_home
ren ER65412 food_out
ren ER65413 food_deliv

ren ER65415 mortgage  
ren ER60051 mortgage_raw1
ren ER60072 mortgage_raw2
gen mortgage_acc = 0
replace mortgage_acc = 1 if mortgage_raw1 >= 99998 | mortgage_raw2 >= 99998

ren ER65416 rent
ren ER60088 rent_raw
gen rent_acc = 0
replace rent_acc = 1 if rent_raw >= 99998
 
ren ER65417 prop_tax 
ren ER60043 prop_tax_raw
gen prop_tax_acc = 0
replace prop_tax_acc = 1 if prop_tax_raw >= 99998 

ren ER65418 homeinsure
ren ER60045 homeinsure_raw
gen homeinsure_acc = 0
replace homeinsure_acc = 1 if homeinsure_raw>=9998

ren ER65419 utilities
ren ER65420 gas_raw
ren ER65421 electric_raw
ren ER65422 water_raw
ren ER65423 otheru_raw
gen utilities_acc = 0
replace utilities_acc = 1 if gas_raw >= 9999999 | electric_raw >= 9999999 | water_raw >= 9999999 | otheru_raw >= 9999999 

ren ER65424 telecom
ren ER60118 telecom_raw
gen telecom_acc = 0
replace telecom_acc = 1 if telecom_raw>=9998

ren ER65426 vehicle_loan
ren ER60821 vehicle_loan_raw1 
ren ER60845 vehicle_loan_raw2
ren ER60869 vehicle_loan_raw3
gen vehicle_loan_acc = 0
replace vehicle_loan_acc = 1 if vehicle_loan_raw1 >=999998 | vehicle_loan_raw2 >=999998 | vehicle_loan_raw3 >=999998 

ren ER65427 vehicle_down
ren ER60818 vehicle_down_raw1
ren ER60842 vehicle_down_raw2
ren ER60866 vehicle_down_raw3
gen vehicle_down_acc = 0
replace vehicle_down_acc = 1 if vehicle_down_raw1 >= 999998 | vehicle_down_raw2 >= 999998 | vehicle_down_raw3 >= 999998 

ren ER65428 vehicle_lease
ren ER60825 vehicle_lease_raw1
ren ER60849 vehicle_lease_raw2
ren ER60873 vehicle_lease_raw3
gen vehicle_lease_acc = 0
replace vehicle_lease_acc = 1 if vehicle_lease_raw1 >=999998 | vehicle_lease_raw2 >=999998 |vehicle_lease_raw3 >=999998 

ren ER65429 auto_ins
ren ER60878 auto_ins_raw
gen auto_ins_acc = 0
replace auto_ins_acc = 1 if auto_ins_raw >= 999998

ren ER65430 other_vehicle
ren ER60882 other_vehicle_raw
gen other_vehicle_acc = 0
replace other_vehicle_acc  = 1 if other_vehicle_raw >=999998

ren ER65432 gasoline
ren ER60883 gasoline_raw
gen gasoline_acc = 0
replace gasoline_acc = 1 if gasoline_raw >= 99998

ren ER65431 auto_repair
ren ER60881 auto_repair_raw
gen auto_repair_acc = 0
replace auto_repair_acc = 1 if auto_repair_raw >=99998

ren ER65433 parking
ren ER60884 parking_raw
gen parking_acc = 0
replace parking_acc = 1 if parking_raw >= 99998 

ren ER65434 bus_train
ren ER60885 bus_train_raw
gen bus_train_acc = 0
replace bus_train_acc = 1 if bus_train_raw >= 99998

ren ER65435 taxi
ren ER60886 taxi_raw
gen taxi_acc = 0
replace taxi_acc = 1 if taxi_raw >= 99998

ren ER65436 other_trans
ren ER60887 other_trans_raw
gen other_trans_acc = 0
replace other_trans_acc = 1 if other_trans_raw >=99998

ren ER65437 education
ren ER60889 education_raw1
ren ER60891 education_raw2
gen education_acc = 0
replace education_acc = 1 if education_raw1 >= 999998 | education_raw2 >= 999998

ren ER65438 childcare
ren ER60697 childcare_raw
gen childcare_acc = 0
replace childcare_raw = 1 if childcare_raw >= 999998

ren ER65440 inst_medical
ren ER64613 inst_medical_raw
gen inst_medical_acc = 0
replace inst_medical_acc = 1 if inst_medical_raw >= 999998

ren ER65441 doctor
ren ER64619 doctor_raw
gen doctor_acc = 0
replace doctor_acc = 1 if doctor_raw >= 999998

ren ER65442 prescription
ren ER64625 prescription_raw
gen prescription_acc = 0
replace prescription_acc = 1 if prescription_raw >= 999998

ren ER65443 health_ins
ren ER64607 health_ins_raw
gen health_ins_acc = 0
replace health_ins_acc = 1 if health_ins_raw >= 99998

ren ER65444 home_repair
ren ER60892 home_repair_raw
gen home_repair_acc = 0
replace home_repair_acc = 1 if home_repair_raw >= 999998

ren ER65445 home_furnish
ren ER60897 home_furnish_raw
gen home_furnish_acc = 0
replace home_furnish_acc = 1 if home_furnish_raw >= 999998

ren ER65446 clothing
ren ER60902 clothing_raw
gen clothing_acc = 0
replace clothing_acc = 1 if clothing_raw >= 999998

ren ER65447 trips
ren ER60907 trips_raw
gen trips_acc = 0
replace trips_acc = 1 if trips_raw >= 999998

ren ER65448 other_rec
ren ER60912 other_rec_raw
gen other_rec_acc = 0
replace other_rec_acc = 1 if other_rec >= 999998


/// INCOME ///

ren ER65349 y
gen truncy=y<1
replace y=1 if y<1


ren ER65253 tyhw
ren ER65314 trhw
ren ER65321 tyoth
ren ER65342 troth
ren ER65343 soc_sec
ren ER65345 soc_secw
ren ER65347 soc_secyoth



* Asset Income *

ren ER65217 renty
ren ER65218 renty_acc
ren ER65245 rentyw
ren ER65246 rentyw_acc

ren ER65219 divi
ren ER65220 divi_acc
ren ER65247 diviw
ren ER65248 diviw_acc

ren ER65221 interest
ren ER65222 interest_acc
ren ER65249 interestw
ren ER65250 interestw_acc

ren ER65223 trustinc
ren ER65224 trustinc_acc
ren ER65251 trustincw
ren ER65252 trustincw_acc

ren ER65198 businc
ren ER65226 busincw


* Generate labor and asset part of farm income *
gen farmlabor=0
replace farmlabor=0.5*ER65195 			if ER65195>0
replace farmlabor=. 					if ER65195>9999997

gen farmasset=0
replace farmasset=0.5*ER65195 			if ER65195>0
replace farmasset=ER65195 				if ER65195<0
replace farmasset=. 					if ER65195>9999997

* Hours Worked *
ren ER65156 hours
ren ER65177 hourw

* Income *
ren ER65216 wages
ren ER65244 wagesw
ren ER65197 y_business
ren ER65225 yw_business

egen ly = rsum(wages farmlabor y_business) 
egen wly= rsum(wagesw yw_business)


/// Prep variables to save ///

#delimit;

keep 

id fsize age sex agew sexw kids marit race* wrace*
outkid house smsa fchg state weight *educ newborn

bus_ass bus_debt cash real_estate real_estate_debt stocks vehicles other_ass
ira other_debt home_equity total_wealth

fstmp food_home food_out food_deliv mortgage rent prop_tax homeinsure utilities telecom
vehicle_loan vehicle_down vehicle_lease auto_ins other_vehicle gasoline auto_repair
parking bus_train taxi other_trans education childcare inst_medical doctor prescription
health_ins home_repair home_furnish clothing trips other_rec

y truncy tyhw trhw tyoth troth soc_sec soc_secw soc_secyoth
farmlabor farmasset hours* wages* y_business yw_business ly wly weeks_emp* hourw
exp_job_years exp_job_months exp_job_weeks empst* wempst* unempl layoff
rent* divi* interest* trust* bus*
; 

#delimit cr

gen year=2015

save "$temp\HH_2015.dta", replace


********************************************************************************************
* 2013
********************************************************************************************

use "$data\2013_H", clear

ren ER53002 id
ren ER53016 fsize
ren ER53017 age
ren ER53018 sex
ren ER53019 agew
ren ER53020 kids
ren ER53021 newborn
ren ER53023 marit
ren ER57659 race1
ren ER57660 race2
ren ER57661 race3
ren ER57662 race4
ren ER57549 wrace1
ren ER57550 wrace2
ren ER57551 wrace3
ren ER57552 wrace4
ren ER54595 outkid
ren ER53030 house
ren ER58218 smsa
ren ER53007 fchg
ren ER53003 state
ren ER58257 weight 
ren ER53148 empst1
ren ER53149 empst2
ren ER53150 empst3
ren ER53411 wempst1
ren ER53412 wempst2
ren ER53413 wempst3
ren ER53704 foodstamps
ren ER58223 educ
ren ER58224 weduc
ren ER57970 weeks_emp
ren ER57991 weeks_empw
ren ER53213 exp_job_years
ren ER53214 exp_job_months
ren ER53215 exp_job_weeks
ren ER53351 unempl 
ren ER53344 layoff

replace weeks_emp =. if weeks_emp == 98 | weeks_emp ==99
replace exp_job_years =. if exp_job_years == 98 | exp_job_years ==99
replace exp_job_months =. if exp_job_months == 98 | exp_job_months ==99
replace exp_job_weeks =. if exp_job_weeks == 98 | exp_job_weeks ==99


replace educ =. if educ == 99
replace weduc =. if weduc == 99

replace unempl =. if unempl == 8 | unempl == 9
replace unempl =0 if unempl == 5

replace layoff =. if layoff == 8 | layoff == 9
replace layoff =0 if layoff == 5

/// ASSETS ///

* Note that these are taken from the imputed values, _acc variables denote imputed entries

ren ER58155 bus_ass
ren ER58156 bus_ass_acc

ren ER58157 bus_debt
ren ER58158 bus_debt_acc

ren ER58161 cash
ren ER58162 cash_acc  

ren ER58165 real_estate
ren ER58166 real_estate_acc

ren ER58167 real_estate_debt
ren ER58168 real_estate_debt_acc

ren ER58171 stocks
ren ER58172 stocks_acc

ren ER58173 vehicles
ren ER58174 vehicles_acc

ren ER58177 other_ass
ren ER58178 other_ass_acc

ren ER58181 ira
ren ER58182 ira_acc

ren ER58185 credit_card
ren ER58186 credit_card_acc

ren ER58189 stu_debt
ren ER58190 stu_debt_acc

ren ER58193 med_debt
ren ER58194 med_debt_acc

ren ER58197 legal_debt
ren ER58198 legal_debt_acc

ren ER58201 fam_debt
ren ER58202 fam_debt_acc

ren ER58205 other_debt 
ren ER58206 other_debt_acc

/* ^^^^ Note this includes student, legal, medical and fam debt */

ren ER58207 home_equity
ren ER58208 home_equity_acc

ren ER58211 total_wealth
ren ER58212 total_wealth_acc

/* ^^^^ Constructed wealth variable, including equity. This imputed variable is constructed as
the sum of values of seven asset types (ER65352, ER65358, ER65362, ER65368, ER65370,
ER65374, ER65378) net of debt value (ER65354, ER65364, ER65382, ER65386, ER65390, ER65394,
ER65398, ER65402) plus value of home equity (ER65404). All missing data were assigned. */

ren ER54836 pension_plan
replace pension_plan=. if pension_plan>999999997
ren ER55052 pension_sp
replace pension_sp=. if pension_sp>999999997

/// EXPENDITURE ///

gen       fstmp=0
cap program drop doit
program def doit
        local i=53708
        while `i' < 53720 {
            replace ER`i'=. if ER`i'>1
        	local i=`i'+1
       }
end
qui doit
egen nm=rsum(ER53708-ER53719)
replace fstmp=ER53705              if ER53706==6 
replace fstmp=ER53705*nm           if ER53706==5 
replace fstmp=ER53705*((26/12)*nm) if ER53706==4  
replace fstmp=ER53705*((52/12)*nm) if ER53706==3
replace fstmp=.					 if ER53705>999997  | ER53706>6 
drop nm


* Note that these are taken from the imputed values, _acc variables denote imputed entries

ren ER58212A2 food_home
ren ER58212A3 food_out
ren ER58212A4 food_deliv

ren ER58212A6 mortgage  
ren ER53050 mortgage_raw1
ren ER53071 mortgage_raw2
gen mortgage_acc = 0
replace mortgage_acc = 1 if mortgage_raw1 >= 99998 | mortgage_raw2 >= 99998

ren ER58212A7 rent
ren ER53087 rent_raw
gen rent_acc = 0
replace rent_acc = 1 if rent_raw >= 99998
 
ren ER58212A8 prop_tax 
ren ER53042 prop_tax_raw
gen prop_tax_acc = 0
replace prop_tax_acc = 1 if prop_tax_raw >= 99998 

ren ER58212A9 homeinsure
ren ER53044 homeinsure_raw
gen homeinsure_acc = 0
replace homeinsure_acc = 1 if homeinsure_raw>=9998

ren ER58212B1 utilities
ren ER58212B2 gas_raw
ren ER58212B3 electric_raw
ren ER58212B4 water_raw
ren ER58212B5 otheru_raw
gen utilities_acc = 0
replace utilities_acc = 1 if gas_raw >= 9999999 | electric_raw >= 9999999 | water_raw >= 9999999 | otheru_raw >= 9999999 

ren ER58212B6 telecom
ren ER53123 telecom_raw
gen telecom_acc = 0
replace telecom_acc = 1 if telecom_raw>=9998

ren ER58212B8 vehicle_loan
ren ER53762 vehicle_loan_raw1 
ren ER53786 vehicle_loan_raw2
ren ER53810 vehicle_loan_raw3
gen vehicle_loan_acc = 0
replace vehicle_loan_acc = 1 if vehicle_loan_raw1 >=999998 | vehicle_loan_raw2 >=999998 | vehicle_loan_raw3 >=999998 

ren ER58212B9 vehicle_down
ren ER53759 vehicle_down_raw1
ren ER53783 vehicle_down_raw2
ren ER53807 vehicle_down_raw3
gen vehicle_down_acc = 0
replace vehicle_down_acc = 1 if vehicle_down_raw1 >= 999998 | vehicle_down_raw2 >= 999998 | vehicle_down_raw3 >= 999998 

ren ER58212C1 vehicle_lease
ren ER53766 vehicle_lease_raw1
ren ER53790 vehicle_lease_raw2
ren ER53814 vehicle_lease_raw3
gen vehicle_lease_acc = 0
replace vehicle_lease_acc = 1 if vehicle_lease_raw1 >=999998 | vehicle_lease_raw2 >=999998 |vehicle_lease_raw3 >=999998 

ren ER58212C2 auto_ins
ren ER53819 auto_ins_raw
gen auto_ins_acc = 0
replace auto_ins_acc = 1 if auto_ins_raw >= 999998

ren ER58212C3 other_vehicle
ren ER53822 other_vehicle_raw
gen other_vehicle_acc = 0
replace other_vehicle_acc  = 1 if other_vehicle_raw >=999998

ren ER58212C5 gasoline
ren ER53824 gasoline_raw
gen gasoline_acc = 0
replace gasoline_acc = 1 if gasoline_raw >= 99998

ren ER58212C4 auto_repair
ren ER53823 auto_repair_raw
gen auto_repair_acc = 0
replace auto_repair_acc = 1 if auto_repair_raw >=99998

ren ER58212C6 parking
ren ER53825 parking_raw
gen parking_acc = 0
replace parking_acc = 1 if parking_raw >= 99998 

ren ER58212C7 bus_train
ren ER53826 bus_train_raw
gen bus_train_acc = 0
replace bus_train_acc = 1 if bus_train_raw >= 99998

ren ER58212C8 taxi
ren ER53827 taxi_raw
gen taxi_acc = 0
replace taxi_acc = 1 if taxi_raw >= 99998

ren ER58212C9 other_trans
ren ER53828 other_trans_raw
gen other_trans_acc = 0
replace other_trans_acc = 1 if other_trans_raw >=99998

ren ER58212D1 education
ren ER53830 education_raw1
ren ER53832 education_raw2
gen education_acc = 0
replace education_acc = 1 if education_raw1 >= 999998 | education_raw2 >= 999998

ren ER58212D2 childcare
ren ER53682 childcare_raw
gen childcare_acc = 0
replace childcare_raw = 1 if childcare_raw >= 999998

ren ER58212D4 inst_medical
ren ER57491 inst_medical_raw
gen inst_medical_acc = 0
replace inst_medical_acc = 1 if inst_medical_raw >= 999998

ren ER58212D5 doctor
ren ER57497 doctor_raw
gen doctor_acc = 0
replace doctor_acc = 1 if doctor_raw >= 999998

ren ER58212D6 prescription
ren ER57503 prescription_raw
gen prescription_acc = 0
replace prescription_acc = 1 if prescription_raw >= 999998

ren ER58212D7 health_ins
ren ER57485 health_ins_raw
gen health_ins_acc = 0
replace health_ins_acc = 1 if health_ins_raw >= 99998

ren ER58212D8 home_repair
ren ER53833 home_repair_raw
gen home_repair_acc = 0
replace home_repair_acc = 1 if home_repair_raw >= 999998

ren ER58212D9 home_furnish
ren ER53838 home_furnish_raw
gen home_furnish_acc = 0
replace home_furnish_acc = 1 if home_furnish_raw >= 999998

ren ER58212E1 clothing
ren ER53843 clothing_raw
gen clothing_acc = 0
replace clothing_acc = 1 if clothing_raw >= 999998

ren ER58212E2 trips
ren ER53848 trips_raw
gen trips_acc = 0
replace trips_acc = 1 if trips_raw >= 999998

ren ER58212E3 other_rec
ren ER53853 other_rec_raw
gen other_rec_acc = 0
replace other_rec_acc = 1 if other_rec >= 999998


/// INCOME ///

ren ER58152 y
gen truncy=y<1
replace y=1 if y<1

ren ER58060 tyhw
ren ER58117 trhw
ren ER58124 tyoth
ren ER58145 troth
ren ER58146 soc_sec
ren ER58148 soc_secw
ren ER58150 soc_secyoth



* Asset Income *

ren ER58039 renty
ren ER58040 renty_acc
ren ER58052 rentyw
ren ER58053 rentyw_acc

ren ER58041 divi
ren ER58042 divi_acc
ren ER58054 diviw
ren ER58055 diviw_acc

ren ER58043 interest
ren ER58044 interest_acc
ren ER58056 interestw
ren ER58057 interestw_acc

ren ER58045 trustinc
ren ER58046 trustinc_acc
ren ER58058 trustincw
ren ER58059 trustincw_acc

ren ER58018 businc
ren ER58048 busincw

* Generate labor and asset part of farm income *
gen farmlabor=0
replace farmlabor=0.5*ER58015 			if ER58015>0
replace farmlabor=. 					if ER58015>9999997

gen farmasset=0
replace farmasset=0.5*ER58015 			if ER58015>0
replace farmasset=ER58015 				if ER58015<0
replace farmasset=. 					if ER58015>9999997

* Hours Worked *
ren ER57976 hours
ren ER57997 hourw

* Income *
ren ER58038 wages
ren ER58050 wagesw
ren ER58017 y_business
ren ER58047 yw_business

egen ly = rsum(wages farmlabor y_business) 
egen wly= rsum(wagesw yw_business)


/// Prep variables to save ///

#delimit;

keep 

id fsize age sex agew kids marit race* wrace*
outkid house smsa fchg state weight foodstamps *educ newborn

bus_ass bus_debt cash real_estate real_estate_debt stocks vehicles other_ass
ira other_debt home_equity total_wealth

fstmp food_home food_out food_deliv mortgage rent prop_tax homeinsure utilities telecom
vehicle_loan vehicle_down vehicle_lease auto_ins other_vehicle gasoline auto_repair
parking bus_train taxi other_trans education childcare inst_medical doctor prescription
health_ins home_repair home_furnish clothing trips other_rec

y truncy tyhw trhw tyoth troth soc_sec soc_secw soc_secyoth
farmlabor farmasset hours* wages* y_business yw_business ly wly weeks_emp* hourw
exp_job_years exp_job_months exp_job_weeks empst* wempst* unempl layoff
rent* divi* interest* trust* bus*
; 

#delimit cr


gen year=2013

save "$temp\HH_2013.dta", replace

********************************************************************************************
* 2011
********************************************************************************************

use "$data\2011_H", clear

ren ER47302 id
ren ER47316 fsize
ren ER47317 age
ren ER47318 sex
ren ER47319 agew
ren ER47320 kids
ren ER47321 newborn
ren ER47323 marit
ren ER51904 race1
ren ER51905 race2
ren ER51906 race3
ren ER51907 race4
ren ER51810 wrace1
ren ER51811 wrace2
ren ER51812 wrace3
ren ER51813 wrace4
ren ER48852 outkid
ren ER47330 house
ren ER52400 smsa
ren ER47307 fchg
ren ER47303 state
ren ER52436 weight 
ren ER47448 empst1
ren ER47449 empst2
ren ER47450 empst3
ren ER47705 wempst1
ren ER47706 wempst2
ren ER47707 wempst3
ren ER48007 foodstamps
ren ER52405 educ
ren ER52406 weduc
ren ER52169 weeks_emp
ren ER52190 weeks_empw
ren ER47513 exp_job_years
ren ER47514 exp_job_months
ren ER47515 exp_job_weeks
ren ER47651 unempl 
ren ER47644 layoff

replace weeks_emp =. if weeks_emp == 98 | weeks_emp ==99
replace exp_job_years =. if exp_job_years == 98 | exp_job_years ==99
replace exp_job_months =. if exp_job_months == 98 | exp_job_months ==99
replace exp_job_weeks =. if exp_job_weeks == 98 | exp_job_weeks ==99

replace unempl =. if unempl == 8 | unempl == 9
replace unempl =0 if unempl == 5

replace layoff =. if layoff == 8 | layoff == 9
replace layoff =0 if layoff == 5

replace educ =. if educ == 99
replace weduc =. if weduc == 99

/// ASSETS ///

* Note that these are taken from the imputed values, _acc variables denote imputed entries

ren ER52346 bus_ass
ren ER52347 bus_ass_acc
gen bus_debt =.
/* ^^^ Note that Prior to 2013 the equity value is given without debt/asset split */

ren ER52350 cash
ren ER52351 cash_acc  

ren ER52354 real_estate
ren ER52355 real_estate_acc
gen real_estate_debt=.

/* ^^^ Note that Prior to 2013 the equity value is given without debt/asset split */

ren ER52358 stocks
ren ER52359 stocks_acc

ren ER52360 vehicles
ren ER52361 vehicles_acc

ren ER52364 other_ass
ren ER52365 other_ass_acc

ren ER52368 ira
ren ER52369 ira_acc

ren ER52372 credit_card
ren ER52373 credit_card_acc

ren ER52376 stu_debt
ren ER52377 stu_debt_acc

ren ER52380 med_debt
ren ER52381 med_debt_acc

ren ER52384 legal_debt
ren ER52385 legal_debt_acc

ren ER52388 fam_debt
ren ER52389 fam_debt_acc

gen other_debt = credit_card + stu_debt + med_debt + legal_debt + fam_debt
/* ^^^^ Note "other debt" to be calculated manually */

ren ER52390 home_equity
ren ER52391 home_equity_acc

ren ER52394 total_wealth
ren ER52395 total_wealth_acc

/* ^^^^ Constructed wealth variable, including equity. This imputed variable is constructed as
the sum of values of seven asset types (ER65352, ER65358, ER65362, ER65368, ER65370,
ER65374, ER65378) net of debt value (ER65354, ER65364, ER65382, ER65386, ER65390, ER65394,
ER65398, ER65402) plus value of home equity (ER65404). All missing data were assigned. */

ren ER49080 pension_plan
replace pension_plan=. if pension_plan>999999997
ren ER49299 pension_sp
replace pension_sp=. if pension_sp>999999997

/// EXPENDITURE ///

gen       fstmp=0
cap program drop doit
program def doit
        local i=48011
        while `i' < 48023 {
            replace ER`i'=. if ER`i'>1
        	local i=`i'+1
       }
end
qui doit
egen nm=rsum(ER48011-ER48022)
replace fstmp=ER48008              if ER48009==6 
replace fstmp=ER48008*nm           if ER48009==5 
replace fstmp=ER48008*((26/12)*nm) if ER48009==4  
replace fstmp=ER48008*((52/12)*nm) if ER48009==3
replace fstmp=.					 if ER48008>999997  | ER48009>6 
drop nm



* Note that these are taken from the imputed values, _acc variables denote imputed entries

ren ER52395A2 food_home
ren ER52395A3 food_out
ren ER52395A4 food_deliv

/* ^^^^ STILL NEED TO DEAL WITH FOOD STAMPS / UTILITY SUBSIDY */

ren ER52395A6 mortgage  
ren ER47350 mortgage_raw1
ren ER47371 mortgage_raw2
gen mortgage_acc = 0
replace mortgage_acc = 1 if mortgage_raw1 >= 99998 | mortgage_raw2 >= 99998

ren ER52395A7 rent
ren ER47387 rent_raw
gen rent_acc = 0
replace rent_acc = 1 if rent_raw >= 99998
 
ren ER52395A8 prop_tax 
ren ER47342 prop_tax_raw
gen prop_tax_acc = 0
replace prop_tax_acc = 1 if prop_tax_raw >= 99998 

ren ER52395A9 homeinsure
ren ER47344 homeinsure_raw
gen homeinsure_acc = 0
replace homeinsure_acc = 1 if homeinsure_raw>=9998

ren ER52395B1 utilities
ren ER52395B2 gas_raw
ren ER52395B3 electric_raw
ren ER52395B4 water_raw
ren ER52395B5 otheru_raw
gen utilities_acc = 0
replace utilities_acc = 1 if gas_raw >= 9999999 | electric_raw >= 9999999 | water_raw >= 9999999 | otheru_raw >= 9999999 

ren ER52395B6 telecom
ren ER47423 telecom_raw
gen telecom_acc = 0
replace telecom_acc = 1 if telecom_raw>=9998

ren ER52395B8 vehicle_loan
ren ER48066 vehicle_loan_raw1 
ren ER48091 vehicle_loan_raw2
ren ER48116 vehicle_loan_raw3
gen vehicle_loan_acc = 0
replace vehicle_loan_acc = 1 if vehicle_loan_raw1 >=999998 | vehicle_loan_raw2 >=999998 | vehicle_loan_raw3 >=999998 

ren ER52395B9 vehicle_down
ren ER48063 vehicle_down_raw1
ren ER48088 vehicle_down_raw2
ren ER48113 vehicle_down_raw3
gen vehicle_down_acc = 0
replace vehicle_down_acc = 1 if vehicle_down_raw1 >= 999998 | vehicle_down_raw2 >= 999998 | vehicle_down_raw3 >= 999998 

ren ER52395C1 vehicle_lease
ren ER48070 vehicle_lease_raw1
ren ER48095 vehicle_lease_raw2
ren ER48120 vehicle_lease_raw3
gen vehicle_lease_acc = 0
replace vehicle_lease_acc = 1 if vehicle_lease_raw1 >=999998 | vehicle_lease_raw2 >=999998 |vehicle_lease_raw3 >=999998 

ren ER52395C2 auto_ins
ren ER48125 auto_ins_raw
gen auto_ins_acc = 0
replace auto_ins_acc = 1 if auto_ins_raw >= 999998

ren ER52395C3 other_vehicle
ren ER48128 other_vehicle_raw
gen other_vehicle_acc = 0
replace other_vehicle_acc  = 1 if other_vehicle_raw >=999998

ren ER52395C5 gasoline
ren ER48130 gasoline_raw
gen gasoline_acc = 0
replace gasoline_acc = 1 if gasoline_raw >= 99998

ren ER52395C4 auto_repair
ren ER48129 auto_repair_raw
gen auto_repair_acc = 0
replace auto_repair_acc = 1 if auto_repair_raw >=99998

ren ER52395C6 parking
ren ER48131 parking_raw
gen parking_acc = 0
replace parking_acc = 1 if parking_raw >= 99998 

ren ER52395C7 bus_train
ren ER48132 bus_train_raw
gen bus_train_acc = 0
replace bus_train_acc = 1 if bus_train_raw >= 99998

ren ER52395C8 taxi
ren ER48133 taxi_raw
gen taxi_acc = 0
replace taxi_acc = 1 if taxi_raw >= 99998

ren ER52395C9 other_trans
ren ER48134 other_trans_raw
gen other_trans_acc = 0
replace other_trans_acc = 1 if other_trans_raw >=99998

ren ER52395D1 education
ren ER48136 education_raw1
ren ER48138 education_raw2
gen education_acc = 0
replace education_acc = 1 if education_raw1 >= 999998 | education_raw2 >= 999998

ren ER52395D2 childcare
ren ER47970 childcare_raw
gen childcare_acc = 0
replace childcare_raw = 1 if childcare_raw >= 999998

ren ER52395D4 inst_medical
ren ER51748 inst_medical_raw
gen inst_medical_acc = 0
replace inst_medical_acc = 1 if inst_medical_raw >= 999998

ren ER52395D5 doctor
ren ER51754 doctor_raw
gen doctor_acc = 0
replace doctor_acc = 1 if doctor_raw >= 999998

ren ER52395D6 prescription
ren ER51760 prescription_raw
gen prescription_acc = 0
replace prescription_acc = 1 if prescription_raw >= 999998

ren ER52395D7 health_ins
ren ER51744 health_ins_raw
gen health_ins_acc = 0
replace health_ins_acc = 1 if health_ins_raw >= 99998

ren ER52395D8 home_repair
ren ER48139 home_repair_raw
gen home_repair_acc = 0
replace home_repair_acc = 1 if home_repair_raw >= 999998

ren ER52395D9 home_furnish
ren ER48144 home_furnish_raw
gen home_furnish_acc = 0
replace home_furnish_acc = 1 if home_furnish_raw >= 999998

ren ER52395E1 clothing
ren ER48149 clothing_raw
gen clothing_acc = 0
replace clothing_acc = 1 if clothing_raw >= 999998

ren ER52395E2 trips
ren ER48154 trips_raw
gen trips_acc = 0
replace trips_acc = 1 if trips_raw >= 999998

ren ER52395E3 other_rec
ren ER48159 other_rec_raw
gen other_rec_acc = 0
replace other_rec_acc = 1 if other_rec >= 999998


/// INCOME ///

ren ER52343 y
gen truncy=y<1
replace y=1 if y<1

ren ER52259 tyhw
ren ER52308 trhw
ren ER52315 tyoth
ren ER52336 troth
ren ER52337 soc_sec
ren ER52339 soc_secw
ren ER52341 soc_secyoth



* Asset Income *

ren ER52238 renty
ren ER52239 renty_acc
ren ER52251 rentyw
ren ER52252 rentyw_acc

ren ER52240 divi
ren ER52241 divi_acc
ren ER52253 diviw
ren ER52254 diviw_acc

ren ER52242 interest
ren ER52243 interest_acc
ren ER52255 interestw
ren ER52256 interestw_acc

ren ER52244 trustinc
ren ER52245 trustinc_acc
ren ER52257 trustincw
ren ER52258 trustincw_acc

ren ER52217 businc
ren ER52247 busincw

* Generate labor and asset part of farm income *
gen farmlabor=0
replace farmlabor=0.5*ER52214 			if ER52214>0
replace farmlabor=. 					if ER52214>9999997

gen farmasset=0
replace farmasset=0.5*ER52214 			if ER52214>0
replace farmasset=ER52214 				if ER52214<0
replace farmasset=. 					if ER52214>9999997

* Hours Worked *
ren ER52175 hours
ren ER52196 hourw

* Income *
ren ER52237 wages
ren ER52249 wagesw
ren ER52216 y_business
ren ER52246 yw_business

egen ly = rsum(wages farmlabor y_business) 
egen wly= rsum(wagesw yw_business)


/// Prep variables to save ///

#delimit;

keep 


id fsize age sex agew kids marit race* wrace*
outkid house smsa fchg state weight foodstamps *educ newborn

bus_ass bus_debt cash real_estate real_estate_debt stocks vehicles other_ass
ira other_debt home_equity total_wealth

fstmp food_home food_out food_deliv mortgage rent prop_tax homeinsure utilities telecom
vehicle_loan vehicle_down vehicle_lease auto_ins other_vehicle gasoline auto_repair
parking bus_train taxi other_trans education childcare inst_medical doctor prescription
health_ins home_repair home_furnish clothing trips other_rec

y truncy tyhw trhw tyoth troth soc_sec soc_secw soc_secyoth
farmlabor farmasset hours* wages* y_business yw_business ly wly weeks_emp* hourw
exp_job_years exp_job_months exp_job_weeks empst* wempst* unempl layoff
rent* divi* interest* trust* bus*
; 

#delimit cr


gen year=2011

save "$temp\HH_2011.dta", replace

********************************************************************************************
* 2009
********************************************************************************************

use "$data\2009_H", clear

ren ER42002 id
ren ER42016 fsize
ren ER42017 age
ren ER42018 sex
ren ER42019 agew
ren ER42020 kids
ren ER42021 newborn
ren ER42023 marit
ren ER46543 race1
ren ER46544 race2
ren ER46545 race3
ren ER46546 race4
ren ER46449 wrace1
ren ER46450 wrace2
ren ER46451 wrace3
ren ER46452 wrace4
ren ER43527 outkid
ren ER42030 house
ren ER46976 smsa  
ren ER42007 fchg
ren ER42003 state
ren ER47012 weight 
ren ER42140 empst1
ren ER42141 empst2
ren ER42142 empst3
ren ER42392 wempst1
ren ER42393 wempst2
ren ER42394 wempst3
ren ER42691 foodstamps
ren ER46981 educ
ren ER46982 weduc
ren ER46761 weeks_emp
ren ER46782 weeks_empw
ren ER42200 exp_job_years
ren ER42201 exp_job_months
ren ER42202 exp_job_weeks
ren ER42338 unempl 
ren ER42331 layoff 

replace weeks_emp =. if weeks_emp == 98 | weeks_emp ==99
replace exp_job_years =. if exp_job_years == 98 | exp_job_years ==99
replace exp_job_months =. if exp_job_months == 98 | exp_job_months ==99
replace exp_job_weeks =. if exp_job_weeks == 98 | exp_job_weeks ==99

replace unempl =. if unempl == 8 | unempl == 9
replace unempl =0 if unempl == 5

replace layoff =. if layoff == 8 | layoff == 9
replace layoff =0 if layoff == 5

replace educ =. if educ == 99
replace weduc =. if weduc == 99

/// ASSETS ///

* Note that these are taken from the imputed values, _acc variables denote imputed entries

ren ER46938 bus_ass
ren ER46939 bus_ass_acc
gen bus_debt =.
/* ^^^ Note that Prior to 2013 the equity value is given without debt/asset split */

ren ER46942 cash
ren ER46943 cash_acc  

ren ER46950 real_estate
ren ER46951 real_estate_acc
gen real_estate_debt=.
/* ^^^ Note that Prior to 2013 the equity value is given without debt/asset split */

ren ER46954 stocks
ren ER46955 stocks_acc

ren ER46956 vehicles
ren ER46957 vehicles_acc

ren ER46960 other_ass
ren ER46961 other_ass_acc

ren ER46964 ira
ren ER46965 ira_acc

ren ER46946 other_debt
ren ER46947 other_debt_acc // Credit card charges, student loans, medical or legal bills, or loans from relatives?

/* ^^^^ Note "other debt" to be calculated manually after 2009*/

ren ER46966 home_equity
ren ER46967 home_equity_acc

ren ER46970 total_wealth
ren ER46971 total_wealth_acc

/* ^^^^ Constructed wealth variable, including equity. This imputed variable is constructed as
sum of values of seven asset types (ER46938, ER46942, ER46950, ER46954, ER46956, ER46960,
ER46964) net of debt value (ER46946) plus value of home equity (ER46966). All missing data
were assigned. */

ren ER43734 pension_plan
replace pension_plan=. if pension_plan>999999997

/* No data on wife pension pre 2011 */ 

/// EXPENDITURE ///
* First need to calculate foodstamps and utility subsidy
gen       fstmp=0
cap program drop doit
program def doit
        local i=42695
        while `i' < 42707 {
            replace ER`i'=. if ER`i'>1
        	local i=`i'+1
       }
end
qui doit
egen nm=rsum(ER42695-ER42706)
replace fstmp=ER42692              if ER42693==6 
replace fstmp=ER42692*nm           if ER42693==5 
replace fstmp=ER42692*((26/12)*nm) if ER42693==4  
replace fstmp=ER42692*((52/12)*nm) if ER42693==3
replace fstmp=.					 if ER42692>999997 | ER42693>6 | ER42693==1 
drop nm

* Note that these are taken from the imputed values, _acc variables denote imputed entries
ren ER46971A2 food_home
ren ER46971A3 food_out
ren ER46971A4 food_deliv

ren ER46971A6 mortgage  
ren ER42045 mortgage_raw1
ren ER42064 mortgage_raw2
gen mortgage_acc = 0
replace mortgage_acc = 1 if mortgage_raw1 >= 99998 | mortgage_raw2 >= 99998

ren ER46971A7 rent
ren ER42080 rent_raw
gen rent_acc = 0
replace rent_acc = 1 if rent_raw >= 99998
 
ren ER46971A8 prop_tax 
ren ER42037 prop_tax_raw
gen prop_tax_acc = 0
replace prop_tax_acc = 1 if prop_tax_raw >= 99998 

ren ER46971A9 homeinsure
ren ER42039 homeinsure_raw
gen homeinsure_acc = 0
replace homeinsure_acc = 1 if homeinsure_raw>=9998

ren ER46971B1 utilities
ren ER46971B2 gas_raw
ren ER46971B3 electric_raw
ren ER46971B4 water_raw
ren ER46971B5 otheru_raw
gen utilities_acc = 0
replace utilities_acc = 1 if gas_raw >= 9999999 | electric_raw >= 9999999 | water_raw >= 9999999 | otheru_raw >= 9999999 

ren ER46971B6 telecom
ren ER42120 telecom_raw
gen telecom_acc = 0
replace telecom_acc = 1 if telecom_raw>=9998

ren ER46971B8 vehicle_loan
ren ER42748 vehicle_loan_raw1 
ren ER42771 vehicle_loan_raw2
ren ER42794 vehicle_loan_raw3
gen vehicle_loan_acc = 0
replace vehicle_loan_acc = 1 if vehicle_loan_raw1 >=999998 | vehicle_loan_raw2 >=999998 | vehicle_loan_raw3 >=999998 

ren ER46971B9 vehicle_down
ren ER42745 vehicle_down_raw1
ren ER42768 vehicle_down_raw2
ren ER42791 vehicle_down_raw3
gen vehicle_down_acc = 0
replace vehicle_down_acc = 1 if vehicle_down_raw1 >= 999998 | vehicle_down_raw2 >= 999998 | vehicle_down_raw3 >= 999998 

ren ER46971C1 vehicle_lease
ren ER42752 vehicle_lease_raw1
ren ER42775 vehicle_lease_raw2
ren ER42798 vehicle_lease_raw3
gen vehicle_lease_acc = 0
replace vehicle_lease_acc = 1 if vehicle_lease_raw1 >=999998 | vehicle_lease_raw2 >=999998 |vehicle_lease_raw3 >=999998 

ren ER46971C2 auto_ins
ren ER42803 auto_ins_raw
gen auto_ins_acc = 0
replace auto_ins_acc = 1 if auto_ins_raw >= 999998

ren ER46971C3 other_vehicle
ren ER42806 other_vehicle_raw
gen other_vehicle_acc = 0
replace other_vehicle_acc  = 1 if other_vehicle_raw >=999998

ren ER46971C5 gasoline
ren ER42808 gasoline_raw
gen gasoline_acc = 0
replace gasoline_acc = 1 if gasoline_raw >= 99998

ren ER46971C4 auto_repair
ren ER42807 auto_repair_raw
gen auto_repair_acc = 0
replace auto_repair_acc = 1 if auto_repair_raw >=99998

ren ER46971C6 parking
ren ER42809 parking_raw
gen parking_acc = 0
replace parking_acc = 1 if parking_raw >= 99998 

ren ER46971C7 bus_train
ren ER42810 bus_train_raw
gen bus_train_acc = 0
replace bus_train_acc = 1 if bus_train_raw >= 99998

ren ER46971C8 taxi
ren ER42811 taxi_raw
gen taxi_acc = 0
replace taxi_acc = 1 if taxi_raw >= 99998

ren ER46971C9 other_trans
ren ER42812 other_trans_raw
gen other_trans_acc = 0
replace other_trans_acc = 1 if other_trans_raw >=99998

ren ER46971D1 education
ren ER42814 education_raw1
ren ER42816 education_raw2
gen education_acc = 0
replace education_acc = 1 if education_raw1 >= 999998 | education_raw2 >= 999998

ren ER46971D2 childcare
ren ER42652 childcare_raw
gen childcare_acc = 0
replace childcare_raw = 1 if childcare_raw >= 999998

ren ER46971D4 inst_medical
ren ER46387 inst_medical_raw
gen inst_medical_acc = 0
replace inst_medical_acc = 1 if inst_medical_raw >= 999998

ren ER46971D5 doctor
ren ER46393 doctor_raw
gen doctor_acc = 0
replace doctor_acc = 1 if doctor_raw >= 999998

ren ER46971D6 prescription
ren ER46399 prescription_raw
gen prescription_acc = 0
replace prescription_acc = 1 if prescription_raw >= 999998

ren ER46971D7 health_ins
ren ER46383 health_ins_raw
gen health_ins_acc = 0
replace health_ins_acc = 1 if health_ins_raw >= 99998

ren ER46971D8 home_repair
ren ER42817 home_repair_raw
gen home_repair_acc = 0
replace home_repair_acc = 1 if home_repair_raw >= 999998

ren ER46971D9 home_furnish
ren ER42822 home_furnish_raw
gen home_furnish_acc = 0
replace home_furnish_acc = 1 if home_furnish_raw >= 999998

ren ER46971E1 clothing
ren ER42827 clothing_raw
gen clothing_acc = 0
replace clothing_acc = 1 if clothing_raw >= 999998

ren ER46971E2 trips
ren ER42832 trips_raw
gen trips_acc = 0
replace trips_acc = 1 if trips_raw >= 999998

ren ER46971E3 other_rec
ren ER42837 other_rec_raw
gen other_rec_acc = 0
replace other_rec_acc = 1 if other_rec >= 999998


/// INCOME ///

ren ER46935 y
gen truncy=y<1
replace y=1 if y<1

ren ER46851 tyhw
ren ER46900 trhw
ren ER46907 tyoth
ren ER46928 troth
ren ER46929 soc_sec
ren ER46931 soc_secw
ren ER46933 soc_secyoth



* Asset Income *

ren ER46830 renty
ren ER46831 renty_acc
ren ER46843 rentyw
ren ER46844 rentyw_acc

ren ER46832 divi
ren ER46833 divi_acc
ren ER46845 diviw
ren ER46846 diviw_acc

ren ER46834 interest
ren ER46835 interest_acc
ren ER46847 interestw
ren ER46848 interestw_acc

ren ER46836 trustinc
ren ER46837 trustinc_acc
ren ER46849 trustincw
ren ER46850 trustincw_acc

ren ER46809 businc
ren ER46839 busincw

* Generate labor and asset part of farm income *
gen farmlabor=0
replace farmlabor=0.5*ER46806 			if ER46806>0
replace farmlabor=. 					if ER46806>9999997

gen farmasset=0
replace farmasset=0.5*ER46806 			if ER46806>0
replace farmasset=ER46806 				if ER46806<0
replace farmasset=. 					if ER46806>9999997

* Hours Worked *
ren ER46767 hours
ren ER46788 hourw

* Income *
ren ER46829 wages
ren ER46841 wagesw
ren ER46808 y_business
ren ER46838 yw_business

egen ly = rsum(wages farmlabor y_business) 
egen wly= rsum(wagesw yw_business)


/// Prep variables to save ///

/// Prep variables to save ///

#delimit;

keep 


id fsize age sex agew kids marit race* wrace*
outkid house smsa fchg state weight foodstamps *educ newborn

bus_ass bus_debt cash real_estate real_estate_debt stocks vehicles other_ass
ira other_debt home_equity total_wealth

fstmp food_home food_out food_deliv mortgage rent prop_tax homeinsure utilities telecom
vehicle_loan vehicle_down vehicle_lease auto_ins other_vehicle gasoline auto_repair
parking bus_train taxi other_trans education childcare inst_medical doctor prescription
health_ins home_repair home_furnish clothing trips other_rec

y truncy tyhw trhw tyoth troth soc_sec soc_secw soc_secyoth
farmlabor farmasset hours* wages* y_business yw_business ly wly weeks_emp* hourw
exp_job_years exp_job_months exp_job_weeks empst* wempst* unempl layoff
rent* divi* interest* trust* bus*
; 

#delimit cr


gen year=2009

save "$temp\HH_2009.dta", replace

********************************************************************************************
* 2007
********************************************************************************************

use "$data\2007_H", clear

* Merge with the supplementary wealth file
ren ER36002 S801
merge 1:1 S801 using "$data\2007_W.dta"

ren S801 id // Corrseponds to wealth file identifier
ren ER36016 fsize
ren ER36017 age
ren ER36018 sex
ren ER36019 agew
ren ER36020 kids
ren ER36021 newborn
ren ER36023 marit
ren ER40565 race1
ren ER40566 race2
ren ER40567 race3
ren ER40568 race4
ren ER40472 wrace1
ren ER40473 wrace2
ren ER40474 wrace3
ren ER40475 wrace4
ren ER37536 outkid
ren ER36029 house
ren ER41034 smsa 
ren ER36007 fchg
ren ER36003 state
ren ER41069 weight 
ren ER36109 empst1
ren ER36110 empst2
ren ER36111 empst3
ren ER36367 wempst1
ren ER36368 wempst2
ren ER36369 wempst3
ren ER36672 foodstamps
ren ER41037 educ
ren ER41038 weduc
ren ER40873 weeks_emp
ren ER40884 weeks_empw
ren ER36165 exp_job_years
ren ER36166 exp_job_months
ren ER36167 exp_job_weeks
ren ER36311 unempl 
ren ER36304 layoff  

replace weeks_emp =. if weeks_emp == 98 | weeks_emp ==99
replace exp_job_years =. if exp_job_years == 98 | exp_job_years ==99
replace exp_job_months =. if exp_job_months == 98 | exp_job_months ==99
replace exp_job_weeks =. if exp_job_weeks == 98 | exp_job_weeks ==99

replace unempl =. if unempl == 8 | unempl == 9
replace unempl =0 if unempl == 5

replace layoff =. if layoff == 8 | layoff == 9
replace layoff =0 if layoff == 5

replace educ =. if educ == 99
replace weduc =. if weduc == 99

/// ASSETS ///

* Note that these are taken from the imputed values, _acc variables denote imputed entries
* Prior to 2009 these taken from separate wealth file

ren S803 bus_ass
ren S803A bus_ass_acc
gen bus_debt=.

/* ^^^ Note that Prior to 2013 the equity value is given without debt/asset split */

ren S805 cash
ren S805A cash_acc  

ren S809 real_estate
ren S809A real_estate_acc
gen real_estate_debt=.

/* ^^^ Note that Prior to 2013 the equity value is given without debt/asset split */

ren S811 stocks
ren S811A stocks_acc

ren S813 vehicles
ren S813A vehicles_acc

ren S815 other_ass
ren S815A other_ass_acc

ren S819 ira
ren S819A ira_acc

ren S807 other_debt
ren S807A other_debt_acc // Credit card charges, student loans, medical or legal bills, or loans from relatives?

/* ^^^^ Note "other debt" to be calculated manually after 2009*/

ren S820 home_equity
ren S820A home_equity_acc

ren S817 total_wealth
ren S817A total_wealth_acc

/* ^^^^ Constructed wealth variable, including equity. This imputed variable is constructed as
sum of values of seven asset types (ER46938, ER46942, ER46950, ER46954, ER46956, ER46960,
ER46964) net of debt value (ER46946) plus value of home equity (ER46966). All missing data
were assigned. */

ren ER37761 pension_plan
replace pension_plan=. if pension_plan>999999997

/* No data on wife pension pre 2011 */ 

/// EXPENDITURE ///

* Foodstamps and utility
gen       fstmp=0
cap program drop doit
program def doit
        local i=36676
        while `i' < 36688 {
            replace ER`i'=. if ER`i'>1
        	local i=`i'+1
       }
end
qui doit
egen nm=rsum(ER36676-ER36687)
replace fstmp=ER36673              if ER36674==6 
replace fstmp=ER36673*nm           if ER36674==5 
replace fstmp=ER36673*((26/12)*nm) if ER36674==4  
replace fstmp=ER36673*((52/12)*nm) if ER36674==3
replace fstmp=.					 if ER36673>999997 | ER36674>6 | ER36674==1 
drop nm


* Note that these are taken from the imputed values, _acc variables denote imputed entries
ren ER41027A2 food_home
ren ER41027A3 food_out
ren ER41027A4 food_deliv

/* ^^^^ STILL NEED TO DEAL WITH FOOD STAMPS / UTILITY SUBSIDY */

ren ER41027A6 mortgage  
ren ER36044 mortgage_raw1
ren ER36056 mortgage_raw2
gen mortgage_acc = 0
replace mortgage_acc = 1 if mortgage_raw1 >= 99998 | mortgage_raw2 >= 99998

ren ER41027A7 rent
ren ER36065 rent_raw
gen rent_acc = 0
replace rent_acc = 1 if rent_raw >= 99998
 
ren ER41027A8 prop_tax 
ren ER36036 prop_tax_raw
gen prop_tax_acc = 0
replace prop_tax_acc = 1 if prop_tax_raw >= 99998 

ren ER41027A9 homeinsure
ren ER36038 homeinsure_raw
gen homeinsure_acc = 0
replace homeinsure_acc = 1 if homeinsure_raw>=9998

ren ER41027B1 utilities
ren ER41027B2 gas_raw
ren ER41027B3 electric_raw
ren ER41027B4 water_raw
ren ER41027B5 otheru_raw
gen utilities_acc = 0
replace utilities_acc = 1 if gas_raw >= 9999999 | electric_raw >= 9999999 | water_raw >= 9999999 | otheru_raw >= 9999999 

ren ER41027B6 telecom
ren ER36091 telecom_raw
gen telecom_acc = 0
replace telecom_acc = 1 if telecom_raw>=9998

ren ER41027B8 vehicle_loan
ren ER36747 vehicle_loan_raw1 
ren ER36775 vehicle_loan_raw2
ren ER36803 vehicle_loan_raw3
gen vehicle_loan_acc = 0
replace vehicle_loan_acc = 1 if vehicle_loan_raw1 >=999998 | vehicle_loan_raw2 >=999998 | vehicle_loan_raw3 >=999998 

ren ER41027B9 vehicle_down
ren ER36744 vehicle_down_raw1
ren ER36772 vehicle_down_raw2
ren ER36800 vehicle_down_raw3
gen vehicle_down_acc = 0
replace vehicle_down_acc = 1 if vehicle_down_raw1 >= 999998 | vehicle_down_raw2 >= 999998 | vehicle_down_raw3 >= 999998 

ren ER41027C1 vehicle_lease
ren ER36751 vehicle_lease_raw1
ren ER36779 vehicle_lease_raw2
ren ER36807 vehicle_lease_raw3
gen vehicle_lease_acc = 0
replace vehicle_lease_acc = 1 if vehicle_lease_raw1 >=999998 | vehicle_lease_raw2 >=999998 |vehicle_lease_raw3 >=999998 

ren ER41027C2 auto_ins
ren ER36812 auto_ins_raw
gen auto_ins_acc = 0
replace auto_ins_acc = 1 if auto_ins_raw >= 999998

ren ER41027C3 other_vehicle
ren ER36815 other_vehicle_raw
gen other_vehicle_acc = 0
replace other_vehicle_acc  = 1 if other_vehicle_raw >=999998

ren ER41027C5 gasoline
ren ER36817 gasoline_raw
gen gasoline_acc = 0
replace gasoline_acc = 1 if gasoline_raw >= 99998

ren ER41027C4 auto_repair
ren ER36816 auto_repair_raw
gen auto_repair_acc = 0
replace auto_repair_acc = 1 if auto_repair_raw >=99998

ren ER41027C6 parking
ren ER36818 parking_raw
gen parking_acc = 0
replace parking_acc = 1 if parking_raw >= 99998 

ren ER41027C7 bus_train
ren ER36819 bus_train_raw
gen bus_train_acc = 0
replace bus_train_acc = 1 if bus_train_raw >= 99998

ren ER41027C8 taxi
ren ER36820 taxi_raw
gen taxi_acc = 0
replace taxi_acc = 1 if taxi_raw >= 99998

ren ER41027C9 other_trans
ren ER36821 other_trans_raw
gen other_trans_acc = 0
replace other_trans_acc = 1 if other_trans_raw >=99998

ren ER41027D1 education
ren ER36823 education_raw1
ren ER36825 education_raw2
gen education_acc = 0
replace education_acc = 1 if education_raw1 >= 999998 | education_raw2 >= 999998

ren ER41027D2 childcare
ren ER36633 childcare_raw
gen childcare_acc = 0
replace childcare_raw = 1 if childcare_raw >= 999998

ren ER41027D4 inst_medical
ren ER40414 inst_medical_raw
gen inst_medical_acc = 0
replace inst_medical_acc = 1 if inst_medical_raw >= 999998

ren ER41027D5 doctor
ren ER40420 doctor_raw
gen doctor_acc = 0
replace doctor_acc = 1 if doctor_raw >= 999998

ren ER41027D6 prescription
ren ER40426 prescription_raw
gen prescription_acc = 0
replace prescription_acc = 1 if prescription_raw >= 999998

ren ER41027D7 health_ins
ren ER40410 health_ins_raw
gen health_ins_acc = 0
replace health_ins_acc = 1 if health_ins_raw >= 99998

ren ER41027D8 home_repair
ren ER36826 home_repair_raw
gen home_repair_acc = 0
replace home_repair_acc = 1 if home_repair_raw >= 999998

ren ER41027D9 home_furnish
ren ER36831 home_furnish_raw
gen home_furnish_acc = 0
replace home_furnish_acc = 1 if home_furnish_raw >= 999998

ren ER41027E1 clothing
ren ER36836 clothing_raw
gen clothing_acc = 0
replace clothing_acc = 1 if clothing_raw >= 999998

ren ER41027E2 trips
ren ER36841 trips_raw
gen trips_acc = 0
replace trips_acc = 1 if trips_raw >= 999998

ren ER41027E3 other_rec
ren ER36846 other_rec_raw
gen other_rec_acc = 0
replace other_rec_acc = 1 if other_rec >= 999998


/// INCOME ///

ren ER41027 y
gen truncy=y<1
replace y=1 if y<1

ren ER40943 tyhw
ren ER40992 trhw
ren ER40999 tyoth
ren ER41020 troth
ren ER41021 soc_sec
ren ER41023 soc_secw
ren ER41025 soc_secyoth



* Asset Income *

ren ER40922 renty
ren ER40923 renty_acc
ren ER40935 rentyw
ren ER40936 rentyw_acc

ren ER40924 divi
ren ER40925 divi_acc
ren ER40937 diviw
ren ER40938 diviw_acc

ren ER40926 interest
ren ER40927 interest_acc
ren ER40939 interestw
ren ER40940 interestw_acc

ren ER40928 trustinc
ren ER40929 trustinc_acc
ren ER40941 trustincw
ren ER40942 trustincw_acc

ren ER40901 businc
ren ER40931 busincw


* Generate labor and asset part of farm income *
gen farmlabor=0
replace farmlabor=0.5*ER40898 			if ER40898>0
replace farmlabor=. 					if ER40898>9999997

gen farmasset=0
replace farmasset=0.5*ER40898 			if ER40898>0
replace farmasset=ER40898 				if ER40898<0
replace farmasset=. 					if ER40898>9999997

* Hours Worked *
ren ER40876 hours
ren ER40887 hourw

* Income *
ren ER40921 wages
ren ER40933 wagesw
ren ER40900 y_business
ren ER40930 yw_business

egen ly = rsum(wages farmlabor y_business) 
egen wly= rsum(wagesw yw_business)

/// Prep variables to save ///

#delimit;

keep 

id fsize age sex agew kids marit race* wrace*
outkid house smsa fchg state weight foodstamps *educ newborn

bus_ass bus_debt cash real_estate real_estate_debt stocks vehicles other_ass
ira other_debt home_equity total_wealth

fstmp food_home food_out food_deliv mortgage rent prop_tax homeinsure utilities telecom
vehicle_loan vehicle_down vehicle_lease auto_ins other_vehicle gasoline auto_repair
parking bus_train taxi other_trans education childcare inst_medical doctor prescription
health_ins home_repair home_furnish clothing trips other_rec

y truncy tyhw trhw tyoth troth soc_sec soc_secw soc_secyoth
farmlabor farmasset hours* wages* y_business yw_business ly wly weeks_emp* hourw
exp_job_years exp_job_months exp_job_weeks empst* wempst* unempl layoff
rent* divi* interest* trust* bus*
; 

#delimit cr


gen year=2007

save "$temp\HH_2007.dta", replace


********************************************************************************************
* 2005
********************************************************************************************

use "$data\2005_H", clear

* Merge with the supplementary wealth file
ren ER25002 S701
merge 1:1 S701 using "$data\2005_W.dta"

ren S701 id // Corrseponds to wealth file identifier
ren ER25016 fsize
ren ER25017 age
ren ER25018 sex
ren ER25019 agew
ren ER25020 kids
ren ER25021 newborn
ren ER25023 marit
ren ER27393 race1
ren ER27394 race2
ren ER27395 race3
ren ER27396 race4
ren ER27297 wrace1
ren ER27298 wrace2
ren ER27299 wrace3
ren ER27300 wrace4
ren ER26518 outkid
ren ER25029 house
ren ER28044 smsa 
ren ER25007 fchg
ren ER25003 state
ren ER28078 weight 
ren ER25104 empst1
ren ER25105 empst2
ren ER25106 empst3
ren ER25362 wempst1
ren ER25363 wempst2
ren ER25364 wempst3
ren ER25654 foodstamps
ren ER28047 educ
ren ER28048 weduc
ren ER27883 weeks_emp
ren ER27894 weeks_empw
ren ER25160 exp_job_years
ren ER25161 exp_job_months
ren ER25162 exp_job_weeks
ren ER25306 unempl 
ren ER25299 layoff  

replace weeks_emp =. if weeks_emp == 98 | weeks_emp ==99
replace exp_job_years =. if exp_job_years == 98 | exp_job_years ==99
replace exp_job_months =. if exp_job_months == 98 | exp_job_months ==99
replace exp_job_weeks =. if exp_job_weeks == 98 | exp_job_weeks ==99

replace educ =. if educ == 99
replace weduc =. if weduc == 99

replace unempl =. if unempl == 8 | unempl == 9
replace unempl =0 if unempl == 5

replace layoff =. if layoff == 8 | layoff == 9
replace layoff =0 if layoff == 5


replace educ =. if educ == 99
replace weduc =. if weduc == 99

/// ASSETS ///

* Note that these are taken from the imputed values, _acc variables denote imputed entries
* Prior to 2009 these taken from separate wealth file

ren S703 bus_ass
ren S703A bus_ass_acc
gen bus_debt=.
/* ^^^ Note that Prior to 2013 the equity value is given without debt/asset split */

ren S705 cash
ren S705A cash_acc  

ren S709 real_estate
ren S709A real_estate_acc
gen real_estate_debt=.
/* ^^^ Note that Prior to 2013 the equity value is given without debt/asset split */

ren S711 stocks
ren S711A stocks_acc

ren S713 vehicles
ren S713A vehicles_acc

ren S715 other_ass
ren S715A other_ass_acc

ren S719 ira
ren S719A ira_acc

ren S707 other_debt
ren S707A other_debt_acc // Credit card charges, student loans, medical or legal bills, or loans from relatives?

/* ^^^^ Note "other debt" to be calculated manually after 2009*/

ren S720 home_equity
ren S720A home_equity_acc

ren S717 total_wealth
ren S717A total_wealth_acc

/* ^^^^ This variable is constructed as sum of values of seven asset types (S703, S705, S709,
S711, S713, S715, S719) net of debt value (S707) plus value of home equity. */

ren ER26725 pension_plan
replace pension_plan=. if pension_plan>999999997

/* No data on wife pension pre 2011 */ 

/// EXPENDITURE ///

gen       fstmp=0
cap program drop doit
program def doit
        local i=25658
        while `i' < 25670 {
            replace ER`i'=. if ER`i'>1
        	local i=`i'+1
       }
end
qui doit
egen nm=rsum(ER25658-ER25669)
replace fstmp=ER25655              if ER25656==6 
replace fstmp=ER25655*nm           if ER25656==5 
replace fstmp=ER25655*((26/12)*nm) if ER25656==4  
replace fstmp=ER25655*((52/12)*nm) if ER25656==3
replace fstmp=.					 if ER25655>999997  | ER25656>6 
drop nm

* Note that these are taken from the imputed values, _acc variables denote imputed entries

ren ER28037A2 food_home
ren ER28037A3 food_out
ren ER28037A4 food_deliv

/* ^^^^ STILL NEED TO DEAL WITH FOOD STAMPS / UTILITY SUBSIDY */

ren ER28037A6 mortgage  
ren ER25044 mortgage_raw1
ren ER25055 mortgage_raw2
gen mortgage_acc = 0
replace mortgage_acc = 1 if mortgage_raw1 >= 99998 | mortgage_raw2 >= 99998

ren ER28037A7 rent
ren ER25063 rent_raw
gen rent_acc = 0
replace rent_acc = 1 if rent_raw >= 99998
 
ren ER28037A8 prop_tax 
ren ER25036 prop_tax_raw
gen prop_tax_acc = 0
replace prop_tax_acc = 1 if prop_tax_raw >= 99998 

ren ER28037A9 homeinsure
ren ER25038 homeinsure_raw
gen homeinsure_acc = 0
replace homeinsure_acc = 1 if homeinsure_raw>=9998

ren ER28037B1 utilities
ren ER28037B2 gas_raw
ren ER28037B3 electric_raw
ren ER28037B4 water_raw
ren ER28037B5 otheru_raw
gen utilities_acc = 0
replace utilities_acc = 1 if gas_raw >= 9999999 | electric_raw >= 9999999 | water_raw >= 9999999 | otheru_raw >= 9999999 

ren ER28037B6 telecom
ren ER25086 telecom_raw
gen telecom_acc = 0
replace telecom_acc = 1 if telecom_raw>=9998

ren ER28037B8 vehicle_loan
ren ER25729 vehicle_loan_raw1 
ren ER25757 vehicle_loan_raw2
ren ER25785 vehicle_loan_raw3
gen vehicle_loan_acc = 0
replace vehicle_loan_acc = 1 if vehicle_loan_raw1 >=999998 | vehicle_loan_raw2 >=999998 | vehicle_loan_raw3 >=999998 

ren ER28037B9 vehicle_down
ren ER25726 vehicle_down_raw1
ren ER25754 vehicle_down_raw2
ren ER25782 vehicle_down_raw3
gen vehicle_down_acc = 0
replace vehicle_down_acc = 1 if vehicle_down_raw1 >= 999998 | vehicle_down_raw2 >= 999998 | vehicle_down_raw3 >= 999998 

ren ER28037C1 vehicle_lease
ren ER25733 vehicle_lease_raw1
ren ER25761 vehicle_lease_raw2
ren ER25789 vehicle_lease_raw3
gen vehicle_lease_acc = 0
replace vehicle_lease_acc = 1 if vehicle_lease_raw1 >=999998 | vehicle_lease_raw2 >=999998 |vehicle_lease_raw3 >=999998 

ren ER28037C2 auto_ins
ren ER25794 auto_ins_raw
gen auto_ins_acc = 0
replace auto_ins_acc = 1 if auto_ins_raw >= 999998

ren ER28037C3 other_vehicle
ren ER25797 other_vehicle_raw
gen other_vehicle_acc = 0
replace other_vehicle_acc  = 1 if other_vehicle_raw >=999998

ren ER28037C5 gasoline
ren ER25799 gasoline_raw
gen gasoline_acc = 0
replace gasoline_acc = 1 if gasoline_raw >= 99998

ren ER28037C4 auto_repair
ren ER25798 auto_repair_raw
gen auto_repair_acc = 0
replace auto_repair_acc = 1 if auto_repair_raw >=99998

ren ER28037C6 parking
ren ER25800 parking_raw
gen parking_acc = 0
replace parking_acc = 1 if parking_raw >= 99998 

ren ER28037C7 bus_train
ren ER25801 bus_train_raw
gen bus_train_acc = 0
replace bus_train_acc = 1 if bus_train_raw >= 99998

ren ER28037C8 taxi
ren ER25802 taxi_raw
gen taxi_acc = 0
replace taxi_acc = 1 if taxi_raw >= 99998

ren ER28037C9 other_trans
ren ER25803 other_trans_raw
gen other_trans_acc = 0
replace other_trans_acc = 1 if other_trans_raw >=99998

ren ER28037D1 education
ren ER25805 education_raw1
ren ER25807 education_raw2
gen education_acc = 0
replace education_acc = 1 if education_raw1 >= 999998 | education_raw2 >= 999998

ren ER28037D2 childcare
ren ER25628 childcare_raw
gen childcare_acc = 0
replace childcare_raw = 1 if childcare_raw >= 999998

ren ER28037D4 inst_medical
ren ER27239 inst_medical_raw
gen inst_medical_acc = 0
replace inst_medical_acc = 1 if inst_medical_raw >= 999998

ren ER28037D5 doctor
ren ER27245 doctor_raw
gen doctor_acc = 0
replace doctor_acc = 1 if doctor_raw >= 999998

ren ER28037D6 prescription
ren ER27251 prescription_raw
gen prescription_acc = 0
replace prescription_acc = 1 if prescription_raw >= 999998

ren ER28037D7 health_ins
ren ER27238 health_ins_raw
gen health_ins_acc = 0
replace health_ins_acc = 1 if health_ins_raw >= 99998

ren ER28037D8 home_repair
ren ER25808 home_repair_raw
gen home_repair_acc = 0
replace home_repair_acc = 1 if home_repair_raw >= 999998

ren ER28037D9 home_furnish
ren ER25813 home_furnish_raw
gen home_furnish_acc = 0
replace home_furnish_acc = 1 if home_furnish_raw >= 999998

ren ER28037E1 clothing
ren ER25818 clothing_raw
gen clothing_acc = 0
replace clothing_acc = 1 if clothing_raw >= 999998

ren ER28037E2 trips
ren ER25823 trips_raw
gen trips_acc = 0
replace trips_acc = 1 if trips_raw >= 999998

ren ER28037E3 other_rec
ren ER25828 other_rec_raw
gen other_rec_acc = 0
replace other_rec_acc = 1 if other_rec >= 999998


/// INCOME ///

ren ER28037 y
gen truncy=y<1
replace y=1 if y<1

ren ER27953 tyhw
ren ER28002 trhw
ren ER28009 tyoth
ren ER28030 troth
ren ER28031 soc_sec
ren ER28033 soc_secw
ren ER28035 soc_secyoth



* Asset Income *

ren ER27932 renty
ren ER27933 renty_acc
ren ER27945 rentyw
ren ER27946 rentyw_acc

ren ER27934 divi
ren ER27935 divi_acc
ren ER27947 diviw
ren ER27948 diviw_acc

ren ER27936 interest
ren ER27937 interest_acc
ren ER27949 interestw
ren ER27950 interestw_acc

ren ER27938 trustinc
ren ER27939 trustinc_acc
ren ER27951 trustincw
ren ER27952 trustincw_acc

ren ER27911 businc
ren ER27941 busincw

* Generate labor and asset part of farm income *
gen farmlabor=0
replace farmlabor=0.5*ER27908 			if ER27908>0
replace farmlabor=. 					if ER27908>9999997

gen farmasset=0
replace farmasset=0.5*ER27908 			if ER27908>0
replace farmasset=ER27908 				if ER27908<0
replace farmasset=. 					if ER27908>9999997

* Hours Worked *
ren ER27886 hours
ren ER27897 hourw

* Income *
ren ER27931 wages
ren ER27943 wagesw
ren ER27910 y_business
ren ER27940 yw_business

egen ly = rsum(wages farmlabor y_business) 
egen wly= rsum(wagesw yw_business)


/// Prep variables to save ///

#delimit;

keep 


id fsize age sex agew kids marit race* wrace*
outkid house smsa fchg state weight foodstamps *educ newborn

bus_ass bus_debt cash real_estate real_estate_debt stocks vehicles other_ass
ira other_debt home_equity total_wealth

fstmp food_home food_out food_deliv mortgage rent prop_tax homeinsure utilities telecom
vehicle_loan vehicle_down vehicle_lease auto_ins other_vehicle gasoline auto_repair
parking bus_train taxi other_trans education childcare inst_medical doctor prescription
health_ins home_repair home_furnish clothing trips other_rec

y truncy tyhw trhw tyoth troth soc_sec soc_secw soc_secyoth
farmlabor farmasset hours* wages* y_business yw_business ly wly weeks_emp* hourw
exp_job_years exp_job_months exp_job_weeks empst* wempst* unempl layoff
rent* divi* interest* trust* bus*
; 

#delimit cr


gen year=2005

save "$temp\HH_2005.dta", replace


********************************************************************************************
* 2003
********************************************************************************************

use "$data\2003_H", clear

* Merge with the supplementary wealth file
ren ER21002 S601
merge 1:1 S601 using "$data\2003_W.dta"

ren S601 id // Corrseponds to wealth file identifier
ren ER21016 fsize
ren ER21017 age
ren ER21018 sex
ren ER21019 agew
ren ER21020 kids
ren ER21021 newborn
ren ER21023 marit
ren ER23426 race1
ren ER23427 race2
ren ER23428 race3
ren ER23429 race4
ren ER23334 wrace1
ren ER23335 wrace2
ren ER23336 wrace3
ren ER23337 wrace4
ren ER22537 outkid
ren ER21043 house
ren ER24145 smsa 
ren ER21007 fchg
ren ER21003 state
ren ER24179 weight 
ren ER21123 empst1
ren ER21124 empst2
ren ER21125 empst3
ren ER21373 wempst1
ren ER21374 wempst2
ren ER21375 wempst3
ren ER21652 foodstamps
ren ER24077 weeks_emp
ren ER24088 weeks_empw
ren ER21171 exp_job_years
ren ER21172 exp_job_months
ren ER21173 exp_job_weeks
ren ER21317 unempl 
ren ER21310 layoff 

replace weeks_emp =. if weeks_emp == 98 | weeks_emp ==99
replace exp_job_years =. if exp_job_years == 98 | exp_job_years ==99
replace exp_job_months =. if exp_job_months == 98 | exp_job_months ==99
replace exp_job_weeks =. if exp_job_weeks == 98 | exp_job_weeks ==99

ren ER24148 educ
ren ER24149 weduc

replace unempl =. if unempl == 8 | unempl == 9
replace unempl =0 if unempl == 5

replace layoff =. if layoff == 8 | layoff == 9
replace layoff =0 if layoff == 5

replace educ =. if educ == 99
replace weduc =. if weduc == 99

/// ASSETS ///

* Note that these are taken from the imputed values, _acc variables denote imputed entries
* Prior to 2009 these taken from separate wealth file

ren S603 bus_ass
ren S603A bus_ass_acc
gen bus_debt=.
/* ^^^ Note that Prior to 2013 the equity value is given without debt/asset split */

ren S605 cash
ren S605A cash_acc  

ren S609 real_estate
ren S609A real_estate_acc
gen real_estate_debt=.
/* ^^^ Note that Prior to 2013 the equity value is given without debt/asset split */

ren S611 stocks
ren S611A stocks_acc

ren S613 vehicles
ren S613A vehicles_acc

ren S615 other_ass
ren S615A other_ass_acc

ren S619 ira
ren S619A ira_acc

ren S607 other_debt
ren S607A other_debt_acc // Credit card charges, student loans, medical or legal bills, or loans from relatives?

/* ^^^^ Note "other debt" to be calculated manually after 2009*/

ren S620 home_equity
ren S620A home_equity_acc

ren S617 total_wealth
ren S617A total_wealth_acc

/* ^^^^ This variable is constructed as sum of values of seven asset types (S703, S705, S709,
S711, S713, S715, S719) net of debt value (S707) plus value of home equity. */

ren ER22744 pension_plan
replace pension_plan=. if pension_plan>999999997

/* No data on wife pension pre 2011 */ 

/// EXPENDITURE ///

gen       fstmp=0
cap program drop doit
program def doit
        local i=21656
        while `i' < 21668 {
            replace ER`i'=. if ER`i'>1
        	local i=`i'+1
       }
end
qui doit
egen nm=rsum(ER21656-ER21668)
replace fstmp=ER21653              if ER21654==6 
replace fstmp=ER21653*nm           if ER21654==5 
replace fstmp=ER21653*((26/12)*nm) if ER21654==4  
replace fstmp=ER21653*((52/12)*nm) if ER21654==3
replace fstmp=.					 if ER21653>999997  | ER21654>6 
drop nm

* Note that these are taken from the imputed values, _acc variables denote imputed entries

ren ER24138A2 food_home
ren ER24138A3 food_out
ren ER24138A4 food_deliv

ren ER24138A6 mortgage  
/* ren ER25044 mortgage_raw1
ren ER25055 mortgage_raw2
gen mortgage_acc = 0
replace mortgage_acc = 1 if mortgage_raw1 >= 99998 | mortgage_raw2 >= 99998 */

ren ER24138A7 rent
/* ren ER25063 rent_raw
gen rent_acc = 0
replace rent_acc = 1 if rent_raw >= 99998 */
 
ren ER24138A8 prop_tax 
/* ren ER25036 prop_tax_raw
gen prop_tax_acc = 0
replace prop_tax_acc = 1 if prop_tax_raw >= 99998 */

ren ER24138A9 homeinsure
/* ren ER25038 homeinsure_raw
gen homeinsure_acc = 0
replace homeinsure_acc = 1 if homeinsure_raw>=9998 */

ren ER24138B1 utilities
ren ER24138B2 gas_raw
ren ER24138B3 electric_raw
ren ER24138B4 water_raw
ren ER24138B5 otheru_raw
gen utilities_acc = 0
replace utilities_acc = 1 if gas_raw >= 9999999 | electric_raw >= 9999999 | water_raw >= 9999999 | otheru_raw >= 9999999 

/*

NO TELECOMS PRE-2005

ren ER28037B6 telecom
ren ER25086 telecom_raw
gen telecom_acc = 0
replace telecom_acc = 1 if telecom_raw>=9998

*/

ren ER24138B7 vehicle_loan
/*
ren ER25729 vehicle_loan_raw1 
ren ER25757 vehicle_loan_raw2
ren ER25785 vehicle_loan_raw3
gen vehicle_loan_acc = 0
replace vehicle_loan_acc = 1 if vehicle_loan_raw1 >=999998 | vehicle_loan_raw2 >=999998 | vehicle_loan_raw3 >=999998 
*/

ren ER24138B8 vehicle_down
/*ren ER25726 vehicle_down_raw1
ren ER25754 vehicle_down_raw2
ren ER25782 vehicle_down_raw3
gen vehicle_down_acc = 0
replace vehicle_down_acc = 1 if vehicle_down_raw1 >= 999998 | vehicle_down_raw2 >= 999998 | vehicle_down_raw3 >= 999998 
*/

ren ER24138B9 vehicle_lease
/*
ren ER25733 vehicle_lease_raw1
ren ER25761 vehicle_lease_raw2
ren ER25789 vehicle_lease_raw3
gen vehicle_lease_acc = 0
replace vehicle_lease_acc = 1 if vehicle_lease_raw1 >=999998 | vehicle_lease_raw2 >=999998 |vehicle_lease_raw3 >=999998 
*/

ren ER24138C1 auto_ins
/*
ren ER25794 auto_ins_raw
gen auto_ins_acc = 0
replace auto_ins_acc = 1 if auto_ins_raw >= 999998
*/

ren ER24138C2 other_vehicle
/* ren ER25797 other_vehicle_raw
gen other_vehicle_acc = 0
replace other_vehicle_acc  = 1 if other_vehicle_raw >=999998
*/

ren ER24138C4 gasoline
/*
ren ER25799 gasoline_raw
gen gasoline_acc = 0
replace gasoline_acc = 1 if gasoline_raw >= 99998
*/

ren ER24138C3 auto_repair
/*
ren ER25798 auto_repair_raw
gen auto_repair_acc = 0
replace auto_repair_acc = 1 if auto_repair_raw >=99998
*/

ren ER24138C5 parking
/*
ren ER25800 parking_raw
gen parking_acc = 0
replace parking_acc = 1 if parking_raw >= 99998 
*/

ren ER24138C6 bus_train
/*ren ER25801 bus_train_raw
gen bus_train_acc = 0
replace bus_train_acc = 1 if bus_train_raw >= 99998
*/

ren ER24138C7 taxi
/*
ren ER25802 taxi_raw
gen taxi_acc = 0
replace taxi_acc = 1 if taxi_raw >= 99998
*/

ren ER24138C8 other_trans
/*
ren ER25803 other_trans_raw
gen other_trans_acc = 0
replace other_trans_acc = 1 if other_trans_raw >=99998
*/

ren ER24138C9 education
/*ren ER25805 education_raw1
ren ER25807 education_raw2
gen education_acc = 0
replace education_acc = 1 if education_raw1 >= 999998 | education_raw2 >= 999998
*/

ren ER24138D1 childcare
/*ren ER25628 childcare_raw
gen childcare_acc = 0
replace childcare_raw = 1 if childcare_raw >= 999998
*/

ren ER24138D3 inst_medical
/*
ren ER27239 inst_medical_raw
gen inst_medical_acc = 0
replace inst_medical_acc = 1 if inst_medical_raw >= 999998
*/

ren ER24138D4 doctor
/*ren ER27245 doctor_raw
gen doctor_acc = 0
replace doctor_acc = 1 if doctor_raw >= 999998
*/

ren ER24138D5 prescription
/*
ren ER27251 prescription_raw
gen prescription_acc = 0
replace prescription_acc = 1 if prescription_raw >= 999998
*/

ren ER24138D6 health_ins
/*
ren ER27238 health_ins_raw
gen health_ins_acc = 0
replace health_ins_acc = 1 if health_ins_raw >= 99998
*/

/*

NOT INCLUDED PRE-2005

ren ER28037D8 home_repair
ren ER25808 home_repair_raw
gen home_repair_acc = 0
replace home_repair_acc = 1 if home_repair_raw >= 999998

ren ER28037D9 home_furnish
ren ER25813 home_furnish_raw
gen home_furnish_acc = 0
replace home_furnish_acc = 1 if home_furnish_raw >= 999998

ren ER28037E1 clothing
ren ER25818 clothing_raw
gen clothing_acc = 0
replace clothing_acc = 1 if clothing_raw >= 999998

ren ER28037E2 trips
ren ER25823 trips_raw
gen trips_acc = 0
replace trips_acc = 1 if trips_raw >= 999998

ren ER28037E3 other_rec
ren ER25828 other_rec_raw
gen other_rec_acc = 0
replace other_rec_acc = 1 if other_rec >= 999998
*/

/// INCOME ///

ren ER24099 y
gen truncy=y<1
replace y=1 if y<1

ren ER24100 tyhw
ren ER24101 trhw
ren ER24102 tyoth
ren ER24103 troth
ren ER24104 soc_sec

/*
ren ER22003 renty
replace renty = renty*52 if ER22004==3
replace renty = renty*26 if ER22004==4
replace renty = renty*12 if ER22004==5
replace renty = . if ER22004>6


ren ER22336 rentyw
replace rentyw = rentyw*52 if ER22337==3
replace rentyw = rentyw*26 if ER22337==4
replace rentyw = rentyw*12 if ER22337==5
replace rentyw = . if ER22337>6
replace rentyw=0if ER22339==5 /* "Is that in addition to your amount, response no" */


ren ER22020 divi
replace divi = divi*52 if ER22021==3
replace divi = divi*26 if ER22021==4
replace divi = divi*12 if ER22021==5
replace divi = . if ER22021>6

ren ER22353 diviw
replace diviw = diviw*52 if ER22354==3
replace diviw = diviw*26 if ER22354==4
replace diviw = diviw*12 if ER22354==5
replace diviw = . if ER22354>6
replace diviw = 0 if ER22356==5

ren ER22037 interest
replace interest = interest*52 if ER22038==3
replace interest = interest*26 if ER22038==4
replace interest = interest*12 if ER22038==5
replace interest = . if ER22038>6

ren ER22370 interestw
replace interestw = interestw*52 if ER22371==3
replace interestw = interestw*26 if ER22371==4
replace interestw = interestw*12 if ER22371==5
replace interestw = . if ER22371>6
replace interestw = 0 if ER22373==5

ren ER22054 trustinc
replace trustinc = trustinc*52 if ER22055==3
replace trustinc = trustinc*26 if ER22055==4
replace trustinc = trustinc*12 if ER22055==5
replace trustinc = . if ER22055>6

ren ER22387 trustincw
replace trustincw = trustincw*52 if ER22388==3
replace trustincw = trustincw*26 if ER22388==4
replace trustincw = trustincw*12 if ER22388==5
replace trustincw = . if ER22388>6
*/

* Generate labor and asset part of farm income *
gen farmlabor=0
replace farmlabor=0.5*ER24105 			if ER24105>0
replace farmlabor=. 					if ER24105>9999997

gen farmasset=0
replace farmasset=0.5*ER24105 			if ER24105>0
replace farmasset=ER24105 				if ER24105<0
replace farmasset=. 					if ER24105>9999997

* Hours Worked *
ren ER24080 hours
ren ER24091 hourw

* Income *
ren ER24116 wages
ren ER24135 wagesw
ren ER24109 y_business
ren ER24111 yw_business

egen ly = rsum(wages farmlabor y_business) 
egen wly= rsum(wagesw yw_business)


/// Prep variables to save ///

#delimit;

keep 


id fsize age sex agew kids marit race* wrace* newborn
outkid house smsa fchg state weight foodstamps *educ

bus_ass bus_debt cash real_estate real_estate_debt stocks vehicles other_ass
ira other_debt home_equity total_wealth

fstmp food_home food_out food_deliv mortgage rent prop_tax homeinsure utilities 
vehicle_loan vehicle_down vehicle_lease auto_ins other_vehicle gasoline auto_repair
parking bus_train taxi other_trans education childcare inst_medical doctor prescription
health_ins 

y truncy tyhw trhw tyoth troth soc_sec
farmlabor farmasset hours* wages* y_business yw_business ly wly empst* wempst*
weeks_emp* hourw
exp_job_years exp_job_months exp_job_weeks unempl layoff
; 

#delimit cr


gen year=2003

save "$temp\HH_2003.dta", replace


********************************************************************************************
* 2001
********************************************************************************************

use "$data\2001_H", clear

* Merge with the supplementary wealth file
ren ER17002 S501
merge 1:1 S501 using "$data\2001_W.dta"

ren S501 id // Corrseponds to wealth file identifier
ren ER17012 fsize
ren ER17013 age
ren ER17014 sex
ren ER17015 agew
ren ER17016 kids
ren ER17017 newborn
ren ER17024 marit
ren ER19989 race1
ren ER19990 race2
ren ER19991 race3
ren ER19992 race4
ren ER19897 wrace1
ren ER19898 wrace2
ren ER19899 wrace3
ren ER19900 wrace4
ren ER19172 outkid
ren ER17044 house
ren ER20378 smsa 
ren ER17007 fchg
ren ER17004 state
ren ER20394 weight 
ren ER17216 empst1
ren ER17217 empst2
ren ER17218 empst3
ren ER17786 wempst1
ren ER17787 wempst2
ren ER17788 wempst3
ren ER18386 foodstamps
ren ER20395 weeks_emp
ren ER20406 weeks_empw
ren ER17254 exp_job_years
ren ER17255 exp_job_months
ren ER17256 exp_job_weeks
ren ER17353 unempl 

replace weeks_emp =. if weeks_emp == 98 | weeks_emp ==99
replace exp_job_years =. if exp_job_years == 98 | exp_job_years ==99
replace exp_job_months =. if exp_job_months == 98 | exp_job_months ==99
replace exp_job_weeks =. if exp_job_weeks == 98 | exp_job_weeks ==99
ren ER20457 educ
ren ER20458 weduc

replace educ =. if educ == 99
replace weduc =. if weduc == 99

replace unempl =. if unempl == 8 | unempl == 9 
replace unempl =0 if unempl == 5 | unempl == 0 



/// ASSETS ///

* Note that these are taken from the imputed values, _acc variables denote imputed entries
* Prior to 2009 these taken from separate wealth file

ren S503 bus_ass
ren S503A bus_ass_acc
gen bus_debt=.
/* ^^^ Note that Prior to 2013 the equity value is given without debt/asset split */

ren S505 cash
ren S505A cash_acc  

ren S509 real_estate
ren S509A real_estate_acc
gen real_estate_debt=.
/* ^^^ Note that Prior to 2013 the equity value is given without debt/asset split */

ren S511 stocks
ren S511A stocks_acc

ren S513 vehicles
ren S513A vehicles_acc

ren S515 other_ass
ren S515A other_ass_acc

ren S519 ira
ren S519A ira_acc

ren S507 other_debt
ren S507A other_debt_acc // Credit card charges, student loans, medical or legal bills, or loans from relatives?

/* ^^^^ Note "other debt" to be calculated manually after 2009*/

ren S520 home_equity
ren S520A home_equity_acc

ren S517 total_wealth
ren S517A total_wealth_acc

/* ^^^^ This variable is constructed as sum of values of seven asset types (S703, S705, S709,
S711, S713, S715, S719) net of debt value (S707) plus value of home equity. */

ren ER19349 pension_plan
replace pension_plan=. if pension_plan>999999997

/* No data on wife pension pre 2011 */ 

/// EXPENDITURE ///

gen       fstmp=0
cap program drop doit
program def doit
        local i=18390
        while `i' < 18402 {
            replace ER`i'=. if ER`i'>1
        	local i=`i'+1
       }
end
qui doit
egen nm=rsum(ER18390-ER18401)
replace fstmp=ER18387              if ER18388==6 
replace fstmp=ER18387*nm           if ER18388==5 
replace fstmp=ER18387*((26/12)*nm) if ER18388==4  
replace fstmp=ER18387*((52/12)*nm) if ER18388==3
replace fstmp=.					 if ER18387>999997  | ER18388>6 
drop nm

* Note that these are taken from the imputed values, _acc variables denote imputed entries

ren ER20456A2 food_home
ren ER20456A3 food_out
ren ER20456A4 food_deliv

ren ER20456A6 mortgage  

/* ren ER25044 mortgage_raw1
ren ER25055 mortgage_raw2
gen mortgage_acc = 0
replace mortgage_acc = 1 if mortgage_raw1 >= 99998 | mortgage_raw2 >= 99998
*/

ren ER20456A7 rent
/* ren ER25063 rent_raw
gen rent_acc = 0
replace rent_acc = 1 if rent_raw >= 99998
*/
 
ren ER20456A8 prop_tax 
/*ren ER25036 prop_tax_raw
gen prop_tax_acc = 0
replace prop_tax_acc = 1 if prop_tax_raw >= 99998 
*/

ren ER20456A9 homeinsure
/* ren ER25038 homeinsure_raw
gen homeinsure_acc = 0
replace homeinsure_acc = 1 if homeinsure_raw>=9998 */

ren ER20456B1 utilities
ren ER20456B2 gas_raw
ren ER20456B3 electric_raw
ren ER20456B4 water_raw
ren ER20456B5 otheru_raw
gen utilities_acc = 0
replace utilities_acc = 1 if gas_raw >= 9999999 | electric_raw >= 9999999 | water_raw >= 9999999 | otheru_raw >= 9999999 

/*

NO TELECOMS PRE-2005

ren ER28037B6 telecom
ren ER25086 telecom_raw
gen telecom_acc = 0
replace telecom_acc = 1 if telecom_raw>=9998

*/

ren ER20456B7 vehicle_loan
/*
ren ER25729 vehicle_loan_raw1 
ren ER25757 vehicle_loan_raw2
ren ER25785 vehicle_loan_raw3
gen vehicle_loan_acc = 0
replace vehicle_loan_acc = 1 if vehicle_loan_raw1 >=999998 | vehicle_loan_raw2 >=999998 | vehicle_loan_raw3 >=999998 
*/

ren ER20456B8 vehicle_down
/*ren ER25726 vehicle_down_raw1
ren ER25754 vehicle_down_raw2
ren ER25782 vehicle_down_raw3
gen vehicle_down_acc = 0
replace vehicle_down_acc = 1 if vehicle_down_raw1 >= 999998 | vehicle_down_raw2 >= 999998 | vehicle_down_raw3 >= 999998 
*/

ren ER20456B9 vehicle_lease
/*
ren ER25733 vehicle_lease_raw1
ren ER25761 vehicle_lease_raw2
ren ER25789 vehicle_lease_raw3
gen vehicle_lease_acc = 0
replace vehicle_lease_acc = 1 if vehicle_lease_raw1 >=999998 | vehicle_lease_raw2 >=999998 |vehicle_lease_raw3 >=999998 
*/

ren ER20456C1 auto_ins
/*
ren ER25794 auto_ins_raw
gen auto_ins_acc = 0
replace auto_ins_acc = 1 if auto_ins_raw >= 999998
*/

ren ER20456C2 other_vehicle
/* ren ER25797 other_vehicle_raw
gen other_vehicle_acc = 0
replace other_vehicle_acc  = 1 if other_vehicle_raw >=999998
*/

ren ER20456C4 gasoline
/*
ren ER25799 gasoline_raw
gen gasoline_acc = 0
replace gasoline_acc = 1 if gasoline_raw >= 99998
*/

ren ER20456C3 auto_repair
/*
ren ER25798 auto_repair_raw
gen auto_repair_acc = 0
replace auto_repair_acc = 1 if auto_repair_raw >=99998
*/

ren ER20456C5 parking
/*
ren ER25800 parking_raw
gen parking_acc = 0
replace parking_acc = 1 if parking_raw >= 99998 
*/

ren ER20456C6 bus_train
/*ren ER25801 bus_train_raw
gen bus_train_acc = 0
replace bus_train_acc = 1 if bus_train_raw >= 99998
*/

ren ER20456C7 taxi
/*
ren ER25802 taxi_raw
gen taxi_acc = 0
replace taxi_acc = 1 if taxi_raw >= 99998
*/

ren ER20456C8 other_trans
/*
ren ER25803 other_trans_raw
gen other_trans_acc = 0
replace other_trans_acc = 1 if other_trans_raw >=99998
*/

ren ER20456C9 education
/*ren ER25805 education_raw1
ren ER25807 education_raw2
gen education_acc = 0
replace education_acc = 1 if education_raw1 >= 999998 | education_raw2 >= 999998
*/

ren ER20456D1 childcare
/*ren ER25628 childcare_raw
gen childcare_acc = 0
replace childcare_raw = 1 if childcare_raw >= 999998
*/

ren ER20456D3 inst_medical
/*
ren ER27239 inst_medical_raw
gen inst_medical_acc = 0
replace inst_medical_acc = 1 if inst_medical_raw >= 999998
*/

ren ER20456D4 doctor
/*ren ER27245 doctor_raw
gen doctor_acc = 0
replace doctor_acc = 1 if doctor_raw >= 999998
*/

ren ER20456D5 prescription
/*
ren ER27251 prescription_raw
gen prescription_acc = 0
replace prescription_acc = 1 if prescription_raw >= 999998
*/

ren ER20456D6 health_ins
/*
ren ER27238 health_ins_raw
gen health_ins_acc = 0
replace health_ins_acc = 1 if health_ins_raw >= 99998
*/

/*

NOT INCLUDED PRE-2005

ren ER28037D8 home_repair
ren ER25808 home_repair_raw
gen home_repair_acc = 0
replace home_repair_acc = 1 if home_repair_raw >= 999998

ren ER28037D9 home_furnish
ren ER25813 home_furnish_raw
gen home_furnish_acc = 0
replace home_furnish_acc = 1 if home_furnish_raw >= 999998

ren ER28037E1 clothing
ren ER25818 clothing_raw
gen clothing_acc = 0
replace clothing_acc = 1 if clothing_raw >= 999998

ren ER28037E2 trips
ren ER25823 trips_raw
gen trips_acc = 0
replace trips_acc = 1 if trips_raw >= 999998

ren ER28037E3 other_rec
ren ER25828 other_rec_raw
gen other_rec_acc = 0
replace other_rec_acc = 1 if other_rec >= 999998
*/

/// INCOME ///

ren ER20456 y
gen truncy=y<1
replace y=1 if y<1

ren ER20449 tyhw
ren ER20450 trhw
ren ER20453 tyoth
ren ER20454 troth
ren ER20455 soc_sec

/* SOCIAL SECURITY INCOME NOT SPLIT BY MEMBER PRE 2003*/
/*
ren ER18634 renty
replace renty = renty*52 if ER18635==3
replace renty = renty*26 if ER18635==4
replace renty = renty*12 if ER18635==5
replace renty = . if ER18635>6


ren ER18650 divi
replace divi = divi*52 if ER18651==3
replace divi = divi*26 if ER18651==4
replace divi = divi*12 if ER18651==5
replace divi = . if ER18651>6

ren ER18966 diviw
replace diviw = diviw*52 if ER18967==3
replace diviw = diviw*26 if ER18967==4
replace diviw = diviw*12 if ER18967==5
replace diviw = . if ER18967>6

ren ER18666 interest
replace interest = interest*52 if ER18667==3
replace interest = interest*26 if ER18667==4
replace interest = interest*12 if ER18667==5
replace interest = . if ER18667>6

ren ER18982 interestw
replace interestw = interestw*52 if ER18983==3
replace interestw = interestw*26 if ER18983==4
replace interestw = interestw*12 if ER18983==5
replace interestw = . if ER18983>6

ren ER18682 trustinc
replace trustinc = trustinc*52 if ER18683==3
replace trustinc = trustinc*26 if ER18683==4
replace trustinc = trustinc*12 if ER18683==5
replace trustinc = . if ER18683>6

ren ER18998 trustincw
replace trustincw = trustincw*52 if ER18999==3
replace trustincw = trustincw*26 if ER18999==4
replace trustincw = trustincw*12 if ER18999==5
*/


* Generate labor and asset part of farm income *
gen farmlabor=0
replace farmlabor=0.5*ER20420 			if ER20420>0
replace farmlabor=. 					if ER20420>9999997

gen farmasset=0
replace farmasset=0.5*ER20420 			if ER20420>0
replace farmasset=ER20420 				if ER20420<0
replace farmasset=. 					if ER20420>9999997

* Hours Worked *
ren ER20399 hours
ren ER20410 hourw

* Income *
ren ER20443 wages
ren ER20447 wagesw
ren ER20422 y_business
ren ER20444 yw_business

egen ly = rsum(wages farmlabor y_business) 
egen wly= rsum(wagesw yw_business)


/// Prep variables to save ///

#delimit;

keep 


id fsize age sex agew kids marit race* wrace* newborn
outkid house smsa fchg state weight foodstamps *educ

bus_ass bus_debt cash real_estate real_estate_debt stocks vehicles other_ass
ira other_debt home_equity total_wealth

fstmp food_home food_out food_deliv mortgage rent prop_tax homeinsure utilities 
vehicle_loan vehicle_down vehicle_lease auto_ins other_vehicle gasoline auto_repair
parking bus_train taxi other_trans education childcare inst_medical doctor prescription
health_ins 

y truncy tyhw trhw tyoth troth soc_sec
farmlabor farmasset hours* wages* y_business yw_business ly wly empst* wempst* hourw
weeks_emp*
exp_job_years exp_job_months exp_job_weeks unempl
; 

#delimit cr


gen year=2001

save "$temp\HH_2001.dta", replace

********************************************************************************************
* 1999
********************************************************************************************

use "$data\1999_H", clear

* Merge with the supplementary wealth file
ren ER13002 S401
merge 1:1 S401 using "$data\1999_W.dta"

ren S401 id // Corrseponds to wealth file identifier
ren ER13009 fsize
ren ER13010 age
ren ER13011 sex
ren ER13012 agew
ren ER13013 kids
ren ER13014 newborn
ren ER16423 marit
ren ER15928 race1
ren ER15929 race2
ren ER15930 race3
ren ER15931 race4
ren ER15836 wrace1
ren ER15837 wrace2
ren ER15838 wrace3
ren ER15839 wrace4
ren ER14976 outkid
ren ER13041 house
ren ER16432 smsa 
ren ER13008A fchg
ren ER13004 state
ren ER16518 weight 
ren ER13205 empst1
ren ER13206 empst2
ren ER13207 empst3
ren ER13717 wempst1
ren ER13718 wempst2
ren ER13719 wempst3
ren ER14255 foodstamps
ren ER16467 weeks_emp
ren ER16478 weeks_empw
ren ER13243 exp_job_years
ren ER13244 exp_job_months
ren ER13245 exp_job_weeks
ren ER13330 unempl 


replace weeks_emp =. if weeks_emp == 98 | weeks_emp ==99
replace exp_job_years =. if exp_job_years == 98 | exp_job_years ==99
replace exp_job_months =. if exp_job_months == 98 | exp_job_months ==99
replace exp_job_weeks =. if exp_job_weeks == 98 | exp_job_weeks ==99

ren ER16516 educ
ren ER16517 weduc

replace unempl =. if unempl == 8 | unempl == 9 
replace unempl =0 if unempl == 5 | unempl == 0 

replace educ =. if educ == 99
replace weduc =. if weduc == 99


/// ASSETS ///

* Note that these are taken from the imputed values, _acc variables denote imputed entries
* Prior to 2009 these taken from separate wealth file

ren S403 bus_ass
ren S403A bus_ass_acc
gen bus_debt=.
/* ^^^ Note that Prior to 2013 the equity value is given without debt/asset split */

ren S405 cash
ren S405A cash_acc  

ren S409 real_estate
ren S409A real_estate_acc
gen real_estate_debt=.
/* ^^^ Note that Prior to 2013 the equity value is given without debt/asset split */

ren S411 stocks
ren S411A stocks_acc

ren S413 vehicles
ren S413A vehicles_acc

ren S415 other_ass
ren S415A other_ass_acc

ren S419 ira
ren S419A ira_acc

ren S407 other_debt
ren S407A other_debt_acc // Credit card charges, student loans, medical or legal bills, or loans from relatives?

/* ^^^^ Note "other debt" to be calculated manually after 2009*/

ren S420 home_equity
ren S420A home_equity_acc

ren S417 total_wealth
ren S417A total_wealth_acc

/* ^^^^ This variable is constructed as sum of values of seven asset types (S703, S705, S709,
S711, S713, S715, S719) net of debt value (S707) plus value of home equity. */

ren ER15181 pension_plan
replace pension_plan=. if pension_plan>999999997

/* No data on wife pension pre 2011 */ 

/// EXPENDITURE ///

gen       fstmp=0
cap program drop doit
program def doit
        local i=14258
        while `i' < 14270 {
            replace ER`i'=. if ER`i'>1
        	local i=`i'+1
       }
end
qui doit
egen nm=rsum(ER14258-ER14269)
replace fstmp=ER14256              if ER14257==6 
replace fstmp=ER14256*nm           if ER14257==5 
replace fstmp=ER14256*((26/12)*nm) if ER14257==4  
replace fstmp=ER14256*((52/12)*nm) if ER14257==3
replace fstmp=.					 if ER14256>999997  | ER14257>6 
drop nm

* Note that these are taken from the imputed values, _acc variables denote imputed entries

ren ER16515A2 food_home
ren ER16515A3 food_out
ren ER16515A4 food_deliv

ren ER16515A6 mortgage  

/* ren ER25044 mortgage_raw1
ren ER25055 mortgage_raw2
gen mortgage_acc = 0
replace mortgage_acc = 1 if mortgage_raw1 >= 99998 | mortgage_raw2 >= 99998
*/

ren ER16515A7 rent
/* ren ER25063 rent_raw
gen rent_acc = 0
replace rent_acc = 1 if rent_raw >= 99998
*/
 
ren ER16515A8 prop_tax 
/*ren ER25036 prop_tax_raw
gen prop_tax_acc = 0
replace prop_tax_acc = 1 if prop_tax_raw >= 99998 
*/

ren ER16515A9 homeinsure
/* ren ER25038 homeinsure_raw
gen homeinsure_acc = 0
replace homeinsure_acc = 1 if homeinsure_raw>=9998 */

ren ER16515B1 utilities
ren ER16515B2 gas_raw
ren ER16515B3 electric_raw
ren ER16515B4 water_raw
ren ER16515B5 otheru_raw
gen utilities_acc = 0
replace utilities_acc = 1 if gas_raw >= 9999999 | electric_raw >= 9999999 | water_raw >= 9999999 | otheru_raw >= 9999999 

/*

NO TELECOMS PRE-2005

ren ER28037B6 telecom
ren ER25086 telecom_raw
gen telecom_acc = 0
replace telecom_acc = 1 if telecom_raw>=9998

*/

ren ER16515B7 vehicle_loan
/*
ren ER25729 vehicle_loan_raw1 
ren ER25757 vehicle_loan_raw2
ren ER25785 vehicle_loan_raw3
gen vehicle_loan_acc = 0
replace vehicle_loan_acc = 1 if vehicle_loan_raw1 >=999998 | vehicle_loan_raw2 >=999998 | vehicle_loan_raw3 >=999998 
*/

ren ER16515B8 vehicle_down
/*ren ER25726 vehicle_down_raw1
ren ER25754 vehicle_down_raw2
ren ER25782 vehicle_down_raw3
gen vehicle_down_acc = 0
replace vehicle_down_acc = 1 if vehicle_down_raw1 >= 999998 | vehicle_down_raw2 >= 999998 | vehicle_down_raw3 >= 999998 
*/

ren ER16515B9 vehicle_lease
/*
ren ER25733 vehicle_lease_raw1
ren ER25761 vehicle_lease_raw2
ren ER25789 vehicle_lease_raw3
gen vehicle_lease_acc = 0
replace vehicle_lease_acc = 1 if vehicle_lease_raw1 >=999998 | vehicle_lease_raw2 >=999998 |vehicle_lease_raw3 >=999998 
*/

ren ER16515C1 auto_ins
/*
ren ER25794 auto_ins_raw
gen auto_ins_acc = 0
replace auto_ins_acc = 1 if auto_ins_raw >= 999998
*/

ren ER16515C2 other_vehicle
/* ren ER25797 other_vehicle_raw
gen other_vehicle_acc = 0
replace other_vehicle_acc  = 1 if other_vehicle_raw >=999998
*/

ren ER16515C4 gasoline
/*
ren ER25799 gasoline_raw
gen gasoline_acc = 0
replace gasoline_acc = 1 if gasoline_raw >= 99998
*/

ren ER16515C3 auto_repair
/*
ren ER25798 auto_repair_raw
gen auto_repair_acc = 0
replace auto_repair_acc = 1 if auto_repair_raw >=99998
*/

ren ER16515C5 parking
/*
ren ER25800 parking_raw
gen parking_acc = 0
replace parking_acc = 1 if parking_raw >= 99998 
*/

ren ER16515C6 bus_train
/*ren ER25801 bus_train_raw
gen bus_train_acc = 0
replace bus_train_acc = 1 if bus_train_raw >= 99998
*/

ren ER16515C7 taxi
/*
ren ER25802 taxi_raw
gen taxi_acc = 0
replace taxi_acc = 1 if taxi_raw >= 99998
*/

ren ER16515C8 other_trans
/*
ren ER25803 other_trans_raw
gen other_trans_acc = 0
replace other_trans_acc = 1 if other_trans_raw >=99998
*/

ren ER16515C9 education
/*ren ER25805 education_raw1
ren ER25807 education_raw2
gen education_acc = 0
replace education_acc = 1 if education_raw1 >= 999998 | education_raw2 >= 999998
*/

ren ER16515D1 childcare
/*ren ER25628 childcare_raw
gen childcare_acc = 0
replace childcare_raw = 1 if childcare_raw >= 999998
*/

ren ER16515D3 inst_medical
/*
ren ER27239 inst_medical_raw
gen inst_medical_acc = 0
replace inst_medical_acc = 1 if inst_medical_raw >= 999998
*/

ren ER16515D4 doctor
/*ren ER27245 doctor_raw
gen doctor_acc = 0
replace doctor_acc = 1 if doctor_raw >= 999998
*/

ren ER16515D5 prescription
/*
ren ER27251 prescription_raw
gen prescription_acc = 0
replace prescription_acc = 1 if prescription_raw >= 999998
*/

ren ER16515D6 health_ins
/*
ren ER27238 health_ins_raw
gen health_ins_acc = 0
replace health_ins_acc = 1 if health_ins_raw >= 99998
*/

/*

NOT INCLUDED PRE-2005

gen home_repair =.
ren ER25808 home_repair_raw
gen home_repair_acc = 0
replace home_repair_acc = 1 if home_repair_raw >= 999998

ren ER28037D9 home_furnish
ren ER25813 home_furnish_raw
gen home_furnish_acc = 0
replace home_furnish_acc = 1 if home_furnish_raw >= 999998

ren ER28037E1 clothing
ren ER25818 clothing_raw
gen clothing_acc = 0
replace clothing_acc = 1 if clothing_raw >= 999998

ren ER28037E2 trips
ren ER25823 trips_raw
gen trips_acc = 0
replace trips_acc = 1 if trips_raw >= 999998

ren ER28037E3 other_rec
ren ER25828 other_rec_raw
gen other_rec_acc = 0
replace other_rec_acc = 1 if other_rec >= 999998
*/

/// INCOME ///

ren ER16462 y
gen truncy=y<1
replace y=1 if y<1

ren ER16452 tyhw
ren ER16454 trhw
ren ER16456 tyoth
ren ER16458 troth
ren ER16460 soc_sec

/* SOCIAL SECURITY INCOME NOT SPLIT BY MEMBER PRE 2003*/

/*
* Asset Income *

ren ER65217 renty
ren ER65218 renty_acc
ren ER65245 rentyw
ren ER65246 renty_acc

ren ER65219 divi
ren ER65220 divi_acc
ren ER65247 diviw
ren ER65248 diviy_acc

ren ER65221 interest
ren ER65222 interest_acc
ren ER65249 interestw
ren ER65250 interestw_acc

ren ER65223 trustinc
ren ER65223 trustinc_acc
ren ER65251 trustincw
ren ER65252 trustincw_acc
*/

* Generate labor and asset part of farm income *
gen farmlabor=0
replace farmlabor=0.5*ER16448 			if ER16448>0
replace farmlabor=. 					if ER16448>9999997

gen farmasset=0
replace farmasset=0.5*ER16448 			if ER16448>0
replace farmasset=ER16448 				if ER16448<0
replace farmasset=. 					if ER16448>9999997

* Hours Worked *
ren ER16471 hours
ren ER16482 hourw

* Income *
ren ER16463 wages
ren ER16465 wagesw
ren ER16490 y_business
ren ER16511 yw_business

egen ly = rsum(wages farmlabor y_business) 
egen wly= rsum(wagesw yw_business)


/// Prep variables to save ///

#delimit;

keep 


id fsize age sex agew kids marit race* wrace* newborn
outkid house smsa fchg state weight foodstamps *educ

bus_ass bus_debt cash real_estate real_estate_debt stocks vehicles other_ass
ira other_debt home_equity total_wealth

fstmp food_home food_out food_deliv mortgage rent prop_tax homeinsure utilities 
vehicle_loan vehicle_down vehicle_lease auto_ins other_vehicle gasoline auto_repair
parking bus_train taxi other_trans education childcare inst_medical doctor prescription
health_ins 

y truncy tyhw trhw tyoth troth soc_sec
farmlabor farmasset hours* wages* y_business yw_business ly wly empst* wempst* 
weeks_emp* hourw
exp_job_years exp_job_months exp_job_weeks unempl
; 

#delimit cr


gen year=1999

save "$temp\HH_1999.dta", replace


/// Collate it all

use "$temp\HH_1999", clear

forvalues i = 2001(2)2017{
append using "$temp\HH_`i'.dta"

}
save "$temp\HH_Panel.dta", replace

forvalues i = 1997(2)2017{
erase "$temp\HH_`i'.dta"
}

cd "$do"

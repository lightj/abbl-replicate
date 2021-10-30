*******************************************************************************
* Codes to replicate the sample in ABBL (2021)
* This version: January 2021
* Questions: jdlight@uchicago.edu
*******************************************************************************

********************************************************************************************
* This file does the basic sample selection (similar to BPP 2008, BPS 2016, ABB 2017)
********************************************************************************************
*
* 1. Clean variables 
* 2. NBER's TAXSIM (can ignore this if using pre-tax vars)
* 
********************************************************************************************

cd "$do"
u "$temp\Full_Panel_1.dta",clear

********************************************************************************************
* 1. Construct Panel with Complete Information, Recode Variables
********************************************************************************************

* Drop if SEO sample 
keep if seo==0
* Drop pre-2005
keep if year >= 1999 & year <=2017


/// I generate indicators for whether in the selected sample ///
gen sample_selectA = 1 /* If it relates to the sample */ 
gen sample_selectB = 1 /* If it relates to missing data */

* Drop if missing state
replace sample_selectB=0 if state==.
replace sample_selectB=0 if state==0
replace sample_selectB=0 if state ==99

* Pre-2005 social security is coded jointly
replace soc_secw =0 if year <2005 & soc_secw==.
replace soc_sec=0 if year<1999 & soc_sec==.

* Drop with no information on other family
replace sample_selectB=0 if kids == . 			  
replace sample_selectB=0 if fsize ==.  			  

* Keep married households
replace sample_selectA=0 if marit!=1 | agew==0

* Drop unstable HH
replace sample_selectA=0 if fchg>1

* Drop female household heads
replace sample_selectA=0 if sex==2

* Recode employment status
gen empst=.
gen tempempst =.
forvalues y=1(1)7{
replace tempempst=`y' if empst1==`y' | empst2==`y' | empst3==`y'
}
replace tempempst=. if empst1==.
replace tempempst = empst if year <1993
replace empst=tempempst

gen wempst=.
gen tempwempst =.
forvalues y=1(1)7{
replace tempwempst=`y' if wempst1==`y' | wempst2==`y' | wempst3==`y'
}
replace tempwempst=. if wempst1==.
replace tempwempst = wempst if year <1993
replace wempst=tempwempst


* Recode education so consistent across waves
gen neweduc=educ
replace neweduc =. if year <=1989
gen newweduc=weduc
replace newweduc=. if year<1989
egen eductemp=max(neweduc), by(person)
replace neweduc = eductemp
drop eductemp
egen eductemp=max(newweduc), by(person)
replace newweduc = eductemp
drop eductemp

gen tempeduc=.
replace tempeduc = 1 if (educ == 1 | educ ==0) & year <= 1989
replace tempeduc = 2 if educ == 2 & year <= 1989
replace tempeduc = 3 if educ == 3 & educ <=11 & year <= 1989
replace tempeduc = 4 if (educ == 4 | educ == 5) & year <= 1989
replace tempeduc = 5 if educ == 6 & year <= 1989
replace tempeduc = 6 if (educ == 7 | educ == 8) & year <= 1989
replace tempeduc =. if educ ==9 & year <= 1989


replace tempeduc = 1 if educ <=5 & year > 1989
replace tempeduc = 2 if educ > 5 & educ <=8 & year > 1989
replace tempeduc = 3 if educ > 8 & educ <=11 & year > 1989
replace tempeduc = 4 if educ == 12 & year > 1989
replace tempeduc = 5 if educ > 12 & educ <=15 & year > 1989
replace tempeduc = 6 if educ > 15 & educ <=17 & year > 1989
replace educ = tempeduc

gen tempweduc=.
replace tempweduc = 1 if (weduc == 1 | weduc ==0) & year <= 1989
replace tempweduc = 2 if weduc == 2 & year <= 1989
replace tempweduc = 3 if weduc == 3 & weduc <=11 & year <= 1989
replace tempweduc = 4 if (weduc == 4 | weduc == 5) & year <= 1989
replace tempweduc = 5 if weduc == 6 & year <= 1989
replace tempweduc = 6 if (weduc == 7 | weduc == 8) & year <= 1989
replace tempweduc =. if weduc ==9 & year <= 1989

replace tempweduc = 1 if weduc <=5 & year > 1989
replace tempweduc = 2 if weduc > 5 & weduc <=8 & year > 1989
replace tempweduc = 3 if weduc > 8 & weduc <=11 & year > 1989
replace tempweduc = 4 if weduc == 12 & year > 1989
replace tempweduc = 5 if weduc > 12 & weduc <=15 & year > 1989
replace tempweduc = 6 if weduc > 15 & weduc <=17 & year > 1989
replace weduc = tempweduc

egen eductemp=max(educ), by(person)
replace educ = eductemp
replace sample_selectB=0 if eductemp==.

egen weductemp=max(weduc), by(person)
replace weduc = weductemp
replace sample_selectB=0 if weductemp==.


replace sample_selectB=0 if educ==.
replace sample_selectB=0 if weduc==.


* Higher education dummy
gen he=.
replace he = 2 if educ==6
replace he = 1 if educ==5
replace he = 0 if inrange(educ,1,4)

* Higher education dummy
gen whe=.
replace whe = 2 if weduc==6
replace whe = 1 if weduc==5
replace whe = 0 if inrange(weduc,1,4)

* Code family composition change
gen fam_change_other=.
replace fam_change_other=0 if fchg==0
replace fam_change_other=1 if fchg==1

gen fam_change_spouse=.
replace fam_change_spouse=0 if fchg==0
replace fam_change_spouse=1 if fchg==2

* Code firm specific experience
gen experience = (52*exp_job_years + 4*exp_job_months + exp_job_weeks)/52

* Code FTE dummy
gen full_time=.
replace full_time =1 if empst==1
replace full_time=0 if empst==2 | empst==3

gen wfull_time=.
replace wfull_time =1 if wempst==1
replace wfull_time=0 if wempst==2 | wempst==3

* Code FTE dummy
gen full_time2=.
replace full_time2 =1 if empst==1  | empst==2 | empst==3
replace full_time2=0 if inrange(empst,4,8)

gen wfull_time2=.
replace wfull_time2 =1 if wempst==1  | wempst==2 | wempst==3
replace wfull_time2=0 if inrange(wempst,4,8)

* Code to combine unempl and layoff
gen unempl2 = .
replace unempl2 = 1 if unempl == 0 & layoff ==0
replace unempl2 = 1 if unempl == 0 & layoff ==.
replace unempl2 = 0 if unempl ==1 | layoff ==1 


* Drop those with always missing info on age 
replace age=. if age ==0 | age >110
egen n=sum(person!=.), by(person)
egen na=sum(age==. | age==0), by(person)
tabstat na, by(na) stat(N)
replace sample_selectB=0 if n==na
drop n na

replace agew=. if agew==0 | agew>110
egen n=sum(person!=.),by(person)
egen na=sum(agew==.),by(person)
tabstat na, by(na) stat(N)
replace sample_selectB=0 if n==na
drop n na

* Recode age so that there is no gap or jump
egen lasty=max(year) if age!=., by(person)
replace sample_selectB=0 if last==.
gen lastage=age if year==lasty
gen b=year-lastage
replace b=0 if b==.
egen yb=sum(b),by(person)
replace age=year-yb
drop lasty lastage b

egen lasty=max(year) if agew!=., by(person)
replace sample_selectB=0 if lasty==.
gen lastage=agew if year==lasty
gen b=year-lastage
replace b=0 if b==.
egen ybw=sum(b),by(person)
replace agew=year-ybw

drop lasty lastage b

* Takes into account the retrospective nature of the data
replace age=age-1
replace agew=agew-1

* Code race consistently across waves
gen race=.
ren race old_race

* Race: 1 is white, 2 black, 3 others
replace race1=old_race if race1==. & old_race!=.

forvalues i=1/4 {
	gen     rc=race`i' if race`i'<=2
	replace rc=3    if race`i'>=3 & race`i'<=7
	drop race`i'
	ren rc race`i'

	gen     wrc=wrace`i' if wrace`i'<=2
	replace wrc=3    if wrace`i'>=3 & wrace`i'<=7
	drop wrace`i'
	ren wrc wrace`i'
}


gen race=race1
replace race=3 if race2!=0 & race2!=race1		/* Any case of more than one reported race is also considered "other" */ 
gen wrace = wrace1
replace wrace = 3 if wrace2!=0 & wrace2!=wrace1

drop *race1 *race2 *race3 *race4 
replace sample_selectB=0 if race==. | race==0
replace sample_selectB=0 if wrace==. | wrace==0
* Recode newborn
gen newborntemp=.
replace newborntemp=0 if inrange(newborn,2,20)
replace newborntemp=1 if newborn==1
replace newborn=newborntemp

gen agebin=.
replace agebin=1 if age>=25 & age <43
replace agebin=0 if age>=43 & age <60

gen agewbin=.
replace agewbin=1 if agew>=25 & agew <43
replace agewbin=0 if agew>=43 & agew <60

gen sample_select=0
replace sample_select=1 if sample_selectA==1 & sample_selectB==1

order sample_selectA sample_selectB sample_select

replace house=. if house == 9999999

********************************************************************************************
* Cohorts
********************************************************************************************
gen coh=.
forvalues i = 1930(10)2000{
replace coh=(`i'-1920)/10 if yb>=`i' & yb<`i'+10
}

gen cohw=.
forvalues i = 1930(10)2000{
replace cohw=(`i'-1920)/10 if ybw>=`i' & ybw<`i'+10
}

********************************************************************************************
* 2. Prepare for NBER's TAXSIM Package:
***  You need entries for ALL TAXSIM inputs or it won't run (i.e. set to 0 if we don't have): 
********************************************************************************************

gen taxsimid=1
replace taxsimid=sum(taxsimid)

preserve
drop if sample_select==0
replace year = year-1 /* Account for retrospective nature of income data */

* Generate TAXSIM compatible state codes from PSID codes
gen state2=.
replace state2=1 if state ==1
replace state2=state+1 if state >= 2 & state <=10
replace state2=state+2 if state >=11 & state <=49
replace state2 = 2 if state ==50
replace state2 = 12 if state ==51
replace state=state2

gen mstat =2 /* Make PSID marital code */

gen page=age /* Head age */

gen sage =agew /* Wife age */

gen depx = kids /* Dependents */

gen pwages = max(ly, wly) /* Primary taxpayer.*/

gen swages = min(ly, wly) /* Secondary taxpayer */

gen dividends = divi + diviw + trustinc  + trustincw + businc + busincw /* Treat: trusts, royalty income, business asset income as divis too */
replace dividends=0 if dividends==. | dividends<0 /* Should be none missing - all imputed post 2005 */
replace dividends=0 /*NO ASSET INCOME*/

gen intrec=interest + interestw /* Should be none missing - all imputed post 2005 */
replace intrec=0 if intrec==. | intrec<0
replace intrec=0 /*NO ASSET INCOME*/

gen stcg =0 /* Not including capital gains */

gen ltcg=0 /* Not including capital gains */

gen otherprop= renty + rentyw /* Should be none missing - all imputed post 2005 */
replace otherprop=0 if otherprop==. | otherprop<0
replace otherprop=0 /*NO ASSET INCOME*/

gen nonprop=0 /* Putting the transfers through together in transfers */

gen pensions=0 /* Looking at in work only */

gen gssi=soc_sec + soc_secw

replace gssi = 0 if gssi==. | gssi <0

gen ui=0 /* Putting all social security throgh gssi */

gen transfers=trhw
replace transfers = 0 if transfers==. | transfers <0

gen rentpaid=rent /* Not including housing */
replace rentpaid=0 if rentpaid ==. | rentpaid<0
replace rentpaid = 0.06*house if house!=. & house!=0

gen proptax = prop_tax 
replace proptax = 0 if proptax==. | proptax<0 /* Should be none missing - all imputed post 2005 */

gen otheritem =0 /* Not using other itemized deductions */

replace childcare =0 if childcare==. | childcare<0 /* Should be none missing - all imputed post 2005 */

replace mortgage = 0  /* Could add mortgage interest and medical deductibles */

* Keep only variables to be uploaded to TAXSIM
keep taxsimid year state mstat page sage depx pwages swages dividends intrec stcg ltcg otherprop nonprop pensions gssi ui transfers rentpaid proptax childcare mortgage
order taxsimid year state mstat page sage depx pwages swages dividends intrec stcg ltcg otherprop nonprop pensions gssi ui transfers rentpaid proptax childcare mortgage 

taxsim27 

* Output to directory
save "$TAXSIM\TAXSIM_main.dta", replace
restore


save "$temp\Full_Panel_2b.dta", replace

cd "$programs"


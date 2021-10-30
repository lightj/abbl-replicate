*******************************************************************************
* Codes to replicate the sample in ABBL (2021)
* This version: January 2021
* Questions: jdlight@uchicago.edu
*******************************************************************************

********************************************************************************************
* This file constructs the panel panel of household heads from the longitudinal file
********************************************************************************************

*\ Structure of File:
*
* 1. Load most recent individual file (contains all current/historic heads)
* 2. Rename relevant variables
* 3. Attach to households
*
********************************************************************************************
* 1. Load Individual file 
********************************************************************************************

clear all

cd "$output"

* Load the individual file
u "$data\2017IND.dta", clear

#delimit;
keep  	ER30001	ER30002 /* 1968 Family ID and Sequence # which can be used for unique id */

		ER30020 ER30067 ER30117 ER30160 ER30217 ER30283 ER30343 ER30399 ER30463
		ER30535 ER30606 ER30689 ER30806 ER33201 ER33401 ER33501 ER33601 ER33701 
		ER33801 ER33901	ER34001 ER34101 ER34201 ER34301 ER34501 /* Interview number */
		
        ER30021 ER30068 ER30118 ER30161 ER30218 ER30284 ER30344 ER30400 ER30464 
		ER30536 ER30607 ER30690 ER30807 ER33202 ER33402 ER33502	ER33602 ER33702 
		ER33802 ER33902 ER34002	ER34102 ER34202 ER34302 ER34502 /* Sequence number */
		
      	ER30022 ER30069 ER30119 ER30162 ER30219 ER30285 ER30345 ER30401 ER30465 
		ER30537 ER30608 ER30691 ER30808 ER33203 ER33403 ER33503 ER33603 ER33703
		ER33803 ER33903	ER34003 ER34103 ER34203 ER34303 ER34503 /* Relationship to Head */
		
		ER33506 ; /* Year individual born + education*/
#delimit cr

rename ER33506 DOB

* Generate uniqe identifier
gen persid = (ER30001*1000) + ER30002

* Drop Latino households
drop if ER30001>=7001 & ER30001<=9308 

* Dummy for SEO sample
gen seo=ER30001>=5000 & ER30001<=7000

* Drop individuals who are never household heads
drop if ER33503!=10 & ER33603!=10 & ER33703!=10 & ER33803!=10 & ER33903!=10 ///
& ER34003!=10 & ER34103!=10	& ER34203!=10 & ER34303!=10	& ER34503!=10 ///
& ER30022!=10 & ER30069!=10 & ER30119!=10 & ER30162!=10 ///
& ER30219!=10 & ER30285!=10 & ER30345!=10 & ER30401!=10 ///
& ER30465!=10 & ER30537!=10 & ER30608!=10 & ER30691!=10 ///
& ER30808!=10 & ER33203!=10 & ER33403!=10 & ER33403!=10

* Rename hh interview # as yearly ID
ren ER30020 id1969
ren ER30067 id1971
ren ER30117 id1973
ren ER30160 id1975
ren ER30217 id1977
ren ER30283 id1979
ren ER30343 id1981
ren ER30399 id1983
ren ER30463 id1985
ren	ER30535 id1987
ren ER30606 id1989
ren ER30689 id1991
ren ER30806 id1993
ren ER33201 id1995
ren ER33401 id1997
ren ER33501 id1999
ren ER33601 id2001
ren ER33701 id2003
ren ER33801 id2005
ren ER33901 id2007
ren ER34001 id2009
ren ER34101 id2011
ren ER34201 id2013
ren ER34301 id2015
ren ER34501 id2017

* Generate dummy for HH each year
gen pid1969=ER30022==1 & ER30021>=1 & ER30021<=20
gen pid1971=ER30069==1 & ER30068>=1 & ER30068<=20
gen pid1973=ER30119==1 & ER30118>=1 & ER30118<=20
gen pid1975=ER30162==1 & ER30161>=1 & ER30161<=20
gen pid1977=ER30219==1 & ER30218>=1 & ER30218<=20
gen pid1979=ER30285==1 & ER30284>=1 & ER30284<=20
gen pid1981=ER30345==1 & ER30344>=1 & ER30344<=20
gen pid1983=ER30401==10 & ER30400>=1 & ER30400<=20
gen pid1985=ER30465==10 & ER30464>=1 & ER30464<=20
gen pid1987=ER30537==10 & ER30536>=1 & ER30536<=20
gen pid1989=ER30608==10 & ER30607>=1 & ER30607<=20
gen pid1991=ER30691==10 & ER30690>=1 & ER30690<=20
gen pid1993=ER30808==10 & ER30807>=1 & ER30807<=20
gen pid1995=ER33203==10 & ER33202>=1 & ER33202<=20
gen pid1997=ER33403==10 & ER33402>=1 & ER33402<=20
gen pid1999=ER33503==10 & ER33502>=1 & ER33502<=20
gen pid2001=ER33603==10 & ER33602>=1 & ER33602<=20
gen pid2003=ER33703==10 & ER33702>=1 & ER33702<=20
gen pid2005=ER33803==10 & ER33802>=1 & ER33802<=20
gen pid2007=ER33903==10 & ER33902>=1 & ER33902<=20
gen pid2009=ER34003==10 & ER34002>=1 & ER34002<=20
gen pid2011=ER34103==10 & ER34102>=1 & ER34102<=20
gen pid2013=ER34203==10 & ER34202>=1 & ER34202<=20
gen pid2015=ER34303==10 & ER34302>=1 & ER34302<=20
gen pid2017=ER34503==10 & ER34502>=1 & ER34502<=20

* Keep relevant variables
keep id* pid* seo persid DOB

* Reshape into panel data format
gen person=1
replace person=sum(person)
reshape long id pid,i(person) j(year)


* Drop any individual in any year for which they are no head of household
drop if pid==0

* Label person id
label var person "UNIQUE PERSON IDENTIFER (HH HEADS ONLY)" 

compress


save "$temp\I_Panel.dta", replace

merge 1:1 id year using "$temp\HH_Panel.dta"
drop if _merge!=3
drop _merge

by person, sort: gen numwav=_N

sort numwav person year 

save "$temp\Full_Panel_1.dta", replace

cd "$programs"

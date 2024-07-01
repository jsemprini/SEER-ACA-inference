***wild meta-replication***

clear all

cd "WORKING DIRECTORY PATH"

set more off
  quietly log
  local logon = r(status)
  if "`logon'" == "on" { 
	log close 
	}
log using wild-test, replace text


foreach a in all  {
	clear all

cd "WORKING DIRECTORY PATH WITH RAW SEER DATA"

	import delimited `a'.txt
		
	**create binary outcomes**
	gen distant=.
	replace distant=0 if v1=="In situ" | v1=="Localized" | v1=="Regional"
	replace distant=1 if v1=="Distant"
	
	gen mo2tx=.
	replace mo2tx=0 if v2=="000" | v2=="001"
	destring(v2), force replace
	replace mo2tx=1 if v2>1 & v2!=.
	
	tab v3, gen(dum_v3)
	gen nosurg=.
	replace nosurg=0 if dum_v38==1
	replace nosurg=1 if dum_v38!=1 & dum_v31!=1 & dum_v37!=1 & dum_v39!=1
	
	gen dead=0
	replace dead=1 if v4=="2020" | v4=="Alive at last contact"
	
	global y_outcome distant mo2tx nosurg dead
	
	**create controls **
	
	rename v5 year
	
	tab v6, gen(dum_state)
	gen statefips=.
	foreach n of numlist 1/16{
		replace statefips=`n' if dum_state`n'==1
	}
	
	rename v6 state
	
	drop dum_state*
	
	tab v7, gen(dum_age)
	
	gen male=0
	replace male=1 if v8=="Male"
	
	tab v9, gen(dum_race)
	
	tab v10, gen(dum_mar)
	
	tab v11, gen(dum_hhinc)
	
	drop dum_hhinc11
	
	tab v12, gen(dum_metro)
	
	rename dum_metro1 metro
	
	drop dum_metro*
	
	global x_controls dum_age* male dum_race* dum_mar* dum_hhinc* metro
	
	
	drop if v13=="Breast" & male==1

	***code treatment (KFF 2023)***
	
gen post=0
replace post=1 if year>=2014

gen medicaid=1
replace medicaid=0 if state=="Georgia" | state=="Idaho" | state=="Texas" |  state=="New York" | state=="Massachusetts"
 
 drop if state=="Louisiana" 
 
 gen expand=post*medicaid
 
 
 ***analysis***

 foreach y in $y_outcome{
 eststo: areg `y' i.expand i.year i.($x_controls)  , vce(ols) absorb(statefips)

  eststo: areg `y' i.expand i.year i.($x_controls)  , vce(robust) absorb(statefips)

 eststo: areg `y' i.expand i.year i.($x_controls)  , vce(cluster statefips) absorb(statefips)
  
eststo: boottest 1.expand , cluster(statefips) rep(999) seed(123) bootcluster(statefips) nograph weight(webb) ptyp(equaltail)

estadd scalar p1=r(p)
matrix b = r(CI)
matrix list b
estadd scalar lowci=b[1,1]
estadd scalar highci=b[1,2]
 }
 
 
 
 

 esttab using `a'-ci.csv, keep(*expand*) sca(p1 lowci highci) b(3) ci(3) replace
 esttab using `a'-p.csv, keep(*expand*) sca(p1 lowci highci) b(3) p(3) replace

 estimates clear 
 
	}
	
	clear all 
	
	
foreach a in topten top_pro top_lb top_mel top_bla top_kid top_uter top_pan top_thy  top_nhl screen sc_fb sc_cvx sc_crc {
	clear all




	import delimited `a'.txt
		
	**create binary outcomes**
	gen distant=.
	replace distant=0 if v1=="In situ" | v1=="Localized" | v1=="Regional"
	replace distant=1 if v1=="Distant"
	
	gen mo2tx=.
	replace mo2tx=0 if v2=="000" | v2=="001"
	destring(v2), force replace
	replace mo2tx=1 if v2>1 & v2!=.
	
	tab v3, gen(dum_v3)
	gen nosurg=.
	replace nosurg=0 if dum_v38==1
	replace nosurg=1 if dum_v38!=1 & dum_v31!=1 & dum_v37!=1 & dum_v39!=1
	
	gen dead=0
	replace dead=1 if v4=="2020" | v4=="Alive at last contact"
	
	global y_outcome distant mo2tx nosurg dead
	
	**create controls **
	
	rename v5 year
	
	tab v6, gen(dum_state)
	gen statefips=.
	foreach n of numlist 1/16{
		replace statefips=`n' if dum_state`n'==1
	}
	
	rename v6 state
	
	drop dum_state*
	
	tab v7, gen(dum_age)
	
	gen male=0
	replace male=1 if v8=="Male"
	
	tab v9, gen(dum_race)
	
	tab v10, gen(dum_mar)
	
	tab v11, gen(dum_hhinc)
	
	drop dum_hhinc11
	
	tab v12, gen(dum_metro)
	
	rename dum_metro1 metro
	
	drop dum_metro*
	
		tab v13, gen(dum_site)

	global x_controls dum_age* male dum_race* dum_mar* dum_hhinc* metro dum_site*
	
	drop if v13=="Breast" & male==1
	

	***code treatment (KFF 2023)***
	
gen post=0
replace post=1 if year>=2014

gen medicaid=1
replace medicaid=0 if state=="Georgia" | state=="Idaho" | state=="Texas" |  state=="New York" | state=="Massachusetts"
 
 drop if state=="Louisiana" 
 
 
 gen expand=post*medicaid
 
 
 ***analysis***

 foreach y in $y_outcome{
 eststo: areg `y' i.expand i.year i.($x_controls)  , vce(ols) absorb(statefips)

  eststo: areg `y' i.expand i.year i.($x_controls)  , vce(robust) absorb(statefips)

 eststo: areg `y' i.expand i.year i.($x_controls)  , vce(cluster statefips) absorb(statefips)
  
eststo: boottest 1.expand , cluster(statefips) rep(999) seed(123) bootcluster(statefips) nograph weight(webb) ptyp(equaltail)

estadd scalar p1=r(p)
matrix b = r(CI)
matrix list b
estadd scalar lowci=b[1,1]
estadd scalar highci=b[1,2]
 }
 
 
 
 

 esttab using `a'-ci.csv, keep(*expand*) sca(p1 lowci highci) b(3) ci(3) replace
 esttab using `a'-p.csv, keep(*expand*) sca(p1 lowci highci) b(3) p(3) replace

 estimates clear 
 
	}
	
	clear all 
	
	
	clear all 
	
	
foreach a in  top_leuk  {
	clear all




	import delimited `a'.txt
	
	
	**create binary outcomes**
	gen distant=.
	replace distant=0 if v1=="In situ" | v1=="Localized" | v1=="Regional"
	replace distant=1 if v1=="Distant"
	
	gen mo2tx=.
	replace mo2tx=0 if v2=="000" | v2=="001"
	destring(v2), force replace
	replace mo2tx=1 if v2>1 & v2!=.
	
	tab v3, gen(dum_v3)
	gen nosurg=.
	replace nosurg=0 if dum_v38==1
	replace nosurg=1 if dum_v38!=1 & dum_v31!=1 & dum_v37!=1 & dum_v39!=1
	
	gen dead=0
	replace dead=1 if v4=="2020" | v4=="Alive at last contact"
	
	global y_outcome mo2tx nosurg dead
	
	**create controls **
	
	rename v5 year
	
	tab v6, gen(dum_state)
	gen statefips=.
	foreach n of numlist 1/16{
		replace statefips=`n' if dum_state`n'==1
	}
	
	rename v6 state
	
	drop dum_state*
	
	tab v7, gen(dum_age)
	
	gen male=0
	replace male=1 if v8=="Male"
	
	tab v9, gen(dum_race)
	
	tab v10, gen(dum_mar)
	
	tab v11, gen(dum_hhinc)
	
	drop dum_hhinc11
	
	tab v12, gen(dum_metro)
	
	rename dum_metro1 metro
	
	drop dum_metro*
	
		tab v13, gen(dum_site)

	global x_controls dum_age* male dum_race* dum_mar* dum_hhinc* metro dum_site*
	
	drop if v13=="Breast" & male==1
	

	***code treatment (KFF 2023)***
	
gen post=0
replace post=1 if year>=2014

gen medicaid=1
replace medicaid=0 if state=="Georgia" | state=="Idaho" | state=="Texas" |  state=="New York" | state=="Massachusetts"
 
 drop if state=="Louisiana" 
 
 
 gen expand=post*medicaid
 
 
 ***analysis***

 foreach y in $y_outcome{
 eststo: areg `y' i.expand i.year i.($x_controls)  , vce(ols) absorb(statefips)

  eststo: areg `y' i.expand i.year i.($x_controls)  , vce(robust) absorb(statefips)

 eststo: areg `y' i.expand i.year i.($x_controls)  , vce(cluster statefips) absorb(statefips)
  
eststo: boottest 1.expand , cluster(statefips) rep(999) seed(123) bootcluster(statefips) nograph weight(webb) ptyp(equaltail)

estadd scalar p1=r(p)
matrix b = r(CI)
matrix list b
estadd scalar lowci=b[1,1]
estadd scalar highci=b[1,2]
 }
 
 
 
 

 esttab using `a'-ci.csv, keep(*expand*) sca(p1 lowci highci) b(3) ci(3) replace
 esttab using `a'-p.csv, keep(*expand*) sca(p1 lowci highci) b(3) p(3) replace

 estimates clear 
 
	}
	
	clear all 
	
	exit, STATA
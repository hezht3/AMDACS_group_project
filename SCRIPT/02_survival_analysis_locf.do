*##############################*
*# AMDACS group project       #*
*# Aim: survival analysis     #*
*# Date: October 5, 2022      #*
*##############################*

cd "/Users/zhengtinghe/Library/CloudStorage/OneDrive-JohnsHopkins/Course/340.728.01 - Advanced Methods for Design and Analysis of Cohort Studies/project/AMDACS_group_project"

*Sevly Code:
*cd /Users/sevly/Library/CloudStorage/OneDrive-JohnsHopkins/AMDACS_group_project

capture log close
log using "02_survival_analysis.log", replace

*######################################*
*# LOCF, CKD diagnosis as time origin #*
*######################################*
* Inport LOCF cleaned dataset
use INPUT/Cleaned/jointfromckd_survival_locf.dta

gen SBPgroup1 = 0
replace SBPgroup1 = 1 if SBPgroup == 1
gen SBPgroup2 = 0
replace SBPgroup2 = 1 if SBPgroup == 2
gen SBPgroup3 = 0
replace SBPgroup3 = 1 if SBPgroup == 3
gen SBPgroup4 = 0
replace SBPgroup4 = 1 if SBPgroup == 4

gen income_recode = 1
replace income_recode = 2 if income == 3 | income == 4 | income == 5 | income == 6
replace income_recode = 3  if income == 7
replace income_recode = 4 if income == 8

gen race_recode = 1
replace race_recode = 2 if race == 2
replace race_recode = 3 if race !=1 & race !=2

* Descriptive analysis
tab renalstatus SBPgroup, col   // CKD progression in 4 SBP groups

* Set Data as Survival Time Data Incorporating Late Entries
stset YearsFromCKDTransition, enter(YearsFromCKD) failure(renalstatus)

* Non-parametric Estimates of Survival Function for SBP groups (Kaplan-Meier method)
set scheme s1color
sts graph, by(SBPgroup) risktable(0 5 10 15 20 25 30) xtitle(Years after CKD diagnosis) ytitle(Probability of no CKD progression)

*********************************
* Cox proportional hazard model *
*********************************
*** Crude model
stcox SBPgroup2 SBPgroup3 SBPgroup4, nohr   // display parameters
stcox SBPgroup2 SBPgroup3 SBPgroup4   // display hazard ratios
estat ic
estimates store cox_crude

stcox SBPgroup, nohr   // trend test

*** Adjusted model
stcox SBPgroup2 SBPgroup3 SBPgroup4 i.race_recode lbw Glomdx1NG0 i.income_recode icu anemia, nohr   // display parameters
stcox SBPgroup2 SBPgroup3 SBPgroup4 i.race_recode lbw Glomdx1NG0 i.income_recode icu anemia   // display hazard ratios
estat ic
estimates store cox_adjust

***********************************
* Concentional Weibull regression *
***********************************
*** Crude model
streg SBPgroup2 SBPgroup3 SBPgroup4, dist(weibull) time   // display parameters
lincom _cons + SBPgroup2   // Get beta
lincom _cons + SBPgroup3
lincom _cons + SBPgroup4

streg SBPgroup2 SBPgroup3 SBPgroup4, dist(weibull) hr   // display hazard ratios
streg SBPgroup2 SBPgroup3 SBPgroup4, dist(weibull) tr   // display time ratios

estat ic
estimates store weibull_crude

*** Adjusted model
streg SBPgroup2 SBPgroup3 SBPgroup4 i.race_recode lbw Glomdx1NG0 i.income_recode icu anemia, dist(weibull) time   // display parameters
lincom _cons + SBPgroup2   // Get beta
lincom _cons + SBPgroup3
lincom _cons + SBPgroup4

streg SBPgroup2 SBPgroup3 SBPgroup4 i.race_recode lbw Glomdx1NG0 i.income_recode icu anemia, dist(weibull) hr   // display hazard ratios
streg SBPgroup2 SBPgroup3 SBPgroup4 i.race_recode lbw Glomdx1NG0 i.income_recode icu anemia, dist(weibull) tr   // display time ratios

estat ic
estimates store weibull_adjust

*********************************
* Non-proportional hazard model *
*********************************
sum _t , detail   // Determine the 75th percentile of the times which is close to 16

*** Crude model
stcox SBPgroup2 SBPgroup3 SBPgroup4, tvc(SBPgroup2 SBPgroup3 SBPgroup4) texp((_t - 16)*(_t<=16)) nohr   // display parameters
stcox SBPgroup2 SBPgroup3 SBPgroup4, tvc(SBPgroup2 SBPgroup3 SBPgroup4) texp((_t - 16)*(_t<=16))   // display hazard ratios
estat ic
estimates store nphazard_crude
lrtest cox_crude nphazard_crude   // No significant improvements

*** Adjusted model
stcox SBPgroup2 SBPgroup3 SBPgroup4 i.race_recode lbw Glomdx1NG0 i.income_recode icu anemia, tvc(SBPgroup2 SBPgroup3 SBPgroup4) texp((_t - 16)*(_t<=16)) nohr   // display parameters
stcox SBPgroup2 SBPgroup3 SBPgroup4 i.race_recode lbw Glomdx1NG0 i.income_recode icu anemia, tvc(SBPgroup2 SBPgroup3 SBPgroup4) texp((_t - 16)*(_t<=16))   // display hazard ratios
estat ic
estimates store nphazard_adjust
lrtest cox_adjust nphazard_adjust   // No significant improvements

******************************
* Conventional GG regression *
******************************
*** Overall
quietly: streg, dist(ggamma) tr
gen convGGSurv = gammaptail(e(kappa)^(-2),(e(kappa)^(-2))*(exp(-[_t]_cons)*YearsFromCKDTransition)^(e(kappa)/e(sigma)))
replace convGGSurv = 1 - convGGSurv if e(kappa)<0   /*(0 real changes made)*/
sts graph, surv xtitle(Years after CKD diagnosis) ytitle(Probability of no CKD progression) addplot(line convGGSurv YearsFromCKDTransition, sort)

*** Crude model
streg SBPgroup2 SBPgroup3 SBPgroup4, dist(ggamma) time   // Display parameters
lincom _cons + SBPgroup2   // Get beta
lincom _cons + SBPgroup3
lincom _cons + SBPgroup4

streg SBPgroup2 SBPgroup3 SBPgroup4, dist(ggamma) tr   // Display time ratio

display chi2tail(1, 2*(300.72771-299.92337))   // conventional GG not significantly better
display 2*(300.72771-299.92337)

ggtaxonomy
estat ic
estimates store ggconven_crude

gen convGGSurv1 = .
gen convGGSurv2 = .
gen convGGSurv3 = .
gen convGGSurv4 = .
replace convGGSurv1 = gammaptail(e(kappa)^(-2),(e(kappa)^(-2))*(exp(-[_t]_cons)*YearsFromCKDTransition)^(e(kappa)/e(sigma))) ///
	if SBPgroup1 == 1
replace convGGSurv2 = gammaptail(e(kappa)^(-2),(e(kappa)^(-2))*(exp(-([_t]_cons + [_t]SBPgroup2))*YearsFromCKDTransition)^(e(kappa)/e(sigma))) if SBPgroup2 == 1
replace convGGSurv3 = gammaptail(e(kappa)^(-2),(e(kappa)^(-2))*(exp(-([_t]_cons + [_t]SBPgroup3))*YearsFromCKDTransition)^(e(kappa)/e(sigma))) if SBPgroup3 == 1
replace convGGSurv4 = gammaptail(e(kappa)^(-2),(e(kappa)^(-2))*(exp(-([_t]_cons + [_t]SBPgroup4))*YearsFromCKDTransition)^(e(kappa)/e(sigma))) if SBPgroup4 == 1

sts graph, by(SBPgroup) ///
	plot1(lcolor(blue)) plot2(lcolor(orange)) plot3(lcolor(green)) plot4(lcolor(red)) ///
	xtitle(Years after CKD diagnosis) ytitle(Probability of survival from CKD progression) ///
	addplot(line convGGSurv1 YearsFromCKDTransition if SBPgroup1 == 1, sort lp(dash) lc(blue) leg(off) || ///
			line convGGSurv2 YearsFromCKDTransition if SBPgroup2 == 1, sort lp(dash) lc(orange) leg(off) || ///
			line convGGSurv3 YearsFromCKDTransition if SBPgroup3 == 1, sort lp(dash) lc(green) leg(off) || ///
			line convGGSurv4 YearsFromCKDTransition if SBPgroup4 == 1, sort lp(dash) lc(red) leg(off) ///
			legend(label(1 "(0, 50)") label(2 "[50, 90)") label(3 "[90, 95)") label(4 "[95, 100]")))

*** Adjusted model
streg SBPgroup2 SBPgroup3 SBPgroup4 i.race_recode lbw Glomdx1NG0 i.income_recode icu anemia, dist(ggamma) time   // Display parameters
lincom _cons + SBPgroup2   // Get beta
lincom _cons + SBPgroup3
lincom _cons + SBPgroup4

streg SBPgroup2 SBPgroup3 SBPgroup4 i.race_recode lbw Glomdx1NG0 i.income_recode icu anemia, dist(ggamma) tr   // Display time ratio

ggtaxonomy
estat ic
estimates store ggconven_adjust

display chi2tail(1, 2*(202.24555-200.65977))   // (p = 0.075)
display 2*(202.24555-200.65977)

***************************
* Saturated GG regression *
***************************
*** Crude model
streg SBPgroup2 SBPgroup3 SBPgroup4, ///
	anc(SBPgroup2 SBPgroup3 SBPgroup4) ///
	anc2(SBPgroup2 SBPgroup3 SBPgroup4) ///
	dist(ggamma) time   // Display parameters

lincom [_t]_cons + [_t]SBPgroup2   // Get beta
lincom [_t]_cons + [_t]SBPgroup3
lincom [_t]_cons + [_t]SBPgroup4

nlcom exp([lnsigma]_cons)   // Get sigma
nlcom exp([lnsigma]_cons + [lnsigma]SBPgroup2)
nlcom exp([lnsigma]_cons + [lnsigma]SBPgroup3)
nlcom exp([lnsigma]_cons + [lnsigma]SBPgroup4)

lincom [kappa]_cons + [kappa]SBPgroup2   // Get kappa
lincom [kappa]_cons + [kappa]SBPgroup3
lincom [kappa]_cons + [kappa]SBPgroup4

estat ic
estimates store ggsatura_crude

display chi2tail(12-6, 2*(299.92337-297.79334))   // P-value = 0.64, saturated not significantly better
display 2*(299.92337-297.79334)

* Bootstrap to obtain CI for differences of quantiles
bootstrap, seed(19210428) reps(500) ///
	saving(./OUTPUT/500saturatedGG.dta, replace): ///
streg SBPgroup2 SBPgroup3 SBPgroup4, ///
	anc(SBPgroup2 SBPgroup3 SBPgroup4) ///
	anc2(SBPgroup2 SBPgroup3 SBPgroup4) ///
	dist(ggamma) time iter(50)

*** Reproduce crude model stratify by SBPgroup to get ggtaxonomy
***** SBP Group 1
quietly: streg if SBPgroup1 == 1, dist(ggamma)
ggtaxonomy
replace convGGSurv1 = gammaptail(e(kappa)^(-2),(e(kappa)^(-2))*(exp(-[_t]_cons)*YearsFromCKDTransition)^(e(kappa)/e(sigma)))
replace convGGSurv1 = 1 - convGGSurv1 if e(kappa)<0   /*(0 real changes made)*/
sts graph if SBPgroup1 == 1, surv xtitle(Years after CKD diagnosis) ytitle(Probability of no CKD progression) addplot(line convGGSurv1 YearsFromCKDTransition, sort)
predict hazardSBPgroup1, hazard
***** SBP Group 2
quietly: streg if SBPgroup2 == 1, dist(ggamma)
ggtaxonomy
replace convGGSurv2 = gammaptail(e(kappa)^(-2),(e(kappa)^(-2))*(exp(-[_t]_cons)*YearsFromCKDTransition)^(e(kappa)/e(sigma)))
replace convGGSurv2 = 1 - convGGSurv2 if e(kappa)<0   /*(0 real changes made)*/
sts graph if SBPgroup2 == 1, surv xtitle(Years after CKD diagnosis) ytitle(Probability of no CKD progression) addplot(line convGGSurv2 YearsFromCKDTransition, sort)
predict hazardSBPgroup2, hazard
***** SBP Group 3
quietly: streg if SBPgroup3 == 1, dist(ggamma)
ggtaxonomy
replace convGGSurv3 = gammaptail(e(kappa)^(-2),(e(kappa)^(-2))*(exp(-[_t]_cons)*YearsFromCKDTransition)^(e(kappa)/e(sigma)))
replace convGGSurv3 = 1 - convGGSurv3 if e(kappa)<0   /*(0 real changes made)*/
sts graph if SBPgroup3 == 1, surv xtitle(Years after CKD diagnosis) ytitle(Probability of no CKD progression) addplot(line convGGSurv3 YearsFromCKDTransition, sort)
predict hazardSBPgroup3, hazard
***** SBP Group 4
quietly: streg if SBPgroup4 == 1, dist(ggamma)
ggtaxonomy
replace convGGSurv4 = gammaptail(e(kappa)^(-2),(e(kappa)^(-2))*(exp(-[_t]_cons)*YearsFromCKDTransition)^(e(kappa)/e(sigma)))
replace convGGSurv4 = 1 - convGGSurv4 if e(kappa)<0   /*(0 real changes made)*/
sts graph if SBPgroup4 == 1, surv xtitle(Years after CKD diagnosis) ytitle(Probability of no CKD progression) addplot(line convGGSurv4 YearsFromCKDTransition, sort)
predict hazardSBPgroup4, hazard

sts graph, by(SBPgroup) ///
	plot1(lcolor(blue)) plot2(lcolor(orange)) plot3(lcolor(green)) plot4(lcolor(red)) ///
	xtitle(Years after CKD diagnosis) ytitle(Probability of survival from CKD progression) ///
	addplot(line convGGSurv1 YearsFromCKDTransition if SBPgroup1 == 1, sort lp(dash) lc(blue) leg(off) || ///
			line convGGSurv2 YearsFromCKDTransition if SBPgroup2 == 1, sort lp(dash) lc(orange) leg(off) || ///
			line convGGSurv3 YearsFromCKDTransition if SBPgroup3 == 1, sort lp(dash) lc(green) leg(off) || ///
			line convGGSurv4 YearsFromCKDTransition if SBPgroup4 == 1, sort lp(dash) lc(red) leg(off) ///
			legend(label(1 "(0, 50)") label(2 "[50, 90)") label(3 "[90, 95)") label(4 "[95, 100]")))


****Revisit this model; Jason says it's wrong the most unacceptable to use race as continous. He suggests that we do a simple model first where we are doing streg with SPBgroup modifying beta sigma and kappa. This should converge, but it doesn't converge. If the values are not different for the Sigma and Kappa don't argument to sigma and kapps. The Crude fully saturated GG to see if the Sigma and Kappa changes, if it doesn't change. We do not have modification to do it. Then we do not need to modify it, then add all the rest. 


*** Adjusted model
streg SBPgroup2 SBPgroup3 SBPgroup4  i.race_recode lbw Glomdx1NG0 i.income_recode icu anemia, ///
	anc(SBPgroup2 SBPgroup3 SBPgroup4) ///
	anc2(SBPgroup2 SBPgroup3 SBPgroup4) ///
	dist(ggamma) time   // Display parameter
	
lincom [_t]_cons	
lincom [_t]_cons + [_t]SBPgroup2   // Get beta
lincom [_t]_cons + [_t]SBPgroup3
lincom [_t]_cons + [_t]SBPgroup4

nlcom exp([lnsigma]_cons)   // Get sigma
nlcom exp([lnsigma]_cons + [lnsigma]SBPgroup2)
nlcom exp([lnsigma]_cons + [lnsigma]SBPgroup3)
nlcom exp([lnsigma]_cons + [lnsigma]SBPgroup4)

lincom [kappa]_cons
lincom [kappa]_cons + [kappa]SBPgroup2   // Get kappa
lincom [kappa]_cons + [kappa]SBPgroup3
lincom [kappa]_cons + [kappa]SBPgroup4

estat ic
estimates store ggsatura_adjust

display chi2tail(12-6, 2*(200.65977-197.9085))   // P-value = 0.48, saturated not significantly better
display 2*(200.65977-197.9085)

* Bootstrap to obtain CI for differences of quantiles
bootstrap, seed(19210428) reps(500) ///
	saving(./OUTPUT/500saturatedGG_adjusted.dta, replace): ///
streg SBPgroup2 SBPgroup3 SBPgroup4  i.race_recode lbw Glomdx1NG0 i.income_recode icu anemia, ///
	dist(ggamma) time iter(50)

*** Reproduce crude model stratify by SBPgroup to get ggtaxonomy
***** SBP Group 1
quietly: streg  SBPgroup2 SBPgroup3 SBPgroup4 i.race_recode lbw Glomdx1NG0 i.income_recode icu anemia if SBPgroup1 == 1, dist(ggamma)
ggtaxonomy
***** SBP Group 2
quietly: streg SBPgroup2 SBPgroup3 SBPgroup4  i.race_recode lbw Glomdx1NG0 i.income_recode icu anemia if SBPgroup2 == 1, dist(ggamma)
ggtaxonomy
***** SBP Group 3
quietly: streg  SBPgroup2 SBPgroup3 SBPgroup4 i.race_recode lbw Glomdx1NG0 i.income_recode icu anemia if SBPgroup3 == 1, dist(ggamma)
ggtaxonomy
***** SBP Group 4
quietly: streg  SBPgroup2 SBPgroup3 SBPgroup4 i.race_recode lbw Glomdx1NG0 i.income_recode icu anemia if SBPgroup4 == 1, dist(ggamma)
ggtaxonomy

****************
* Final models *
****************

/* Crude model */
*Weibull for SBPgroup 1
streg if SBPgroup1 == 1, d(weibull) time
replace convGGSurv1 = gammaptail(1^(-2),(1^(-2))*(exp(-[_t]_cons)*YearsFromCKDTransition)^(1/.5965208))
estat ic

*Saturated for SBPgroup 2
streg if SBPgroup2 == 1, dist(ggamma)
replace convGGSurv2 = gammaptail(e(kappa)^(-2),(e(kappa)^(-2))*(exp(-[_t]_cons)*YearsFromCKDTransition)^(e(kappa)/e(sigma)))
replace convGGSurv2 = 1 - convGGSurv2 if e(kappa)<0   /*(0 real changes made)*/
estat ic

*ammag for SBPgroup 3
	
gen event= renalstatus /* the program assumes the variable event is 1 for uncensored times and 0 for censored times*/
capture program drop ammag_loglik 
program ammag_loglik 
	args lnlik beta lnsigma /*lnsigma=log(sigma) */
	tempvar kappa kappatom2 xentry xexit Sentry fexit Sexit  /*temporary variables*/
	quietly {
		 gen double `kappa' = exp(-`lnsigma')     // ammag whereby shape = 1/scale
				 
         gen double `kappatom2' = (`kappa')^(-2)
         gen double `xentry' = `kappatom2'*(exp(-(`beta'))*YearsFromCKD)^(`kappa'/exp(`lnsigma'))
		 gen double `xexit'  = `kappatom2'*(exp(-(`beta'))*YearsFromCKDTransition)^(`kappa'/exp(`lnsigma'))
		 
	     gen double `fexit' = ( abs(`kappa')*`xexit'*gammaden(`kappatom2',1,0,`xexit')/(exp(`lnsigma')*YearsFromCKDTransition) ) if event==1
	 
	     gen double `Sexit' = ( 1 - gammaptail(`kappatom2',`xexit') ) if (event==0 & `kappa'<0)
         replace    `Sexit' = ( gammaptail(`kappatom2',`xexit') )     if (event==0 & `kappa'>0)
	     gen double `Sentry' = ( 1 - gammaptail(`kappatom2',`xentry') ) if (`kappa'<0)
         replace    `Sentry' = ( gammaptail(`kappatom2',`xentry') )     if (`kappa'>0)
         
         replace `lnlik'   = ln( `fexit'/`Sentry')              if event==1  /* event at exit */
	     replace `lnlik'   = ln( `Sexit'/`Sentry')              if event==0  /*censored at exit */
		 }
end

ml model lf ammag_loglik  (beta: ) (lnsigma: )  if SBPgroup3 == 1   // ammag for SBPgroup 3
ml maximize
nlcom (sigma: exp([lnsigma]_cons))
replace convGGSurv3 = gammaptail((1/.4477858)^(-2),((1/.4477858)^(-2))*(exp(-2.818067)*YearsFromCKDTransition)^((1/.4477858)/.4477858))
estat ic

*standard gamma for SBPgroup 4

capture program drop gamma_loglik
program gamma_loglik 
	args lnlik beta lnsigma /*lnsigma=log(sigma) */
	tempvar kappa kappatom2 xentry xexit Sentry fexit Sexit  /*temporary variables*/
	quietly {
		 gen double `kappa' = exp(`lnsigma')     // gamma whereby shape = scale
				 
         gen double `kappatom2' = (`kappa')^(-2)
         gen double `xentry' = `kappatom2'*(exp(-(`beta'))*YearsFromCKD)^(`kappa'/exp(`lnsigma'))
		 gen double `xexit'  = `kappatom2'*(exp(-(`beta'))*YearsFromCKDTransition)^(`kappa'/exp(`lnsigma'))
		 
	     gen double `fexit' = ( abs(`kappa')*`xexit'*gammaden(`kappatom2',1,0,`xexit')/(exp(`lnsigma')*YearsFromCKDTransition) ) if event==1
	 
	     gen double `Sexit' = ( 1 - gammaptail(`kappatom2',`xexit') ) if (event==0 & `kappa'<0)
         replace    `Sexit' = ( gammaptail(`kappatom2',`xexit') )     if (event==0 & `kappa'>0)
	     gen double `Sentry' = ( 1 - gammaptail(`kappatom2',`xentry') ) if (`kappa'<0)
         replace    `Sentry' = ( gammaptail(`kappatom2',`xentry') )     if (`kappa'>0)
         
         replace `lnlik'   = ln( `fexit'/`Sentry')              if event==1  /* event at exit */
	     replace `lnlik'   = ln( `Sexit'/`Sentry')              if event==0  /*censored at exit */
		 }
end

ml model lf gamma_loglik  (beta: ) (lnsigma: )  if SBPgroup4 == 1
set more off
ml maximize
nlcom (sigma: exp([lnsigma]_cons))
replace convGGSurv4 = gammaptail((.5708262)^(-2),((.5708262)^(-2))*(exp(-2.39798)*YearsFromCKDTransition)^((.5708262)/.5708262))
estat ic

sts graph, by(SBPgroup) ///
	plot1(lcolor(blue)) plot2(lcolor(orange)) plot3(lcolor(green)) plot4(lcolor(red)) ///
	xtitle(Years after CKD diagnosis) ytitle(Probability of survival from CKD progression) ///
	addplot(line convGGSurv1 YearsFromCKDTransition if SBPgroup1 == 1, sort lp(dash) lc(blue) leg(off) || ///
			line convGGSurv2 YearsFromCKDTransition if SBPgroup2 == 1, sort lp(dash) lc(orange) leg(off) || ///
			line convGGSurv3 YearsFromCKDTransition if SBPgroup3 == 1, sort lp(dash) lc(green) leg(off) || ///
			line convGGSurv4 YearsFromCKDTransition if SBPgroup4 == 1, sort lp(dash) lc(red) leg(off) ///
			legend(label(1 "(0, 50)") label(2 "[50, 90)") label(3 "[90, 95)") label(4 "[95, 100]")))
			
/* Adjusted model */
*Weibull for SBPgroup 1
streg i.race_recode lbw Glomdx1NG0 i.income_recode icu anemia if SBPgroup1 == 1, d(weibull) time
replace convGGSurv1 = gammaptail(1^(-2),(1^(-2))*(exp(-[_t]_cons)*YearsFromCKDTransition)^(1/.5965208))
estat ic

*Saturated for SBPgroup 2
streg i.race_recode lbw Glomdx1NG0 i.income_recode icu anemia if SBPgroup2 == 1, dist(ggamma)
replace convGGSurv2 = gammaptail(e(kappa)^(-2),(e(kappa)^(-2))*(exp(-[_t]_cons)*YearsFromCKDTransition)^(e(kappa)/e(sigma)))
replace convGGSurv2 = 1 - convGGSurv2 if e(kappa)<0   /*(0 real changes made)*/
estat ic

*ammag for SBPgroup 3
ml model lf ammag_loglik  (beta: ) (lnsigma: )  i.race_recode lbw Glomdx1NG0 i.income_recode icu anemia if SBPgroup3 == 1   // ammag for SBPgroup 3
ml maximize
nlcom (sigma: exp([lnsigma]_cons))

*standard gamma for SBPgroup 4
ml model lf gamma_loglik  (beta: ) (lnsigma: )  i.race_recode lbw Glomdx1NG0 i.income_recode icu anemia if SBPgroup4 == 1
set more off
ml maximize


use ./OUTPUT/500saturatedGG.dta, clear

gen beta_SBP1 = _t_b_cons
gen sigma_SBP1 = exp(lnsigma_b_cons)
gen kappa_SBP1 = kappa_b_cons
gen beta_SBP2 = _t_b_cons + _t_b_SBPgroup2
gen sigma_SBP2 = exp(lnsigma_b_cons + lnsigma_b_SBPgroup2)
gen kappa_SBP2 = kappa_b_cons + kappa_b_SBPgroup2
gen beta_SBP3 = _t_b_cons + _t_b_SBPgroup3
gen sigma_SBP3 = exp(lnsigma_b_cons + lnsigma_b_SBPgroup3)
gen kappa_SBP3 = kappa_b_cons + kappa_b_SBPgroup3
gen beta_SBP4 = _t_b_cons + _t_b_SBPgroup4
gen sigma_SBP4 = exp(lnsigma_b_cons + lnsigma_b_SBPgroup4)
gen kappa_SBP4 = kappa_b_cons + kappa_b_SBPgroup4

save ./OUTPUT/500saturatedGG.dta, replace


use ./OUTPUT/500saturatedGG_adjusted.dta, clear

gen beta_SBP1 = _t_b_cons
gen sigma_SBP1 = exp(_eq2_b_lnsigma)
gen kappa_SBP1 = _eq2_b_kappa
gen beta_SBP2 = _t_b_cons + _t_b_SBPgroup2
gen sigma_SBP2 = exp(_eq2_b_lnsigma)
gen kappa_SBP2 = _eq2_b_kappa
gen beta_SBP3 = _t_b_cons + _t_b_SBPgroup3
gen sigma_SBP3 = exp(_eq2_b_lnsigma)
gen kappa_SBP3 = _eq2_b_kappa
gen beta_SBP4 = _t_b_cons + _t_b_SBPgroup4
gen sigma_SBP4 = exp(_eq2_b_lnsigma)
gen kappa_SBP4 = _eq2_b_kappa

save ./OUTPUT/500saturatedGG_adjusted.dta, replace
			
			
******************
* Summary tables *
******************
esttab cox_crude weibull_crude ggconven_crude ggsatura_crude, cells(b(fmt(%4.2f)) ci(fmt(%4.2f)) p(fmt(%4.3f)))
esttab nphazard_crude, cells(b(fmt(%4.2f)) ci(fmt(%4.2f)) p(fmt(%4.3f)))
esttab cox_adjust weibull_adjust ggconven_adjust ggsatura_adjust, cells(b(fmt(%4.2f)) ci(fmt(%4.2f)) p(fmt(%4.3f)))
esttab nphazard_adjust, cells(b(fmt(%4.2f)) ci(fmt(%4.2f)) p(fmt(%4.3f)))

save ./OUTPUT/data_for_plot.dta

log close



***************************
* Secondary Data Analysis *
***************************
/* By Sevly: To conceptualize race in our model we can use interaction terms. 
We can't stratify our model separately and we get the Betas, the confidence intervals we cannot use. The correct inference would be to use Race and the Exposure* The confidence interval is bloated because we stratified by the model." The confidence interval would not be the same. It becomes a problem when a group is not powered enough.* 
Adding interaction term from the full model and doing the lin_coms and etc. Strafiying is good because it helps us get the main effect, but it doesn't allow us to do statistical inference to help justify why would you do it*/ 

streg SBPgroup2 SBPgroup3 SBPgroup4  i.race_recode lbw Glomdx1NG0 i.income, anc(SBPgroup2 SBPgroup3 SBPgroup4) ///
	anc2(SBPgroup2 SBPgroup3 SBPgroup4) ///
	dist(ggamma) time  
	

***Subpop 	
svyset, subpop(race_recode = 1): streg BPgroup2 SBPgroup3 SBPgroup4  lbw Glomdx1NG0 i.income, anc(SBPgroup2 SBPgroup3 SBPgroup4) ///
	anc2(SBPgroup2 SBPgroup3 SBPgroup4) ///
	dist(ggamma) time  

gen SBPgroup2xRace = SBPgroup2*i_race_recode 
	
streg SBPgroup2 SBPgroup3 SBPgroup4 i.race_recode lbw Glomdx1NG0 i.income, ///
	dist(ggamma) time   // gg model with no specified sigma and kappa. 

	
	



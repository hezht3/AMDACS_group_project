*##############################*
*# AMDACS group project       #*
*# Aim: survival analysis     #*
*# Date: October 5, 2022      #*
*##############################*

cd "/Users/zhengtinghe/Library/CloudStorage/OneDrive-JohnsHopkins/Course/340.728.01 - Advanced Methods for Design and Analysis of Cohort Studies/project/AMDACS_group_project"

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

*** Adjusted model
stcox SBPgroup2 SBPgroup3 SBPgroup4 i.race lbw Glomdx1NG0 i.income, nohr   // display parameters
stcox SBPgroup2 SBPgroup3 SBPgroup4 i.race lbw Glomdx1NG0 i.income   // display hazard ratios
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
streg SBPgroup2 SBPgroup3 SBPgroup4 i.race lbw Glomdx1NG0 i.income, dist(weibull) time   // display parameters
lincom _cons + SBPgroup2   // Get beta
lincom _cons + SBPgroup3
lincom _cons + SBPgroup4

streg SBPgroup2 SBPgroup3 SBPgroup4 i.race lbw Glomdx1NG0 i.income, dist(weibull) hr   // display hazard ratios
streg SBPgroup2 SBPgroup3 SBPgroup4 i.race lbw Glomdx1NG0 i.income, dist(weibull) tr   // display time ratios

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

*** Adjusted model
stcox SBPgroup2 SBPgroup3 SBPgroup4 i.race lbw Glomdx1NG0 i.income, tvc(SBPgroup2 SBPgroup3 SBPgroup4) texp((_t - 16)*(_t<=16)) nohr   // display parameters
stcox SBPgroup2 SBPgroup3 SBPgroup4 i.race lbw Glomdx1NG0 i.income, tvc(SBPgroup2 SBPgroup3 SBPgroup4) texp((_t - 16)*(_t<=16))   // display hazard ratios
estat ic
estimates store nphazard_adjust

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

ggtaxonomy
estat ic
estimates store ggconven_crude

*** Adjusted model
streg SBPgroup2 SBPgroup3 SBPgroup4 i.race lbw Glomdx1NG0 i.income, dist(ggamma) time   // Display parameters
lincom _cons + SBPgroup2   // Get beta
lincom _cons + SBPgroup3
lincom _cons + SBPgroup4

streg SBPgroup2 SBPgroup3 SBPgroup4 i.race lbw Glomdx1NG0 i.income, dist(ggamma) tr   // Display time ratio

ggtaxonomy
estat ic
estimates store ggconven_adjust

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
nlcom exp([lnsigma]_cons + [lnsigma]SBPgroup2)
nlcom exp([lnsigma]_cons + [lnsigma]SBPgroup2)

lincom [kappa]_cons + [kappa]SBPgroup2   // Get kappa
lincom [kappa]_cons + [kappa]SBPgroup3
lincom [kappa]_cons + [kappa]SBPgroup4

estat ic
estimates store ggsatura_crude

*** Reproduce crude model stratify by SBPgroup to get ggtaxonomy
***** SBP Group 1
quietly: streg if SBPgroup1 == 1, dist(ggamma)
ggtaxonomy
replace convGGSurv = gammaptail(e(kappa)^(-2),(e(kappa)^(-2))*(exp(-[_t]_cons)*YearsFromCKDTransition)^(e(kappa)/e(sigma)))
replace convGGSurv = 1 - convGGSurv if e(kappa)<0   /*(0 real changes made)*/
sts graph if SBPgroup1 == 1, surv xtitle(Years after CKD diagnosis) ytitle(Probability of no CKD progression) addplot(line convGGSurv YearsFromCKDTransition, sort)
predict hazardSBPgroup1, hazard
***** SBP Group 2
quietly: streg if SBPgroup2 == 1, dist(ggamma)
ggtaxonomy
replace convGGSurv = gammaptail(e(kappa)^(-2),(e(kappa)^(-2))*(exp(-[_t]_cons)*YearsFromCKDTransition)^(e(kappa)/e(sigma)))
replace convGGSurv = 1 - convGGSurv if e(kappa)<0   /*(0 real changes made)*/
sts graph if SBPgroup2 == 1, surv xtitle(Years after CKD diagnosis) ytitle(Probability of no CKD progression) addplot(line convGGSurv YearsFromCKDTransition, sort)
predict hazardSBPgroup2, hazard
***** SBP Group 3
quietly: streg if SBPgroup3 == 1, dist(ggamma)
ggtaxonomy
replace convGGSurv = gammaptail(e(kappa)^(-2),(e(kappa)^(-2))*(exp(-[_t]_cons)*YearsFromCKDTransition)^(e(kappa)/e(sigma)))
replace convGGSurv = 1 - convGGSurv if e(kappa)<0   /*(0 real changes made)*/
sts graph if SBPgroup3 == 1, surv xtitle(Years after CKD diagnosis) ytitle(Probability of no CKD progression) addplot(line convGGSurv YearsFromCKDTransition, sort)
predict hazardSBPgroup3, hazard
***** SBP Group 4
quietly: streg if SBPgroup4 == 1, dist(ggamma)
ggtaxonomy
replace convGGSurv = gammaptail(e(kappa)^(-2),(e(kappa)^(-2))*(exp(-[_t]_cons)*YearsFromCKDTransition)^(e(kappa)/e(sigma)))
replace convGGSurv = 1 - convGGSurv if e(kappa)<0   /*(0 real changes made)*/
sts graph if SBPgroup4 == 1, surv xtitle(Years after CKD diagnosis) ytitle(Probability of no CKD progression) addplot(line convGGSurv YearsFromCKDTransition, sort)
predict hazardSBPgroup4, hazard

*** Adjusted model
streg SBPgroup2 SBPgroup3 SBPgroup4  race lbw Glomdx1NG0 income, ///
	anc(SBPgroup2 SBPgroup3 SBPgroup4) ///
	anc2(SBPgroup2 SBPgroup3 SBPgroup4) ///
	dist(ggamma) time   // Display parameters, did not code race and income as nominal disjoint due to converge issues

lincom [_t]_cons + [_t]SBPgroup2   // Get beta
lincom [_t]_cons + [_t]SBPgroup3
lincom [_t]_cons + [_t]SBPgroup4

nlcom exp([lnsigma]_cons)   // Get sigma
nlcom exp([lnsigma]_cons + [lnsigma]SBPgroup2)
nlcom exp([lnsigma]_cons + [lnsigma]SBPgroup2)
nlcom exp([lnsigma]_cons + [lnsigma]SBPgroup2)

lincom [kappa]_cons + [kappa]SBPgroup2   // Get kappa
lincom [kappa]_cons + [kappa]SBPgroup3
lincom [kappa]_cons + [kappa]SBPgroup4

estat ic
estimates store ggsatura_adjust

******************
* Summary tables *
******************
estimates table cox_crude weibull_crude nphazard_crude ggconven_crude ggsatura_crude, stats(aic)
estimates table cox_adjust weibull_adjust nphazard_adjust ggconven_adjust ggsatura_adjust, stats(aic)

log close

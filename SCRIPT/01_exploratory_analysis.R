##############################
# AMDACS group project       #
# Aim: exploratory analysis  #
# Date: September 9, 2022    #
##############################

setwd("/Users/zhengtinghe/Library/CloudStorage/OneDrive-JohnsHopkins/Course/340.728.01 - Advanced Methods for Design and Analysis of Cohort Studies/project/AMDACS_group_project")

require(tidyverse)
require(survival)
require(survminer)
require(ldatools)
require(gtsummary)

jointfromckd <- read_rds("./INPUT/Cleaned/jointfromckd.rds")


# Generate composite outcome
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Ref: doi: 10.1001/jamanetworkopen.2019.21213
# Composite renal outcome defined as any of:
# 1. Initiation of RRT (dialysis or transplant)
# 2. eGFR < 15 mL/min/1.73 m2
# 3. 50% eGFR reduction from baseline
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Missingness:
# RRTstatusAtTransition: no missing
# gfr: 112 missing
# Composite outcome missingness definition: no RRT & missing GFR, 110 missing
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
jointfromckd <- jointfromckd %>% 
    group_by(id) %>% 
    mutate(renalstatus = case_when(RRTstatusAtTransition == 1 | RRTstatusAtTransition == 2 ~ 1,
                                   gfr < 15 ~ 1,
                                   gfr < 0.5 * gfr[1] ~ 1,
                                   RRTstatusAtTransition == 0 & is.na(gfr) ~ as.numeric(NA),
                                   TRUE ~ 0))   # 304 cases (351 person-periods)


# Exposure group
jointfromckd <- jointfromckd %>% 
    mutate(SBPgroup = cut(SBPpercentile, c(0, 50, 90, 95, 101), right = FALSE))


# Set up time origin & time metric

### Generate visit indicator variable
jointfromckd <- jointfromckd %>% group_by(id) %>% mutate(visit = 1:n())

### Exclude participants not at risk (prevalent cases)
prevalent_cases <- jointfromckd %>% filter(YearsFromBaseline == 0 & renalstatus == 1) %>% pull(id)
jointfromckd <- jointfromckd %>% filter(!(id %in% prevalent_cases))

### Drop participants with missing data
jointfromckd <- jointfromckd %>% drop_na(renalstatus, SBPgroup)

### Drop participants with problematic follow-up time (end time < start time)
jointfromckd <- jointfromckd %>% 
    filter(YearsFromBaselineTransition > YearsFromBaseline |
           YearsFromCKDTransition > YearsFromCKD)   # 1 participant

### Examine "gaps"
jointfromckd %>% group_by(id) %>% mutate(visit_gap = visit - lag(visit)) %>% ungroup() %>% pull(visit_gap) %>% table()
#    1    2    3    4    8 
# 3572  104    9    1    1

### Construct risk sets
jointfromckd <- jointfromckd %>% 
    mutate(SurvObj_FromCKD = Surv(YearsFromCKD, YearsFromCKDTransition, event = renalstatus)) %>% 
    mutate(SurvObj_FromBaseline = Surv(YearsFromBaseline, YearsFromBaselineTransition, event = renalstatus))


# Kaplan-Meier estimates of survival function

### From CKD
fit_FromCKD <- survfit(SurvObj_FromCKD ~ SBPgroup, data = jointfromckd)
summary(fit_FromCKD)$table
tbl_survfit(fit_FromCKD, times = seq(5, 20, 5))

tiff("./OUTPUT/Exploratory_analysis/Figure 1. survival curve (from CKD).tiff",
     width = 3000, height = 2000, pointsize = 10, res = 300)
ggsurvplot(fit_FromCKD, data = jointfromckd,
           palette = "jama", linetype = 1, censor.size = 0.8,
           xlab = "Duration of follow-up (in years)", ylab = "Probability of survival",
           break.x.by = 5, #conf.int = TRUE,
           legend.title = "SBP percentile", legend.labs = c("[0,50)", "[50,90)", "[90,95)", "[95,100)"),
           risk.table = TRUE, risk.table.title = "No. at risk",
           risk.table.height = 0.15, cumevents.height = 0.15, risk.table.x.text = FALSE, cumevents.x.text = FALSE,
           font.main = c(13, "bold", "black"), font.x = c(13, "bold", "black"), font.y = c(13, "bold", "black"),
           font.tickslab = c(13, "bold", "black"), font.legend = c(13, "bold", "black"),
           tables.theme = theme_cleantable())
dev.off()

tiff("./OUTPUT/Exploratory_analysis/Figure 2. cumulative hazard (from CKD).tiff",
     width = 3000, height = 2000, pointsize = 10, res = 300)
ggsurvplot(fit_FromCKD, data = jointfromckd, fun = "cumhaz", #conf.int = TRUE,
           palette = "jama", linetype = 1, censor.size = 0.8,
           xlab = "Duration of follow-up (in years)", ylab = "Cumulative hazard",
           legend.title = "SBP percentile", legend.labs = c("[0,50)", "[50,90)", "[90,95)", "[95,100)"))
dev.off()

tiff("./OUTPUT/Exploratory_analysis/Figure 3. cloglog (from CKD).tiff",
     width = 3000, height = 2000, pointsize = 10, res = 300)
ggsurvplot(fit_FromCKD, data = jointfromckd, fun = "cloglog", #conf.int = TRUE,
           palette = "jama", linetype = 1, censor.size = 0.8,
           xlab = "Duration of follow-up (in years)",
           legend.title = "SBP percentile", legend.labs = c("[0,50)", "[50,90)", "[90,95)", "[95,100)"))
dev.off()

cox_FromCKD <- coxph(SurvObj_FromCKD ~ SBPgroup, data = jointfromckd)
cox_FromCKD %>% tbl_regression(exponentiate = TRUE)
tiff("./OUTPUT/Exploratory_analysis/Figure 4. cox-snell residuals (from CKD).tiff",
     width = 3000, height = 2000, pointsize = 10, res = 300)
gg_coxsnell(cox_FromCKD) + geom_abline(slope = 1, intercept = 0) + theme_minimal()
dev.off()


### From baseline
fit_FromBaseline <- survfit(SurvObj_FromBaseline ~ SBPgroup, data = jointfromckd)
summary(fit_FromBaseline)$table
tbl_survfit(fit_FromBaseline, times = seq(5, 20, 5))

tiff("./OUTPUT/Exploratory_analysis/Figure 1. survival curve (from Baseline).tiff",
     width = 3000, height = 2000, pointsize = 10, res = 300)
ggsurvplot(fit_FromBaseline, data = jointfromckd,
           palette = "jama", linetype = 1, censor.size = 0.8,
           xlab = "Duration of follow-up (in years)", ylab = "Probability of survival",
           break.x.by = 2.5, #conf.int = TRUE,
           legend.title = "SBP percentile", legend.labs = c("[0,50)", "[50,90)", "[90,95)", "[95,100)"),
           risk.table = TRUE, risk.table.title = "No. at risk",
           risk.table.height = 0.15, cumevents.height = 0.15, risk.table.x.text = FALSE, cumevents.x.text = FALSE,
           font.main = c(13, "bold", "black"), font.x = c(13, "bold", "black"), font.y = c(13, "bold", "black"),
           font.tickslab = c(13, "bold", "black"), font.legend = c(13, "bold", "black"),
           tables.theme = theme_cleantable())
dev.off()

tiff("./OUTPUT/Exploratory_analysis/Figure 2. cumulative hazard (from Baseline).tiff",
     width = 3000, height = 2000, pointsize = 10, res = 300)
ggsurvplot(fit_FromBaseline, data = jointfromckd, fun = "cumhaz", #conf.int = TRUE,
           palette = "jama", linetype = 1, censor.size = 0.8,
           xlab = "Duration of follow-up (in years)", ylab = "Cumulative hazard",
           legend.title = "SBP percentile", legend.labs = c("[0,50)", "[50,90)", "[90,95)", "[95,100)"))
dev.off()

tiff("./OUTPUT/Exploratory_analysis/Figure 3. cloglog (from Baseline).tiff",
     width = 3000, height = 2000, pointsize = 10, res = 300)
ggsurvplot(fit_FromBaseline, data = jointfromckd, fun = "cloglog", #conf.int = TRUE,
           palette = "jama", linetype = 1, censor.size = 0.8,
           xlab = "Duration of follow-up (in years)",
           legend.title = "SBP percentile", legend.labs = c("[0,50)", "[50,90)", "[90,95)", "[95,100)"))
dev.off()

cox_FromBaseline <- coxph(SurvObj_FromBaseline ~ SBPgroup, data = jointfromckd)
cox_FromBaseline %>% tbl_regression(exponentiate = TRUE)
tiff("./OUTPUT/Exploratory_analysis/Figure 4. cox-snell residuals (from Baseline).tiff",
     width = 3000, height = 2000, pointsize = 10, res = 300)
gg_coxsnell(cox_FromBaseline) + geom_abline(slope = 1, intercept = 0) + theme_minimal()
dev.off()

##############################
# AMDACS group project       #
# Aim: exploratory analysis  #
# Date: September 9, 2022    #
##############################

setwd("/Users/zhengtinghe/Library/CloudStorage/OneDrive-JohnsHopkins/Course/340.728.01 - Advanced Methods for Design and Analysis of Cohort Studies/project/AMDACS_group_project")

require(tidyverse)
require(haven)
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


# Set up time origin & time metric

### Generate visit variable
jointfromckd <- jointfromckd %>% group_by(id) %>% mutate(visit = 1:n())

### Exclude participants not at risk (prevalent cases)
prevalent_cases <- jointfromckd %>% filter(YearsFromBaseline == 0 & renalstatus == 1) %>% pull(id)
jointfromckd <- jointfromckd %>% filter(!(id %in% prevalent_cases))
rm(prevalent_cases)

### Resolve person-periods with problematic follow-up time (end time < start time)
### Based on suggestions from office hours
jointfromckd %>% 
    filter(id == "4148" | id == "5653") %>% 
    select(id, contains("Years"), VisitAtTransition)
jointfromckd <- jointfromckd %>% 
    mutate(YearsFromCKD = ifelse(id == "4148" & VisitAtTransition == "4",
                                 (15.3 + 17.1)/2, YearsFromCKD)) %>% 
    mutate(YearsFromBaseline = ifelse(id == "4148" & VisitAtTransition == "4",
                                      (1.2 + 2.97)/2, YearsFromBaseline)) %>% 
    mutate(YearsFromCKD = ifelse(id == "5653" & VisitAtTransition == "5",
                                 (2.60 + 4.48)/2, YearsFromCKD)) %>% 
    mutate(YearsFromBaseline = ifelse(id == "5653" & VisitAtTransition == "5",
                                      (1.93 + 3.81)/2, YearsFromBaseline)) %>% 
    mutate(YearsFromCKDTransition = ifelse(id == "4148" & VisitAtTransition == "3",
                                           (15.3 + 17.1)/2, YearsFromCKDTransition)) %>% 
    mutate(YearsFromBaselineTransition = ifelse(id == "4148" & VisitAtTransition == "3",
                                                (1.2 + 2.97)/2, YearsFromBaselineTransition)) %>% 
    mutate(YearsFromCKDTransition = ifelse(id == "5653" & VisitAtTransition == "4",
                                           (2.60 + 4.48)/2, YearsFromCKDTransition)) %>% 
    mutate(YearsFromBaselineTransition = ifelse(id == "5653" & VisitAtTransition == "4",
                                                (1.93 + 3.81)/2, YearsFromBaselineTransition))

### Drop person-periods after event
jointfromckd <- jointfromckd %>% 
    group_by(id) %>% 
    mutate(firstevent = min(which(renalstatus == 1 | row_number() == n()))) %>%
    filter(row_number() <= firstevent) %>% 
    mutate(firstevent = NULL)


#############################
# Data set 1: complete case #
#############################

# Complete case in exposure, covariates, outcome
jointfromckd_complete <- jointfromckd %>% 
    drop_na(renalstatus, SBPpercentile, DBPpercentile,
            creatinine:income, Heightpercentile:BMIpercentile, AgeAtCKD:Glomdx1NG0)

### Examine "gaps": 1 - no gaps, 2 - 1 gap
jointfromckd_complete %>% group_by(id) %>% mutate(visit_gap = visit - lag(visit)) %>% ungroup() %>% pull(visit_gap) %>% table()
#    1    2    3    4    5    6 
# 1766  346   60    7    4    1 

### Exposure group
jointfromckd_complete <- jointfromckd_complete %>% 
    mutate(SBPgroup = cut(SBPpercentile, c(0, 50, 90, 95, 101), right = FALSE)) %>% 
    mutate(DBPgroup = cut(DBPpercentile, c(0, 50, 90, 95, 101), right = FALSE)) %>% 
    mutate(across(contains("group"),
                  ~ case_when(.x == "[0,50)" ~ 1,
                              .x == "[50,90)" ~ 2,
                              .x == "[90,95)" ~ 3,
                              .x == "[95,101)" ~ 4))) %>% 
    mutate(BPgroup = max(SBPgroup, DBPgroup)) %>% 
    mutate(visit = 1:n())

write_rds(jointfromckd_complete, file = "./INPUT/Cleaned/jointfromckd_survival_complete.rds")
write_dta(jointfromckd_complete, path = "./INPUT/Cleaned/jointfromckd_survival_complete.dta")


####################################
# Data set 2: LOCF & outcome rules #
####################################

### Missing values in event status
###>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
### Summary statistics
table(is.na(jointfromckd$renalstatus))   # 105 person-periods have missing event status
jointfromckd %>% filter(is.na(renalstatus)) %>% distinct(id) %>% nrow()   # 85 participants with missing event status
jointfromckd %>% 
    group_by(id) %>% 
    filter(row_number() == n()) %>% 
    filter(is.na(renalstatus)) %>% 
    nrow()   # 53 participants' missing event status is in last person-periods
###<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
###>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
### Aims: For person-periods of each participant
### 1. Drop person-periods after last person-period with non missing event status (treat as right censored)
### 2. Fill event status of person-periods with missing event status before last person-period with non missing event status as 0
### Rules:
### 1. Group by `id`
### 2. Generate `indi` as:
###    1.1 If event status is missing for a person-period, `indi` = 0
###    1.2 If event status is not missing for a person-period, `indi` = 1
### 3. Generate `rev_cumsum_indi` as cumulative sum of `indi` from last row to first row
###    3.1 For person-periods before last person-period with non missing event status
###     & last person-period with non missing event status, `rev_cumsum_indi` > 0
###    3.2 For person-period after last person-period with non missing event status, `rev_cumsum_indi` = 0
### 4. Filter person-period with `rev_cumsum_indi` > 0
### 5. Change missing `renalstatus` for rest person-periods to 0
###<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

jointfromckd_locf <- jointfromckd %>% 
    group_by(id) %>% 
    mutate(indi = ifelse(is.na(renalstatus), 0, 1)) %>% 
    mutate(rev_cumsum_indi = rev(cumsum(rev(indi)))) %>% 
    filter(rev_cumsum_indi > 0) %>% 
    mutate(renalstatus = ifelse(is.na(renalstatus), 0, renalstatus)) %>% 
    mutate(indi = NULL, rev_cumsum_indi = NULL)

### Last observation carried 1 cycle forward for missing exposures & covariates
jointfromckd_locf <- jointfromckd_locf %>% 
    group_by(id) %>% 
    fill(- contains("Years"), - gfr, - RRTstatusAtTransition, - VisitAtTransition, - id,
         .direction = "downup")

### Examine "gaps": 1 - no gaps, 2 - 1 gap
jointfromckd_locf %>% group_by(id) %>% mutate(visit_gap = visit - lag(visit)) %>% ungroup() %>% pull(visit_gap) %>% table()
###### no gaps

### Exposure group
jointfromckd_locf <- jointfromckd_locf %>% 
    mutate(SBPgroup = cut(SBPpercentile, c(0, 50, 90, 95, 101), right = FALSE)) %>% 
    mutate(DBPgroup = cut(DBPpercentile, c(0, 50, 90, 95, 101), right = FALSE)) %>% 
    mutate(across(contains("group"),
                  ~ case_when(.x == "[0,50)" ~ 1,
                              .x == "[50,90)" ~ 2,
                              .x == "[90,95)" ~ 3,
                              .x == "[95,101)" ~ 4))) %>% 
    mutate(BPgroup = max(SBPgroup, DBPgroup)) %>% 
    mutate(visit = 1:n())

jointfromckd_locf <- jointfromckd_locf %>% 
    drop_na(BPgroup,
            creatinine:income, Heightpercentile:BMIpercentile, AgeAtCKD:Glomdx1NG0)

write_rds(jointfromckd_locf, file = "./INPUT/Cleaned/jointfromckd_survival_locf.rds")
write_dta(jointfromckd_locf, path = "./INPUT/Cleaned/jointfromckd_survival_locf.dta")



#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Kaplan-Meier estimates of survival function

### From CKD
tiff("./OUTPUT/Exploratory_analysis/Figure 0. spaghetti plot (from CKD).tiff",
     width = 3000, height = 2000, pointsize = 10, res = 300)
jointfromckd %>% 
    ungroup() %>% 
    ggplot(aes(x = YearsFromCKDTransition, y = SBPpercentile)) +
    geom_line(aes(group = id), size = 0.1) +
    xlab("Duration of follow-up (in years)") +
    ylab("SBP percentile") +
    theme_minimal()
dev.off()

tiff("./OUTPUT/Exploratory_analysis/Figure 0. spaghetti plot (from CKD)_mean.tiff",
     width = 3000, height = 2000, pointsize = 10, res = 300)
jointfromckd %>% 
    ungroup() %>% 
    ggplot(aes(x = YearsFromCKDTransition, y = SBPpercentile)) +
    geom_line(aes(group = id), size = 0.1, alpha = 0.3) +
    geom_smooth() +
    xlab("Duration of follow-up (in years)") +
    ylab("SBP percentile") +
    theme_minimal()
dev.off()

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
tiff("./OUTPUT/Exploratory_analysis/Figure 0. spaghetti plot (from Baseline).tiff",
     width = 3000, height = 2000, pointsize = 10, res = 300)
jointfromckd %>% 
    ungroup() %>% 
    ggplot(aes(x = YearsFromBaselineTransition, y = SBPpercentile)) +
    geom_line(aes(group = id), size = 0.1) +
    xlab("Duration of follow-up (in years)") +
    ylab("SBP percentile") +
    theme_minimal()
dev.off()

tiff("./OUTPUT/Exploratory_analysis/Figure 0. spaghetti plot (from Baseline)_mean.tiff",
     width = 3000, height = 2000, pointsize = 10, res = 300)
jointfromckd %>% 
    ungroup() %>% 
    ggplot(aes(x = YearsFromBaselineTransition, y = SBPpercentile)) +
    geom_line(aes(group = id), size = 0.1, alpha = 0.3) +
    geom_smooth() +
    xlab("Duration of follow-up (in years)") +
    ylab("SBP percentile") +
    theme_minimal()
dev.off()

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

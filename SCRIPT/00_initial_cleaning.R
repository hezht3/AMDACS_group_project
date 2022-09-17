###########################
# AMDACS group project    #
# Aim: initial cleaning   #
# Date: September 9, 2022 #
###########################

setwd("/Users/zhengtinghe/Library/CloudStorage/OneDrive-JohnsHopkins/Course/340.728.01 - Advanced Methods for Design and Analysis of Cohort Studies/project/AMDACS_group_project")

require(tidyverse)
require(haven)
require(data.table)


##############
# ckdtoevent #
##############

# Read in dataset
ckdtoevent <- read.table("./INPUT/Origin/ckdtoevent.dat",
                         col.names = c("id", "male1female0", "race", "lbw", "premature", "sga", "icu", "Glomdx1NG0",
                                       "AgeAtCKD", "AgeAtStudyEntry", "gfr", "creatinine", "cystatin", "PrCrRatio",
                                       "anemia", "income", "MaternalEDU", "SBPpercentile", "DBPpercentile",
                                       "Heightpercentile", "Weightpercentile", "BMIpercentile", "QofLbyParent",
                                       "QofLbyChild", "IQ", "StudyEntrytoRRTorExit", "Trans2Dial1noRRT0atExit"))

# Change `.` to `NA`, variable class
ckdtoevent <- ckdtoevent %>%
    mutate(across(everything(), ~ ifelse(.x == ".", as.numeric(NA), .x))) %>% 
    mutate(across(everything(), ~ as.numeric(.x))) %>% 
    mutate(id = as.character(id)) %>% 
    mutate(across(c("male1female0", "race", "lbw", "premature", "sga", "icu", "Glomdx1NG0", "anemia", "income", "MaternalEDU", "Trans2Dial1noRRT0atExit"),
                  ~ factor(.x)))

# Add variable labels
ckdtoevent <- ckdtoevent %>% 
    mutate_at(vars(id), funs(setattr(., "label", "ID number of study participant"))) %>% 
    mutate_at(vars(male1female0), funs(setattr(., "label", "Sex"))) %>% 
    mutate_at(vars(race), funs(setattr(., "label", "race"))) %>% 
    mutate_at(vars(lbw), funs(setattr(., "label", "low birth weight"))) %>% 
    mutate_at(vars(premature), funs(setattr(., "label", "Prematurity"))) %>% 
    mutate_at(vars(sga), funs(setattr(., "label", "small for gestational age"))) %>% 
    mutate_at(vars(icu), funs(setattr(., "label", "intensive care unit at birth"))) %>% 
    mutate_at(vars(Glomdx1NG0), funs(setattr(., "label", "primary chronic kidney disease diagnosis"))) %>% 
    mutate_at(vars(AgeAtCKD), funs(setattr(., "label", "age at the time of CKD onset"))) %>% 
    mutate_at(vars(AgeAtStudyEntry), funs(setattr(., "label", "age at the time of the baseline study visit"))) %>% 
    mutate_at(vars(gfr), funs(setattr(., "label", "baseline glomerular filtration rate (ml/min per 1.73m2)"))) %>% 
    mutate_at(vars(creatinine), funs(setattr(., "label", "baseline serum creatinine (mg/dL)"))) %>% 
    mutate_at(vars(cystatin), funs(setattr(., "label", "baseline cystatin C (mg/L)"))) %>% 
    mutate_at(vars(PrCrRatio), funs(setattr(., "label", "baseline urine protein / urine creatinine"))) %>% 
    mutate_at(vars(anemia), funs(setattr(., "label", "baseline haemoglobin < age-sex-race specific 5th percentile"))) %>% 
    mutate_at(vars(income), funs(setattr(., "label", "baseline household income"))) %>% 
    mutate_at(vars(MaternalEDU), funs(setattr(., "label", "baseline level of education of participant’s mom"))) %>% 
    mutate_at(vars(SBPpercentile), funs(setattr(., "label", "baseline age‐sex‐height specific systolic BP percentile"))) %>% 
    mutate_at(vars(DBPpercentile), funs(setattr(., "label", "baseline age‐sex‐height specific diastolic BP percentile"))) %>% 
    mutate_at(vars(Heightpercentile), funs(setattr(., "label", "baseline age‐sex specific height percentile"))) %>% 
    mutate_at(vars(Heightpercentile), funs(setattr(., "label", "baseline age‐sex specific height percentile"))) %>% 
    mutate_at(vars(Weightpercentile), funs(setattr(., "label", "baseline age‐sex specific weight percentile"))) %>% 
    mutate_at(vars(BMIpercentile), funs(setattr(., "label", "baseline age‐sex specific body mass index percentile"))) %>% 
    mutate_at(vars(QofLbyParent), funs(setattr(., "label", "baseline quality of life score – parent report"))) %>% 
    mutate_at(vars(QofLbyChild), funs(setattr(., "label", "baseline quality of life score – participant report"))) %>% 
    mutate_at(vars(IQ), funs(setattr(., "label", "baseline intelligent quotient score"))) %>% 
    mutate_at(vars(StudyEntrytoRRTorExit), funs(setattr(., "label", "total years as RRT‐free after baseline study visit"))) %>% 
    mutate_at(vars(Trans2Dial1noRRT0atExit), funs(setattr(., "label", "renal replacement therapy (RRT) status at exit")))

# Write out dataset
write_rds(ckdtoevent, file = "./INPUT/Cleaned/ckdtoevent.rds")
write_dta(ckdtoevent, path = "./INPUT/Cleaned/ckdtoevent.dta")


################
# jointfromckd #
################

# Read in dataset
jointfromckd <- read.table("./INPUT/Origin/jointfromckd.dat",
                          col.names = c("YearsFromCKD", "YearsFromBaseline", "gfr", "creatinine", "cystatin", "PrCrRatio", "anemia", "income",
                                        "SBPpercentile", "DBPpercentile", "Heightpercentile", "Weightpercentile", "BMIpercentile",
                                        "YearsFromCKDTransition", "YearsFromBaselineTransition", "RRTstatusAtTransition", "VisitAtTransition",
                                        "id", "AgeAtCKD", "male1female0", "race", "lbw", "premature", "sga", "icu", "Glomdx1NG0"))

# Change `.` to `NA`, variable class
jointfromckd <- jointfromckd %>%
    mutate(across(everything(), ~ ifelse(.x == ".", as.numeric(NA), .x))) %>% 
    mutate(across(everything(), ~ as.numeric(.x))) %>% 
    mutate(id = as.character(id)) %>% 
    mutate(across(c("anemia", "income", "RRTstatusAtTransition", "VisitAtTransition", "male1female0", "race", "lbw", "premature", "sga", "icu",
                    "Glomdx1NG0"),
                  ~ factor(.x)))

# Add variable labels
jointfromckd <- jointfromckd %>% 
    mutate_at(vars(YearsFromCKD), funs(setattr(., "label", "years from CKD"))) %>% 
    mutate_at(vars(YearsFromBaseline), funs(setattr(., "label", "years from baseline visit"))) %>% 
    mutate_at(vars(gfr), funs(setattr(., "label", "glomerular filtration rate (ml/min per 1.73m2)"))) %>% 
    mutate_at(vars(creatinine), funs(setattr(., "label", "serum creatinine (mg/dL)"))) %>% 
    mutate_at(vars(cystatin), funs(setattr(., "label", "cystatin C (mg/L)"))) %>% 
    mutate_at(vars(PrCrRatio), funs(setattr(., "label", "urine protein:creatinine"))) %>% 
    mutate_at(vars(anemia), funs(setattr(., "label", "haemoglobin < age-sex-race specific 5th percentile"))) %>% 
    mutate_at(vars(income), funs(setattr(., "label", "household income"))) %>% 
    mutate_at(vars(SBPpercentile), funs(setattr(., "label", "systolic blood pressure percentile (age‐gender‐height specific)"))) %>% 
    mutate_at(vars(DBPpercentile), funs(setattr(., "label", "diastolic blood pressure percentile (age‐gender‐height specific)"))) %>% 
    mutate_at(vars(Heightpercentile), funs(setattr(., "label", "height percentile (age‐gender specific)"))) %>% 
    mutate_at(vars(Weightpercentile), funs(setattr(., "label", "weight percentile (age‐gender specific)"))) %>% 
    mutate_at(vars(BMIpercentile), funs(setattr(., "label", "body mass index percentile (age‐gender specific)"))) %>% 
    mutate_at(vars(YearsFromCKDTransition), funs(setattr(., "label", "Years from CKD at the end of the visit window (transition)"))) %>% 
    mutate_at(vars(YearsFromBaselineTransition), funs(setattr(., "label", "Years from baseline visit at the end of the visit window (transition)"))) %>% 
    mutate_at(vars(RRTstatusAtTransition), funs(setattr(., "label", "RRT status at the end of the visit window (transition)"))) %>% 
    mutate_at(vars(VisitAtTransition), funs(setattr(., "label", "visit number of next visit"))) %>% 
    mutate_at(vars(id), funs(setattr(., "label", "ID number of study participant"))) %>% 
    mutate_at(vars(AgeAtCKD), funs(setattr(., "label", "age at the time of CKD onset"))) %>% 
    mutate_at(vars(male1female0), funs(setattr(., "label", "Sex of partipant"))) %>% 
    mutate_at(vars(race), funs(setattr(., "label", "Race of participant"))) %>% 
    mutate_at(vars(lbw), funs(setattr(., "label", "low birth weight"))) %>% 
    mutate_at(vars(premature), funs(setattr(., "label", "prematurity"))) %>% 
    mutate_at(vars(sga), funs(setattr(., "label", "small for gestational age"))) %>% 
    mutate_at(vars(icu), funs(setattr(., "label", "intensive care unit at birth"))) %>% 
    mutate_at(vars(Glomdx1NG0), funs(setattr(., "label", "primary chronic kidney disease diagnosis")))

# Write out dataset
write_rds(jointfromckd, file = "./INPUT/Cleaned/jointfromckd.rds")
write_dta(jointfromckd, path = "./INPUT/Cleaned/jointfromckd.dta")

##############################
# AMDACS group project       #
# Aim: descriptive analysis  #
# Date: October 23, 2022     #
##############################

setwd("/Users/zhengtinghe/Library/CloudStorage/OneDrive-JohnsHopkins/Course/340.728.01 - Advanced Methods for Design and Analysis of Cohort Studies/project/AMDACS_group_project")

require(tidyverse)
require(haven)
require(survival)
require(survminer)
require(gtsummary)

jointfromckd <- read_dta("./OUTPUT/data_for_plot.dta")

jointfromckd %>%
    group_by(id) %>%
    mutate(sur_diff = YearsFromCKDTransition - YearsFromCKD) %>%
    mutate(sur_time = sum(sur_diff)) %>%
    ungroup() %>%
    distinct(id, .keep_all = TRUE) %>%
    summarise(min = min(sur_time),
              median = median(sur_time),
              `25th` = quantile(sur_time, 0.25),
              `75th` = quantile(sur_time, 0.75),
              max = max(sur_time))

pyears(Surv(YearsFromCKD, YearsFromCKDTransition, renalstatus) ~ 1, data = jointfromckd, scale = 1)

jointfromckd %>%
    group_by(id) %>%
    mutate(id_num = 1:n()) %>% 
    filter(id_num == n()) %>%
    ungroup() %>%
    summarise(event = sum(renalstatus)/n())

tbl_survfit(survfit(Surv(YearsFromCKD, YearsFromCKDTransition, renalstatus) ~ 1,
                    data = jointfromckd), times = c(5, 10, 15, 20, 25, 30), reverse = TRUE)
tbl_survfit(survfit(Surv(YearsFromCKD, YearsFromCKDTransition, renalstatus) ~ SBPgroup,
                    data = jointfromckd), times = c(5, 10, 15, 20, 25, 30), reverse = TRUE)

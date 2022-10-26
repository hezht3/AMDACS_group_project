##############################
# AMDACS group project       #
# Aim: plotting              #
# Date: October 23, 2022     #
##############################

setwd("/Users/zhengtinghe/Library/CloudStorage/OneDrive-JohnsHopkins/Course/340.728.01 - Advanced Methods for Design and Analysis of Cohort Studies/project/AMDACS_group_project")

require(tidyverse)
require(haven)
require(survival)
require(survminer)
require(flexsurv)
require(patchwork)

jointfromckd <- read_dta("./OUTPUT/data_for_plot.dta")
bootstrap <- read_dta("./OUTPUT/500saturatedGG.dta")
bootstrap_adjusted <- read_dta("./OUTPUT/500saturatedGG_adjusted.dta")


################
# Crude models #
################

###########################
# Figure 1. Survival plot #
###########################

fit <- survfit(Surv(YearsFromCKD, YearsFromCKDTransition, renalstatus) ~ SBPgroup,
               data = jointfromckd)
summary(fit)$table

tiff("./OUTPUT/Figures/SFigure 4. cloglog.tiff",
     width = 3000, height = 2300, pointsize = 10, res = 300)
ggsurvplot(fit, data = jointfromckd, fun = "event",
           palette = "jama", linetype = 1, censor.size = 0.8,
           xlab = "Years after CKD diagnosis", ylab = "cloglog transformation of survival from CKD progression",
           break.x.by = 5, #conf.int = TRUE,
           legend.title = "SBP percentile", legend.labs = c("<50th", "50th-90th", "90th-95th", ">95th"),
           # risk.table = TRUE, risk.table.title = "No. at risk",
           # risk.table.height = 0.15, cumevents.height = 0.15, risk.table.x.text = FALSE, cumevents.x.text = FALSE,
           font.main = c(13, "bold", "black"), font.x = c(13, "bold", "black"), font.y = c(13, "bold", "black"),
           font.tickslab = c(13, "bold", "black"), font.legend = c(13, "bold", "black"),
           tables.theme = theme_cleantable())
dev.off()

tiff("./OUTPUT/Figures/SFigure 4. cumhaz.tiff",
     width = 3000, height = 2300, pointsize = 10, res = 300)
ggsurvplot(fit, data = jointfromckd, fun = "cumhaz",
           palette = "jama", linetype = 1, censor.size = 0.8,
           xlab = "Years after CKD diagnosis", ylab = "cloglog transformation of survival from CKD progression",
           break.x.by = 5, #conf.int = TRUE,
           legend.title = "SBP percentile", legend.labs = c("<50th", "50th-90th", "90th-95th", ">95th"),
           # risk.table = TRUE, risk.table.title = "No. at risk",
           # risk.table.height = 0.15, cumevents.height = 0.15, risk.table.x.text = FALSE, cumevents.x.text = FALSE,
           font.main = c(13, "bold", "black"), font.x = c(13, "bold", "black"), font.y = c(13, "bold", "black"),
           font.tickslab = c(13, "bold", "black"), font.legend = c(13, "bold", "black"),
           tables.theme = theme_cleantable())
dev.off()

tiff("./OUTPUT/Figures/SFigure 4. schoenfeld residuals.tiff",
     width = 3000, height = 2300, pointsize = 10, res = 300)
ggcoxdiagnostics(coxph(Surv(YearsFromCKD, YearsFromCKDTransition, renalstatus) ~ SBPgroup,
                       data = jointfromckd), type = "schoenfeld")
dev.off()

fit <- survfit(Surv(YearsFromCKD, YearsFromCKDTransition, renalstatus) ~ 1,
               data = jointfromckd)
tiff("./OUTPUT/Figures/SFigure 5. size of risk sets.tiff",
     width = 1500, height = 1000, pointsize = 10, res = 300)
tibble(n.size = summary(fit)$n.risk, time = summary(fit)$time) %>% 
    ggplot(aes(x = time, y = n.size)) +
    geom_point() +
    xlab("Years from CKD diagnosis") +
    ylab("Size of risk sets") +
    theme_minimal()
dev.off()

#############################
# Figure 2. Relative hazard #
#############################

ttx = seq(0, 30, 1)
hazardSBP1 = hgengamma(ttx,3.262489,.5965208,1)
hazardSBP2 = hgengamma(ttx,3.1287,.410369,1.687639)
hazardSBP3 = hgengamma(ttx,2.818067,.4477858,1/.4477858)
hazardSBP4 = hgengamma(ttx,2.39798,.5708262,.5708262)
hazardratio <- tibble(hazardSBP1 = hazardSBP1, hazardSBP2 = hazardSBP2,
                      hazardSBP3 = hazardSBP3, hazardSBP4 = hazardSBP4) %>% 
    mutate(hazardSBP2to1 = hazardSBP2/hazardSBP1,
           hazardSBP3to1 = hazardSBP3/hazardSBP1,
           hazardSBP4to1 = hazardSBP4/hazardSBP1) %>% 
    mutate(time = seq(0, 30, 1)) %>% 
    select(time, contains("to1")) %>% 
    pivot_longer(cols = -time, names_to = "HR", values_to = "value") %>% 
    mutate(HR = factor(HR, levels = c("hazardSBP2to1", "hazardSBP3to1", "hazardSBP4to1"),
                       labels = c("SBP percentile 50-90 vs 0-50",
                                  "SBP percentile 90-95 vs 0-50",
                                  "SBP percentile 95-100 vs 0-50")))

# SBPgroup2 vs SBPgroup1
ninetyCISBP2toSBP1at123yes = bootstrap %>%
    mutate(mu1 = beta_SBP1, sigma1 = sigma_SBP1, Q1 = kappa_SBP1,
           mu2 = beta_SBP2, sigma2 = sigma_SBP2, Q2 = kappa_SBP2) %>%
    mutate(rh5 = hgengamma(5,mu = mu2,sigma = sigma2, Q = Q2) / hgengamma(5,mu = mu1,sigma = sigma1, Q = Q1), # get the ratio
           rh10 = hgengamma(10,mu = mu2,sigma = sigma2, Q = Q2) / hgengamma(10,mu = mu1,sigma = sigma1, Q = Q1),
           rh15 = hgengamma(15,mu = mu2,sigma = sigma2, Q = Q2) / hgengamma(15,mu = mu1,sigma = sigma1, Q = Q1),
           rh20 = hgengamma(20,mu = mu2,sigma = sigma2, Q = Q2) / hgengamma(20,mu = mu1,sigma = sigma1, Q = Q1),
           rh25 = hgengamma(25,mu = mu2,sigma = sigma2, Q = Q2) / hgengamma(25,mu = mu1,sigma = sigma1, Q = Q1),
           rh30 = hgengamma(30,mu = mu2,sigma = sigma2, Q = Q2) / hgengamma(30,mu = mu1,sigma = sigma1, Q = Q1)) %>%
    select(rh5,rh10,rh15,rh20,rh25,rh30)
hazardratio %>% 
    filter(HR == "SBP percentile 50-90 vs 0-50") %>% 
    ggplot(aes(x = time, y = value)) +
    geom_line(color = "#ef8632") +
    geom_point(color = "#ef8632") +
    geom_errorbar(aes(x = 5,
                      ymin = quantile(ninetyCISBP2toSBP1at123yes[,1], probs=c(.05,.95), na.rm = TRUE)[[1]],
                      ymax = quantile(ninetyCISBP2toSBP1at123yes[,1], probs=c(.05,.95), na.rm = TRUE)[[2]]),
                  color = "#ef8632") +
    geom_errorbar(aes(x = 10,
                      ymin = quantile(ninetyCISBP2toSBP1at123yes[,2], probs=c(.05,.95), na.rm = TRUE)[[1]],
                      ymax = quantile(ninetyCISBP2toSBP1at123yes[,2], probs=c(.05,.95), na.rm = TRUE)[[2]]),
                  color = "#ef8632") +
    geom_errorbar(aes(x = 15,
                      ymin = quantile(ninetyCISBP2toSBP1at123yes[,3], probs=c(.05,.95), na.rm = TRUE)[[1]],
                      ymax = quantile(ninetyCISBP2toSBP1at123yes[,3], probs=c(.05,.95), na.rm = TRUE)[[2]]),
                  color = "#ef8632") +
    geom_errorbar(aes(x = 20,
                      ymin = quantile(ninetyCISBP2toSBP1at123yes[,4], probs=c(.05,.95), na.rm = TRUE)[[1]],
                      ymax = quantile(ninetyCISBP2toSBP1at123yes[,4], probs=c(.05,.95), na.rm = TRUE)[[2]]),
                  color = "#ef8632") +
    geom_errorbar(aes(x = 25,
                      ymin = quantile(ninetyCISBP2toSBP1at123yes[,5], probs=c(.05,.95), na.rm = TRUE)[[1]],
                      ymax = quantile(ninetyCISBP2toSBP1at123yes[,5], probs=c(.05,.95), na.rm = TRUE)[[2]]),
                  color = "#ef8632") +
    geom_hline(yintercept = 1, color = "#0000f5") +
    scale_color_manual(values = c("#0389ff", "#0dd3ff", "#1c79c0"), guide = "none") +
    ggtitle("SBP percentile 50-90 vs 0-50") +
    xlab("Years after CKD diagnosis") +
    ylab("Relative hazard") +
    theme_minimal() -> p1

# SBPgroup3 vs SBPgroup1
ninetyCISBP3toSBP1at123yes = bootstrap %>%
    mutate(mu1 = beta_SBP1, sigma1 = sigma_SBP1, Q1 = kappa_SBP1,
           mu2 = beta_SBP3, sigma2 = sigma_SBP3, Q2 = kappa_SBP3) %>%
    mutate(rh5 = hgengamma(5,mu = mu2,sigma = sigma2, Q = Q2) / hgengamma(5,mu = mu1,sigma = sigma1, Q = Q1), # get the ratio
           rh10 = hgengamma(10,mu = mu2,sigma = sigma2, Q = Q2) / hgengamma(10,mu = mu1,sigma = sigma1, Q = Q1),
           rh15 = hgengamma(15,mu = mu2,sigma = sigma2, Q = Q2) / hgengamma(15,mu = mu1,sigma = sigma1, Q = Q1),
           rh20 = hgengamma(20,mu = mu2,sigma = sigma2, Q = Q2) / hgengamma(20,mu = mu1,sigma = sigma1, Q = Q1),
           rh25 = hgengamma(25,mu = mu2,sigma = sigma2, Q = Q2) / hgengamma(25,mu = mu1,sigma = sigma1, Q = Q1),
           rh30 = hgengamma(30,mu = mu2,sigma = sigma2, Q = Q2) / hgengamma(30,mu = mu1,sigma = sigma1, Q = Q1)) %>%
    select(rh5,rh10,rh15,rh20,rh25,rh30)
hazardratio %>% 
    filter(HR == "SBP percentile 90-95 vs 0-50") %>% 
    ggplot(aes(x = time, y = value)) +
    geom_line(color = "#377e22") +
    geom_point(color = "#377e22") +
    geom_errorbar(aes(x = 5,
                      ymin = quantile(ninetyCISBP3toSBP1at123yes[,1], probs=c(.05,.95), na.rm = TRUE)[[1]],
                      ymax = quantile(ninetyCISBP3toSBP1at123yes[,1], probs=c(.05,.95), na.rm = TRUE)[[2]]),
                  color = "#377e22") +
    geom_errorbar(aes(x = 10,
                      ymin = quantile(ninetyCISBP3toSBP1at123yes[,2], probs=c(.05,.95), na.rm = TRUE)[[1]],
                      ymax = quantile(ninetyCISBP3toSBP1at123yes[,2], probs=c(.05,.95), na.rm = TRUE)[[2]]),
                  color = "#377e22") +
    geom_errorbar(aes(x = 15,
                      ymin = quantile(ninetyCISBP3toSBP1at123yes[,3], probs=c(.05,.95), na.rm = TRUE)[[1]],
                      ymax = quantile(ninetyCISBP3toSBP1at123yes[,3], probs=c(.05,.95), na.rm = TRUE)[[2]]),
                  color = "#377e22") +
    geom_errorbar(aes(x = 20,
                      ymin = quantile(ninetyCISBP3toSBP1at123yes[,4], probs=c(.05,.95), na.rm = TRUE)[[1]],
                      ymax = quantile(ninetyCISBP3toSBP1at123yes[,4], probs=c(.05,.95), na.rm = TRUE)[[2]]),
                  color = "#377e22") +
    geom_errorbar(aes(x = 25,
                      ymin = quantile(ninetyCISBP3toSBP1at123yes[,5], probs=c(.05,.95), na.rm = TRUE)[[1]],
                      ymax = quantile(ninetyCISBP3toSBP1at123yes[,5], probs=c(.05,.95), na.rm = TRUE)[[2]]),
                  color = "#377e22") +
    geom_hline(yintercept = 1, color = "#0000f5") +
    scale_color_manual(values = c("#0389ff", "#0dd3ff", "#1c79c0"), guide = "none") +
    ggtitle("SBP percentile 90-95 vs 0-50") +
    xlab("Years after CKD diagnosis") +
    ylab("Relative hazard") +
    theme_minimal() -> p2

# SBPgroup4 vs SBPgroup1
ninetyCISBP4toSBP1at123yes = bootstrap %>%
    mutate(mu1 = beta_SBP1, sigma1 = sigma_SBP1, Q1 = kappa_SBP1,
           mu2 = beta_SBP4, sigma2 = sigma_SBP4, Q2 = kappa_SBP4) %>%
    mutate(rh5 = hgengamma(5,mu = mu2,sigma = sigma2, Q = Q2) / hgengamma(5,mu = mu1,sigma = sigma1, Q = Q1), # get the ratio
           rh10 = hgengamma(10,mu = mu2,sigma = sigma2, Q = Q2) / hgengamma(10,mu = mu1,sigma = sigma1, Q = Q1),
           rh15 = hgengamma(15,mu = mu2,sigma = sigma2, Q = Q2) / hgengamma(15,mu = mu1,sigma = sigma1, Q = Q1),
           rh20 = hgengamma(20,mu = mu2,sigma = sigma2, Q = Q2) / hgengamma(20,mu = mu1,sigma = sigma1, Q = Q1),
           rh25 = hgengamma(25,mu = mu2,sigma = sigma2, Q = Q2) / hgengamma(25,mu = mu1,sigma = sigma1, Q = Q1),
           rh30 = hgengamma(30,mu = mu2,sigma = sigma2, Q = Q2) / hgengamma(30,mu = mu1,sigma = sigma1, Q = Q1)) %>%
    select(rh5,rh10,rh15,rh20,rh25,rh30)
hazardratio %>% 
    filter(HR == "SBP percentile 95-100 vs 0-50") %>% 
    ggplot(aes(x = time, y = value)) +
    geom_line(color = "#ea3323") +
    geom_point(color = "#ea3323") +
    geom_errorbar(aes(x = 5,
                      ymin = quantile(ninetyCISBP4toSBP1at123yes[,1], probs=c(.05,.95), na.rm = TRUE)[[1]],
                      ymax = quantile(ninetyCISBP4toSBP1at123yes[,1], probs=c(.05,.95), na.rm = TRUE)[[2]]),
                  color = "#ea3323") +
    geom_errorbar(aes(x = 10,
                      ymin = quantile(ninetyCISBP4toSBP1at123yes[,2], probs=c(.05,.95), na.rm = TRUE)[[1]],
                      ymax = quantile(ninetyCISBP4toSBP1at123yes[,2], probs=c(.05,.95), na.rm = TRUE)[[2]]),
                  color = "#ea3323") +
    geom_errorbar(aes(x = 15,
                      ymin = quantile(ninetyCISBP4toSBP1at123yes[,3], probs=c(.05,.95), na.rm = TRUE)[[1]],
                      ymax = quantile(ninetyCISBP4toSBP1at123yes[,3], probs=c(.05,.95), na.rm = TRUE)[[2]]),
                  color = "#ea3323") +
    geom_errorbar(aes(x = 20,
                      ymin = quantile(ninetyCISBP4toSBP1at123yes[,4], probs=c(.05,.95), na.rm = TRUE)[[1]],
                      ymax = quantile(ninetyCISBP4toSBP1at123yes[,4], probs=c(.05,.95), na.rm = TRUE)[[2]]),
                  color = "#ea3323") +
    geom_errorbar(aes(x = 25,
                      ymin = quantile(ninetyCISBP4toSBP1at123yes[,5], probs=c(.05,.95), na.rm = TRUE)[[1]],
                      ymax = quantile(ninetyCISBP4toSBP1at123yes[,5], probs=c(.05,.95), na.rm = TRUE)[[2]]),
                  color = "#ea3323") +
    geom_hline(yintercept = 1, color = "#0000f5") +
    scale_color_manual(values = c("#0389ff", "#0dd3ff", "#1c79c0"), guide = "none") +
    ggtitle("SBP percentile 95-100 vs 0-50") +
    xlab("Years after CKD diagnosis") +
    ylab("Relative hazard") +
    theme_minimal() -> p3

tiff("./OUTPUT/Figures/Figure 2. relative hazard.tiff",
     width = 3000, height = 1000, pointsize = 10, res = 300)
p1 + p2 + p3
dev.off()

hazardratio %>% 
    filter(HR == "SBP percentile 95-100 vs 0-50") %>% 
    filter(time %in% seq(5, 30, 5))
c(1:6) %>% 
    map(~ print(quantile(ninetyCISBP4toSBP1at123yes[,.x], probs=c(.05,.95), na.rm = TRUE)))

###########################
# Figure 3. Relative time #
###########################

ttp = seq(0.05,0.70,.01)
timesSBP1 = qgengamma(ttp,3.262489,.5965208,1)
timesSBP2 = qgengamma(ttp,3.1287,.410369,1.687639)
timesSBP3 = qgengamma(ttp,2.818067,.4477858,1/.4477858)
timesSBP4 = qgengamma(ttp,2.39798,.5708262,.5708262)
timeratio <- tibble(timesSBP1 = timesSBP1, timesSBP2 = timesSBP2,
                    timesSBP3 = timesSBP3, timesSBP4 = timesSBP4) %>% 
    mutate(timesSBP2to1 = timesSBP2/timesSBP1,
           timesSBP3to1 = timesSBP3/timesSBP1,
           timesSBP4to1 = timesSBP4/timesSBP1) %>% 
    mutate(percentile = seq(0.05,0.70,.01)*100) %>% 
    select(percentile, contains("to1")) %>% 
    pivot_longer(cols = -percentile, names_to = "TR", values_to = "value") %>% 
    mutate(TR = factor(TR, levels = c("timesSBP2to1", "timesSBP3to1", "timesSBP4to1"),
                       labels = c("SBP percentile 50-90 vs 0-50",
                                  "SBP percentile 90-95 vs 0-50",
                                  "SBP percentile 95-100 vs 0-50")))

# SBPgroup2 vs SBPgroup1
ninetyCISBP2toSBP1at102030percent = bootstrap %>%
    mutate(mu1 = beta_SBP1, sigma1 = sigma_SBP1, Q1 = kappa_SBP1,
           mu2 = beta_SBP2, sigma2 = sigma_SBP2, Q2 = kappa_SBP2) %>%
    mutate(pt10 = qgengamma(0.1,mu = mu2,sigma = sigma2, Q = Q2)/ qgengamma(0.1,mu = mu1,sigma = sigma1, Q = Q1), # get the ratio
           pt20 = qgengamma(0.2,mu = mu2,sigma = sigma2, Q = Q2)/ qgengamma(0.2,mu = mu1,sigma = sigma1, Q = Q1),
           pt30 = qgengamma(0.3,mu = mu2,sigma = sigma2, Q = Q2)/ qgengamma(0.3,mu = mu1,sigma = sigma1, Q = Q1),
           pt40 = qgengamma(0.4,mu = mu2,sigma = sigma2, Q = Q2)/ qgengamma(0.4,mu = mu1,sigma = sigma1, Q = Q1),
           pt50 = qgengamma(0.5,mu = mu2,sigma = sigma2, Q = Q2)/ qgengamma(0.5,mu = mu1,sigma = sigma1, Q = Q1),
           pt60 = qgengamma(0.6,mu = mu2,sigma = sigma2, Q = Q2)/ qgengamma(0.6,mu = mu1,sigma = sigma1, Q = Q1)) %>%
    select(pt10,pt20,pt30,pt40,pt50,pt60)
timeratio %>% 
    filter(TR == "SBP percentile 50-90 vs 0-50") %>% 
    ggplot(aes(x = percentile, y = value)) +
    geom_line(color = "#ef8632") +
    geom_point(color = "#ef8632") +
    geom_errorbar(aes(x = 10,
                      ymin = quantile(ninetyCISBP2toSBP1at102030percent[,1], probs=c(.05,.95), na.rm = TRUE)[[1]],
                      ymax = quantile(ninetyCISBP2toSBP1at102030percent[,1], probs=c(.05,.95), na.rm = TRUE)[[2]]),
                  color = "#ef8632") +
    geom_errorbar(aes(x = 20,
                      ymin = quantile(ninetyCISBP2toSBP1at102030percent[,2], probs=c(.05,.95), na.rm = TRUE)[[1]],
                      ymax = quantile(ninetyCISBP2toSBP1at102030percent[,2], probs=c(.05,.95), na.rm = TRUE)[[2]]),
                  color = "#ef8632") +
    geom_errorbar(aes(x = 30,
                      ymin = quantile(ninetyCISBP2toSBP1at102030percent[,3], probs=c(.05,.95), na.rm = TRUE)[[1]],
                      ymax = quantile(ninetyCISBP2toSBP1at102030percent[,3], probs=c(.05,.95), na.rm = TRUE)[[2]]),
                  color = "#ef8632") +
    geom_errorbar(aes(x = 40,
                      ymin = quantile(ninetyCISBP2toSBP1at102030percent[,4], probs=c(.05,.95), na.rm = TRUE)[[1]],
                      ymax = quantile(ninetyCISBP2toSBP1at102030percent[,4], probs=c(.05,.95), na.rm = TRUE)[[2]]),
                  color = "#ef8632") +
    geom_errorbar(aes(x = 50,
                      ymin = quantile(ninetyCISBP2toSBP1at102030percent[,5], probs=c(.05,.95), na.rm = TRUE)[[1]],
                      ymax = quantile(ninetyCISBP2toSBP1at102030percent[,5], probs=c(.05,.95), na.rm = TRUE)[[2]]),
                  color = "#ef8632") +
    geom_errorbar(aes(x = 60,
                      ymin = quantile(ninetyCISBP2toSBP1at102030percent[,6], probs=c(.05,.95), na.rm = TRUE)[[1]],
                      ymax = quantile(ninetyCISBP2toSBP1at102030percent[,6], probs=c(.05,.95), na.rm = TRUE)[[2]]),
                  color = "#ef8632") +
    geom_hline(yintercept = 1, color = "#0000f5") +
    scale_color_manual(values = c("#0389ff", "#0dd3ff", "#1c79c0"), guide = "none") +
    ggtitle("SBP percentile 50-90 vs 0-50") +
    xlab("Percent of CKD progression") +
    ylab("Relative time") +
    theme_minimal() -> p1

# SBPgroup3 vs SBPgroup1
ninetyCISBP3toSBP1at102030percent = bootstrap %>%
    mutate(mu1 = beta_SBP1, sigma1 = sigma_SBP1, Q1 = kappa_SBP1,
           mu2 = beta_SBP3, sigma2 = sigma_SBP3, Q2 = kappa_SBP3) %>%
    mutate(pt10 = qgengamma(0.1,mu = mu2,sigma = sigma2, Q = Q2)/ qgengamma(0.1,mu = mu1,sigma = sigma1, Q = Q1), # get the ratio
           pt20 = qgengamma(0.2,mu = mu2,sigma = sigma2, Q = Q2)/ qgengamma(0.2,mu = mu1,sigma = sigma1, Q = Q1),
           pt30 = qgengamma(0.3,mu = mu2,sigma = sigma2, Q = Q2)/ qgengamma(0.3,mu = mu1,sigma = sigma1, Q = Q1),
           pt40 = qgengamma(0.4,mu = mu2,sigma = sigma2, Q = Q2)/ qgengamma(0.4,mu = mu1,sigma = sigma1, Q = Q1),
           pt50 = qgengamma(0.5,mu = mu2,sigma = sigma2, Q = Q2)/ qgengamma(0.5,mu = mu1,sigma = sigma1, Q = Q1),
           pt60 = qgengamma(0.6,mu = mu2,sigma = sigma2, Q = Q2)/ qgengamma(0.6,mu = mu1,sigma = sigma1, Q = Q1)) %>%
    select(pt10,pt20,pt30,pt40,pt50,pt60)
timeratio %>% 
    filter(TR == "SBP percentile 90-95 vs 0-50") %>% 
    ggplot(aes(x = percentile, y = value, group = TR, color = TR)) +
    geom_line(color = "#377e22") +
    geom_point(color = "#377e22") +
    geom_errorbar(aes(x = 10,
                      ymin = quantile(ninetyCISBP3toSBP1at102030percent[,1], probs=c(.05,.95), na.rm = TRUE)[[1]],
                      ymax = quantile(ninetyCISBP3toSBP1at102030percent[,1], probs=c(.05,.95), na.rm = TRUE)[[2]]),
                  color = "#377e22") +
    geom_errorbar(aes(x = 20,
                      ymin = quantile(ninetyCISBP3toSBP1at102030percent[,2], probs=c(.05,.95), na.rm = TRUE)[[1]],
                      ymax = quantile(ninetyCISBP3toSBP1at102030percent[,2], probs=c(.05,.95), na.rm = TRUE)[[2]]),
                  color = "#377e22") +
    geom_errorbar(aes(x = 30,
                      ymin = quantile(ninetyCISBP3toSBP1at102030percent[,3], probs=c(.05,.95), na.rm = TRUE)[[1]],
                      ymax = quantile(ninetyCISBP3toSBP1at102030percent[,3], probs=c(.05,.95), na.rm = TRUE)[[2]]),
                  color = "#377e22") +
    geom_errorbar(aes(x = 40,
                      ymin = quantile(ninetyCISBP3toSBP1at102030percent[,4], probs=c(.05,.95), na.rm = TRUE)[[1]],
                      ymax = quantile(ninetyCISBP3toSBP1at102030percent[,4], probs=c(.05,.95), na.rm = TRUE)[[2]]),
                  color = "#377e22") +
    geom_errorbar(aes(x = 50,
                      ymin = quantile(ninetyCISBP3toSBP1at102030percent[,5], probs=c(.05,.95), na.rm = TRUE)[[1]],
                      ymax = quantile(ninetyCISBP3toSBP1at102030percent[,5], probs=c(.05,.95), na.rm = TRUE)[[2]]),
                  color = "#377e22") +
    geom_errorbar(aes(x = 60,
                      ymin = quantile(ninetyCISBP3toSBP1at102030percent[,6], probs=c(.05,.95), na.rm = TRUE)[[1]],
                      ymax = quantile(ninetyCISBP3toSBP1at102030percent[,6], probs=c(.05,.95), na.rm = TRUE)[[2]]),
                  color = "#377e22") +
    geom_hline(yintercept = 1, color = "#0000f5") +
    scale_color_manual(values = c("#0389ff", "#0dd3ff", "#1c79c0"), guide = "none") +
    ggtitle("SBP percentile 90-95 vs 0-50") +
    xlab("Percent of CKD progression") +
    ylab("Relative time") +
    theme_minimal() -> p2

# SBPgroup4 vs SBPgroup1
ninetyCISBP4toSBP1at102030percent = bootstrap %>%
    mutate(mu1 = beta_SBP1, sigma1 = sigma_SBP1, Q1 = kappa_SBP1,
           mu2 = beta_SBP4, sigma2 = sigma_SBP4, Q2 = kappa_SBP4) %>%
    mutate(pt10 = qgengamma(0.1,mu = mu2,sigma = sigma2, Q = Q2)/ qgengamma(0.1,mu = mu1,sigma = sigma1, Q = Q1), # get the ratio
           pt20 = qgengamma(0.2,mu = mu2,sigma = sigma2, Q = Q2)/ qgengamma(0.2,mu = mu1,sigma = sigma1, Q = Q1),
           pt30 = qgengamma(0.3,mu = mu2,sigma = sigma2, Q = Q2)/ qgengamma(0.3,mu = mu1,sigma = sigma1, Q = Q1),
           pt40 = qgengamma(0.4,mu = mu2,sigma = sigma2, Q = Q2)/ qgengamma(0.4,mu = mu1,sigma = sigma1, Q = Q1),
           pt50 = qgengamma(0.5,mu = mu2,sigma = sigma2, Q = Q2)/ qgengamma(0.5,mu = mu1,sigma = sigma1, Q = Q1),
           pt60 = qgengamma(0.6,mu = mu2,sigma = sigma2, Q = Q2)/ qgengamma(0.6,mu = mu1,sigma = sigma1, Q = Q1)) %>%
    select(pt10,pt20,pt30,pt40,pt50,pt60)
timeratio %>% 
    filter(TR == "SBP percentile 95-100 vs 0-50") %>% 
    ggplot(aes(x = percentile, y = value, group = TR, color = TR)) +
    geom_line(color = "#ea3323") +
    geom_point(color = "#ea3323") +
    geom_errorbar(aes(x = 10,
                      ymin = quantile(ninetyCISBP4toSBP1at102030percent[,1], probs=c(.05,.95), na.rm = TRUE)[[1]],
                      ymax = quantile(ninetyCISBP4toSBP1at102030percent[,1], probs=c(.05,.95), na.rm = TRUE)[[2]]),
                  color = "#ea3323") +
    geom_errorbar(aes(x = 20,
                      ymin = quantile(ninetyCISBP4toSBP1at102030percent[,2], probs=c(.05,.95), na.rm = TRUE)[[1]],
                      ymax = quantile(ninetyCISBP4toSBP1at102030percent[,2], probs=c(.05,.95), na.rm = TRUE)[[2]]),
                  color = "#ea3323") +
    geom_errorbar(aes(x = 30,
                      ymin = quantile(ninetyCISBP4toSBP1at102030percent[,3], probs=c(.05,.95), na.rm = TRUE)[[1]],
                      ymax = quantile(ninetyCISBP4toSBP1at102030percent[,3], probs=c(.05,.95), na.rm = TRUE)[[2]]),
                  color = "#ea3323") +
    geom_errorbar(aes(x = 40,
                      ymin = quantile(ninetyCISBP4toSBP1at102030percent[,4], probs=c(.05,.95), na.rm = TRUE)[[1]],
                      ymax = quantile(ninetyCISBP4toSBP1at102030percent[,4], probs=c(.05,.95), na.rm = TRUE)[[2]]),
                  color = "#ea3323") +
    geom_errorbar(aes(x = 50,
                      ymin = quantile(ninetyCISBP4toSBP1at102030percent[,5], probs=c(.05,.95), na.rm = TRUE)[[1]],
                      ymax = quantile(ninetyCISBP4toSBP1at102030percent[,5], probs=c(.05,.95), na.rm = TRUE)[[2]]),
                  color = "#ea3323") +
    geom_errorbar(aes(x = 60,
                      ymin = quantile(ninetyCISBP4toSBP1at102030percent[,6], probs=c(.05,.95), na.rm = TRUE)[[1]],
                      ymax = quantile(ninetyCISBP4toSBP1at102030percent[,6], probs=c(.05,.95), na.rm = TRUE)[[2]]),
                  color = "#ea3323") +
    geom_hline(yintercept = 1, color = "#0000f5") +
    scale_color_manual(values = c("#0389ff", "#0dd3ff", "#1c79c0"), guide = "none") +
    ggtitle("SBP percentile 95-100 vs 0-50") +
    xlab("Percent of CKD progression") +
    ylab("Relative time") +
    theme_minimal() -> p3

tiff("./OUTPUT/Figures/Figure 3. relative time.tiff",
     width = 3000, height = 1000, pointsize = 10, res = 300)
p1 + p2 + p3
dev.off()

timeratio %>% 
    filter(TR == "SBP percentile 95-100 vs 0-50") %>% 
    filter(percentile %in% seq(10, 60, 10))
c(1:6) %>% 
    map(~ print(quantile(ninetyCISBP4toSBP1at102030percent[,.x], probs=c(.05,.95), na.rm = TRUE)))

##################
# Adjsuted model #
##################

#############################
# Figure 2. Relative hazard #
#############################

ttx = seq(0, 30, 1)
hazardSBP1 = hgengamma(ttx,5.114027,.5840102,.592973)
hazardSBP2 = hgengamma(ttx,5.114027-.2755161,.5840102,.592973)
hazardSBP3 = hgengamma(ttx,5.114027-.6535616,.5840102,.592973)
hazardSBP4 = hgengamma(ttx,5.114027-.7762131,.5840102,.592973)
hazardratio <- tibble(hazardSBP1 = hazardSBP1, hazardSBP2 = hazardSBP2,
                      hazardSBP3 = hazardSBP3, hazardSBP4 = hazardSBP4) %>% 
    mutate(hazardSBP2to1 = hazardSBP2/hazardSBP1,
           hazardSBP3to1 = hazardSBP3/hazardSBP1,
           hazardSBP4to1 = hazardSBP4/hazardSBP1) %>% 
    mutate(time = seq(0, 30, 1)) %>% 
    select(time, contains("to1")) %>% 
    pivot_longer(cols = -time, names_to = "HR", values_to = "value") %>% 
    mutate(HR = factor(HR, levels = c("hazardSBP2to1", "hazardSBP3to1", "hazardSBP4to1"),
                       labels = c("SBP percentile 50-90 vs 0-50",
                                  "SBP percentile 90-95 vs 0-50",
                                  "SBP percentile 95-100 vs 0-50")))

hazardratio %>% 
    ggplot(aes(x = time, y = value, group = HR, color = HR)) +
    geom_line() +
    geom_hline(yintercept = 1, color = "#0000f5") +
    scale_color_manual("Time ratio", values = c("#ef8632", "#377e22", "#ea3323")) +
    xlab("Years after CKD diagnosis") +
    ylab("Relative hazard") +
    theme_minimal() -> p1


###########################
# Figure 3. Relative time #
###########################

ttp = seq(0.05,0.70,.01)
timesSBP1 = qgengamma(ttp,5.114027,.5840102,.592973)
timesSBP2 = qgengamma(ttp,5.114027-.2755161,.5840102,.592973)
timesSBP3 = qgengamma(ttp,5.114027-.6535616,.5840102,.592973)
timesSBP4 = qgengamma(ttp,5.114027-.7762131,.5840102,.592973)
timeratio <- tibble(timesSBP1 = timesSBP1, timesSBP2 = timesSBP2,
                    timesSBP3 = timesSBP3, timesSBP4 = timesSBP4) %>% 
    mutate(timesSBP2to1 = timesSBP2/timesSBP1,
           timesSBP3to1 = timesSBP3/timesSBP1,
           timesSBP4to1 = timesSBP4/timesSBP1) %>% 
    mutate(percentile = seq(0.05,0.70,.01)*100) %>% 
    select(percentile, contains("to1")) %>% 
    pivot_longer(cols = -percentile, names_to = "TR", values_to = "value") %>% 
    mutate(TR = factor(TR, levels = c("timesSBP2to1", "timesSBP3to1", "timesSBP4to1"),
                       labels = c("SBP percentile 50-90 vs 0-50",
                                  "SBP percentile 90-95 vs 0-50",
                                  "SBP percentile 95-100 vs 0-50")))
timeratio %>% 
    ggplot(aes(x = percentile, y = value, group = TR, color = TR)) +
    geom_line() +
    geom_hline(yintercept = 1, color = "#0000f5") +
    scale_color_manual("Time ratio", values = c("#ef8632", "#377e22", "#ea3323")) +
    xlab("Percent of CKD progression") +
    ylab("Relative time") +
    theme_minimal() -> p2

tiff("./OUTPUT/Figures/Sfigure 5. HR TR_adjust.tiff",
     width = 3000, height = 1000, pointsize = 10, res = 300)
p1 + p2
dev.off()

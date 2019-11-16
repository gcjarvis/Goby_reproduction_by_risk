## Description: Reanalyzing recollection data for MS
# Author: George C Jarvis
# Date: Mon Oct 7 12:24:00 2019
# Notes: After emailing with M. Steele, he thought it was best
#         to reanalyze my data with mixed models
#         with two key differences: 1) including "trial" as a random effect
#                                   2) coding T6 uncaged the same as high-risk
#       3) for recollections, M. Steele suggested just noting it in the figure for 
#           visual counts, and noting it in the figure caption
#           -- I agree with that, so I went with count data for these revised analyses
#             -not with the proportion of fish recollected
#
#       Q: Use Poisson dist. for mixed model?
# Overall findings: A.) With all trials pooled, accounting for trial as a random
#     factor in a mixed-effect model, and using a Poisson distribution (count data),
#     there were no differences in the number of fish recollected.
#     Suggests changes in behavior (less seen in HR reefs) without changes
#     in mortality (recollections) or reproduction (egg counts similar)
# --------------

#loading packages####
rm(list=ls())

library(sciplot)
library(lme4)
library(car)
library(lmerTest)
library(dplyr)
library(ggplot2)
library(MASS)
library(nlme)
library(pwr)
library(HH)#for ancova and plots
library(vegan)
library(agricolae)#for tukey post-hoc test

#importing dataset, adding number of gobies on each reef, ordering treatments####
#includes a column ("Treatment") where uncaged and HR are coded as "High"
#also includes a column ("T6.comparison") where uncaged and high are separated
reco<-read.csv("Data/2019.10.8.recollection.data.csv", na.strings = "")
reco<-read.csv("Data/2019.10.8.recollection.data.csv")
#adding a column for survivorship (proportion of initial population recollected)
reco$Survivorship<-reco$Count/20
#subsetting T6 only for comparison
reco.T6<-reco[c(91:110),c(1:6)]

#A. running mixed models with combined HR factor and trial as both a fixed and random factor####
# 1) it does seem like a Poisson distribution is a better fit than the normal dist.
#   when using the untransformed data. No diff between Hi and Un in T6, but skewed data
#   because of low ss for that trial

#not including trial
mod.1n<-lm(Count ~ Treatment, data=reco)
hist(resid(mod.1n))
qqnorm(resid(mod.1n))
qqline(resid(mod.1n))
anova(mod.1n)
Anova(mod.1n)

#trial as a fixed factor, normal distribution
mod.1f<-lm(Count ~ Treatment*Trial, data=reco)
hist(resid(mod.1f))
qqnorm(resid(mod.1f))
qqline(resid(mod.1f))
anova(mod.1f)
Anova(mod.1f)

#fixed, poiss distribution
mod.1fp<-glm(Count ~ Treatment*Trial,family=poisson, data=reco)
hist(resid(mod.1fp))
qqnorm(resid(mod.1fp))
qqline(resid(mod.1fp))
anova(mod.1fp)
Anova(mod.1fp)
#looks better for assumptions, but the results are the same: no diff among trt's

#mixed model, normal distribution, used in MS
mod.1<-lmer(Count ~ Treatment + (1|Trial), data=reco)
hist(resid(mod.1))
qqnorm(resid(mod.1))
qqline(resid(mod.1))
anova(mod.1)
Anova(mod.1)
#looks pretty good to me, no difference in the number of fish recollected, when
# high and uncaged are combined and comapred among all trials

#same model, but now with survivorship, instead of raw counts
mod.1s<-lmer(Survivorship ~ Treatment + (1|Trial), data=reco)
hist(resid(mod.1s))
qqnorm(resid(mod.1s))
qqline(resid(mod.1s))
anova(mod.1s)
Anova(mod.1s)

#results for survivorship data

#Type III Analysis of Variance Table with Satterthwaite's method
#            Sum Sq   Mean Sq NumDF  DenDF F value Pr(>F)
#Treatment 0.011194 0.0055972     2 102.39  0.2163 0.8059

#mixed model, poisson distribution, think I should go with this one

#poisson distribution looks much better for recollections, in terms of assumptions
#think I should go with this model, but will see what M. Steele has to say
mod.1p<-glmer(Count ~ Treatment + (1|Trial),family=poisson, data=reco)
hist(resid(mod.1p))
qqnorm(resid(mod.1p))
qqline(resid(mod.1p))
anova(mod.1p)
Anova(mod.1p, type = "III") 
summary(mod.1p)
fixed.effects(mod.1p)
ranef(mod.1p)

# trial 6 only, normal distribution
mod.t6<-lm(Count ~ T6.comparison, data=reco.T6)
hist(resid(mod.t6))
qqnorm(resid(mod.t6))
qqline(resid(mod.t6))#pretty skewed, but low ss will cause that? Had 5 reps of each trt
anova(mod.t6)
#no sig. difference among treatments
summary(mod.t6)

#trial 6 only, poisson dist, no need to use mixed model
mod.t6p<-glm(Count ~ T6.comparison,family = poisson, data=reco.T6)
hist(resid(mod.t6p))
qqnorm(resid(mod.t6p))
qqline(resid(mod.t6p))
anova(mod.t6p)
Anova(mod.t6p)
summary(mod.t6p)
fixed.effects(mod.t6p)
ranef(mod.t6p)

#plots

#recollections with control and high risk pooled together
bargraph.CI(x.factor = Treatment, response = Count, 
            legend=TRUE, main="recollections, all trials pooled", 
            data = reco)

#recollections including distinction between control and high risk, just for trial 6
bargraph.CI(x.factor = T6.comparison, response = Count, 
            legend=TRUE, main="recollections, all trials pooled", 
            data = reco.T6)
#citing the lmerTest package, used to get stats on mixed models
#example cite("R, boot-package", refs, textual = TRUE)

citation("lmerTest")

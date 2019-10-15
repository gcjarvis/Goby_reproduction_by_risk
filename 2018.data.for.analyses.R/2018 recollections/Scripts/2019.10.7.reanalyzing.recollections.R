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
#     Suggests changes in behavior (more seen in high-risk viz. surveys) without changes
#     -in mortality (recollections) or reproduction (egg counts similar)
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

#subsetting T6 only for comparison
reco.T6<-reco[c(91:110),c(1:6)]

#A. running mixed models with combined HR factor and trial as both a fixed and random factor####
# 1) it does seem like a Poisson distribution is a better fit than the normal dist.
#   when using the untransformed data. No diff between Hi and Un in T6, but skewed data
#   because of low ss for that trial

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

#mixed model, normal distribution
mod.1<-lmer(Count ~ Treatment + (1|Trial), data=reco)
hist(resid(mod.1))
qqnorm(resid(mod.1))
qqline(resid(mod.1))
anova(mod.1)
Anova(mod.1)
#looks pretty good to me, no difference in the number of fish recollected, when
# high and uncaged are combined and comapred among all trials

#mixed model, poisson distribution, think I should go with this one

#poisson distribution looks much better for recollections, in terms of assumptions
#think I should go with this model, but will see what M. Steele has to say
mod.1p<-glmer(Count ~ Treatment + (1|Trial),family=poisson, data=reco)
hist(resid(mod.1p))
qqnorm(resid(mod.1p))
qqline(resid(mod.1p))
anova(mod.1p)
Anova(mod.1p) 
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

#counts
lineplot.CI(Day,den.max,group=Treatment,legend = TRUE,
            main="Visual surveys, combined HR and Uncaged", 
            xlab="Day", ylab="May number fo fish seen", 
            data=viz.surv)

#square-root transformed data for den max
lineplot.CI(Day,sqrt.dm,group=Treatment,legend = TRUE,
            main="Visual surveys, combined HR and Uncaged", 
            xlab="Day", ylab="May number fo fish seen (sqrt)", 
            data=viz.surv)

#wanting to see how it looks as a bargraph by day
bargraph.CI(x.factor = Day, response = den.max, 
            group= Treatment, legend=TRUE, main="all trials, HR combined grouped by trial", 
            data = viz.surv)

#no grouping factor, seems like H>L=M, strange...with den.max
#goes against what Mark found...
bargraph.CI(x.factor = Treatment, response = den.max, 
            legend=TRUE, main="den.max, HR combined grouped by trial", 
            data = viz.surv)

#with Density
bargraph.CI(x.factor = Treatment, response = Density, 
            legend=TRUE, main="Density, HR combined grouped by trial", 
            data = viz.surv)

#going to run the same model again but with distiction between
# HR and uncaged treatments
#results: U=H>L=M
#seems like uncaged is driving trend upward when combined
# with HR, but in general, doesn't bring the average up that much,
# likely because of low replication
bargraph.CI(x.factor = T6.comparison, response = den.max, 
            legend=TRUE, main="den.max, HR separated caged.uncaged", 
            data = viz.surv)

#now just with density
bargraph.CI(x.factor = T6.comparison, response = Density, 
            legend=TRUE, main="density, HR separated caged.uncaged", 
            data = viz.surv)
#seems like the average values are slightly lower when we use den instead of den.max
# that makes sense

#no grouping factor
bargraph.CI(x.factor = Treatment, response = Egg.count, 
            legend=TRUE, main="all trials, HR combined w uncaged, no trial", 
            data = repro)
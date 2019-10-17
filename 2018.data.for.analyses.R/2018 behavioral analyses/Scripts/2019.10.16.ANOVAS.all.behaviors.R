# Description: reanalyzing behavioral observations
# Author: George C Jarvis
# Date: Wed Oct 16 14:39:14 2019
# Notes: 1) After discussion with M. Steele, decided it might be better to analyze as 
#       separate anova's for each behavior, instead of doing a mixed-model MANCOVA.
#       I will code it the same way as I have done all the others, as a mixed model
#       with trial included as a random factor
#       2) I didn't do observations for uncaged reefs, so don't need to worry about
#       changing the labels in the .csv file
# --------------

#clear workspace
rm(list=ls())

#loading packages####
library(sciplot)
library(lme4)
library(lmerTest)
library(dplyr)
library(ggplot2)
library(plyr)
#packages needed for MANCOVA
library(car)#type III SS
#library(psych)#descriptive statistics
#library(effects)#adjusted means
#easiest mancova package
#library(jmv)

#loading df####
behave<-read.csv("Data/2019.10.8.behavior.ss.csv", na.strings = "")
head(behave)
#going to test and see if it makes a difference in models if trial is coded 
#as a factor instead of an integer, it's the same
behave$tfactor<-as.factor(behave$Trial)

#behaviors and order of script:####
# 1) proportion of time exposed
# 2) movements per minute
# 3) bites per minute
# 4) total distance moved
# 5) courtship displays per minute
#note: going to (try to) run the same model for all behaviors, with trial as a 
# random factor as mixed-model ANOVAS
# I may include density as a covariate, as it may have affected behavior
# Re: density, I'll try it both ways, without density in the model as a 
# mixed-model ANOVA, and with is in the model as as a mixed-model ANCOVA
#

# 1) proportion of time exposed, no effect, but I want to try different ANOVA's####
# without density in model
#trial different outcome if coded as a factor??? NO
mod.ex<-lmer(proportion.exposed~ Treatment  + (1|Trial), data=behave)
hist(resid(mod.ex))
qqnorm(resid(mod.ex))
qqline(resid(mod.ex))
anova(mod.ex) # sig. effect
Anova(mod.ex) #sig.effect

#a. sig effect of treatment on proportion of time exposed, not accounting for
#   the number of fish seen on the reef at the time of the observation

#rerunning analysis as a mixed model ancova (density as covariate)
mod.ex.anc<-lmer(proportion.exposed~ Treatment * density  + (1|Trial), data=behave)
hist(resid(mod.ex.anc))#looks like a slightly better model in terms of assumptions
qqnorm(resid(mod.ex.anc))
qqline(resid(mod.ex.anc))
anova(mod.ex.anc)#type III anova = no treatment effect, no interactive effect
Anova(mod.ex.anc)#Analysis of deviance, type II wald test = sig. trt. effect,
#                                                     no density, no interactive

#plotting
#grouped by trial
bargraph.CI(x.factor = Treatment, response = proportion.exposed, 
            group= Trial, legend=TRUE, main="proportion of time exposed, grouped by trial", 
            data = behave)
#no grouping factor
bargraph.CI(x.factor = Treatment, response = proportion.exposed, 
            legend=TRUE, main="proportion of time exposed, all trials combined", 
            data = behave)

# seems to be less time exposed in high-risk treatment, need to run tukey test

# 2) movements per minute, sig.effect####
mod.mo<-lmer(movements.min ~ Treatment  + (1|Trial), data=behave)
hist(resid(mod.mo))
qqnorm(resid(mod.mo))
qqline(resid(mod.mo))
anova(mod.mo) # no effect of treatment
Anova(mod.mo) # no effect of treatment

#rerunning analysis as a mixed model ancova (density as covariate)
mod.mo.anc<-lmer(movements.min~ Treatment * density  + (1|Trial), data=behave)
hist(resid(mod.mo.anc))#looks like a slightly better model in terms of assumptions
qqnorm(resid(mod.mo.anc))
qqline(resid(mod.mo.anc))
anova(mod.mo.anc)#type III anova = sig. treatment effect, no interactive effect
Anova(mod.mo.anc)#Analysis of deviance, type II wald test =  no trt. effect,
#                                                     no density, no interactive

#plotting
#grouped by trial
bargraph.CI(x.factor = Treatment, response = movements.min, 
            group= Trial, legend=TRUE, main="movement rate, grouped by trial", 
            data = behave)
#no grouping factor
bargraph.CI(x.factor = Treatment, response = movements.min, 
            legend=TRUE, main="movement rate, all trials combined", 
            data = behave)

, # looks like gobies in low risk spent more time exposed than in high and medium
# need to run tukey test

# 3) bites per minute, no effect####
mod.bi<-lmer(bites.min~ Treatment  + (1|Trial), data=behave)
hist(resid(mod.bi))
qqnorm(resid(mod.bi))
qqline(resid(mod.bi))
anova(mod.bi) # no effect
Anova(mod.bi) # no effect

#rerunning analysis as a mixed model ancova (density as covariate)
mod.bi.anc<-lmer(bites.min ~ Treatment * density  + (1|Trial), data=behave)
hist(resid(mod.bi.anc))#looks about the same as ANOVA in terms of assumptions
qqnorm(resid(mod.bi.anc))
qqline(resid(mod.bi.anc))
anova(mod.bi.anc)#type III anova = no treatment effect, no interactive effect
Anova(mod.bi.anc)#Analysis of deviance, type II wald test = no trt. effect,
#                                                     no density, no interactive

#plotting
#grouped by trial
bargraph.CI(x.factor = Treatment, response = bites.min, 
            group= Trial, legend=TRUE, main="bite rate, grouped by trial", 
            data = behave)
#no grouping factor
bargraph.CI(x.factor = Treatment, response = bites.min, 
            legend=TRUE, main="bite rate, all trials combined", 
            data = behave)

# seem to be higher bite rates in HR, but not significantly different than med and low



# 4) total distance moved, no effect####
mod.di<-lmer(total.dist.moved~ Treatment  + (1|Trial), data=behave)
hist(resid(mod.di))
qqnorm(resid(mod.di))
qqline(resid(mod.di))
anova(mod.di) # no effect
Anova(mod.di) # no effect

#rerunning analysis as a mixed model ancova (density as covariate)
mod.di.anc<-lmer(total.dist.moved ~ Treatment * density  + (1|Trial), data=behave)
hist(resid(mod.di.anc))#better fit than ANOVA in terms of assumptions
qqnorm(resid(mod.di.anc))
qqline(resid(mod.di.anc))
anova(mod.di.anc)#type III anova = no treatment effect, no interactive effect
Anova(mod.di.anc)#Analysis of deviance, type II wald test = no trt. effect,
#                                                     no density, no interactive

#plotting
#grouped by trial
bargraph.CI(x.factor = Treatment, response = total.dist.moved, 
            group= Trial, legend=TRUE, main="linear distance moved (cm), grouped by trial", 
            data = behave)
#no grouping factor
bargraph.CI(x.factor = Treatment, response = total.dist.moved, 
            legend=TRUE, main="linear distance moved (cm), all trials combined", 
            data = behave)

# seemed to move the same distance, regardless of treatment

# 5) courtship displays per minute, no effect####
mod.co<-lmer(courtship.min ~ Treatment  + (1|Trial), data=behave)
hist(resid(mod.co))
qqnorm(resid(mod.co))
qqline(resid(mod.co))
anova(mod.co) # no effect
Anova(mod.di) # no effect

#rerunning analysis as a mixed model ancova (density as covariate)
mod.co.anc<-lmer(courtship.min ~ Treatment * density  + (1|Trial), data=behave)
hist(resid(mod.co.anc))#better fit than ANOVA in terms of assumptions
qqnorm(resid(mod.co.anc))
qqline(resid(mod.co.anc))
anova(mod.co.anc)#type III anova = no treatment effect, no interactive effect
Anova(mod.co.anc)#Analysis of deviance, type II wald test = no trt. effect,
#                                                     no density, no interactive

#plotting
#grouped by trial
bargraph.CI(x.factor = Treatment, response = courtship.min, 
            group= Trial, legend=TRUE, main="courtship rate, grouped by trial", 
            data = behave)
#no grouping factor
bargraph.CI(x.factor = Treatment, response = courtship.min, 
            legend=TRUE, main="courship rate, all trials combined", 
            data = behave)

# no diff, but seems to trend as M>L>H



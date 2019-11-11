# Description: reanalyzing behavioral observations
# Author: George C Jarvis
# Date: Wed Oct 16 14:39:14 2019
# Notes: 1) After discussion with M. Steele, decided it might be better to analyze as 
#       separate anova's for each behavior, instead of doing a mixed-model MANCOVA.
#       I will code it the same way as I have done all the others, as a mixed model
#       with trial included as a random factor
#       2) I didn't do observations for uncaged reefs, so don't need to worry about
#       changing the labels in the .csv file
#       3) 2019.10.25 - decided to use the number of fish recollected from each reef
#         as the covariate, because density is confounded by treatent
#         - In the models where recollections weren't significant, decided to drop from the model,
#             but I should probably check with Mark on that one. It seems like a more 
#             powerful analysis if I remove that factor from the model
#         -- I also did posthoc tests with multiple comparisons for the mixed models
#         - found that exposure time was sig. lower in high-risk treatment, which is
#             consistent with the results from my density counts
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
library(agricolae)#can't use standard tukey test for mixed models
library(multcomp)# post-hoc for fixed effects in mixed model

#loading df####
behave<-read.csv("Data/2019.10.8.behavior.ss.csv")
behave<-na.omit(behave) #no NA's to remove
#with recollection data
behave<-read.csv("Data/2019.10.25.behavior.includes.recollections.csv")
behave<-na.omit(behave)
head(behave)
#going to test and see if it makes a difference in models if trial is coded 
#as a factor instead of an integer, it's the same
#behave$tfactor<-as.factor(behave$Trial)
#making the variable "avg.inhab" ((20+reco)/2), rounded to the nearest whole fish
behave$avg.inhab<-(ceiling((behave$Recollection+20)/2))
#going to include the average number of inhabitants as the covariate,
# -not the number of fish recollected

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

# 1) proportion of time exposed, sig. effect, ran post-hoc test for mult comp: L(a),M(ab),H(b)####
# without density in model
#trial different outcome if coded as a factor??? NO
mod.ex<-lmer(proportion.exposed~ Treatment + (1|Trial), data=behave)
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

#rerunning analyses, with number of fish recollected from each reef as covariate
#if there is no correlation between number fo fish recollected and frequency of behavior
# I would say remove it from model
mod.ex.anc.reco<-lmer(proportion.exposed~ Treatment * Recollection  + (1|Trial), data=behave)
hist(resid(mod.ex.anc.reco))#looks like a slightly better model in terms of assumptions
qqnorm(resid(mod.ex.anc.reco))
qqline(resid(mod.ex.anc.reco))
anova(mod.ex.anc.reco)#type III anova = no treatment effect, no interactive effect
Anova(mod.ex.anc.reco)#Analysis of deviance, type II wald test = sig. trt. effect,
#                                                     no density, no interactive

#plotting relationship between recollections and exposure time
plot(behave$Recollection, behave$proportion.exposed)
abline(lm(behave$proportion.exposed~behave$Recollection))
#hardly any relationship there, consider dropping from model

mod.ex.anc.no.reco<-lmer(proportion.exposed~ Treatment + (1|Trial), data=behave)
hist(resid(mod.ex.anc.no.reco))#looks like a slightly better model in terms of assumptions
qqnorm(resid(mod.ex.anc.no.reco))
qqline(resid(mod.ex.anc.no.reco))
anova(mod.ex.anc.no.reco)#type III anova = treatment effect
Anova(mod.ex.anc.no.reco)#Analysis of deviance, type II wald test = sig. trt. effect

#mixed model with avg.inhab as covariate
mod.avg.inhab.mm<-lmer(proportion.exposed~ Treatment * avg.inhab  + (1|Trial), data=behave)
hist(resid(mod.avg.inhab.mm))#looks like a slightly better model in terms of assumptions
qqnorm(resid(mod.avg.inhab.mm))
qqline(resid(mod.avg.inhab.mm))
anova(mod.avg.inhab.mm)#type III anova = no treatment effect, no interactive effect
Anova(mod.avg.inhab.mm)#Analysis of deviance, type II wald test = sig. trt. effect,
#                                                     no density, no interactive

#plotting relationship between avg.inhab and exposure time
plot(behave$avg.inhab, behave$proportion.exposed)
abline(lm(behave$proportion.exposed~behave$avg.inhab))
#hardly any relationship there, dropped from the model, found sig.
# -treatment effect

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
# with multcomp package
#summary(glht(YOUR MODEL, linfct=mcp(YOUR FIXED FACTOR="Tukey")))
summary(glht(mod.ex, linfct=mcp(Treatment="Tukey")))
#seems like the sig. difference is between low and high, but not med and high, 
# or low and med, so L(a)>M(ab)>H(b), if you were going to graph it
#I'm not sure if you can do post-hoc with mixed models, but this multcomp
# package (multiple comparisons) may be a workaround

# 2) movements per minute, no sig.effect, no post-hoc####
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

#rerunning as mixed model with recollections as covariate
#plotting relationship between recollections and movements
plot(behave$Recollection, behave$movements.min)
abline(lm(behave$movements.min~behave$Recollection))
#hardly any relationship there, consider dropping from model

mod.mo.anc.reco<-lmer(movements.min~ Treatment *Recollection + (1|Trial), data=behave)
hist(resid(mod.mo.anc.reco))#looks like a slightly better model in terms of assumptions
qqnorm(resid(mod.mo.anc.reco))
qqline(resid(mod.mo.anc.reco))
anova(mod.mo.anc.reco)#type III anova = no treatment effect
Anova(mod.mo.anc.reco)#Analysis of deviance, type II wald test = sig. trt. effect

#see first model for stats when recollection dropped from model

mod.mo.avg.inhab<-lmer(movements.min~ Treatment *avg.inhab + (1|Trial), data=behave)
hist(resid(mod.mo.avg.inhab))#looks like a slightly better model in terms of assumptions
qqnorm(resid(mod.mo.avg.inhab))
qqline(resid(mod.mo.avg.inhab))
anova(mod.mo.avg.inhab)#type III anova = no treatment effect
Anova(mod.mo.avg.inhab)#Analysis of deviance, type II wald test = sig. trt. effect
 
#no difference, even when using avg.inhab, so dropping from model

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

# 3) bites per minute, no effect, no post-hoc####
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

#rerunning analysis as a mixed model ancova (recollections as covariate)
mod.bi.anc.reco<-lmer(bites.min ~ Treatment * Recollection  + (1|Trial), data=behave)
hist(resid(mod.bi.anc.reco))#looks about the same as ANOVA in terms of assumptions
qqnorm(resid(mod.bi.anc.reco))
qqline(resid(mod.bi.anc.reco))
anova(mod.bi.anc.reco)#type III anova = no treatment effect, no interactive effect
Anova(mod.bi.anc.reco)#Analysis of deviance, type II wald test = no trt. effect,
#                                                     no density, no interactive

#plotting relationship between recollections and bite rate
plot(behave$Recollection, behave$bites.min)
abline(lm(behave$bites.min~behave$Recollection))

#using avg.inhab as covariate
mod.bi.avg.inhab<-lmer(bites.min ~ Treatment * avg.inhab  + (1|Trial), data=behave)
hist(resid(mod.bi.avg.inhab))#looks about the same as ANOVA in terms of assumptions
qqnorm(resid(mod.bi.avg.inhab))
qqline(resid(mod.bi.avg.inhab))
anova(mod.bi.avg.inhab)#type III anova = no treatment effect, no interactive effect
Anova(mod.bi.avg.inhab)#Analysis of deviance, type II wald test = no trt. effect,
#                                                     no density, no interactive

#plotting relationship between recollections and bite rate
plot(behave$avg.inhab, behave$bites.min)
abline(lm(behave$bites.min~behave$avg.inhab))

#no relationship, see first model for stats when term is dropped

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



# 4) total distance moved, no effect of treatment, but positive effect of recollection, keep recollection in the model? I included what I think is an explanation for the behavior####
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


#rerunning analysis as a mixed model ancova (recollection as covariate)
mod.di.anc.reco<-lmer(total.dist.moved ~ Treatment * Recollection + (1|Trial), data=behave)
hist(resid(mod.di.anc.reco))#better fit than ANOVA in terms of assumptions
qqnorm(resid(mod.di.anc.reco))
qqline(resid(mod.di.anc.reco))
anova(mod.di.anc.reco)#type III anova = no treatment effect, higher recollection=more movements, no interactive effect
Anova(mod.di.anc.reco)#Analysis of deviance, type II wald test = no trt. effect,
#                                                     recollection effect, no interactive

#plotting relationship between recollections and distance moved
plot(behave$Recollection, behave$total.dist.moved)
abline(lm(behave$total.dist.moved~behave$Recollection))
#positive relationship

mod.di.rec<-lm(total.dist.moved~Recollection, data=behave)
hist(resid(mod.di.rec))
qqnorm(resid(mod.di.rec))
qqline(resid(mod.di.rec))
anova(mod.di.rec)
Anova(mod.di.rec)
summary(mod.di.rec)
#shows clear relationship, but only has an R-squared of 0.078, so not a great fit
#rationale: 1) more fish on reefs may have been associated with more interactions?
# - that will show in the courtship analyses
# 2) may have been a proxy for safety? more fish on reef overall = more boldness?

#maybe slight positive relationship?

#will keep recollection in the model, I suppose, because here it was meaningful
#I'm not sure what it means, but I'll keep it in the model

#using avg.inhab as covariate, no interaction, so will drop interaction and rerun model
mod.di.avg.inhab<-lmer(total.dist.moved ~ Treatment + avg.inhab + (1|Trial), data=behave)
hist(resid(mod.di.avg.inhab))#better fit than ANOVA in terms of assumptions
qqnorm(resid(mod.di.avg.inhab))
qqline(resid(mod.di.avg.inhab))
anova(mod.di.avg.inhab)#type III anova = no treatment effect, higher avg.inhab=more movements, no interactive effect
Anova(mod.di.avg.inhab)#Analysis of deviance, type II wald test = no trt. effect,
#                                                     avg.inhab effect, no interactive

#plotting relationship between recollections and distance moved
plot(behave$avg.inhab, behave$total.dist.moved)
abline(lm(behave$total.dist.moved~behave$avg.inhab))
#positive relationship

#wonder if I should plot this like an ANCOVA? Not sure

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
Anova(mod.co) # no effect

#rerunning analysis as a mixed model ancova (density as covariate)
mod.co.anc<-lmer(courtship.min ~ Treatment * density  + (1|Trial), data=behave)
hist(resid(mod.co.anc))#better fit than ANOVA in terms of assumptions
qqnorm(resid(mod.co.anc))
qqline(resid(mod.co.anc))
anova(mod.co.anc)#type III anova = no treatment effect, no interactive effect
Anova(mod.co.anc)#Analysis of deviance, type II wald test = no trt. effect,
#                                                     no density, no interactive

#rerunning with recollections as covariate
mod.co.anc.reco<-lmer(courtship.min ~ Treatment * Recollection  + (1|Trial), data=behave)
hist(resid(mod.co.anc.reco))#better fit than ANOVA in terms of assumptions
qqnorm(resid(mod.co.anc.reco))
qqline(resid(mod.co.anc.reco))
anova(mod.co.anc.reco)#type III anova = no treatment effect, no interactive effect
Anova(mod.co.anc.reco)#Analysis of deviance, type II wald test = no trt. effect,
#                                                     no density, no interactive

#no effect of recollection

#rerunning with avg.inhab as covariate
mod.co.anv.inhab<-lmer(courtship.min ~ Treatment * avg.inhab  + (1|Trial), data=behave)
hist(resid(mod.co.anv.inhab))#better fit than ANOVA in terms of assumptions
qqnorm(resid(mod.co.anv.inhab))
qqline(resid(mod.co.anv.inhab))
anova(mod.co.anv.inhab)#type III anova = no treatment effect, no interactive effect
Anova(mod.co.anv.inhab)#Analysis of deviance, type II wald test = no trt. effect,
#                                                     no density, no interactive

#removed from model and reran

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




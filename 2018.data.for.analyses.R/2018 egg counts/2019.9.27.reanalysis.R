# Description: Reanalyzing data for MS
# Author: George C Jarvis
# Date: Mon Sep 30 08:48:33 2019
# Notes: After emailing with M. Steele, he thought it was best
#         to reanalyze my reproduction and behavior data with mixed models
#         with two key differences: 1) including "trial" as a random effect
#                                   2) coding T6 uncaged the same as high-risk
#       I will still include average # of fish on reef as a covariate 
#       (using average of initial and final number of gobies on the reefs as
#       covariate in ANCOVAS) for T6 I will keep that treatment as uncaged so 
#       I can compare it to the high-risk treatment to show that there were
#       no statistical differences. But that will only be that way for 
#       the single, short, low-rep trial.
#
# Overall findings: A.) With all trials pooled, and accounting for trial as a random
#     factor in a mixed-effect model, there was no diff in reproduction among
#     treatments. Less fish = less reproduction, but reproduction in each treatment
#     did not depend on the number of fish on the reefs (i.e. no interactive effects)
#     B.) no difference in output between high-risk and ungacaged treatments,
#     which justifies pooling them as "High-risk" in the overall analyses
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

#importing dataset, adding number of gobies on each reef, ordering treatments####
#includes a cloumn ("Treatment") where uncaged and HR are coded as "High"
#also includes a column ("T6.comparison") where uncaged and high are separated
repro<-read.csv("Data/new.data.2019.9.30.csv", na.strings = "")

#adding column for average density, rounded to nearest whole number of fish
repro$avg.inhab<-(ceiling((repro$Recollection+20)/2))

#ordering "Treatment" and "T.6.comparison"
repro$Treatment<-ordered(repro$Treatment, levels=c("Low","Medium","High"))
repro$T6.comparison<-ordered(repro$T6.comparison, levels=c("Low","Medium","High","Uncaged"))

#A. running models with combined HR factor and trial as both a fixed and random factor####
# 1) trial as a fixed factor (don't think this is justified)
mod.1<-lm(Egg.count ~ Treatment * avg.inhab *Trial, data=repro)
hist(resid(mod.1))
qqnorm(resid(mod.1))
qqline(resid(mod.1))
anova(mod.1)
#a. no sig. effect of treatment, with HR and uncaged combined
#b. L>M>H, but it seems like there is too much variation see differences
#c. sig. differences in reproduction based on the number of gobies, the trial, 
# and interaction between the two
#d. no three-way interaction

#plots

#grouped by trial
bargraph.CI(x.factor = Treatment, response = Egg.count, 
            group= Trial, legend=TRUE, main="all trials, HR combined grouped by trial", 
            data = repro)
#no grouping factor
bargraph.CI(x.factor = Treatment, response = Egg.count, 
            legend=TRUE, main="all trials, HR combined w uncaged, no trial", 
            data = repro)

# 2) trial as a random factor
mod.2<-lmer(Egg.count~ Treatment * avg.inhab + (1|Trial), data=repro)
hist(resid(mod.2))
qqnorm(resid(mod.2))
qqline(resid(mod.2))
anova(mod.2)
#a. no effect of treatment on output
#b. there was an effect of number of gobies on output
#bi. general trend seems to be, more gobies, higher reproduction, 
#   regardless of trt, and this occurs independent of trial effects;
#   meaning that this is the trend even in Part A, w/trial as a fixed effect
#c. no interactive effects of treatment and number of gobies

lineplot.CI(avg.inhab,Egg.count,group=Treatment,legend = TRUE,
            main="Reproduction vs. number of gobies", 
            xlab="Number of gobies", ylab="Reproduction (mean +/- se)", 
            data=repro)

#B. running models with HR and uncaged separately, with only trial 6 data####
#Rationale: want to see if we can justify pooling high-risk and uncaged as a single treatment
#   From goby's perspective, they are the same because predator access is the same
#Note: no longer need to include trial as a factor, only one trial with uncaged
#     because of that, only need to run one model

#subsetting raw data for only Trial 6 (rows 91-110, in raw data)
repro.t6 <- repro[c(91:110),c(1:10)]

# 1) trial as a fixed factor
mod.3<-lm(Egg.count ~ T6.comparison * avg.inhab, data=repro.t6)
hist(resid(mod.3))
qqnorm(resid(mod.3))
qqline(resid(mod.3))
anova(mod.3)
#a. no sig. effect of treatment, with HR and uncaged analyzed separately
#b. U>H>M>L, but it seems like there is too much variation see differences
#c. no sig. differences of any factors on output, including interaction 

#plots

#reproduction by treatment
bargraph.CI(x.factor = T6.comparison, response = Egg.count, 
            legend=TRUE, main="Trial 6 only, HR separated", 
            data = repro.t6)

#no effect of the number fo gobies on reproductive output 
lineplot.CI(avg.inhab,Egg.count,group=T6.comparison,legend = TRUE,
            main="Trial 6 only: Reproduction vs. number of gobies", 
            xlab="Number of gobies", ylab="Reproduction (mean +/- se)", 
            data=repro.t6)

# bii. general trend seems to be, more gobies, higher reproduction, 
#   regardless of trt, and this occurs independent of trial effects;
#   BUT, only surveyed three times, so only have data for high densities 
#   (11-13 gobies/plot), so that could be why we don't see a big difference

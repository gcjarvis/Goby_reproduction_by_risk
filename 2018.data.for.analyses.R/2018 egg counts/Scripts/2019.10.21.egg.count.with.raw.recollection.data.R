# Description: reanalyzing, using raw number of gobies recollected as covariate
# Author: George C Jarvis
# Date: Mon Oct 21 16:34:41 2019
# Notes: I think this is more representative than using an average (20+ # reco/2),
#       so I want to go with this over a calculated metric
# --------------

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
#repro$avg.inhab<-(ceiling((repro$Recollection+20)/2))

#ordering "Treatment" and "T.6.comparison"
repro$Treatment<-ordered(repro$Treatment, levels=c("Low","Medium","High"))
repro$T6.comparison<-ordered(repro$T6.comparison, levels=c("Low","Medium","High","Uncaged"))

#model
mod.2<-lmer(Egg.count~ Treatment * Recollection + (1|Trial), data=repro)
hist(resid(mod.2))
qqnorm(resid(mod.2))
qqline(resid(mod.2))
anova(mod.2)
Anova(mod.2) 

#plotting relationship between recolelctions and egg count
plot(repro$Recollection, repro$Egg.count)
abline(lm(repro$Egg.count~repro$Recollection))
#positive relationship, but much more data for lower numbers of fish recollected

mod.ER<-lm(Egg.count~Recollection, data=repro)
hist(resid(mod.ER))
qqnorm(resid(mod.ER))
qqline(resid(mod.ER))
anova(mod.ER)
Anova(mod.ER)
summary(mod.ER)
#shows clear relationship, but only has an R-squared of 0.385
#this makes sense becasue of the large spread of repro for each recollection value
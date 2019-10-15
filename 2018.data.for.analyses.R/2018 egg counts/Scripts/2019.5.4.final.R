# Description: final analyses for egg counts
# Author: George C Jarvis
# Date: Sat May 04 12:23:57 2019
# --------------

##NOTE: using average of initial and final number of gobies on the reefs as
##      covariate in ANCOVAS

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

#broken down by trial for ANCOVA with jmv package (hopefully)
#adding extra column for average number of fish on the reef throughout the trials
#might eventually consider doing ceiling(average), to get numbers of while fish
#but not entirely necessary

#doing with non-integers:

egg.t1.3<-read.csv("Data/new.data.2019.4.23.t1.3.csv", na.strings = "") #uses adjusted counts for density
egg.t1.3$avg.inhab<-(ceiling((egg.t1.3$Recollection+20)/2))
#egg.t1.3$Trial<-as.factor(egg.t1.3$Trial)
#for t4.5
egg.t4.5<-read.csv("Data/new.data.2019.4.23.no.t6.csv", na.strings = "") #uses adjusted counts for density
egg.t4.5$avg.inhab<-(ceiling((egg.t4.5$Recollection+20)/2))
#for t6
egg.t6<-read.csv("Data/new.data.2019.4.23.csv", na.strings = "") #uses adjusted counts for density
egg.t6$avg.inhab<-(ceiling((egg.t6$Recollection+20)/2))
egg.t6<-egg.t6[egg.t6$Trial==6,]#subsetting only t6
mapvalues(egg.t6$Treatment, from = "Control", to = "Uncaged")


#models and plots (maybe)
#can run linear models, don't need to call them ancovas
#2017
ancova(Egg.count ~ Trial * avg.inhab, data=egg.t1.3) #interaction between inhab and trial
ancova(Egg.count ~ Treatment * Trial, data=egg.t1.3) #interaction between inhab and trial

#trial as a fixed factor
mod.2017<-lm(Egg.count ~ Treatment * avg.inhab *Trial, data=egg.t1.3)
hist(resid(mod.2017))
qqnorm(resid(mod.2017))
qqline(resid(mod.2017))
anova(mod.2017)

bargraph.CI(x.factor = avg.inhab, response = Egg.count, 
            group= Trial, legend=TRUE, main="trials 1-3", 
            data = egg.t1.3)

#trial as a random factor
mod.2017a<-lmer(Egg.count~ Treatment * avg.inhab + (1|Trial), data=egg.t1.3)
hist(resid(mod.2017a))
qqnorm(resid(mod.2017a))
qqline(resid(mod.2017a))
anova(mod.2017a)

#for exp. 1, there seemed to be more of a positive effect of the number 
# of gobies inhabiting the reef on egg counts in trial 1
# those effects were not present in trials 2 or 3
#could be due to seasonal effect (earliest trial might have been different)
#but I'm not sure

#include as a fixed or random factor based on that?\
# my inclination, given that it doesn't affect treatment effects, is
#to include it as a random factor
#but not sure whether to run analysis of deviance or analysis of variance?

#2018.t4.5
#trial as a fixed factor
mod.2018.t.4.5<-lm(Egg.count ~ Treatment * avg.inhab *Trial, data=egg.t4.5)
hist(resid(mod.2018.t.4.5))
qqnorm(resid(mod.2018.t.4.5))
qqline(resid(mod.2018.t.4.5))
anova(mod.2018.t.4.5)

#trial as a random factor
mod.2018.t.4.5a<-lmer(Egg.count~ Treatment * avg.inhab + (1|Trial), data=egg.t4.5)
hist(resid(mod.2018.t.4.5a))
qqnorm(resid(mod.2018.t.4.5a))
qqline(resid(mod.2018.t.4.5a))
anova(mod.2018.t.4.5a)
Anova(mod.2018.t.4.5a)

ancova(Egg.count ~ Treatment * avg.inhab, data=egg.t4.5) #interaction between inhab and trial
ancova(Egg.count ~ Treatment * Trial, data=egg.t4.5) #interaction between inhab and trial

#experiment 2 shows no effects of trial, so would include as random factor

#2018.t6 --> only ran one trial, so no need to include as factor
#trial as a fixed factor
mod.2018.t.6<-lm(Egg.count ~ Treatment * avg.inhab, data=egg.t6)
hist(resid(mod.2018.t.6))
qqnorm(resid(mod.2018.t.6))
qqline(resid(mod.2018.t.6))
anova(mod.2018.t.6)

ancova(Egg.count ~ Treatment * avg.inhab, data=egg.t6) #interaction between inhab and trial


# still see a positive trend in egg counts with increased numbers of fish
#on the reefs, but not a significant trend here

##using experiment as a factor to see if I can see a treatment effect

#doesn't seem like there is a trt effect, p = 0.15, even with best model

egg.exp<-read.csv("Data/new.data.2019.5.9.exp1.2.factor.csv",na.strings = "")
egg.exp$avg.inhab<-((egg.exp$Recollection+20)/2)

mod1<-lm(Egg.count ~ Treatment * avg.inhab, data=egg.exp)
hist(resid(mod1))
qqnorm(resid(mod1))
qqline(resid(mod1))
anova(mod1)

bargraph.CI(x.factor = Treatment, response = Egg.count, 
            legend=FALSE, main="all exp", xlab="Treatment", 
            ylab="Eggs per reef (mean +/- se)", data = egg.exp)

mod2<-lm(Egg.count ~ Treatment * avg.inhab*Trial*Experiment, data=egg.exp)
hist(resid(mod2))
qqnorm(resid(mod2))
qqline(resid(mod2))
anova(mod2)

mod2<-lm(Egg.count ~ Treatment * avg.inhab*Experiment, data=egg.exp)
hist(resid(mod2))
qqnorm(resid(mod2))
qqline(resid(mod2))
anova(mod2)

mod2<-lmer(Egg.count ~ Treatment * avg.inhab*Experiment + (1|Trial), data=egg.exp)
hist(resid(mod2))
qqnorm(resid(mod2))
qqline(resid(mod2))
anova(mod2)         

#####rerun without trial for written thesis
#2017
mod1<-lm(Egg.count~Treatment*avg.inhab,data=egg.t1.3)
hist(resid(mod1))
anova(mod1)

tapply(egg.t1.3$Egg.count,egg.t1.3$Treatment,mean)
#Low   Medium     High 
#3642.833 3023.722 3070.778 

#2018t4.5
mod2<-lm(Egg.count~Treatment*avg.inhab,data=egg.t4.5)
hist(resid(mod2))
anova(mod2)

tapply(egg.t4.5$Egg.count,egg.t4.5$Treatment,mean)
#Low   Medium     High 
#17138.58 14902.17 14002.67 

#2018.t6
mod3<-lm(Egg.count~Treatment*avg.inhab,data=egg.t6)
hist(resid(mod3))
anova(mod3)

tapply(egg.t6$Egg.count,egg.t6$Treatment,mean)
#Low  Medium    High Control 
#2293.6  2306.4  2698.0  3530.6 
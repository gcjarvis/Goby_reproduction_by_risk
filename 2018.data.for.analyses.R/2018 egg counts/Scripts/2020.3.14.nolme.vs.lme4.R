# Description: rerunning models for egg counts based off of mariana and giulia's comments 
# Author: George C Jarvis
# Date: Tue Mar 10 20:25:02 2020
# Notes: Trying to figure out how to go about analyzing my data...still
#         Mainly, 1) do I include trial as a fixed vs. random effect in my model
#                 2) also wondering if I include year at all? Trial is coded in a way that distinguishes it among years
#
#
# --------------

rm(list=ls())

#packages

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
library(multcomp)

#maybe don't run this until models have been constructed
options(contrasts = c("contr.sum","contr.poly")) #this is important, run before ANOVA, will set SS to type III

repro<-read.csv("Data/new.data.2019.9.30.csv")
repro<-na.omit(repro) # no NA's to omit

#data manipulation####
#adding column for average density, rounded to nearest whole number of fish
repro$avg.inhab<-(ceiling((repro$Recollection+20)/2))

#adding a column for year (as a proxy for tagging procedure), where trials 1-3 = 2017, and 4-6 = 2018
repro$Year <- ifelse(repro$Trial <=3, 2017, 2018)
#want it as a factor? Going to make a variable with it as a factor, run the model, and see if I get different results
repro$Year.fact<- as.factor(repro$Year)

#adding egg/week variable to the dataset
repro<-repro %>%
  mutate(egg.week = ifelse(Trial<4, Egg.count/1,
                           ifelse(Trial == 4| Trial == 5, (ceiling(Egg.count/4)),
                                  ifelse(Trial == 6, (ceiling(Egg.count/2)), NA))))

#making trial a factor
repro$Trial.fact<-as.factor(repro$Trial)

#models####

#using lme4 package, not nlme
mod1<-lm(egg.week~Treatment*Year.fact*Trial.fact+avg.inhab,repro)
summary(mod1)
anova(mod1)

mod2<-lm(egg.week~Treatment*Year.fact*Trial.fact*avg.inhab,repro)
summary(mod2)
anova(mod2)

#there are some three-way interactions happening here, so not really sure
# - what the best course of action is
# I do think that I need to make sure that Trial is coded as a factor, 
# - not as an integer, because I think it makes a difference

#now looking at mixed model, and then going to compare to linear model

mod.1me<-lmer(egg.week~Treatment*Year.fact*avg.inhab+(1|Year.fact:Trial.fact),repro, method=ML)
summary(mod.1me)
anova(mod.1me)

#comparing models
anova(mod2,mod.1me)


mod2.2.luk<-lme(egg.week~(Treatment*Year.fact)+ Treatment+
                  avg.inhab+Year.fact,random=~1|Trial,repro)
summary(mod2.2.luk)
anova(mod2.2.luk, type='marginal')






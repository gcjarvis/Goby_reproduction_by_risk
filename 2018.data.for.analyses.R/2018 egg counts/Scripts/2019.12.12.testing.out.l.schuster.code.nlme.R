# Description: testing lukas's code for nlme nested design vs. my previous models
# Author: George C Jarvis
# Date: Thu Dec 12 23:12:13 2019
# Notes: Lukas seems to think that this will be better for my nested ANCOVA models
#       for egg counts. We'll see.  I want to compare results for both models
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


#modeling####
#original model with lme4, with weekly output per reef as response variable, shown is the reduced model 
# - (from "2019.7.12.adding.year.factor" script)
mod2.2<-lmer(egg.week~Treatment*Year.fact+avg.inhab+(1|Year.fact:Trial),
             data=repro)
hist(resid(mod2.2))
qqnorm(resid(mod2.2))
qqline(resid(mod2.2))
anova(mod2.2)
Anova(mod2.2)
summary(mod2.2) # same results, will check out a model comparison?

#Lukas's model with nlme, let's do the full model and then reduce if needed
mod2.luk<-lme(egg.week~Treatment*Year.fact*avg.inhab,random=~1|Trial,repro,method="REML")
summary(mod2.luk)
anova(mod2.luk, type='marginal')

#running reduced model
mod2.1.luk<-lme(egg.week~(Treatment*Year.fact)+(avg.inhab*Year.fact)+(Treatment*avg.inhab)+ Treatment+
               avg.inhab+Year.fact,random=~1|Trial,repro,method="REML")
summary(mod2.1.luk)
anova(mod2.1.luk, type='marginal')

#running further reduced model
mod2.2.luk<-lme(egg.week~(Treatment*Year.fact)+ Treatment+
                  avg.inhab+Year.fact,random=~1|Trial,repro,method="REML")
summary(mod2.2.luk)
anova(mod2.2.luk, type='marginal')

anova(mod2.luk,mod2.1.luk,mod2.2.luk) #aic value is highest for the least-reduced model, but 
# also might not be meaningful because I'm comparing models with different fixed effects


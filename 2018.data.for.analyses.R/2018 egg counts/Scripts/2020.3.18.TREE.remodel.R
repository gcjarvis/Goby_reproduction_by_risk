# Description: Reanalyzing mixed models with data from Arnqvist 2020 TREE supplementary material
# Author: George C Jarvis
# Date: Wed Mar 18 14:00:53 2020
# Notes: I think my analyses are pseudoreplicated, because the way that they are currently
# - being analyzed assumes that reefs are treated as independent replicates, when they really are not.
# - Instead, I need to take into consideration that reefs, and treatments, nested within trial, which
# - are nested within year, are the correct level of replication
#
# As such, I'm going to model my data after the more suitable models that are outlined in this TREE
# - paper, specifically, models D.1 (‘random slopes’ model, will have to do in SYSTAT), 
# - E (RM-ANOVA form), and F (full ‘random slopes’ model) from Box 2 in this paper
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
library(multcomp)

#options(contrasts = c("contr.sum","contr.poly")) #this is important, run before ANOVA, will set SS to type III

#working with dummy data from paper first####

#loading data
mmc1<-read.csv("Data/mmc1.csv")
mmc1a<-na.omit(mmc1)# for E design

#modeling dummy data####

#E (RM-ANOVA)

mod1 <- aov(RESPONSE ~ TREAT * factor(AC) + Error(POP), data=mmc1)
summary(mod1)

#reformatted data and took cell means (means for each response for each )
mod1a <- aov(cell.means ~ TREAT * factor(AC) + Error(POP), data=mmc1a)
summary(mod1a)

#F (full ‘random slopes’ model)

library (lme4)
mod<-lmer(RESPONSE~TREAT*AC+(1+AC|POP),data=mmc1)
Anova(mod, type=c("II","III", 2, 3), test.statistic=c("F"))

#importing goby dataset, adding number of gobies on each reef, ordering treatments####
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

# square-root transforming reprodcution per week, not sure how to analyze this quite yet,
# - but will figure that out later

repro$sqrt.egg<-sqrt(repro$egg.week)

#subsetting data to do t6 comparison between HR caged and uncaged treatments
#subset trial
repro.t6<-repro[repro$Trial==6,]
#subset treatments, caged and uncaged only (T6.comparison)

#subsetting by multiple character factors within a variable. Awesome code!
repro.t6.comp<-repro.t6[repro.t6$T6.comparison == "High" | repro.t6$T6.comparison == "Uncaged", ]
View(repro.t6.comp)



#trying the new models with my data####

#dummy data
mod<-lmer(RESP~T*AC+(1+AC|POP),data=xxxx)
Anova(mod, type=c("II","III", 2, 3), test.statistic=c("F"))

modg<-lmer(egg.week~Treatment*Year*avg.inhab+(1+Year|Trial),data=repro)
Anova(modg, type=c("II","III", 2, 3), test.statistic=c("F"))

#looking at model assumptions using egg.week vs. sqrt.egg

#Note: I modeled this as fixed, just to test assumptions
mod.e<-lm(egg.week~Treatment*Year*Trial*avg.inhab,data=repro)
hist(resid(mod.e))
plot(mod.e)
summary(mod.e)
anova(mod.e)

mod.se<-lm(sqrt.egg~Treatment*Year*Trial*avg.inhab,data=repro)
hist(resid(mod.se))
plot(mod.se)
summary(mod.se)
anova(mod.se)




# Description: reanalyzing recollection data, with new model framework
# Author: George C Jarvis
# Date: Mon Dec 09 21:22:41 2019
# Notes: going to include year, and trial nested within year
# then run the models. Though I don't think I will reduce it this time
# because there aren't any interaction terms that don't involve the categorical factors (year.fact and treatment)
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

#read data
reco<-read.csv("Data/2019.10.8.recollection.data.csv")
reco<-na.omit(reco)
#adding column for survivorship, dividing counts by 
#- 20 (initial number of fish on reefs)
reco$Survivorship<-reco$Count/20
#already have a column for year, going to make it a factor now
reco$Year.fact<- as.factor(reco$Year)

# - and I think it's okay to leave it as a proportion, and not whole numbers
# I think I should anlayze it this way, instead of doing the raw recollections? Will
# ask Mark
#for now, I'll analyze it with the raw data, and also with the proportional data to 
# see how the results change

#models####
 
#using raw counts, including year, and trial nested within year, year as numeric
mod1<-lmer(Count~Treatment*Year.fact+(1|Year.fact:Trial), data=reco)
hist(resid(mod1))
qqnorm(resid(mod1))
qqline(resid(mod1))
anova(mod1)
Anova(mod1)
summary(mod1) #no differences in either categorical factor, or interaction
ranef(mod1) # it seems like there may have been some differences by trial 
# - fewer recollected in trials 1,2, and 6

#using proportional data for survival
mod2<-lmer(Survivorship~Treatment*Year.fact+(1|Year.fact:Trial), data=reco)
hist(resid(mod2))
qqnorm(resid(mod2))
qqline(resid(mod2))
anova(mod2)
Anova(mod2)
summary(mod2) #no differences in either categorical factor, or interaction
ranef(mod2)

#plotting####
#ordering treatments
reco$Treatment.ord<-ordered(reco$Treatment, levels=c("Low","Medium","High"))
reco$T6.comparison.ord<-ordered(reco$T6.comparison, levels=c("Low","Medium","High","Uncaged"))

#recollection by treatment and year, raw counts
bargraph.CI(x.factor = Treatment.ord, response = Count, 
            group= Year.fact, legend=TRUE, main="Recollection by treatment and year", 
            data = reco, ylab="Gobies Recollected")
#no year factor
bargraph.CI(x.factor = Treatment.ord, response = Count, 
            main="Recollection by treatment", 
            data = reco, ylab="Gobies Recollected")

#same thing, but with proportion surviving
bargraph.CI(x.factor = Treatment.ord, response = Survivorship, 
            group= Year.fact, legend=TRUE, main="Survivorship by treatment and year", 
            data = reco, ylab="Survivorship")
#no year factor
bargraph.CI(x.factor = Treatment.ord, response = Survivorship, 
            main="Survivorship by treatment", 
            data = reco, ylab="Survivorship")





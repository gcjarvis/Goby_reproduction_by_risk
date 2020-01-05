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

options(contrasts = c("contr.sum","contr.poly")) #this is important, run before ANOVA, will set SS to type III

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

#subsetting data to do t6 comparison between HR caged and uncaged treatments
#subset trial
reco.t6<-reco[reco$Trial==6,]
#subset treatments, caged and uncaged only (T6.comparison)

#subsetting by multiple character factors within a variable. Awesome code!
reco.t6.comp<-reco.t6[reco.t6$T6.comparison == "High" | reco.t6$T6.comparison == "Control", ]
View(reco.t6.comp)

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

#using nlme package to model recollections
mod2.2.luk<-lme(Survivorship~Treatment*Year.fact,random=~1|Trial,reco,method="REML")
hist(resid(mod2.2.luk))
qqnorm(resid(mod2.2.luk))
qqline(resid(mod2.2.luk))
summary(mod2.2.luk)
anova(mod2.2.luk, type='marginal')

#arithmetic means for comparison and ss
library(FSA)

Summarize(Survivorship~Treatment,
          data=reco,
          digits=3)

#testing for differences in output between HR caged and uncaged treatments####

#NOTES: 1) trial 6 only, using "repro.t6.comp" df, and "T6.comparison" for treatments
# 2) it's no longer a nested ANCOVA, it's just an ANCOVA with trt and avg.inhab,
# - so it's just a linear model
# 3) will add the reduced model results to the table for egg counts,
# - and the full model results to the supplementary table for egg counts

#full model
modt6.luk<-lm(Survivorship~T6.comparison,data=reco.t6.comp)
summary(modt6.luk)
anova(modt6.luk)

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





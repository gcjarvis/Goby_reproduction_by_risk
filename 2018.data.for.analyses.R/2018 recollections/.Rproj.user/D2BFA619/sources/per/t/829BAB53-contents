##############################################################################
# recollection data                                                         #
#   George Jarvis                                                           #
#    10/29/2018                                                             #
##############################################################################

rm(list=ls())

#load packages
library(sciplot)
library(lme4)
library(lmerTest)
library(car)
library(dplyr)
library(ggplot2)
library(extrafont)
#library(ggplot)

getwd()
reco<-read.csv("Data/recollections.cleaned.10.29.18.csv")
#sc$Trial<-as.factor(sc$Trial)#made trial a factor
#only including adults that we tagged initially (immatures were listed as NA in datasheet)
#sc<- sc[complete.cases(sc), ]#122 rows of data

head(reco)
tail(reco)

#looking at the proportion of total fish recollected from each reef
reco.mod<-lm(prop.total.recollected~treatment, data=reco)
hist(resid(reco.mod))
qqnorm(resid(reco.mod))
boxplot(resid(reco.mod))
aov(reco.mod)
anova(reco.mod)
#no diff in the proportion of fish recollected on each reef

bargraph.CI(x.factor = treatment, response = prop.total.recollected, xlab="Risk treatment", ylab="Proportion of population recollected", data = reco)
#no diff, looks like high had lowest proportion of fish recollected

#number of immatures/treatment (proxy for recruitment??)
bargraph.CI(x.factor = treatment, response = num.immature, main="number of immatures", data = reco)

#now want to break it down by sex
#proportion of males
bargraph.CI(x.factor = treatment, response = prop.male.recollected, main="proportion of males recollected", data = reco)

#prop. female
bargraph.CI(x.factor = treatment, response = prop.female.recollected, main="proportion of females recollected", data = reco)
#equal proportions of males and females recollected

#now want to look at biomass
#total female biomass
bargraph.CI(x.factor = treatment, response = total.female.biomass, main="biomass of females recollected", data = reco)

#total male biomass
bargraph.CI(x.factor = treatment, response = total.male.biomass, main="biomass of males recollected", data = reco)
#INTERESTING: looking at this, it seems like I collected smaller (lighter) males in the
# high-risk treatment, despite the fact that I recollected the same prop of males 
# among all three risk treatments
# I doubt any of these are statistically significant, however.

#######model selection for results######
mod1.reco<-lm(total.recollection~treatment, data=reco)
hist(resid(mod1.reco))
qqnorm(resid(mod1.reco))
boxplot(resid(mod1.reco))
aov(mod1.reco)
anova(mod1.reco)

mod2.reco<-lmer(total.recollection~treatment+(1|deployment.day), data=reco)
hist(resid(mod2.reco))
qqnorm(resid(mod2.reco))
qqline(resid(mod2.reco))

mod3.reco.poiss<-glmer(total.recollection~treatment+(1|deployment.day),family=poisson,data=reco)
hist(resid(mod3.reco.poiss))

#going to try gamma with the proportion recollected data
mod4.reco.gamma<-glmer(prop.total.recollected~treatment+(1|deployment.day),family=Gamma,data=reco)
hist(resid(mod4.reco.gamma))

AIC(mod1.reco)
AIC(mod3.reco.poiss)#better AIC value

#########this is the model tht I used for my initial analyses##########
mod3.reco.poiss<-glmer(total.recollection~treatment+(1|deployment.day),family=poisson,data=reco)
hist(resid(mod3.reco.poiss))
Anova(mod3.reco.poiss,type="II")

bargraph.CI(x.factor = treatment, response = total.recollection, main="proportion of total pop recollected", data = reco)
#no diff, looks like high had lowest proportion of fish recollected

bargraph.CI(x.factor = treatment, response = prop.total.recollected, main="proportion of total pop recollected", data = reco)
#no diff, looks like high had lowest proportion of fish recollected

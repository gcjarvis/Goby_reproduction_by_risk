# Description: Script for behavioral analyses from Jarvis and Steele 
# Author: George C Jarvis
# Date: Sat Apr 11 17:25:38 2020
# Notes: No behavioral observations for Trial 6, so no comparisons of behaviors between high-risk
## caged and uncaged treatments
# --------------

rm(list=ls())

library(lme4)
library(lmerTest)
library(dplyr)
library(ggplot2)
library(emmeans) #for generating least-squares adjusted means from models 

#importing data####
behave<-read.csv("Data/2019.10.25.behavior.includes.recollections.csv")

#data manipulation####
#adding a column for year (as a proxy for tagging procedure), where trials 1-3 = 2017, and 4-6 = 2018
# NOTE: there were no behavioral observations for trial 6
behave$Year <- ifelse(behave$Trial <=3, 2017, 2018)
#making Year and Trial factors
behave$Year<- as.factor(behave$Year)
behave$Trial<- as.factor(behave$Trial)

#making the variable "avg.inhab" ((20+reco)/2), rounded to the nearest whole fish
#using the average number of inhabitants per reef as the covariate in mixed models
behave$avg.inhab<-(ceiling((behave$Recollection+20)/2))
pairs(behave)# data viz, but also seems to mess up the plot parameters?
dev.off() #turns off the plot parameters

#exporting data for non-R users
#write.csv(behave,"Data\\behavioral_analyses_after_data_wrangling.csv", row.names = FALSE)

#analyses####
# analyzing mixed models with log-likelihood estimates and chi-square tests.
# start with full model, then reduce model first by non-significant random effects, then by NS. fixed effects
# want to definitely leave in Trial term as random effect, and Treatement, Year, T x Y, and avg.inhab as fixed effects.
# remove all NS. interactions with the covariate (fixed and random)

options(contrasts = c("contr.sum","contr.poly")) #this is important, run before ANOVA, will set SS to type III

# 1. proportion of time exposed ####
pe<-lmer(proportion.exposed~Treatment*Year*avg.inhab+(1|Trial) + (1|Treatment:Trial) +
           (1|Trial:avg.inhab)+(1|avg.inhab:Treatment:Trial), REML=F, behave)
hist(resid(pe))
qqnorm(resid(pe))
qqline(resid(pe))
plot(pe)
summary(pe) #don't seem to be very many differences based on trial
anova(pe)

#removing three-way interaction of random effect
pe2<-update(pe, .~. -(1|avg.inhab:Treatment:Trial))
summary(pe2)
anova(pe2)

anova(pe,pe2) #no difference in models when three-way interaction with random intercept is removed

#next logical removal is the 1|treatment:trial term

pe3<-update(pe2, .~. -(1|Treatment:Trial))
summary(pe3)
anova(pe3)

anova(pe2,pe3) #no difference

#next logical removal is the 1|Trial:avg.inhab term

pe4<-update(pe3, .~. -(1|Trial:avg.inhab)) #all interactions with covariate that were N.S. were removed, including random effects
## in fact, it may put more variance in the overall model, including variance for fixed effects
summary(pe4)
anova(pe4)

anova(pe3,pe4)

#trial effect
pe5<-update(pe2,.~. -(1|Trial))
summary(pe5)
anova(pe5)

anova(pe2,pe5) # no sig. effect of trial (p=1, so maybe pool by trials? or does it even matter to keep it in?)

#might leave trial in just because I want it to soak up some variation, although it is very small

#testing effects of fixed factors, going to start with highest-order interactions with the covariate

pe6<-update(pe4, .~. -(Treatment:Year:avg.inhab))
summary(pe6) 
anova(pe6)

anova(pe4,pe6)

#now removing Year.fact:avg.inhab (the order of removal seems arbitrary, but that's the next one)
pe7<-update(pe6, .~. -(Year:avg.inhab))
summary(pe7) 
anova(pe7)

anova(pe7,pe6)

#now want to take out Treatment:avg.inhab

pe8<-update(pe7, .~. -(Treatment:avg.inhab)) 
summary(pe8) 
anova(pe8)

anova(pe7,pe8)

#me8 is the final model that I will likely end up with (those are all of the fixed and random factors that I'm interested in)
hist(resid(pe8))#not bad
#LS-adjusted means from model
emmeans(pe8, pairwise~Treatment) #warning message re: interactions, but I think it's okay
boxplot(proportion.exposed~Treatment,data=behave)# variances don't look too bad, and medians look pretty good to me (i.e. match up with LS-means for the most part)

#dropping the fixed effects that I have reduced the model to do the log-likelihood estimates
#NOTE: m8 is the fully-reduced model, so just have to jeep iterating that model +/- individual fixed factors of interest

# - avg.inhab, but keeping in other fixed factors

pe9<-update(pe8, .~. -(avg.inhab))
summary(pe9) anova(pe9)

anova(pe8,pe9)

# - Treatment*Year.fact

pe10<-update(pe8, .~. -(Treatment:Year))
summary(pe10)
anova(pe10)

anova(pe8,pe10)

# - Treatment

pe11<-update(pe10,.~. -(Treatment))
summary(pe11)
anova(pe11)

anova(pe10,pe11) #significant effect of treatment (P = 0.01)

# - Year
#again, dropping from the model that does not include the Treatment:Year interaction

pe12<-update(pe10,.~. -(Year))
summary(pe12)
anova(pe12)

anova(pe10,pe12)

#not sure how best to display results for chi-square estimates with multiple levels of a factor

# 1b. linear distance traveled ####
dm<-lmer(total.dist.moved~Treatment*Year*avg.inhab+(1|Trial) + (1|Treatment:Trial) +
           (1|Trial:avg.inhab)+(1|avg.inhab:Treatment:Trial), REML=F, behave)
hist(resid(dm))
qqnorm(resid(dm))
qqline(resid(dm))
plot(dm)
summary(dm) #seems to be a bit more variation by trial
anova(dm)

#removing three-way interaction of random effect
dm2<-update(dm, .~. -(1|avg.inhab:Treatment:Trial))
summary(dm2)
anova(dm2)

anova(dm,dm2) #no difference in models when three-way interaction with random intercept is removed

#next logical removal is the 1|treatment:trial term

dm3<-update(dm2, .~. -(1|Treatment:Trial))
summary(dm3)
anova(dm3)

anova(dm2,dm3) #no difference

#next logical removal is the 1|Trial:avg.inhab term

dm4<-update(dm3, .~. -(1|Trial:avg.inhab)) #all interactions with covariate that were N.S. were removed, including random effects
## in fact, it may put more variance in the overall model, including variance for fixed effects
summary(dm4)
anova(dm4)

anova(dm3,dm4)

#trial effect
dm5<-update(dm2,.~. -(1|Trial))
summary(dm5)
anova(dm5)

anova(dm2,dm5) #  sig. effect of trial (P = 0.01), so keeping in the model

#testing effects of fixed factors, going to start with highest-order interactions with the covariate

dm6<-update(dm4, .~. -(Treatment:Year:avg.inhab))
summary(dm6) 
anova(dm6)

anova(dm4,dm6) # P = 0.07, so dropped from model 

#now removing Year.fact:avg.inhab (the order of removal seems arbitrary, but that's the next one)
dm7<-update(dm6, .~. -(Year:avg.inhab))
summary(dm7) 
anova(dm7)

anova(dm7,dm6)

#now want to take out Treatment:avg.inhab

dm8<-update(dm7, .~. -(Treatment:avg.inhab)) 
summary(dm8) 
anova(dm8) #significant effect of treatment at P = 0.047565, and sig. effect of avg.inhab,
## it seems that the difference is when I include Treatment*Year in the model

#If I have trt*Year in the model, then there is a sig.effect of trrt.by itself
#when I take trt*year out of the model, then there is no effect of Trt. by itself

#this is a tough call, because my final model will have Treatment:Year in it, but it will be 
## analyzed as a chi-square test, which won't run if I include trt*year in the model...

anova(dm7,dm8)

#me8 is the final model that I will likely end up with (those are all of the fixed and random factors that I'm interested in)
hist(resid(dm8))#not bad
#LS-adjusted means from model
emmeans(dm8, pairwise~Treatment) #warning message re: interactions, but I think it's okay
boxplot(total.dist.moved~Treatment,data=behave)# variances don't look too bad, and medians look pretty good to me (i.e. match up with LS-means for the most part)

#dropping the fixed effects that I have reduced the model to do the log-likelihood estimates
#NOTE: m8 is the fully-reduced model, so just have to jeep iterating that model +/- individual fixed factors of interest

# - avg.inhab, but keeping in other fixed factors

dm9<-update(dm8, .~. -(avg.inhab))
summary(dm9) anova(dm9)

anova(dm8,dm9)

# - Treatment*Year.fact

dm10<-update(dm8, .~. -(Treatment:Year))
summary(dm10)
anova(dm10)

anova(dm8,dm10)

# - Treatment

dm11<-update(dm10,.~. -(Treatment))
summary(dm11)
anova(dm11)

anova(dm10,dm11) # interesting. Chi-square test shows no effect of treatment in the model, but an 
## ANOVA of the recuded model (dm8) shows a sig. effect of treatment. I'm assuming it takes into account the
## effect of trial as well? This is an impass. If I analyze with chi-square test, then treatment has no effect
## but if I analyze as ANOVA, there is an effect. I'm inclined to go with chi-square results

#at least that's what I'll go with for now

# - Year
#again, dropping from the model that does not include the Treatment:Year interaction

dm12<-update(dm10,.~. -(Year))
summary(dm12)
anova(dm12)

anova(dm10,dm12)

#seeing if I can remove treatment from the model while keeping trt*year in as well

dm13<-update(dm8,.~. -(Treatment))
summary(dm13)
anova(dm13)

anova(dm8,dm13) #well, this stinks. I can't have it both ways

#the issue is that there is a big effect of trial for linear distance moved, so there is 
## good evidence to keep it in the model and analyze with chi-square tests, but there is 
## also a significant effect of treatment when the full model (including trt*year) is analyzed
## as an ANOVA

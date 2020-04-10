# Description: Script for reproductive output from Jarvis and Steele 
# Author: George C Jarvis
# Date: Thu Apr 09 15:14:09 2020
# Notes:
# --------------

rm(list=ls())

library(nlme)
library(lme4)
library(lmerTest)
library(dplyr)
library(ggplot2)
library(visreg) #visualizing linear models
library(emmeans) #for generating least-squares adjusted means from models (will likely help when I have to
## back-transform means if data transformation is needed)

#importing dataset####
repro<-read.csv("Data/goby_reproduction.2020.4.9.csv")

#intial data viz
pairs(repro)


#data manipulation####

#adding column for average number of inhabitants throughout the trial, rounded to nearest whole number of fish
repro$avg.inhab<-(ceiling((repro$Recollection+20)/2))

#adding egg/week variable to the dataset for comparisons among trials of differing length
repro<-repro %>%
  mutate(egg.week = ifelse(Trial<4, Egg.count/1,
                           ifelse(Trial == 4| Trial == 5, (ceiling(Egg.count/4)),
                                  ifelse(Trial == 6, (ceiling(Egg.count/2)), NA))))

#adding a column for year (as a proxy for tagging procedure), where trials 1-3 = 2017, and 4-6 = 2018
repro$Year <- ifelse(repro$Trial <=3, 2017, 2018)
#making Year and Trial factors
repro$Year<- as.factor(repro$Year)
repro$Trial<-as.factor(repro$Trial)

#sqrt-transforming egg counts to better satisfy assumptions
repro$sqrt.egg.week<-sqrt(repro$egg.week)

#subsetting data from Trial 6 for comparison of caged vs. uncaged treatments
#subset trial
repro.t6<-repro[repro$Trial==6,]
#subset treatments, caged and uncaged only to test for cage effects
## Note that Treatment factor for trial 6 models should be replaced with "T6.comparison"
repro.t6<-repro.t6[repro.t6$T6.comparison == "High" | repro.t6$T6.comparison == "Uncaged", ]

#exporting wrangled data for those that are not working with these data in R
#Data from Trials 1-5
#write.csv(repro,"Data\\egg_counts_after_data_wrangling.2020.4.7.csv", row.names = FALSE)
#Data from Trial 6 only
#write.csv(repro.t6,"Data\\Cage_effects_egg_counts_after_data_wrangling.2020.4.7.csv", row.names = FALSE)


#analyses####
# analyzing mixed models with log-likelihood estimates and chi-square tests.
# start with full model, then reduce model first by non-significant random effects, then by NS. fixed effects
# want to definitely leave in Trial term as random effect, and Treatement, Year, T x Y, and avg.inhab as fixed effects.
# remove all NS. interactions with the covariate (fixed and random)

options(contrasts = c("contr.sum","contr.poly")) #this is important, run before ANOVA, will set SS to type III

#model structure for random effects: all slopes are the same (effcets are the same among trials), but intercepts are random (magnitude differs)
## Note: use "REML=F" (maximum likelihood estimates) to compare models with log-likelihood estimates (Pinheiro & Bates, 2000; Bolker et al., 2009)

# 1. all trials, raw data for reproduction per week####
# full model
me<-lmer(egg.week ~ Treatment*Year*avg.inhab + (1|Trial) + (1|Treatment:Trial) +
           (1|Trial:avg.inhab)+(1|avg.inhab:Treatment:Trial), REML=F, repro)
hist(resid(me))
qqnorm(resid(me))
qqline(resid(me))
plot(me)
summary(me)
anova(me)

#removing three-way interaction of random effect
me2<-update(me, .~. -(1|avg.inhab:Treatment:Trial))
summary(me2)
anova(me2)

anova(me,me2) #no difference in models when three-way interaction with random intercept is removed

#next logical removal is the 1|treatment:trial term

me3<-update(me2, .~. -(1|Treatment:Trial))
summary(me3)
anova(me3)

anova(me2,me3) #no difference

#next logical removal is the 1|Trial:avg.inhab term

me4<-update(me3, .~. -(1|Trial:avg.inhab)) #all interactions with covariate that were N.S. were removed, including random effects
## in fact, it may put more variance in the overall model, including variance for fixed effects
summary(me4) # trial accounts for 6% of the residual variance
anova(me4)

anova(me3,me4) # P = 0.3939, doesn't suggest any difference due to the dropping of the random effect of Trial:avg.inhab

#trial effect
me5<-update(me2,.~. -(1|Trial))
summary(me5)
anova(me5)

anova(me2,me5) # no sig. effect of trial

# decided to leave the random factor of 1|Trial in the model, after removing all other nonsignificant random effects at P>0.25
## as tested for with log-likelihood tests

# me4 is the model where all random effects aside from 1|Trial have been removed

#testing effects of fixed factors, going to start with highest-order interactions with the covariate

me6<-update(me4, .~. -(Treatment:Year:avg.inhab))
summary(me6) # trial accounts for 6% of the residual variance
anova(me6)

anova(me4,me6)# P = 0.5013

#now removing Year.fact:avg.inhab (the order of removal seems arbitrary, but that's the next one)
me7<-update(me6, .~. -(Year:avg.inhab))
summary(me7) # trial variance is still 6% of residual variance for random effects
anova(me7)

anova(me7,me6)# P = 0.779

#now want to take out Treatment:avg.inhab

me8<-update(me7, .~. -(Treatment:avg.inhab)) 
summary(me8) # trial variance is still 6% of residual variance for random effects
anova(me8)

anova(me7,me8)# P = 0.7064

#me8 is the final model that I will likely end up with (those are all of the fixed and random factors that I'm interested in)
hist(resid(me8))#not bad
#LS-adjusted means from model
emmeans(me8, pairwise~Treatment) #warning message re: interactions, but I think it's okay
boxplot(egg.week~Treatment,data=repro)# variances don't look too bad, and medians look pretty good to me (i.e. match up with LS-means for the most part)

#dropping the fixed effects that I have reduced the model to do the log-likelihood estimates
#NOTE: m8 is the fully-reduced model, so just have to jeep iterating that model +/- individual fixed factors of interest

# - avg.inhab, but keeping in other fixed factors

## should see a sig. log-likelihood result
me9<-update(me8, .~. -(avg.inhab))
summary(me9) # residual variance for random effects just shot up a ton (soaked up all variance from avg.inhab)
anova(me9)

anova(me8,me9)#chisq = 38.021      df = 1      p = 7e-10, super significant, so don't want to remove that term
#also going to include the estimate for avg.inhab from the previous model (m8)

# - Treatment*Year.fact

me10<-update(me8, .~. -(Treatment:Year))
summary(me10)
anova(me10)

anova(me8,me10) #chisq = 0.5296     df= 2     p = 0.7674

# - Treatment

me11<-update(me10,.~. -(Treatment)) #was getting 0 df when dropping treatment, 
# and it is because Treatment:Year was still in the model
summary(me11)
anova(me11)

anova(me10,me11) #chisq = 2.26     df = 2         p = 0.32 #this is acceptable, amanda has this in her paper for a factor

# - Year
#again, dropping from the model that does not include the Treatment:Year interaction

me12<-update(me10,.~. -(Year))
summary(me12)
anova(me12)

anova(me10,me12)

# 1a. Trial 6 data only, comparing high-risk caged to uncaged treatments in Trial 6 with ANCOVA####
#note: Year is not needed anymore, and Trial is removed as well

me13<-lm(egg.week~T6.comparison*avg.inhab,data=repro.t6)
hist(resid(me13)) # not so great, but might look better after data transformation
qqnorm(resid(me13))
qqline(resid(me13))

anova(me13)

# removing non-significant interaction of covariate

me14<-lm(egg.week~T6.comparison+avg.inhab,data=repro.t6)
hist(resid(me14)) #not a super great fit, but I think it's okay
qqnorm(resid(me14))
qqline(resid(me14))

anova(me13,me14) # no difference in model after interaction was removed

# LS-means for treatment from model:
emmeans(me14, pairwise~T6.comparison)
boxplot(egg.week~T6.comparison,data=repro.t6)# variances don't look too bad

# 1b. same analysis with log-likelihood estimates, but with sqrt.eggs as response####
## all I have to do is change the orignial model
mes<-lmer(sqrt.egg.week ~ Treatment*Year*avg.inhab + (1|Trial) + (1|Treatment:Trial) +
           (1|Trial:avg.inhab)+(1|avg.inhab:Treatment:Trial), REML=F, repro)
hist(resid(mes))
qqnorm(resid(mes))
qqline(resid(mes))
plot(mes)
summary(mes)
anova(mes)

anova(me,mes)# seems like the log-transformed data might be a better model, based on lower AIC

#removing three-way interaction of random effect
mes2<-update(mes, .~. -(1|avg.inhab:Treatment:Trial))
summary(mes2)
anova(mes2)

anova(mes,mes2) #no difference in models when three-way interaction with random intercept is removed

#next logical removal is the 1|treatment:trial term

mes3<-update(mes2, .~. -(1|Treatment:Trial))
summary(mes3)
anova(mes3)

anova(mes2,mes3) #no difference

#next logical removal is the 1|Trial:avg.inhab term

mes4<-update(mes3, .~. -(1|Trial:avg.inhab)) #all interactions with covariate that were N.S. were removed, including random effects
## in fact, it may put more variance in the overall model, including variance for fixed effects
summary(mes4) # trial accounts for 4% of the residual variance
anova(mes4)

anova(mes3,mes4) # P = 1, which is different from the analyses for raw data

#trial effect
mes5<-update(mes2,.~. -(1|Trial))
summary(mes5)
anova(mes5)

anova(mes2,mes5) # no sig. effect of trial

# decided to leave the random factor of 1|Trial in the model, after removing all other nonsignificant random effects at P>0.25
## as tested for with log-likelihood tests

# me4 is the model where all random effects aside from 1|Trial have been removed

#testing effects of fixed factors, going to start with highest-order interactions with the covariate

mes6<-update(mes4, .~. -(Treatment:Year:avg.inhab))
summary(mes6) # trial accounts for 4% of the residual variance
anova(mes6)

anova(mes4,mes6)# P = 0.42

#now removing Year.fact:avg.inhab (the order of removal seems arbitrary, but that's the next one)
mes7<-update(mes6, .~. -(Year:avg.inhab))
summary(mes7) # trial variance is still 6% of residual variance for random effects
anova(mes7)

anova(mes7,mes6)# P = 0.91

#now want to take out Treatment:avg.inhab

mes8<-update(mes7, .~. -(Treatment:avg.inhab)) 
summary(mes8) # trial variance is still 6% of residual variance for random effects
anova(mes8)

anova(mes7,mes8)# P = 0.86

#me8 is the final model that I will likely end up with (those are all of the fixed and random factors that I'm interested in)
hist(resid(mes8))#not bad
#LS-adjusted means from model
emmeans(mes8, pairwise~Treatment) #warning message re: interactions, but I think it's okay
boxplot(sqrt.egg.week~Treatment,data=repro)# variances don't look too bad, and medians look pretty good to me (i.e. match up with LS-means for the most part)

#dropping the fixed effects that I have reduced the model to do the log-likelihood estimates
#NOTE: m8 is the fully-reduced model, so just have to jeep iterating that model +/- individual fixed factors of interest

# - avg.inhab, but keeping in other fixed factors

## should see a sig. log-likelihood result
me9<-update(me8, .~. -(avg.inhab))
summary(me9) # residual variance for random effects just shot up a ton (soaked up all variance from avg.inhab)
anova(me9)

anova(me8,me9)#chisq = 38.021      df = 1      p = 7e-10, super significant, so don't want to remove that term
#also going to include the estimate for avg.inhab from the previous model (m8)

# - Treatment*Year.fact

me10<-update(me8, .~. -(Treatment:Year))
summary(me10)
anova(me10)

anova(me8,me10) #chisq = 0.5296     df= 2     p = 0.7674

# - Treatment

me11<-update(me10,.~. -(Treatment)) #was getting 0 df when dropping treatment, 
# and it is because Treatment:Year was still in the model
summary(me11)
anova(me11)

anova(me10,me11) #chisq = 2.26     df = 2         p = 0.32 #this is acceptable, amanda has this in her paper for a factor

# - Year
#again, dropping from the model that does not include the Treatment:Year interaction

me12<-update(me10,.~. -(Year))
summary(me12)
anova(me12)

anova(me10,me12)

# 1a. Trial 6 data only, comparing high-risk caged to uncaged treatments in Trial 6 with ANCOVA####
#note: Year is not needed anymore, and Trial is removed as well

me13<-lm(egg.week~T6.comparison*avg.inhab,data=repro.t6)
hist(resid(me13)) # not so great, but might look better after data transformation
qqnorm(resid(me13))
qqline(resid(me13))

anova(me13)

# removing non-significant interaction of covariate

me14<-lm(egg.week~T6.comparison+avg.inhab,data=repro.t6)
hist(resid(me14)) #not a super great fit, but I think it's okay
qqnorm(resid(me14))
qqline(resid(me14))

anova(me13,me14) # no difference in model after interaction was removed

# LS-means for treatment from model:
emmeans(me14, pairwise~T6.comparison)
boxplot(egg.week~T6.comparison,data=repro.t6)# variances don't look too bad


me<-lmer(sqrt(egg.week) ~ Treatment*Year.fact*avg.inhab + (1|Trial.fact) + (1|Treatment:Trial.fact) +
           (1|Trial.fact:avg.inhab)+(1|avg.inhab:Treatment:Trial.fact), REML=F, repro)
hist(resid(me)) #I honestly think that 
qqnorm(resid(me))
qqline(resid(me))
plot(me)

summary(me)
anova(me) #Okay, this seems to be estimating denDF for year, avg.inhab, and avg.inhab*year sort of correctly
#but there are tons of error messages
coef(me)

me2<-update(me, .~. -(1|avg.inhab:Treatment:Trial.fact))
summary(me2)
anova(me2) #Okay, this seems to be estimating denDF for year, avg.inhab, and avg.inhab*year sort of correctly
#but there are tons of error messages
coef(me2)

anova(me,me2) #no difference in models when three-waty interaction with random intercept is removed

#next logical removal is the 1|treatment:trial term

me3<-update(me2, .~. -(1|Treatment:Trial.fact))
summary(me3)
anova(me3)

anova(me2,me3) #no difference

#next logical removal is the 1|Trial:avg.inhab term, but I'm skeptical to take that out b/c it accounts for 12% of the
## residual variance in random effects (in fact, there's more variance attributed to this term than trial alone)
## fortunately, can test whwether it makes sense to remove that term with log-likelihood test comapring m3 with m4

#the key question will be whether to evaluate the effect of avg.inhab with random slope (different effect of covariate among trials),
## and random intercept (same effect of covariate, but different magnitude among trials)

#should compare model with avg.inhab|Trial.fact vs. 1|Trial:avg.inhab and see what log-likelihood shows

me4<-update(me3, .~. -(1|Trial:avg.inhab)) #all interactions with covariate that were N.S. were removed, including random effects
## in fact, it may put more variance in the overall model, including variance for fixed effects
summary(me4) # trial accounts for 6% of the residual variance
anova(me4)

anova(me3,me4) # P = 0.3939, doesn't suggest any differnce due to the dropping of the random effect of Trial:avg.inhab

#so, decided to leave the random factor of 1|Trial in the model, after removing all other nonsignificant randm effects at P>0.25
## as tested for with log-likelihood tests
#NOTE: this model structure assumes all random effects had random intercepts (1|Random...), but not random slopes among trials

#okay, now doing the same thing with fixed factors, going to start with highest-order interactions with the covariate

me5<-update(me4, .~. -(Treatment:Year.fact:avg.inhab))
summary(me5) # trial accounts for 6% of the residual variance
anova(me5)

anova(me4,me5)# P = 0.5013

#now removing Year.fact:avg.inhab (the order of removal seems arbitrary, but that's the next one)
me6<-update(me5, .~. -(Year.fact:avg.inhab))
summary(me6) # trial variance is still 6% of residual variance for random effects
anova(me6)

anova(me5,me6)# P = 0.779

#now want to take out Treatment:avg.inhab

me7<-update(me6, .~. -(Treatment:avg.inhab))
summary(me7) # trial variance is still 6% of residual variance for random effects
#NOTE: this summary shows me the estimate for the effect of avg inhab (741.827),
## indicating that for every goby added, there's an increase of about 740 eggs on average

#the trouble that I'm anticipating is that there might be a case where one of the years or the interactions
## is significant, not just the term itself (e.g. avg.inhab)

# you'll see that Treatment1 and Treatment2 are in the model, so I guess I can show what treatment si driving the
## trend? We'll cross that bridge when we come to it (likely behavior analyses)

anova(me7)

anova(me6,me7)# P = 0.7064

#dropping the fixed effects that I have reduced the model to do the log-likelihood estimates
#NOTE: m7 is the fully-reduced model, so just have to jeep iterating that model +/- individual fixed factors

# - avg.inhab, but keeping in other fixed factors

## should see a sig. log-likelihood result
me8<-update(me7, .~. -(avg.inhab))
summary(me8) # residual variance for random effects just shot up a ton (soaked up all variance from avg.inhab)
anova(me8)

anova(me7,me8)#chisq = 38.021      df = 1      p = 7e-10, super significant, so don't want to remove that term

# - Treatment*Year.fact

me9<-update(me7, .~. -(Treatment:Year.fact))
summary(me9)
anova(me9)

anova(me7,me9) #chisq = 0.5296     df= 2     p = 0.7674

# - Treatment

me10<-update(me7,.~. -(Treatment))
summary(me10)
anova(me10)

anova(me7,me10)

#not done yet...

#next step will be to go back and rewrite the model where the slope for avg.inhab AND the intercept for avg.inhab is random

#the point that I'm not sure about is whether to code the interactions as avg.inhab|Trial, or as avg.inhab|avg.inhab:Trial

# it seems to me that it is redundant, and that I might want to include avg.inhab|Trial in addition to 1|avg.inhab:Trial

#stock code to be applied to my data
lme(y ~ time * tx, 
    random = list(therapist = ~time * tx, 
                  subjects = ~time),
    data=df)

newmod<-lme(egg.week~Treatment*Year.fact*avg.inhab, 
            random=list(Trial.fact=~1,Trial.fact:avg.inhab=~1),data=repro, method="ML")

#comparing lme to lmer models ####
# have to put all of the interactions into single terms so that I can list them in 
repro$trial_inhab_treatment<-paste0(repro$Trial,repro$avg.inhab,repro$Treatment)
repro$trial_treatment<-paste0(repro$Trial,repro$Treatment)
repro$trial_inhab<-paste0(repro$Trial,repro$avg.inhab)
#View(repro$trial_inhab)

#maybe if avg.inhab were a factor?
repro$trial_inhab_factor_treatment<-paste0(repro$Trial.fact,as.factor(repro$avg.inhab),repro$Treatment)
repro$trial_treatment<-paste0(repro$Trial.fact,repro$Treatment)
repro$trial_inhab_factor<-paste0(repro$Trial.fact,as.factor(repro$avg.inhab))

#full model
newmod<-lme(egg.week~Treatment*Year*avg.inhab, 
            random=list(Trial=~1,trial_treatment=~1,trial_inhab=~1,trial_inhab_treatment=~1),data=repro, method="ML")
summary(newmod)
anova(newmod)
#seems like there's something going on with the trial_inhab term, and that's why the models aren't the same, and also why all terms
## that have avg.inhab included in them do not have lower dendf

newmod1<-lme(egg.week~Treatment*Year.fact*avg.inhab, 
             random=list(Trial.fact=~1,trial_treatment=~1,trial_inhab_factor=~1,trial_inhab_factor_treatment=~1),data=repro, method="ML")
summary(newmod1)
anova(newmod1)

#now want to go back and chack the full model with all of the random effects to see if denom df are correct!

#lmer full model
me<-lmer(egg.week ~ Treatment*Year.fact*avg.inhab + (1|Trial.fact) + (1|Treatment:Trial.fact) +
           (1|Trial.fact:avg.inhab)+(1|avg.inhab:Treatment:Trial.fact), REML=F, repro)
summary(me)
anova(me) #Okay, this seems to be estimating denDF for year, avg.inhab, and avg.inhab*year sort of correctly
#but there are tons of error messages
#coef(me)

#SUMMARY: not entirely the same, although it does have a closer df 

boxplot(repro$egg.week~repro$avg.inhab)

n<-1883/(sqrt(379033+135203+3245737))
n

View(repro)

#testing out visreg package with my model to try and figure out why the interept is acting funny

#anyway, not going to put too much thought into this, but it seems weird
library(visreg)

model <- lmer(sqrt(egg.week) ~ Treatment + avg.inhab + Year.fact + Treatment:Year.fact +(1|Trial.fact), data=repro)
visreg(model, "Treatment", by="Year.fact")
visreg(model, "Treatment")

#not sure why, when I include avg.inhab in the model,
## that it shows values for reproduction that are less than 0?
## It might have something to do with the structure of the model?
## It's only noticeable in the plot for 2017, which shows negative values

or (if no interaction):
  visreg(model, ‘Treatment’)



# --------------------------------------
# Description: quick analyses before meeting with Mark
# Author: George C. Jarvis
# Date: Tue Jul 28 06:54:38 2020
# Notes: rerun original analysis now that I've corrected Trial 1 eggs per week (was 2 weeks long, not 1)
# --------------------------------------

rm(list=ls())

library(lme4)
library(lmerTest)
library(dplyr)
library(ggplot2)
library(emmeans) #for generating least-squares adjusted means from models (will likely help when I have to
## back-transform means if data transformation is needed)

#importing dataset####
repro<-read.csv("Data/egg.counts_correced_trial_1.csv")
repro1<-read.csv("Data/15.7.20.Jarvis_egg_counts-per_capita_biomass_fem.csv")
repro1<-na.omit(repro1) #remove cases where no females were recollected

repro1$Year<- as.factor(repro1$Year)
repro1$Trial<-as.factor(repro1$Trial)

#intial data viz
pairs(repro)

#data manipulation####

#adding column for average number of inhabitants throughout the trial, rounded to nearest whole number of fish
#repro$avg.inhab<-(ceiling((repro$Recollection+20)/2))

#adding egg/week variable to the dataset for comparisons among trials of differing length
#repro<-repro %>%
#  mutate(egg.week = ifelse(Trial<4, Egg.count/1,
#                           ifelse(Trial == 4| Trial == 5, (ceiling(Egg.count/4)),
#                                  ifelse(Trial == 6, (ceiling(Egg.count/2)), NA))))

#adding a column for year (as a proxy for tagging procedure), where trials 1-3 = 2017, and 4-6 = 2018
repro$Year <- ifelse(repro$Trial <=3, 2017, 2018)
#making Year and Trial factors
repro$Year<- as.factor(repro$Year)
repro$Trial<-as.factor(repro$Trial)
repro$Treatment.fact<-as.factor(repro$Treatment)

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

# 1a. all trials, raw data for reproduction per week####
# full model
me<-lmer(egg.week ~ Treatment*Year*avg.inhab + (1|Trial) + (1|Treatment:Trial) +
           (1|Trial:avg.inhab)+(1|avg.inhab:Treatment:Trial), REML=F, repro)
hist(resid(me))
qqnorm(resid(me))
qqline(resid(me))
plot(me)
summary(me)
anova(me)

boxplot(egg.week~Treatment,repro)

me$coeff

ranef(me)$Trial

coef(me)

rand(me)

levels(repro$Treatment.fact)

#same model, but with separate slopes for each trial ((avg.inhab|Treatment:Trial))
#doesn't appear to make much of a difference, and the two models are not 
#mer<-lmer(egg.week ~ Treatment*Year*avg.inhab + (1|Trial) + (1|Treatment:Trial) +
#           (1|Trial:avg.inhab)+(avg.inhab|Treatment:Trial), REML=F, repro)
#hist(resid(mer))
#qqnorm(resid(mer))
#qqline(resid(mer))
#plot(mer)
#summary(mer)
#anova(mer)

#anova(me,mer) #not going to use this model in any of the the further analyses


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

# 1b. Trial 6 data only, comparing high-risk caged to uncaged treatments in Trial 6 with ANCOVA####
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

# 1c. same analysis with log-likelihood estimates, but with sqrt.eggs as response####
## all I have to do is change the orignial model
mes<-lmer(sqrt.egg.week ~ Treatment*Year*avg.inhab + (1|Trial) + (1|Treatment:Trial) +
            (1|Trial:avg.inhab)+(1|avg.inhab:Treatment:Trial), REML=F, repro)
hist(resid(mes))
qqnorm(resid(mes))
qqline(resid(mes))
plot(mes)
summary(mes)
anova(mes)

rand(mes)# trial is the only random effect that seems to exaplain any of the variance

#note: this is only the case with sqrt-transformed data

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
mes9<-update(mes8, .~. -(avg.inhab))
summary(mes9) # residual variance for random effects just shot up a ton (soaked up all variance from avg.inhab)
anova(mes9)

anova(mes8,mes9)#chisq = 27.78      df = 1      p = 1.36e-10, super significant, so don't want to remove that term
#also going to include the estimate for avg.inhab from the previous model (mes8)

# after doing log-likelihood tests, only random effect that explained any variance was 
# trial, so I'm going to rerun the reduced model with that, and the non-significant
# interactions between fixed effects



# - Treatment*Year.fact

mes10<-update(mes8, .~. -(Treatment:Year))
summary(mes10)
anova(mes10)

anova(mes8,mes10) #chisq = 0.193     df= 2     p = 0.908

# - Treatment

mes11<-update(mes10,.~. -(Treatment)) #was getting 0 df when dropping treatment, 
# and it is because Treatment:Year was still in the model
summary(mes11)
anova(mes11)

anova(mes10,mes11) #chisq = 1.79    df = 2         p = 0.41 #this is acceptable, amanda has this in her paper for a factor

# - Year
#again, dropping from the model that does not include the Treatment:Year interaction

mes12<-update(mes10,.~. -(Year))
summary(mes12)
anova(mes12)

anova(mes10,mes12) #chisq = 0.002    df = 2         p = 0.41

# 1d. Trial 6 data only, comparing square root data for high-risk caged to uncaged treatments in Trial 6 with ANCOVA####
#note: Year is not needed anymore, and Trial is removed as well

mes13<-lm(sqrt.egg.week~T6.comparison*avg.inhab,data=repro.t6)
hist(resid(mes13)) # not so great, but might look better after data transformation
qqnorm(resid(mes13))
qqline(resid(mes13))

anova(mes13)

# removing non-significant interaction of covariate

mes14<-lm(sqrt.egg.week~T6.comparison+avg.inhab,data=repro.t6)
hist(resid(mes14)) #not a super great fit, but I think it's okay
qqnorm(resid(mes14))
qqline(resid(mes14))

anova(mes13,mes14) # no difference in model after interaction was removed

anova(mes14)

# LS-means for treatment from model:
emmeans(mes14, pairwise~T6.comparison)
boxplot(sqrt.egg.week~T6.comparison,data=repro.t6)# variances don't look too bad

#plotting ANCOVA figure for RAW data; this is what I present in the manuscript####

#Order the treatments: Low --> Med --> High
repro$Treatment.ord<-ordered(repro$Treatment,levels=c("Low","Medium","High"))

#png("Output/2020.4.11.gobies.per.reef.9.5x5.5.300dpi.png", width = 9.5, height = 5.7, units = 'in', res = 300)

anc1<-ggplot(repro, aes(avg.inhab, egg.week, shape=Treatment.ord, linetype=Treatment.ord, col=Treatment.ord)) +
  geom_smooth(method="lm", se=FALSE, show.legend = TRUE)  +
  geom_point(size=3)+
  theme_classic()+
  labs(x=(bquote('Gobies '~Reef^ -1*'')),y=(bquote('Reproducton '~Reef^ -1~ Week^-1*'')))+
  expand_limits(y=0)+
  scale_color_manual(values=c("black", "#666666", "grey"))+
  scale_linetype_manual(values=c("solid", "dashed", "twodash"))+
  theme(axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"), 
        axis.title=element_text(size=20))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), 
        axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0)), 
        axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(legend.text=element_text(size=16)) +
  theme(legend.title =element_text(size=17))+
  scale_x_continuous(breaks=c(10,12,14,16,18,20)) + scale_y_continuous(limits = c(0,12000))+
  labs(color  = "Perceived Risk", linetype = "Perceived Risk", shape = "Perceived Risk")
anc1

#dev.off()

# rerunning same model, but with trial as fixed factor ####

fixed_T<-aov(egg.week ~ Treatment*Year*avg.inhab*(Trial%in%Year), repro)
hist(resid(fixed_T))
qqnorm(resid(fixed_T))
qqline(resid(fixed_T))
plot(fixed_T)
summary(fixed_T)
anova(fixed_T)

boxplot(egg.week~Treatment,repro)

me2<-update(fixed_T, .~. -(Treatment:Year:avg.inhab:Trial))
summary(me2)
anova(me2)

anova(fixed_T,me2)

me3<-update(me3, .~. -(Year:avg.inhab:Trial))
summary(me3)
anova(me3)

anova(me2,me3) #error?

me4<-update(me3, .~. -(Treatment:Year:avg.inhab))
summary(me4)
anova(me4)

anova(me3,me4)

me5<-update(me4, .~. -(Year:avg.inhab))
summary(me5)
anova(me5)

anova(me4,me5)

m1<-aov(egg.week ~ Treatment*Year+avg.inhab+(Trial%in%Year), repro)
hist(resid(m1))
qqnorm(resid(m1))
qqline(resid(m1))
plot(m1)
summary(m1)
anova(m1)


m1a<-aov(egg.week ~ Treatment*Year+avg.inhab+(Trial%in%Year)+Treatment*(Trial%in%Year), repro)
hist(resid(m1a))
qqnorm(resid(m1a))
qqline(resid(m1a))
plot(m1a)
summary(m1a)
anova(m1a)

anova(m1,m1a)

boxplot(egg.week~(Trial%in%Year), repro)

boxplot(egg.week~Treatment*(Trial%in%Year), repro)

# runnning per capita output as response variable with same model
head(repro1)

mcap<-aov(reproduction_per_capita_female_biomass_per_day_per_gram ~ Treatment*Year+recollection_female+(Trial%in%Year), repro1)
hist(resid(mcap))
qqnorm(resid(mcap))
qqline(resid(mcap))
plot(mcap)
summary(mcap)
anova(mcap)

boxplot(reproduction_per_capita_female_biomass_per_day_per_gram~(Trial%in%Year), repro1)

mcap1<-aov(reproduction_per_capita_female_biomass_per_day_per_gram ~ Treatment*Year+recollection_female+Treatment*(Trial%in%Year), repro1)
hist(resid(mcap1))
qqnorm(resid(mcap1))
qqline(resid(mcap1))
plot(mcap1)
summary(mcap1)
anova(mcap1) #no trial x treatment effect

boxplot(recollection_female~Treatment, repro1)

mcap2<-aov(reproduction_per_capita_female_biomass_per_day_per_gram ~ Treatment*Year+recollection_female
           +Treatment*(Trial%in%Year)+Treatment*recollection_female, repro1)
hist(resid(mcap2))
qqnorm(resid(mcap2))
qqline(resid(mcap2))
plot(mcap2)
summary(mcap2)
anova(mcap2)

# findings: 1) sig. differences by trial when included as a nested fixed factor. 
# 2) Same qualitative results for per capita reproduction per female biomass per day (reproduction per female per gram per day)
# 3) 

#trial as random
mcap3<-lmer(reproduction_per_capita_female_biomass_per_day_per_gram ~ Treatment*Year+recollection_female+
           Treatment*recollection_female+(1|Trial)+(1|Treatment:Trial)+(1|recollection_female:Trial), REML = F, repro1)
hist(resid(mcap3))
qqnorm(resid(mcap3))
qqline(resid(mcap3))
plot(mcap3)
summary(mcap3)
anova(mcap3)

### biomass recollected

mcap<-aov(per_capita_female_biomass ~ Treatment*Year+recollection_female+(Trial%in%Year), repro1)
hist(resid(mcap))
qqnorm(resid(mcap))
qqline(resid(mcap))
plot(mcap)
summary(mcap)
anova(mcap)

mcap<-aov(per_capita_female_biomass ~ Treatment*Year+recollection_female+Treatment*(Trial%in%Year), repro1)
hist(resid(mcap))
qqnorm(resid(mcap))
qqline(resid(mcap))
plot(mcap)
summary(mcap)
anova(mcap)

boxplot(per_capita_female_biomass~Treatment, repro1)

mcap<-lmer(per_capita_female_biomass ~ Treatment*Year+recollection_female+(1|Trial), repro1)
hist(resid(mcap))
qqnorm(resid(mcap))
qqline(resid(mcap))
plot(mcap)
summary(mcap)
anova(mcap)


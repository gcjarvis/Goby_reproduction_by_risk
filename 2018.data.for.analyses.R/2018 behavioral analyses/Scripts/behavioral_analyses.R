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
# NOTE: there were no behavioral observations for high-risk uncaged reefs in trial 6
behave$Year <- ifelse(behave$Trial <=3, 2017, 2018)
#making Year and Trial factors
behave$Year<- as.factor(behave$Year)
behave$Trial<- as.factor(behave$Trial)
#View(behave)

#making the variable "avg.inhab" ((20+reco)/2), rounded to the nearest whole fish
#using the average number of inhabitants per reef as the covariate in mixed models
behave$avg.inhab<-(ceiling((behave$Recollection+20)/2))

#data viz
#pairs(behave)# data viz, but also seems to mess up the plot parameters?
#dev.off() #turns off the plot parameters

#changing levels of Treatment so the reference level is Low, 
## - then Medium ("Treatment1"), and High ("Treatment2"), this will make sense for the summary output
behave$Treatment <- factor(behave$Treatment, levels = c("Low", "Medium", "High"))
levels(behave$Treatment)

#exporting data for non-R users
#write.csv(behave,"Data\\behavioral_analyses_after_data_wrangling.csv", row.names = FALSE)

# 1. analyses####
# analyzing mixed models with log-likelihood estimates and chi-square tests.
# start with full model, then reduce model first by non-significant random effects, then by NS. fixed effects
# At the very least, I will include Trial term as random effect, and Treatement, Year, T x Y, and avg.inhab as fixed effects.
# remove all NS. interactions with the covariate (fixed and random)

#NOTE: the values for estimates that I included in the results (supplementary table for behaviors)
## - came from the summary estimates in the fully-reduced model for fixed factors; 
## - (random factors don't have estimates)

options(contrasts = c("contr.sum","contr.poly")) #this is important, run before ANOVA, will set SS to type III

# 1a. proportion of time exposed ####
pe<-lmer(proportion.exposed~Treatment*Year*avg.inhab+(1|Trial) + (1|Treatment:Trial) +
           (1|Trial:avg.inhab)+(1|avg.inhab:Treatment:Trial), REML=F, behave)
hist(resid(pe))
qqnorm(resid(pe))
qqline(resid(pe))
plot(pe)
summary(pe) #don't seem to be very many differences based on trial
anova(pe)

rand(pe) # exposure did not change based on any random effects, removing trial from the
# model and running as a standard ANCOVA

#testing full model first, then reducing if necessary
pe.1<-lm(proportion.exposed~Treatment*Year*avg.inhab, data = behave)
hist(resid(pe.1))
qqnorm(resid(pe.1))
qqline(resid(pe.1))
plot(pe.1)
summary(pe.1)
anova(pe.1)

#no sig. effect of 3-way interaction or between fixed effects and covariate; removing
# those terms from the model

pe.2<-lm(proportion.exposed~Treatment*Year+avg.inhab, data = behave)
hist(resid(pe.2))
qqnorm(resid(pe.2))
qqline(resid(pe.2))
plot(pe.2)
summary(pe.2)
anova(pe.2)

#calculating LS means based on this linear model
emmeans(pe.2, pairwise~Treatment)

#removing covariate b/c not significant
pe.3<-lm(proportion.exposed~Treatment*Year, data = behave)
anova(pe.3)


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
ran<-ranef(pe8)$Trial

coef(pe8)$Trial

anova(pe7,pe8)

pe8a<-lm(proportion.exposed~Treatment*Year*avg.inhab, data=behave) 
summary(pe8a) 
anova(pe8a)

anova(pe8,pe8a)

#me8 is the final model that I will likely end up with (those are all of the fixed and random factors that I'm interested in)
hist(resid(pe8))#not bad
#LS-adjusted means from model
emmeans(pe8, pairwise~Treatment) #warning message re: interactions, but I think it's okay
boxplot(proportion.exposed~Treatment,data=behave)# variances don't look too bad, and medians look pretty good to me (i.e. match up with LS-means for the most part)

boxplot(proportion.exposed~Trial,data=behave)

#dropping the fixed effects that I have reduced the model to do the log-likelihood estimates
#NOTE: m8 is the fully-reduced model, so just have to jeep iterating that model +/- individual fixed factors of interest

# - avg.inhab, but keeping in other fixed factors

pe9<-update(pe8, .~. -(avg.inhab))
summary(pe9) 
anova(pe9)

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

# 1b. linear distance traveled ####
dm<-lmer(total.dist.moved~Treatment*Year*avg.inhab+(1|Trial) + (1|Treatment:Trial) +
           (1|Trial:avg.inhab)+(1|avg.inhab:Treatment:Trial), REML=F, behave)
hist(resid(dm))
qqnorm(resid(dm))
qqline(resid(dm))
plot(dm)
boxplot(total.dist.moved~Treatment,data=behave) #slightly higher variance in Low treatment
summary(dm) #seems to be a bit more variation by trial
anova(dm)

rand(dm) #trial is significant, all others are not, so going to remove them from model

#checking full model, then reducing fixed effects if necessary
dm.1<-lmer(total.dist.moved~Treatment*Year*avg.inhab+(1|Trial), REML=F, data = behave)
anova(dm.1)

#removing non-significant interactions (3-way and NS int with covariate)
dm.2<-lmer(total.dist.moved~Treatment*Year+avg.inhab+(1|Trial), REML=F, data = behave)
anova(dm.2)

#seeing a treatment effect now? Treatment: P = 0.047, so I would say that P<0.05,
# but barely

#checking out LS means

emmeans(dm.2, pairwise~Treatment)
boxplot(total.dist.moved~avg.inhab,data=behave) #total distance moved increased with
# an increased number of gobies on the reef, regardless of risk level

#removing number of inhabitants
dm.3<-lmer(total.dist.moved~Treatment*Year+(1|Trial), REML=F, data = behave)
anova(dm.3)

emmeans(dm.3, pairwise~Treatment)
boxplot(total.dist.moved~avg.inhab,data=behave)

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
anova(dm8) 

anova(dm7,dm8)

#me8 is the final model that I will likely end up with (those are all of the fixed and random factors that I'm interested in)
hist(resid(dm8))#not bad
#LS-adjusted means from model
emmeans(dm8, pairwise~Treatment) #warning message re: interactions, but I think it's okay
boxplot(total.dist.moved~Treatment,data=behave)# variances don't look too bad, and medians look pretty good to me (i.e. match up with LS-means for the most part)
library(visreg)
visreg(dm8)
visreg(dm8, "Treatment",by="Year")

#dropping the fixed effects that I have reduced the model to do the log-likelihood estimates
#NOTE: m8 is the fully-reduced model, so just have to jeep iterating that model +/- individual fixed factors of interest

# - avg.inhab, but keeping in other fixed factors

dm9<-update(dm8, .~. -(avg.inhab))
summary(dm9) 
anova(dm9)

anova(dm8,dm9)

# - Treatment:Year

dm10<-update(dm8, .~. -(Treatment:Year))
summary(dm10)
anova(dm10)

anova(dm8,dm10)

# - Treatment

dm11<-update(dm10,.~. -(Treatment))
summary(dm11)
anova(dm11)

anova(dm10,dm11) # chi = 5.017      df = 2    P = 0.08139 .

# - Year
#again, dropping from the model that does not include the Treatment:Year interaction

dm12<-update(dm10,.~. -(Year))
summary(dm12)
anova(dm12)

anova(dm10,dm12)


#want to see how avg.inhab may have differed among treatments

anc1<-ggplot(behave, aes(avg.inhab, total.dist.moved, shape=Treatment.o, linetype=Treatment.o, col=Treatment.o)) +
  geom_smooth(method="lm", se=FALSE, show.legend = TRUE)  +
  geom_point(size=2)+
  theme_classic()+
  labs(x=(bquote('Gobies '~Reef^ -1*'')),y=(bquote('Total dist. moved '~Reef^ -1~ Week^-1*'')))+
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
  scale_x_continuous(breaks=c(10,12,14,16,18,20)) + scale_y_continuous(limits = c(0,500))+
  labs(color  = "Perceived Risk", linetype = "Perceived Risk", shape = "Perceived Risk")
anc1

# 1c. foraging rate ####
fr<-lmer(bites.min~Treatment*Year*avg.inhab+(1|Trial) + (1|Treatment:Trial) +
           (1|Trial:avg.inhab)+(1|avg.inhab:Treatment:Trial), REML=F, behave)
hist(resid(fr))
qqnorm(resid(fr))
qqline(resid(fr))
plot(fr)
boxplot(bites.min~Treatment,data=behave) #slightly higher variance in Low treatment
summary(fr) #seems to be a bit more variation by trial
anova(fr)

rand(fr) #trial explains some of the variance, will keep it in, all other random
# effects removed

fr.1<-lmer(bites.min~Treatment*Year*avg.inhab+(1|Trial), REML=F, data=behave)
rand(fr.1)
anova(fr.1)

#removing higher-level interactions

fr.2<-lmer(bites.min~Treatment*Year+avg.inhab+(1|Trial), REML=F, data=behave)
summary(fr.2)
anova(fr.2)

#seems that the variance explained by trial goes away when I take away the higher-level
# interactions between fixed effects, going to remove it

fr.3<-lm(bites.min~Treatment*Year+avg.inhab, data=behave)
anova(fr.3)

#removing avg.inhab
fr.4<-lm(bites.min~Treatment*Year, data=behave)
anova(fr.4)

#removing three-way interaction of random effect
fr2<-update(fr, .~. -(1|avg.inhab:Treatment:Trial))
summary(fr2)
anova(fr2)

anova(fr,fr2) #no difference in models when three-way interaction with random intercept is removed

#next logical removal is the 1|treatment:trial term

fr3<-update(fr2, .~. -(1|Treatment:Trial)) # error message, model failed to converge
summary(fr3)
anova(fr3)

anova(fr2,fr3) #no difference

#next logical removal is the 1|Trial:avg.inhab term

fr4<-update(fr3, .~. -(1|Trial:avg.inhab))
summary(fr4)
anova(fr4)

anova(fr3,fr4)

#trial effect
fr5<-update(fr2,.~. -(1|Trial))
summary(fr5)
anova(fr5)

anova(fr2,fr5)

#testing effects of fixed factors, going to start with highest-order interactions with the covariate

fr6<-update(fr4, .~. -(Treatment:Year:avg.inhab))
summary(fr6) 
anova(fr6)

anova(fr4,fr6) 

#now removing Year.fact:avg.inhab (the order of removal seems arbitrary, but that's the next one)
fr7<-update(fr6, .~. -(Year:avg.inhab))
summary(fr7) 
anova(fr7)

anova(fr7,fr6)

#now want to take out Treatment:avg.inhab

fr8<-update(fr7, .~. -(Treatment:avg.inhab)) 
summary(fr8) 
anova(fr8) 

anova(fr7,fr8)

#me8 is the final model that I will likely end up with (those are all of the fixed and random factors that I'm interested in)
hist(resid(fr8))#not bad
#LS-adjusted means from model
emmeans(fr8, pairwise~Treatment) #warning message re: interactions, but I think it's okay
boxplot(bites.min~Treatment,data=behave)# variances don't look too bad, and medians look pretty good to me (i.e. match up with LS-means for the most part)
library(visreg)
visreg(fr8)
visreg(fr8, "Treatment",by="Year")

#dropping the fixed effects that I have reduced the model to do the log-likelihood estimates
#NOTE: m8 is the fully-reduced model, so just have to jeep iterating that model +/- individual fixed factors of interest

# - avg.inhab, but keeping in other fixed factors

fr9<-update(fr8, .~. -(avg.inhab))
summary(fr9) 
anova(fr9)

anova(fr8,fr9)

# - Treatment:Year

fr10<-update(fr8, .~. -(Treatment:Year))
summary(fr10)
anova(fr10)

anova(fr8,fr10)

# - Treatment

fr11<-update(fr10,.~. -(Treatment))
summary(fr11)
anova(fr11)

anova(fr10,fr11) # chi = 5.017      df = 2    P = 0.08139 .

# - Year
#again, dropping from the model that does not include the Treatment:Year interaction

fr12<-update(fr10,.~. -(Year))
summary(fr12)
anova(fr12)

anova(fr10,fr12)


# 1d. interactins with conspecifics ####
ci<-lmer(courtship.min~Treatment*Year*avg.inhab+(1|Trial) + (1|Treatment:Trial) +
           (1|Trial:avg.inhab)+(1|avg.inhab:Treatment:Trial), REML=F, behave)
hist(resid(ci))
qqnorm(resid(ci))
qqline(resid(ci))
plot(ci)
boxplot(courtship.min~Treatment,data=behave) #slightly higher variance in Low treatment
summary(ci) #seems to be a bit more variation by trial
anova(ci)

rand(ci) #will keep trial, explains some of the variance

ci.1<-lmer(courtship.min~Treatment*Year*avg.inhab+(1|Trial), REML=F, behave)
anova(ci.1)

#removing NS interactions
ci.2<-lmer(courtship.min~Treatment*Year+avg.inhab+(1|Trial), REML=F, behave)
anova(ci.2)

# removing covariate
ci.3<-lmer(courtship.min~Treatment*Year+(1|Trial), REML=F, behave)
anova(ci.3)

#removing three-way interaction of random effect
ci2<-update(ci, .~. -(1|avg.inhab:Treatment:Trial))
summary(ci2)
anova(ci2)

anova(ci,ci2) #no difference in models when three-way interaction with random intercept is removed

#next logical removal is the 1|treatment:trial term

ci3<-update(ci2, .~. -(1|Treatment:Trial)) # error message, model failed to converge
summary(ci3)
anova(ci3)

anova(ci2,ci3) #no difference

#next logical removal is the 1|Trial:avg.inhab term

ci4<-update(ci3, .~. -(1|Trial:avg.inhab))
summary(ci4)
anova(ci4)

anova(ci3,ci4)

#trial effect
ci5<-update(ci2,.~. -(1|Trial))
summary(ci5)
anova(ci5)

anova(ci2,ci5)

#testing effects of fixed factors, going to start with highest-order interactions with the covariate

ci6<-update(ci4, .~. -(Treatment:Year:avg.inhab))
summary(ci6) 
anova(ci6)

anova(ci4,ci6) 

#now removing Year.fact:avg.inhab (the order of removal seems arbitrary, but that's the next one)
ci7<-update(ci6, .~. -(Year:avg.inhab))
summary(ci7) 
anova(ci7)

anova(ci7,ci6)

#now want to take out Treatment:avg.inhab

ci8<-update(ci7, .~. -(Treatment:avg.inhab)) 
summary(ci8) 
anova(ci8) 

anova(ci7,ci8)

#me8 is the final model that I will likely end up with (those are all of the fixed and random factors that I'm interested in)
hist(resid(ci8))#not bad
#LS-adjusted means ciom model
emmeans(ci8, pairwise~Treatment) #warning message re: interactions, but I think it's okay
boxplot(courtship.min~Treatment,data=behave)# variances don't look too bad, and medians look pretty good to me (i.e. match up with LS-means for the most part)
library(visreg)
visreg(ci8)
visreg(ci8, "Treatment",by="Year")

#dropping the fixed effects that I have reduced the model to do the log-likelihood estimates
#NOTE: m8 is the fully-reduced model, so just have to jeep iterating that model +/- individual fixed factors of interest

# - avg.inhab, but keeping in other fixed factors

ci9<-update(ci8, .~. -(avg.inhab))
summary(ci9) 
anova(ci9)

anova(ci8,ci9)

# - Treatment:Year

ci10<-update(ci8, .~. -(Treatment:Year))
summary(ci10)
anova(ci10)

anova(ci8,ci10)

# - Treatment

ci11<-update(ci10,.~. -(Treatment))
summary(ci11)
anova(ci11)

anova(ci10,ci11) 

# - Year
#again, dropping ciom the model that does not include the Treatment:Year interaction

ci12<-update(ci10,.~. -(Year))
summary(ci12)
anova(ci12)

anova(ci10,ci12)

# 1e. movement rate ####
mr<-lmer(movements.min~Treatment*Year*avg.inhab+(1|Trial) + (1|Treatment:Trial) +
           (1|Trial:avg.inhab)+(1|avg.inhab:Treatment:Trial), REML=F, behave)
hist(resid(mr))
qqnorm(resid(mr))
qqline(resid(mr))
plot(mr)
boxplot(movements.min~Treatment,data=behave) #slightly higher variance in Low treatment
summary(mr) #seems to be a bit more variation by trial
anova(mr)

rand(mr) #going to remove trial, seems like it doesn't explain much variation
# P = 0.94

mr.1<-lm(movements.min~Treatment*Year*avg.inhab, data = behave)
anova(mr.1)

# removing NS interactions
mr.2<-lm(movements.min~Treatment*Year+avg.inhab, data = behave)
anova(mr.2)

#removing covariate
mr.3<-lm(movements.min~Treatment*Year, data = behave)
anova(mr.3)

#removing three-way interaction of random effect
mr2<-update(mr, .~. -(1|avg.inhab:Treatment:Trial))
summary(mr2)
anova(mr2)

anova(mr,mr2) #no difference in models when three-way interaction with random intercept is removed

#next logical removal is the 1|treatment:trial term

mr3<-update(mr2, .~. -(1|Treatment:Trial)) # error message, model failed to converge
summary(mr3)
anova(mr3)

anova(mr2,mr3) #no difference

#next logical removal is the 1|Trial:avg.inhab term

mr4<-update(mr3, .~. -(1|Trial:avg.inhab))
summary(mr4)
anova(mr4)

anova(mr3,mr4)

#trial effect
mr5<-update(mr2,.~. -(1|Trial))
summary(mr5)
anova(mr5)

anova(mr2,mr5)

#testing effects of fixed factors, going to start with highest-order interactions with the covariate

mr6<-update(mr4, .~. -(Treatment:Year:avg.inhab))
summary(mr6) 
anova(mr6)

anova(mr4,mr6) 

#now removing Year.fact:avg.inhab (the order of removal seems arbitrary, but that's the next one)
mr7<-update(mr6, .~. -(Year:avg.inhab))
summary(mr7) 
anova(mr7)

anova(mr7,mr6)

#now want to take out Treatment:avg.inhab

mr8<-update(mr7, .~. -(Treatment:avg.inhab)) 
summary(mr8) 
anova(mr8) 

anova(mr7,mr8)

#me8 is the final model that I will likely end up with (those are all of the fixed and random factors that I'm interested in)
hist(resid(mr8))#not bad
#LS-adjusted means mrom model
emmeans(mr8, pairwise~Treatment) #warning message re: interactions, but I think it's okay
boxplot(movements.min~Treatment,data=behave)# variances don't look too bad, and medians look pretty good to me (i.e. match up with LS-means for the most part)
library(visreg)
visreg(mr8)
visreg(mr8, "Treatment",by="Year")

#dropping the fixed effects that I have reduced the model to do the log-likelihood estimates
#NOTE: m8 is the fully-reduced model, so just have to jeep iterating that model +/- individual fixed factors of interest

# - avg.inhab, but keeping in other fixed factors

mr9<-update(mr8, .~. -(avg.inhab))
summary(mr9) 
anova(mr9)

anova(mr8,mr9)

# - Treatment:Year

mr10<-update(mr8, .~. -(Treatment:Year))
summary(mr10)
anova(mr10)

anova(mr8,mr10)

# - Treatment

mr11<-update(mr10,.~. -(Treatment))
summary(mr11)
anova(mr11)

anova(mr10,mr11) 

# - Year
#again, dropping mrom the model that does not include the Treatment:Year interaction

mr12<-update(mr10,.~. -(Year))
summary(mr12)
anova(mr12)

anova(mr10,mr12)


# 2. plotting ####
# proportion of time exposed ####
exp<-with(behave, aggregate((proportion.exposed), list(Treatment=Treatment), mean))
exp$se<-with(behave, aggregate((proportion.exposed), list(Treatment=Treatment), 
                               function(x) sd(x)/sqrt(length(x))))[,2]

exp.plot<- ggplot(exp, aes(x=Treatment, y=x, fill=Treatment)) +
  geom_bar(stat="identity", colour= "black", width = 0.5, position="dodge")+ 
  scale_x_discrete(limits=c("Low","Medium","High"))+
  theme_classic() + 
  labs(x="Risk Treatment", y="Proportion of Time Exposed") +
  theme(legend.position="none") + 
  scale_fill_manual(values=c("grey", "grey", "grey")) +
  theme(axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"), 
        axis.title=element_text(size=20))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), 
        axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0)), 
        axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(legend.text=element_text(size=18)) +
  theme(legend.title =element_text(size=20))+
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0,0),limits = c(0,0.807))
exp.plot + geom_linerange(aes(ymin=x-se, ymax=x+se), size=0.5,   
                           position=position_dodge(.85)) + theme(text = element_text(family="Arial"))

# linear distance traveled####
td<-with(behave, aggregate((total.dist.moved), list(Treatment=Treatment), mean))
td$se<-with(behave, aggregate((total.dist.moved), list(Treatment=Treatment), 
                              function(x) sd(x)/sqrt(length(x))))[,2]

td.plot<- ggplot(td, aes(x=Treatment, y=x, fill=Treatment)) +
  geom_bar(stat="identity", colour= "black", width = 0.5, position="dodge")+ 
  scale_x_discrete(limits=c("Low","Medium","High"))+
  theme_classic() + 
  labs(x="Risk Treatment", y="Linear Distance Traveled (mm)") +
  theme(legend.position="none") + 
  scale_fill_manual(values=c("grey", "grey", "grey")) +
  theme(axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"), 
        axis.title=element_text(size=20))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), 
        axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0)), 
        axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(legend.text=element_text(size=18)) +
  theme(legend.title =element_text(size=20))+
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0,0),limits = c(0,205))
td.plot + geom_linerange(aes(ymin=x-se, ymax=x+se), size=0.5,   
                           position=position_dodge(.85)) + theme(text = element_text(family="Arial"))


# foraging rate (bites per minute)####
fr<-with(behave, aggregate((bites.min), list(Treatment=Treatment), mean))
fr$se<-with(behave, aggregate((bites.min), list(Treatment=Treatment), 
                              function(x) sd(x)/sqrt(length(x))))[,2]

fr.plot<- ggplot(fr, aes(x=Treatment, y=x, fill=Treatment)) +
  geom_bar(stat="identity", colour= "black", width = 0.5, position="dodge")+ 
  scale_x_discrete(limits=c("Low","Medium","High"))+
  theme_classic() + 
  labs(x="Risk Treatment",y=(expression(atop("Foraging", 
                                             paste((bites~min^-1))))))+
  theme(legend.position="none") + 
  scale_fill_manual(values=c("grey", "grey", "grey")) +
  theme(axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"), 
        axis.title=element_text(size=20))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), 
        axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0)), 
        axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(legend.text=element_text(size=18)) +
  theme(legend.title =element_text(size=20))+
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0,0),limits = c(0,1.01))
fr.plot + geom_linerange(aes(ymin=x-se, ymax=x+se), size=0.5,   
                           position=position_dodge(.85)) + theme(text = element_text(family="Arial"))

# interactions with conspecifics (displays per minute)####
cr<-with(behave, aggregate((courtship.min), list(Treatment=Treatment), mean))
cr$se<-with(behave, aggregate((courtship.min), list(Treatment=Treatment), 
                              function(x) sd(x)/sqrt(length(x))))[,2]

cr.plot<- ggplot(cr, aes(x=Treatment, y=x, fill=Treatment)) +
  geom_bar(stat="identity", colour= "black", width = 0.5, position="dodge")+ 
  scale_x_discrete(limits=c("Low","Medium","High"))+
  theme_classic() + 
  labs(x="Risk Treatment", y = "Courtship Displays per Minute")+
  theme(legend.position="none") + 
  scale_fill_manual(values=c("grey", "grey", "grey")) +
  theme(axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"), 
        axis.title=element_text(size=20))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), 
        axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0)), 
        axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(legend.text=element_text(size=18)) +
  theme(legend.title =element_text(size=20))+
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0,0),limits = c(0,0.151))
cr.plot + geom_linerange(aes(ymin=x-se, ymax=x+se), size=0.5,   
                           position=position_dodge(.85)) + theme(text = element_text(family="Arial"))

#movements per minute with rate in parentheses

png("Output/2019.11.14.courtship.rate.with.parentheses.9.5x5.5.300dpi.png", width = 6.5, height = 5.5, units = 'in', res = 300)
png("Output/2019.2.1.interactions.with.parentheses.9.5x5.5.300dpi.png", width = 6.5, height = 5.5, units = 'in', res = 300)

reco.plot<- ggplot(cr, aes(x=Treatment, y=x, fill=Treatment)) +
  geom_bar(stat="identity", colour= "black", width = 0.5, position="dodge")+ 
  scale_x_discrete(limits=c("Low","Medium","High"))+
  theme_classic() + 
  labs(x="Risk Treatment",y=(expression(atop("Interactions with Conspecifics", 
                                             paste((displays~min^-1))))))+
  theme(legend.position="none") + 
  scale_fill_manual(values=c("grey", "grey", "grey")) +
  theme(axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"), 
        axis.title=element_text(size=20))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), 
        axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0)), 
        axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(legend.text=element_text(size=18)) +
  theme(legend.title =element_text(size=20))+
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0,0),limits = c(0,0.151))
reco.plot + geom_linerange(aes(ymin=x-se, ymax=x+se), size=0.5,   
                           position=position_dodge(.85)) + theme(text = element_text(family="Arial"))

dev.off()

# movement rate (movements per minute)####
mm<-with(behave, aggregate((movements.min), list(Treatment=Treatment), mean))
mm$se<-with(behave, aggregate((movements.min), list(Treatment=Treatment), 
                              function(x) sd(x)/sqrt(length(x))))[,2]

mm.plot<- ggplot(mm, aes(x=Treatment, y=x, fill=Treatment)) +
  geom_bar(stat="identity", colour= "black", width = 0.5, position="dodge")+ 
  scale_x_discrete(limits=c("Low","Medium","High"))+
  theme_classic() +
  labs(x="Risk Treatment",y=(expression(atop("Movements min"^-1))))+
  theme(legend.position="none") + 
  scale_fill_manual(values=c("grey", "grey", "grey")) +
  theme(axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"), 
        axis.title=element_text(size=20))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), 
        axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0)), 
        axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(legend.text=element_text(size=18)) +
  theme(legend.title =element_text(size=20))+
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0,0),limits = c(0,1.51),
                                                             labels = scales::number_format(accuracy = 0.01)) #changed to 2 decimal places for movement rate)
mm.plot + geom_linerange(aes(ymin=x-se, ymax=x+se), size=0.5,   
                           position=position_dodge(.85)) + theme(text = element_text(family="Arial"))

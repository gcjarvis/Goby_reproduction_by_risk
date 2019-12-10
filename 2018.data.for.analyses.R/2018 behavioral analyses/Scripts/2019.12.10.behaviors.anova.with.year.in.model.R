# Description: reanalyzing behaviors with new model 
# Author: George C Jarvis
# Date: Tue Dec 10 20:18:19 2019
# Notes: Going to use the same model as before, because there's reason to think
# that the number of inhabitants would affect behaviors:
# -----behavior~treatment*avg.inhab*Year.factor+(1|Year:Trial)---
# Because of that, I will drop nonsignificant interactions with the covariate
# NOTE: these are technically mixed-model ANCOVAS, b/c of covariate 
# - and nested random factor of trial
# --------------

#clear workspace
rm(list=ls())

#loading packages####
library(sciplot)
library(lme4)
library(lmerTest)
library(dplyr)
library(ggplot2)
library(plyr)
#packages needed for MANCOVA
library(car)#type III SS
#library(psych)#descriptive statistics
#library(effects)#adjusted means
#easiest mancova package
#library(jmv)
library(agricolae)#can't use standard tukey test for mixed models
library(multcomp)# post-hoc for fixed effects in mixed model

#loading df####
behave<-read.csv("Data/2019.10.25.behavior.includes.recollections.csv")
behave<-na.omit(behave)

#adding a column for year (as a proxy for tagging procedure), where trials 1-3 = 2017, and 4-6 = 2018
# NOTE: there were no behavioral observations for trial 6
behave$Year <- ifelse(behave$Trial <=3, 2017, 2018)
#Going to make a variable with it as a factor, though I don't think it matters a whole lot
behave$Year.fact<- as.factor(behave$Year)

#making the variable "avg.inhab" ((20+reco)/2), rounded to the nearest whole fish
#using the average number of inhabitants per reef as the covariate
behave$avg.inhab<-(ceiling((behave$Recollection+20)/2))

#modeling for each behavior####

#proportion of time exposed####
pe1<-lmer(proportion.exposed~Treatment*avg.inhab*Year.fact+(1|Year.fact:Trial), data=behave)
hist(resid(pe1))#looks like a slightly better model in terms of assumptions
qqnorm(resid(pe1))
qqline(resid(pe1))
anova(pe1)#type III anova = no treatment effect, no interactive effects
Anova(pe1)#Analysis of deviance, type II wald test = sig. trt. effect,
#                                                     no density, no interactives
summary(pe1)

#reducing the model, removing nonsignificant interactions with covariate
pe1.1<-lmer(proportion.exposed~(Treatment*avg.inhab)+(Treatment*Year.fact)+(avg.inhab*Year.fact)+
            Treatment+avg.inhab+Year.fact+(1|Year.fact:Trial), data=behave)
hist(resid(pe1.1))#looks like a slightly better model in terms of assumptions
qqnorm(resid(pe1.1))
qqline(resid(pe1.1))
anova(pe1.1)#type III anova = no treatment effect, no interactive effect
Anova(pe1.1)#Analysis of deviance, type II wald test = sig. trt. effect,
#                                                     no density, no interactive
summary(pe1.1)

#reducing further, removing the NS interactions with covariate
pe1.2<-lmer(proportion.exposed~Treatment+(Treatment*Year.fact)+
            avg.inhab+Year.fact+(1|Year.fact:Trial), data=behave)
hist(resid(pe1.2))#looks like a slightly better model in terms of assumptions
qqnorm(resid(pe1.2))
qqline(resid(pe1.2))
anova(pe1.2)#type III anova = sig. treatment effect with ANOVA (III)
Anova(pe1.2)#Analysis of deviance, type II wald test = sig. trt. effect,
#                                                     no density, no interactive
summary(pe1.2)

#posthoc tests, I'm not sure this works with covariate
summary(glht(pe1, linfct=mcp(Treatment="Tukey")))
summary(glht(pe1.1, linfct=mcp(Treatment="Tukey")))
summary(glht(pe1.2, linfct=mcp(Treatment="Tukey")))

#plotting
#ordering "Treatment" and "T.6.comparison"
behave$Treatment.ord<-ordered(behave$Treatment, levels=c("Low","Medium","High"))

#no grouping factor
bargraph.CI(x.factor = Treatment.ord, response = proportion.exposed, 
            legend=TRUE, main="proportion of time exposed", 
            data = behave)

#grouped by trial
bargraph.CI(x.factor = Treatment.ord, response = proportion.exposed, 
            group= Trial, legend=TRUE, main="proportion of time exposed, grouped by trial", 
            data = behave)

#grouped by year
bargraph.CI(x.factor = Treatment.ord, response = proportion.exposed, 
            group= Year, legend=TRUE, main="proportion of time exposed, grouped by year", 
            data = behave)

#movements/minute (movement rate)####
mm1<-lmer(movements.min~Treatment*avg.inhab*Year.fact+(1|Year.fact:Trial), data=behave)
hist(resid(mm1))#looks like a slightly better model in terms of assumptions
qqnorm(resid(mm1))
qqline(resid(mm1))
anova(mm1)#type III anova = no treatment effect, no interactive effects
Anova(mm1)#Analysis of deviance, type II wald test = NS trt. effect,
#                                                     no density, no interactives
summary(mm1)

#reducing the model, removing nonsignificant interactions with covariate
mm1.1<-lmer(movements.min~(Treatment*avg.inhab)+(Treatment*Year.fact)+(avg.inhab*Year.fact)+
              Treatment+avg.inhab+Year.fact+(1|Year.fact:Trial), data=behave)
hist(resid(mm1.1))#looks like a slightly better model in terms of assumptions
qqnorm(resid(mm1.1))
qqline(resid(mm1.1))
anova(mm1.1)#type III anova = no treatment effect, no interactive effect
Anova(mm1.1)#Analysis of deviance, type II wald test = NS trt. effect,
#                                                     no density, no interactive
summary(mm1.1)

#reducing further, removing the NS interactions with covariate
mm1.2<-lmer(movements.min~Treatment+(Treatment*Year.fact)+
              avg.inhab+Year.fact+(1|Year.fact:Trial), data=behave)
hist(resid(mm1.2))#looks like a slightly better model in terms of assumptions
qqnorm(resid(mm1.2))
qqline(resid(mm1.2))
anova(mm1.2)#type III anova = sig. treatment effect with ANOVA (III)
Anova(mm1.2)#Analysis of deviance, type II wald test = sig. trt. effect,
#                                                     no density, no interactive
summary(mm1.2)

#plotting
#ordering "Treatment" and "T.6.comparison"
behave$Treatment.ord<-ordered(behave$Treatment, levels=c("Low","Medium","High"))

#no grouping factor
bargraph.CI(x.factor = Treatment.ord, response = movements.min, 
            legend=TRUE, main="movements per minute", 
            data = behave)

#grouped by trial
bargraph.CI(x.factor = Treatment.ord, response = movements.min, 
            group= Trial, legend=TRUE, main="movements.min, grouped by trial", 
            data = behave)

#grouped by year
bargraph.CI(x.factor = Treatment.ord, response = proportion.exposed, 
            group= Year, legend=TRUE, main="movements.min, grouped by year", 
            data = behave)

#bites/min (feeding rate)####
bm1<-lmer(bites.min~Treatment*avg.inhab*Year.fact+(1|Year.fact:Trial), data=behave)
hist(resid(bm1))
qqnorm(resid(bm1))
qqline(resid(bm1))
anova(bm1)#type III anova = sig. average.inhab*year interaction, all else NS
Anova(bm1)#Analysis of deviance, type II wald test = 
#                   sig. average.inhab*year interaction, all else NS
summary(bm1)

#reducing the model, removing the NS three-way interaction
#unlike the two previous models, I'm going to keep in the avg.inhab*year term
br1.1<-lmer(bites.min~(Treatment*avg.inhab)+(Treatment*Year.fact)+(avg.inhab*Year.fact)+
              Treatment+avg.inhab+Year.fact+(1|Year.fact:Trial), data=behave)
hist(resid(br1.1))
qqnorm(resid(br1.1))
qqline(resid(br1.1))
anova(br1.1)#type III anova = sig. avg.inhab*year effect, also a nearly sig. 
#                             year effect (p = 0.05548), still not sig.,
#                             but closer than Chi square result for year
Anova(br1.1)#Analysis of deviance, type II wald test = sig. avg.inhab*year,
#                           NS year effect (p = 0.88739), no others sig.
#                           
# NOTE: this is the first time that ANOVA and Chi square have been notably different in result
# - I feel a bit better about the chi square result (see plot for bites per minute, by year)
# - It doesn't seem like there were differences by year, so I think R is calculating F and p values incorrectly
# - with the mixed model
summary(br1.1)

#reducing further, removing the NS interactions with covariate (-treatment*avg.inhab)
br1.2<-lmer(bites.min~Treatment+(Treatment*Year.fact)+(avg.inhab*Year.fact)+
              avg.inhab+Year.fact+(1|Year.fact:Trial), data=behave)
hist(resid(br1.2))
qqnorm(resid(br1.2))
qqline(resid(br1.2))
anova(br1.2)#type III anova = sig. treatment effect with ANOVA (III)
Anova(br1.2)#Analysis of deviance, type II wald test = sig. trt. effect,
#                                                     no density, no interactive
summary(br1.2)

#plotting
#ordering "Treatment" and "T.6.comparison"
behave$Treatment.ord<-ordered(behave$Treatment, levels=c("Low","Medium","High"))

#no grouping factor
bargraph.CI(x.factor = Treatment.ord, response = bites.min, 
            legend=TRUE, main="bites per minute", 
            data = behave)

#just looking at bite rates by year, doesn't seem correct? I think R 
# is calculating F and p values incorrectly because of the mixed model
bargraph.CI(x.factor = Year.fact, response = bites.min, 
            legend=TRUE, main="bites per minute, by year only", 
            data = behave)

#grouped by trial
bargraph.CI(x.factor = Treatment.ord, response = bites.min, 
            group= Trial, legend=TRUE, main="bites per min, grouped by trial", 
            data = behave)

#grouped by year
bargraph.CI(x.factor = Treatment.ord, response = bites.min, 
            group= Year, legend=TRUE, main="bites per minute, grouped by year", 
            data = behave)

#looking at average inhab*year result
bargraph.CI(x.factor = avg.inhab, response = bites.min, 
            group= Year.fact, legend=TRUE, main="bites per min, by inhab and year", 
            xlab= "avg.inhab",data = behave)

#total distance moved####
dm1<-lmer(total.dist.moved~Treatment*avg.inhab*Year.fact+(1|Year.fact:Trial), data=behave)
hist(resid(dm1))#looks like a slightly better model in terms of assumptions
qqnorm(resid(dm1))
qqline(resid(dm1))
anova(dm1)#type III anova = sig. avg.inhab, all else NS, no 3-way int.
Anova(dm1)#Chi square (same as above) = sig. avg inhab, no 3-way int (p = 0.084)
summary(dm1)

#reducing the model, removing nonsignificant interactions with covariate
dm1.1<-lmer(total.dist.moved~(Treatment*avg.inhab)+(Treatment*Year.fact)+(avg.inhab*Year.fact)+
              Treatment+avg.inhab+Year.fact+(1|Year.fact:Trial), data=behave)
hist(resid(dm1.1))#looks like a slightly better model in terms of assumptions
qqnorm(resid(dm1.1))
qqline(resid(dm1.1))
anova(dm1.1)#type III anova = sig. avg.inhab, all else NS
Anova(dm1.1)#Chi square (same as above) = sig. avg inhab, all else NS
summary(dm1.1)

#reducing further, removing the NS interactions with covariate
dm1.2<-lmer(total.dist.moved~Treatment+(Treatment*Year.fact)+
              avg.inhab+Year.fact+(1|Year.fact:Trial), data=behave)
hist(resid(dm1.2))#looks like a slightly better model in terms of assumptions
qqnorm(resid(dm1.2))
qqline(resid(dm1.2))
anova(dm1.2)#type III anova = sig. inhab effect, nearly a treatment effect (p = 0.058)
Anova(dm1.2)#Chi square (same as above) = sig. avg inhab, no sig. trt. effect (p = 0.083)
# I suppose I'd be more comfortable with the Chi square result, based on the previous models
summary(dm1.2)

#plotting
#ordering "Treatment" and "T.6.comparison"
behave$Treatment.ord<-ordered(behave$Treatment, levels=c("Low","Medium","High"))

#no grouping factor
bargraph.CI(x.factor = Treatment.ord, response = total.dist.moved, 
            legend=TRUE, main="total distance moved", 
            data = behave)

#grouped by trial
bargraph.CI(x.factor = Treatment.ord, response = total.dist.moved, 
            group= Trial, legend=TRUE, main="total.dist.moved, grouped by trial", 
            data = behave)

#grouped by year
bargraph.CI(x.factor = Treatment.ord, response = total.dist.moved, 
            group= Year, legend=TRUE, main="total.dist.moved, grouped by year", 
            data = behave)

#looking at average inhab*year result
bargraph.CI(x.factor = avg.inhab, response = total.dist.moved, 
            group= Year.fact, legend=TRUE, main="distance moved, by inhab and year", 
            xlab= "avg.inhab",data = behave)

#courtship rate####
cr1<-lmer(courtship.min~Treatment*avg.inhab*Year.fact+(1|Year.fact:Trial), data=behave)
hist(resid(cr1))#looks like a slightly better model in terms of assumptions
qqnorm(resid(cr1))
qqline(resid(cr1))
anova(cr1)#type III anova = nothing sig. 
Anova(cr1)#Chi square (same as above) = nothing sig.
summary(cr1)

#reducing the model, removing nonsignificant interactions with covariate
cr1.1<-lmer(courtship.min~(Treatment*avg.inhab)+(Treatment*Year.fact)+(avg.inhab*Year.fact)+
              Treatment+avg.inhab+Year.fact+(1|Year.fact:Trial), data=behave)
hist(resid(cr1.1))#looks like a slightly better model in terms of assumptions
qqnorm(resid(cr1.1))
qqline(resid(cr1.1))
anova(cr1.1)#type III anova = nothing sig.
Anova(cr1.1)#Chi square (same as above) = nothing sig.
summary(cr1.1)

#reducing further, removing the NS interactions with covariate
cr1.2<-lmer(courtship.min~Treatment+(Treatment*Year.fact)+
              avg.inhab+Year.fact+(1|Year.fact:Trial), data=behave)
hist(resid(cr1.2))#looks like a slightly better model in terms of assumptions
qqnorm(resid(cr1.2))
qqline(resid(cr1.2))
anova(cr1.2)#type III anova = sig. inhab effect, nearly a treatment effect (p = 0.058)
Anova(cr1.2)#Chi square (same as above) = sig. avg inhab, no sig. trt. effect (p = 0.083)
# I suppose I'd be more comfortable with the Chi square result, based on the previous models
summary(cr1.2)

#plotting
#ordering "Treatment" and "T.6.comparison"
behave$Treatment.ord<-ordered(behave$Treatment, levels=c("Low","Medium","High"))

#no grouping factor
bargraph.CI(x.factor = Treatment.ord, response = courtship.min, 
            legend=TRUE, main="courtship by treatment", 
            data = behave)

#grouped by trial
bargraph.CI(x.factor = Treatment.ord, response = total.dist.moved, 
            group= Trial, legend=TRUE, main="total.dist.moved, grouped by trial", 
            data = behave)

#grouped by year
bargraph.CI(x.factor = Treatment.ord, response = total.dist.moved, 
            group= Year, legend=TRUE, main="total.dist.moved, grouped by year", 
            data = behave)

#looking at average inhab*year result
bargraph.CI(x.factor = avg.inhab, response = total.dist.moved, 
            group= Year.fact, legend=TRUE, main="distance moved, by inhab and year", 
            xlab= "avg.inhab",data = behave)


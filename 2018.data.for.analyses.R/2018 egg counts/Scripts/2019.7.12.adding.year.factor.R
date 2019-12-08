# Description: changing model to reflect differences in tagging procedures (year added as fixed factor)
# Author: George C Jarvis
# Date: Sat Dec 07 10:24:59 2019
# Notes: I'm going to add year as a fixed factor to the model and code the full model as
#       egg count~Treatment*Year*avg.inhab+(1|Trial:Year)
#       I have to figure out if the model is running it correctly (i.e. all of the interactions, and also the error df)
#
# 2019.8.12 update: go to "reproduction per week" section
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

#modeling####
# including year, and trial nested within year, year as numeric
#mod1<-lmer(Egg.count~Treatment*avg.inhab*Year*(1|Year:Trial), data=repro)
#hist(resid(mod1))
#qqnorm(resid(mod1))
#qqline(resid(mod1))
#anova(mod1, type = "III")
#Anova(mod1)
#summary(mod1) #none of it is sig. can reduce model? also want to make sure I get the same result when year is run as factor

# including year, and trial nested within year, year as factor --> this is the right way to do it, b/c year = tagging method
#mod1.1<-lmer(Egg.count~Treatment*avg.inhab*Year.fact*(1|Year.fact:Trial), data=repro)
#hist(resid(mod1.1))
#qqnorm(resid(mod1.1))
#qqline(resid(mod1.1))
#anova(mod1.1, type = "III")
#Anova(mod1.1)
#summary(mod1.1) #see sig. effect of avg. inhab. (makes sense), year (makes sense), and avg.inhab*year (makes sense)

# SEE PLOT FOR MOD1.1 BELOW TO SEE HOW AVG.INHAB*YEAR AFFECTED REPRODUCTIVE OUTPUT (not super surprising results)

#the interaction makes sense - the number of eggs laid by number of inhabitants changed, 
# because the durations of the experiments were much longer. All this interaction is telling us is that the tagging methods
# allowed for greater retention of gobies over time - which is useful information for someone trying to do this exp. again

#this is the full model --> (Egg.count~Treatment*avg.inhab*Year.fact+(1|Year.fact:Trial)

#now running without the interaction for the random effect to see how the results change
#I think that this is the best model, it tests the fixed effect of year, plus the random variation in intercept among
# trials within year. I don't think you build the models with interactions with the random variable (i.e. * vs. +)
mod1.1.1<-lmer(Egg.count~Treatment*avg.inhab*Year.fact+(1|Year.fact:Trial), data=repro)
hist(resid(mod1.1.1))
qqnorm(resid(mod1.1.1))
qqline(resid(mod1.1.1))
anova(mod1.1.1)
Anova(mod1.1.1)
summary(mod1.1.1) #see sig. effect of avg. inhab. (makes sense), year (makes sense), and avg.inhab*year (makes sense)

#QUESTION: go with chi-squared test or with ANOVA test? Seems to be that I can do either? Doesn't change the result, just the stats

#reduced model 1 --> removed three-way interaction, nonsignificant

#running reduced model taking out the three-way interaction between treatment*avg.inhab*year
mod1.2<-lmer(Egg.count~(Treatment*avg.inhab)+(Treatment*Year.fact)+(avg.inhab*Year.fact)+
               Treatment+avg.inhab+Year.fact+(1|Year.fact:Trial), data=repro)
hist(resid(mod1.2))
qqnorm(resid(mod1.2))
qqline(resid(mod1.2))
anova(mod1.2)
Anova(mod1.2)
summary(mod1.2) #no sig. interaction between treatment*year (2 categorical factors), will remove them as per Mark's suggestion

#further reduced model 1 --> removed nonsignificant interaction between treatment and year

#running a further reduced model without treatment*year interaction, seeing how it affects results
mod1.3<-lmer(Egg.count~(Treatment*avg.inhab)+(avg.inhab*Year.fact)+Treatment+
               avg.inhab+Year.fact+(1|Year.fact:Trial), data=repro)
hist(resid(mod1.3))
qqnorm(resid(mod1.3))
qqline(resid(mod1.3))
anova(mod1.3)
Anova(mod1.3)
summary(mod1.3) # same results, will check out a model comparison?

#comparing all possible models
anova(mod1,mod1.1,mod1.1.1,mod1.2,mod1.3) #seems like model 1.3 (most reduced model) is best? Will ask M. Steele

#comparing models that I sent to M. Steele for review (mod1.1.1, 1.2, and 1.3)
anova(mod1.1.1,mod1.2,mod1.3) #seems like model 1.3 (most reduced model) is best? Will ask M. Steele


#plotting####

#ordering "Treatment" and "T.6.comparison"
repro$Treatment.ord<-ordered(repro$Treatment, levels=c("Low","Medium","High"))
repro$T6.comparison.ord<-ordered(repro$T6.comparison, levels=c("Low","Medium","High","Uncaged"))

#grouped by trial
bargraph.CI(x.factor = Treatment.ord, response = Egg.count, 
            group= Year.fact, legend=TRUE, main="Reproduction by risk and year", 
            data = repro, ylab="egg count")
bargraph.CI(x.factor = Treatment.ord, response = Egg.count, 
            group= Trial, legend=TRUE, main="all trials, HR combined grouped by trial", 
            data = repro)
bargraph.CI(x.factor = avg.inhab, response = Egg.count, 
            group= Year.fact, legend=TRUE, main="reproduction by number of gobies", 
            xlab="avg.inhab", ylab="egg count",
            data = repro)
lineplot.CI(avg.inhab,Egg.count,group=Year.fact,legend = TRUE,main="reproduction by number of gobies", 
            xlab="avg.inhab", ylab="egg count", data=repro) 
# avg.inhab*year stats show that there was a greater effect of avg.inhab on reproductive
#       output in 2018 than in 2017. In fact, it doesn't seem like there was much of a difference in output by avg.inhab,
#       which suggests that those effects manifest over time (i.e. can't see them in short (weeklong) trials)

#reproduction by week (2019.8.12)####
#converting total output to output per week to better reflect the differences
# in the data based on the duration of the trial

#CALCULATING OUTPUT/WEEK: take the output from my raw data and dividing it 
# by the number of weeks that the trial lasted:
#-trials 1-3 = 1 week
#-trials 4-5 = 4 weeks
#-trial 6 = 2 weeks (I emailed M.steele about this on 2019.8.12 to see what he though
# about the difference in trial duration within year, but hopefully this is okay)

#going to try and do this in dplyr with the conditional mutation 
# NOTE: I rounded the egg counts to the nearest whole number, b/c non-while eggs didn't 
# seem like a great variable

repro<-repro %>%
  mutate(egg.week = ifelse(Trial<4, Egg.count/1,
                           ifelse(Trial == 4| Trial == 5, (ceiling(Egg.count/4)),
                                  ifelse(Trial == 6, (ceiling(Egg.count/2)), NA))))
#View(repro)
#base code for future reference, NOTE, you have to name the df, or else the new variable
# won't show up in the 
#df1<-df %>%
#  mutate(g = ifelse(a == 2 | a == 5 | a == 7 | (a == 1 & b == 4), 2,
#                    ifelse(a == 0 | a == 1 | a == 4 | a == 3 |  c == 4, 3, NA)))

#now want to set up the new model, including new response ("egg.year") and year

#full model first, jus to check and see if stats hold up

mod2<-lmer(egg.week~Treatment*avg.inhab*Year.fact+(1|Year.fact:Trial), data=repro)
hist(resid(mod2))
qqnorm(resid(mod2))
qqline(resid(mod2))
anova(mod2)
Anova(mod2)
summary(mod2) #only see sig. effect of avg.inhab, going to run reduced models

#notes for model reduction:
# 1. remove all nonsignificant (p>0.05) interactions with covariate
#   -unless there is an interaction with a higher order interaction involving it.
#   -In the case of mod2, there aren't any higher-order interactions
# 2. keep all categorical factors and their interactions (trt., year)

#going to run reduced model, minus three-way interaction
mod2.1<-lmer(egg.week~(Treatment*avg.inhab)+(avg.inhab*Year.fact)+Treatment+
               avg.inhab+Year.fact+(1|Year.fact:Trial), data=repro)
hist(resid(mod2.1))
qqnorm(resid(mod2.1))
qqline(resid(mod2.1))
anova(mod2.1)
Anova(mod2.1)
summary(mod2.1) # only factor that is sig. is the avg.inhab, will 
# reduce model further to reflect fixed effects without any interactions,
# will also include the random effect of trial nested within year

#maybe running this reduced model
mod2.2<-lmer(egg.week~Treatment+avg.inhab+Year.fact+(1|Year.fact:Trial),
             data=repro)
hist(resid(mod2.2))
qqnorm(resid(mod2.2))
qqline(resid(mod2.2))
anova(mod2.2)
Anova(mod2.2)
summary(mod2.2) # same results, will check out a model comparison?

#comparing new models
anova(mod2,mod2.1,mod2.2) #mod 2.2 seems to be the best mdoel in terms of AIC

#it's interesting that there was no effect of year on reproduction per week
# going to look at that now

#plotting with eggs per week####

#ordering "Treatment" and "T.6.comparison"
repro$Treatment.ord<-ordered(repro$Treatment, levels=c("Low","Medium","High"))
repro$T6.comparison.ord<-ordered(repro$T6.comparison, levels=c("Low","Medium","High","Uncaged"))

bargraph.CI(x.factor = Treatment.ord, response = egg.week, 
            group= Year.fact, legend=TRUE, main="Reproduction per week between years", 
            data = repro, ylab="egg count per reef per week")
bargraph.CI(x.factor = Treatment.ord, response = egg.week, 
            group= Trial, legend=TRUE, main="all trials, HR combined grouped by trial", 
            data = repro)
bargraph.CI(x.factor = avg.inhab, response = egg.week, 
            group= Year.fact, legend=TRUE, main="weekly reproduction by number of gobies", 
            xlab="avg.inhab", ylab="egg count per reef per week",
            data = repro)
lineplot.CI(avg.inhab,egg.week,group=Year.fact,legend = TRUE,main="reproduction by number of gobies", 
            xlab="avg.inhab", ylab="egg count", data=repro) 
#NOTE: re: lineplot, I went back and checked the raw data and saw that there was one 
# reef (reef 7 in trial 3) where I recollected 11 gobies (avg.inhab=15.5, which rounds up to 16)
# but there were no eggs laid. That's why there's a 0 value at 16 for 2017

View(repro)

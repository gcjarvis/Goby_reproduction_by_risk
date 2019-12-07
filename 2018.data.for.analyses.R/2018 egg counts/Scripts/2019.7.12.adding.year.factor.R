# Description: changing model to reflect differences in tagging procedures (year added as fixed factor)
# Author: George C Jarvis
# Date: Sat Dec 07 10:24:59 2019
# Notes: I'm going to add year as a fixed factor to the model and code the full model as
#       egg count~Treatment*Year*avg.inhab+(1|Trial:Year)
#       I have to figure out if the model is running it correctly (i.e. all of the interactions, and also the error df)
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

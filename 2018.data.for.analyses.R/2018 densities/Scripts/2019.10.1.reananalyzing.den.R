## Description: Reanalyzing visual counts for MS
# Author: George C Jarvis
# Date: Mon Sep 30 08:48:33 2019
# Notes: After emailing with M. Steele, he thought it was best
#         to reanalyze my data with mixed models
#         with two key differences: 1) including "trial" as a random effect
#                                   2) coding T6 uncaged the same as high-risk
#       Q: Use "den.max", or "Density"? Will do both, I'm leaning towards den.max
# Overall findings: A.) With all trials pooled, accounting for trial as a random
#     factor in a mixed-effect model, and using a Poisson distribution (count data),
#     there were more fish seen on the high-risk/uncaged reefs than low and med
#     treatments. I think it makes sense to do den.max because it's more representative.
#     Not sure why there were more seen in high-risk, other than maybe they needed to be more
#     vigilant?
# --------------

#loading packages####
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
library(agricolae)#for tukey post-hoc test

#importing dataset, adding number of gobies on each reef, ordering treatments####
#includes a cloumn ("Treatment") where uncaged and HR are coded as "High"
#also includes a column ("T6.comparison") where uncaged and high are separated
viz.surv<-read.csv("Data/density.2019.10.1.csv")
viz.surv$Density<-as.numeric(viz.surv$Density)
viz.surv$den.max<-as.numeric(viz.surv$den.max)

#log-transforming den.max and density for assumptions of normality
viz.surv$log.d<-log(viz.surv$Density + 1)
viz.surv$log.dm<-log(viz.surv$den.max + 1)

#sqrt-transforming den.max and density for assumptions of normality
viz.surv$sqrt.d<-sqrt(viz.surv$Density)
viz.surv$sqrt.dm<-sqrt(viz.surv$den.max)

#ordering "Treatment" and "T.6.comparison"
viz.surv$Treatment<-ordered(viz.surv$Treatment, levels=c("Low","Medium","High"))
viz.surv$T6.comparison<-ordered(viz.surv$T6.comparison, levels=c("Low","Medium","High","Uncaged"))

#subsetting T6 only for comparison
viz.surv.T6<-viz.surv[c(742:821),c(1:15)]

#A. running mixed models with combined HR factor and trial as both a fixed and random factor####
# 1)it does seem like a Poisson distribution is a better fit than the normal dist.
#   when using the untransformed data

#think I should go with this model, but will see what M. Steele has to say
mod.1<-glmer(den.max ~ Treatment + (1|Trial),family=poisson, data=viz.surv)
hist(resid(mod.1))
qqnorm(resid(mod.1))
qqline(resid(mod.1))
anova(mod.1)
Anova(mod.1) 
summary(mod.1)
#chi-squared shows an effect of treatment, showing more fish seen in HR/uncaged treament
fixed.effects(mod.1)
ranef(mod.1)

mod.1a<-glmer(Density ~ Treatment + (1|Trial),family = poisson, data=viz.surv)
hist(resid(mod.1a))
qqnorm(resid(mod.1a))
qqline(resid(mod.1a))
anova(mod.1a)
Anova(mod.1a)
#chi-squared shows p=0.07 for treatment

# looking at log-transformed data, no need to glm
mod.1.log<-lmer(log.dm ~ Treatment + (1|Trial), data=viz.surv)
hist(resid(mod.1.log))
qqnorm(resid(mod.1.log))
qqline(resid(mod.1.log))#looks skewed
anova(mod.1.log)
Anova(mod.1.log)
#no sig. difference among treatments

mod.1a.log<-lmer(log.d ~ Treatment + (1|Trial), data=viz.surv)
hist(resid(mod.1a.log))
qqnorm(resid(mod.1a.log))
qqline(resid(mod.1a.log))
anova(mod.1a.log)
Anova(mod.1a.log)
#no sig. difference among treatments

#going to look at sqrt transformation
mod.1.sqrt<-lmer(sqrt.dm ~ Treatment + (1|Trial), data=viz.surv)
hist(resid(mod.1.sqrt))
qqnorm(resid(mod.1.sqrt))
qqline(resid(mod.1.sqrt))#looks like the best one in terms of normality
anova(mod.1.sqrt)
Anova(mod.1.sqrt)
#no sig. difference among treatments

mod.1a.sqrt<-lmer(sqrt.d ~ Treatment + (1|Trial), data=viz.surv)
hist(resid(mod.1a.sqrt))
qqnorm(resid(mod.1a.sqrt))
qqline(resid(mod.1a.sqrt))
anova(mod.1a.sqrt)
Anova(mod.1a.sqrt)
#no sig. difference among treatments

#plots

#den.max
lineplot.CI(Day,den.max,group=Treatment,legend = TRUE,
            main="Visual surveys, combined HR and Uncaged", 
            xlab="Day", ylab="May number fo fish seen", 
            data=viz.surv)

#square-root transformed data for den max
lineplot.CI(Day,sqrt.dm,group=Treatment,legend = TRUE,
            main="Visual surveys, combined HR and Uncaged", 
            xlab="Day", ylab="May number fo fish seen (sqrt)", 
            data=viz.surv)

#wanting to see how it looks as a bargraph by day
bargraph.CI(x.factor = Day, response = den.max, 
            group= Treatment, legend=TRUE, main="all trials, HR combined grouped by trial", 
            data = viz.surv)

#no grouping factor, seems like H>L=M, strange...with den.max
#goes against what Mark found...
bargraph.CI(x.factor = Treatment, response = den.max, 
            legend=TRUE, main="den.max, HR combined grouped by trial", 
            data = viz.surv)

#with Density
bargraph.CI(x.factor = Treatment, response = Density, 
            legend=TRUE, main="Density, HR combined grouped by trial", 
            data = viz.surv)

#going to run the same model again but with distiction between
# HR and uncaged treatments
#results: U=H>L=M
#seems like uncaged is driving trend upward when combined
# with HR, but in general, doesn't bring the average up that much,
# likely because of low replication
bargraph.CI(x.factor = T6.comparison, response = den.max, 
            legend=TRUE, main="den.max, HR separated caged.uncaged", 
            data = viz.surv)

#now just with density
bargraph.CI(x.factor = T6.comparison, response = Density, 
            legend=TRUE, main="density, HR separated caged.uncaged", 
            data = viz.surv)
#seems like the average values are slightly lower when we use den instead of den.max
# that makes sense

#no grouping factor
bargraph.CI(x.factor = Treatment, response = Density, 
            legend=TRUE, main="all trials, HR combined w uncaged, no trial", 
            data = viz.surv)

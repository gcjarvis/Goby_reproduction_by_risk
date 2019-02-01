##########################################
# 11/05/2018                             #  
# egg counts for WSN                     #
# George Jarvis                          #
##########################################

rm(list=ls())

library(sciplot)
library(lme4)
library(car)
library(lmerTest)

gob.sub<-read.csv("Data/10.29.18.eggs.csv")
#gob.sub$week<-as.factor(gob.sub$week)#need to make week a factor???
View(gob.sub)
head(gob.sub)
gob.sub$week<-as.factor(gob.sub$week)
gob.sub$Treatment<-ordered(gob.sub$Treatment,levels=c("Low","Medium","High"))
gob.sub$Treatment<-as.factor(gob.sub$Treatment)
#gob.egg.back$treatment
bargraph.CI(x.factor = treatment, response = egg.count, main="back eggs", data = gob.egg.back)


#subsetting the data by reef deployment order
#reefs 1-9
gob.egg.front<-gob.sub[gob.sub$reef<10, c("week", "deployment.day", "reef", "treatment", "egg.count", "initial.biomass","final.biomass","change.in.biomass","average.biomass","eggs.initial.biomass","eggs.avg.biomass","nest.count","eggs.per.nest","avg.density.week","average.density.samp.days","avg.egg.density.week","avg.egg.samp.days","round")]
#reefs 10-18
gob.egg.back<-gob.sub[gob.sub$reef>9, c("week", "deployment.day", "reef", "treatment", "egg.count", "initial.biomass","final.biomass","change.in.biomass","average.biomass","eggs.initial.biomass","eggs.avg.biomass","nest.count","eggs.per.nest","avg.density.week","average.density.samp.days","avg.egg.density.week","avg.egg.samp.days","round")]
#wanted to just look at week 1
gob.egg.wk.1<-gob.sub[gob.sub$week==1, c("week", "deployment.day", "reef", "treatment", "egg.count", "initial.biomass","final.biomass","change.in.biomass","average.biomass","eggs.initial.biomass","eggs.avg.biomass","nest.count","eggs.per.nest","avg.density.week","average.density.samp.days","avg.egg.density.week","avg.egg.samp.days","round")]


#model building
#lm fixed factors only
mod1.egg<-lm(egg.count~treatment,data=gob.sub)
hist(resid(mod1.egg))
qqnorm(resid(mod1.egg))
qqline(resid(mod1.egg))

mod2.egg<-lm(egg.count~treatment*week, data = gob.sub)
hist(resid(mod1.egg))
plot(mod2.egg)

mod3.egg<-lm(egg.count~treatment*week*deployment.day, data=gob.sub)
hist(resid(mod3.egg))
anova(mod3.egg)
summary(mod3.egg)
#mixed models
#lmer

mod4.egg<-lmer(egg.count~treatment*week+(1|deployment.day), data = gob.sub)
hist(resid(mod4.egg))
qqnorm(resid(mod4.egg))
qqline(resid(mod4.egg))
anova(mod4.egg)
summary(mod4.egg)

mod4b.egg<-lm(egg.count~treatment*week, data = gob.egg.front)
anova(mod4b.egg)
summary(mod4b.egg)

mod6.egg<-lmer(egg.count~treatment*deployment.day+(1|week), data = gob.sub)
hist(resid(mod6.egg))
anova(mod6.egg)
summary(mod6.egg)
Anova(mod6.egg,type="II")
#not really what I'm interested in. I want to see how reproduction changes over time

mod4a.egg<-lmer(egg.count~treatment+(1|week)+(1|deployment.day), data = gob.sub)#covariate of final biomass
mod4b.egg<-lmer(egg.count~treatment*week+(1|deployment.day), data = gob.sub)
mod4c.egg<-lmer(egg.count~treatment+week+(1|deployment.day), data = gob.sub)
mod4d.egg<-lmer(egg.count~week+(1|deployment.day), data = gob.sub)

#transforming data
#log + 1 transformation

#adding a column to the data to log-transform +1 (for nests that had 0 eggs) the egg counts
gob.sub$log.egg<-(log(gob.sub$egg.count + 1))
View(gob.sub)

#model building with log-transformed data

log.mod1<-lm(log.egg~treatment, data=gob.sub)
hist(resid(log.mod1))#not so good
log.mod2<-lm(log.egg~treatment+week, data=gob.sub)
hist(resid(log.mod2))
log.mod3<-lm(log.egg~treatment*week, data=gob.sub)
hist(resid(log.mod3))
log.mod4<-lmer(log.egg~treatment+week+(1|deployment.day),data=gob.sub)
hist(resid(log.mod4))
log.mod5<-lmer(log.egg~treatment*week+(1|deployment.day), data=gob.sub)
hist(resid(log.mod5))
log.mod.6.poiss<-glmer(log.egg~treatment*week+(1|deployment.day),family=poisson, data=gob.sub)
hist(resid(log.mod.6.poiss))

#square root + 0.5 transformation

gob.sub$sqrt<-(sqrt(gob.sub$egg.count + 0.5))
View(gob.sub)

log.sqrt1<-lm(sqrt~treatment, data=gob.sub)
hist(resid(log.sqrt1))#not so good
log.sqrt2<-lm(sqrt~treatment+week, data=gob.sub)
hist(resid(log.sqrt2))
log.sqrt3<-lm(sqrt~treatment*week, data=gob.sub)
hist(resid(log.sqrt3))
log.sqrt4<-lmer(sqrt~treatment+week+(1|deployment.day),data=gob.sub)
hist(resid(log.sqrt4))
log.sqrt5<-lmer(sqrt~treatment*week+(1|deployment.day), data=gob.sub)
hist(resid(log.sqrt5))
log.sqrt6.poiss<-glmer(sqrt~treatment*week+(1|deployment.day), family=poisson, data=gob.sub)
hist(resid(log.sqrt6.poiss))

###########now want to run analyses accounting for density######

#1. want to look at same model while accounting for density
#2. want to look at models with different response variable that has density built into it

bargraph.CI(x.factor = week, response = avg.egg.density.week, group=treatment, Legend=TRUE, main="eggs per capita weekly density", data = gob.sub)
bargraph.CI(x.factor = week, response = avg.egg.samp.days, group=treatment, Legend=TRUE, main="eggs per capita sampling density", data = gob.sub)
#interesting, let's run some stats on it

#1. model building with total reproductive output and density factored into the model
mod1.egg<-lm(egg.count~treatment*week*average.density.samp.days,data=gob.sub)
hist(resid(mod1.egg))
qqnorm(resid(mod1.egg))
qqline(resid(mod1.egg))
anova(mod1.egg)
summary(mod1.egg)

mod2.egg<-lm(log.egg~treatment*week*average.density.samp.days, data=gob.sub)
hist(resid(mod2.egg))
qqnorm(resid(mod2.egg))
qqline(resid(mod2.egg))
anova(mod2.egg)
summary(mod2.egg)

#lineplot.CI(x.factor = average.density.samp.days, response = egg.count, group=treatment, Legend=TRUE, main="eggs per capita weekly density", data = gob.sub)

#2. 
#model building with different response variables accounting for density
#lm fixed factors only
mod1.egg<-lmer(avg.egg.samp.days~treatment*week+(1|deployment.day),data=gob.sub)
hist(resid(mod1.egg))
qqnorm(resid(mod1.egg))
qqline(resid(mod1.egg))
anova(mod1.egg)
summary(mod1.egg)

mod3.egg<-lm(avg.egg.samp.days~treatment*week*deployment.day, data=gob.sub)
hist(resid(mod3.egg))
anova(mod3.egg)
summary(mod3.egg)
#mixed models
#lmer

mod4.egg<-lmer(egg.count~treatment*week+(1|deployment.day), data = gob.sub)
hist(resid(mod4.egg))
qqnorm(resid(mod4.egg))
qqline(resid(mod4.egg))
anova(mod4.egg)
summary(mod4.egg)

mod4b.egg<-lm(egg.count~treatment*week, data = gob.egg.front)
anova(mod4b.egg)
summary(mod4b.egg)

mod6.egg<-lmer(egg.count~treatment*deployment.day+(1|week), data = gob.sub)
hist(resid(mod6.egg))
anova(mod6.egg)
summary(mod6.egg)
Anova(mod6.egg,type="II")
#not really what I'm interested in. I want to see how reproduction changes over time

mod4a.egg<-lm(egg.count~Treatment*week+avg.density.week, data = gob.sub)#covariate of final biomass
hist(resid(mod4a.egg))
qqnorm(resid(mod4a.egg))
qqline(resid(mod4a.egg))
anova(mod4a.egg)
summary(mod4a.egg)
AIC(mod4a.egg)


#this is the analysis I went with for WSN
mod4b.egg<-lmer(egg.count~Treatment*week+(1|avg.density.week), data = gob.sub)#covariate of final biomass
hist(resid(mod4b.egg))
qqnorm(resid(mod4b.egg))
qqline(resid(mod4b.egg))
anova(mod4b.egg)
summary(mod4b.egg)
AIC(mod4b.egg)
ranef(mod4b.egg)

mod4b.egg<-lmer(egg.count~treatment*week+(1|deployment.day), data = gob.sub)
mod4c.egg<-lmer(egg.count~treatment+week+(1|deployment.day), data = gob.sub)
mod4d.egg<-lmer(egg.count~week+(1|deployment.day), data = gob.sub)

#rounding the data out to a nearest whole number to do a poiss dist.
gob.sub$round<-round(gob.sub$avg.egg.samp.days)
mod1.egg<-glmer(round~treatment*week+(1|deployment.day),family=poisson,data=gob.sub)
hist(resid(mod1.egg))
qqnorm(resid(mod1.egg))
qqline(resid(mod1.egg))
Anova(mod1.egg,type="III")

#####looking at subsetted data from front and back reefs#####

#front
mod1.egg<-lm(avg.egg.samp.days~treatment*week,data=gob.egg.front)
hist(resid(mod1.egg))
qqnorm(resid(mod1.egg))
qqline(resid(mod1.egg))
anova(mod1.egg)
summary(mod1.egg)

bargraph.CI(x.factor = week, response = avg.egg.density.week, group=treatment, Legend=TRUE, main="eggs per capita weekly density", data = gob.egg.front)
bargraph.CI(x.factor = week, response = avg.egg.samp.days, group=treatment, Legend=TRUE, main="eggs per capita sampling density", data = gob.egg.front)


#back
mod1.egg<-lm(avg.egg.samp.days~treatment*week,data=gob.egg.back)
hist(resid(mod1.egg))
qqnorm(resid(mod1.egg))
qqline(resid(mod1.egg))
anova(mod1.egg)
summary(mod1.egg)

#bargraph.CI(x.factor = week, response = eggs.prop.population.remaining, group=treatment, Legend=TRUE, main="eggs per capita weekly density", data = gob.sub)
bargraph.CI(x.factor = week, response = avg.egg.samp.days, group=treatment, Legend=TRUE, main="eggs per capita sampling density", data = gob.egg.back)

mod3.egg<-lm(avg.egg.samp.days~treatment*week*deployment.day, data=gob.sub)
hist(resid(mod3.egg))
anova(mod3.egg)
summary(mod3.egg)

###log-transforming the adjusted count data to account for density
gob.sub$log.egg.den<-log(gob.sub$avg.egg.samp.days+1)
gob.sub$sqrt<-sqrt(gob.sub$avg.egg.samp.days+0.5)

mod1.egg<-lmer(sqrt~treatment*week+(1|deployment.day),data=gob.sub)
hist(resid(mod1.egg))
qqnorm(resid(mod1.egg))
qqline(resid(mod1.egg))
anova(mod1.egg)
summary(mod1.egg)

bargraph.CI(x.factor = week, response = sqrt, group=treatment, Legend=TRUE, main="eggs per capita weekly density", data = gob.sub)
bargraph.CI(x.factor = treatment, response = sqrt, main="eggs per capita sampling density", data = gob.sub)

#####figuring out the best model to meet the assumptions, while accounting for density

mod1.egg<-lm(avg.egg.samp.days~treatment*week, data=gob.sub)#pretty bad
hist(resid(mod1.egg))
qqnorm(resid(mod1.egg))
qqline(resid(mod1.egg))

mod1.egg.dep<-lm(avg.egg.samp.days~treatment*week*deployment.day, data=gob.sub)#worse
hist(resid(mod1.egg.dep))
qqnorm(resid(mod1.egg.dep))
qqline(resid(mod1.egg.dep))

mod1.egg.dep.rand<-lmer(avg.egg.samp.days~treatment*week+(1|deployment.day),data=gob.sub)
hist(resid(mod1.egg.dep.rand))
qqnorm(resid(mod1.egg.dep.rand))
qqline(resid(mod1.egg.dep.rand))
anova(mod1.egg.dep.rand)
summary(mod1.egg.dep.rand)
#slightly better, but still not great, shows that the effect of deployment day doesn't
#really account for a whole lot of the total variance

anova(mod1.egg,mod1.egg.dep)
AIC(mod1.egg)
AIC(mod1.egg.dep)
AIC(mod1.egg.dep.rand)

#seems like this model (mod1.egg.dep.rand) is the best one

# I think that because I am normalizing by the average density on  the sampling days,
# including deployment day doesn't have as big of an impact (they're correlated)

#let's rerun accouting for density
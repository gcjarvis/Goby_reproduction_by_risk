##########################################
# 10/30/2018                             #  
# subsetting egg counts for WSN analyses #
# George Jarvis                          #
##########################################

#want to break down egg counts into two groups: R1-9 and R10-18

#clear workspace
rm(list=ls())

library(sciplot)
library(lme4)
library(car)

#quick density analyses, with just density data from 6/16,6/17, and 6/20
getwd()
gob.sub<-read.csv("Data/10.29.18.eggs.csv")
gob.sub$week<-as.factor(gob.sub$week)#need to make week a factor???
View(gob.sub)
head(gob.sub)
gob.sub<-read.csv("Data/10.29.18.eggs.outliers.a.csv")

#subsetting the data by reef deployment order
#reefs 1-9
gob.egg.front<-gob.sub[gob.sub$reef<10, c("week", "deployment.day", "reef", "treatment", "egg.count", "initial.biomass","final.biomass","change.in.biomass","average.biomass","eggs.initial.biomass","eggs.avg.biomass","nest.count","eggs.per.nest")]
#reefs 10-18
gob.egg.back<-gob.sub[gob.sub$reef>9, c("week", "deployment.day", "reef", "treatment", "egg.count", "initial.biomass","final.biomass","change.in.biomass","average.biomass","eggs.initial.biomass","eggs.avg.biomass","nest.count","eggs.per.nest")]
#wanted to just look at week 1
gob.egg.wk.1<-gob.sub[gob.sub$week==1, c("week", "deployment.day", "reef", "treatment", "egg.count", "initial.biomass","final.biomass","change.in.biomass","average.biomass","eggs.initial.biomass","eggs.avg.biomass","nest.count","eggs.per.nest","avg.egg.density.week")]
gob.egg.wk.2<-gob.sub[gob.sub$week==2, c("week", "deployment.day", "reef", "treatment", "egg.count", "initial.biomass","final.biomass","change.in.biomass","average.biomass","eggs.initial.biomass","eggs.avg.biomass","nest.count","eggs.per.nest","avg.egg.density.week")]
gob.egg.wk.3<-gob.sub[gob.sub$week==3, c("week", "deployment.day", "reef", "treatment", "egg.count", "initial.biomass","final.biomass","change.in.biomass","average.biomass","eggs.initial.biomass","eggs.avg.biomass","nest.count","eggs.per.nest","avg.egg.density.week")]
gob.egg.wk.4<-gob.sub[gob.sub$week==4, c("week", "deployment.day", "reef", "treatment", "egg.count", "initial.biomass","final.biomass","change.in.biomass","average.biomass","eggs.initial.biomass","eggs.avg.biomass","nest.count","eggs.per.nest","avg.egg.density.week")]




#visualization

#FRONT REEFS

#front reefs only by total reproduction (i.e., not including week)
bargraph.CI(x.factor = treatment, response = egg.count, main="front eggs", data = gob.egg.front)
#including week as a grouping factor
bargraph.CI(x.factor = week, response = egg.count, group = treatment, legend = TRUE, main="front eggs by week", xlab="week",ylab="average reproduction per reef", data = gob.egg.front)
bargraph.CI(x.factor = treatment, response = egg.count, group = week, legend = TRUE, main="front eggs by week", xlab="week",ylab="average reproduction per reef", data = gob.egg.front)


#BACK REEFS

#back reefs only by total reproduction (i.e., not including week)
bargraph.CI(x.factor = treatment, response = egg.count, main="back eggs", data = gob.egg.back)
#including week as a grouping factor
bargraph.CI(x.factor = week, response = egg.count, group = treatment, legend = TRUE, main="back eggs by week", xlab="week",ylab="average reproduction per reef", data= gob.egg.back)
bargraph.CI(x.factor = treatment, response = egg.count, group = week, legend = TRUE, main="back eggs by week", xlab="week",ylab="average reproduction per reef", data = gob.egg.back)

#total egg production per week
bargraph.CI(x.factor = treatment, response = egg.count, main="week 1 eggs", xlab="treatment",ylab="average reproduction per reef", data = gob.egg.wk.1)
bargraph.CI(x.factor = treatment, response = egg.count, main="week 2 eggs", xlab="treatment",ylab="average reproduction per reef", data = gob.egg.wk.2)
bargraph.CI(x.factor = treatment, response = egg.count, main="week 3 eggs", xlab="treatment",ylab="average reproduction per reef", data = gob.egg.wk.3)
bargraph.CI(x.factor = treatment, response = egg.count, main="week 4 eggs", xlab="treatment",ylab="average reproduction per reef", data = gob.egg.wk.4)

#average egg production per week
bargraph.CI(x.factor = treatment, response = avg.egg.density.week, main="week 1 eggs", xlab="treatment",ylab="average reproduction per reef", data = gob.egg.wk.1)
bargraph.CI(x.factor = treatment, response = avg.egg.density.week, main="week 2 eggs", xlab="treatment",ylab="average reproduction per reef", data = gob.egg.wk.2)
bargraph.CI(x.factor = treatment, response = avg.egg.density.week, main="week 3 eggs", xlab="treatment",ylab="average reproduction per reef", data = gob.egg.wk.3)
bargraph.CI(x.factor = treatment, response = avg.egg.density.week, main="week 4 eggs", xlab="treatment",ylab="average reproduction per reef", data = gob.egg.wk.4)


#stats for week 1 #only sig. diff is the poiss.dist. model, and that's not a great model fit...
mod.1a.wk.1.egg<-lm(egg.count~treatment, data = gob.egg.wk.1)
hist(resid(mod.1a.wk.1.egg))
qqnorm(resid(mod.1a.wk.1.egg))
qqline(resid(mod.1a.wk.1.egg))
anova(mod.1a.wk.1.egg)
summary(mod.1a.wk.1.egg)

mod.wk1.egg<-lmer(egg.count~treatment+(1|deployment.day), data = gob.egg.wk.1)#random effect of deployment day
hist(resid(mod.wk1.egg))
qqnorm(resid(mod.wk1.egg))
qqline(resid(mod.wk1.egg))
anova(mod.wk1.egg)
summary(mod.wk1.egg)

mod.avg.wk.1<-lmer(avg.egg.density.week~treatment+(1|deployment.day),data=gob.egg.wk.1)
hist(resid(mod.avg.wk.1))
anova(mod.avg.wk.1)
summary(mod.avg.wk.1)

mod.wk1.egg.poiss<-glmer(egg.count~treatment+(1|deployment.day),family=poisson, data = gob.egg.wk.1)

anova(mod.1a.wk.1.egg)
anova(mod.wk1.egg)
summary(mod.wk1.egg)

Anova(mod.wk1.egg,type="II")
Anova(mod.wk1.egg.poiss,type="II")

hist(resid(mod.1a.wk.1.egg))
hist(resid(mod.wk1.egg))
hist(resid(mod.wk1.egg.poiss))

AIC(mod.1a.wk.1.egg)
AIC(mod.wk1.egg)
AIC(mod.wk1.egg.poiss)

#THOUGHT: PLOT DEN. BY WEEK
#might consider grouping by week for analyses, because that is what I present in the repro data

###model selection for best fit

#pooled data models
tot.den<-lm(density~treatment, data = gob.den)
hist(resid(tot.den))
qqnorm(resid(tot.den))
qqline(resid(tot.den))
boxplot(density~treatment, data = gob.den)
#front reefs  only
tot.den.f<-lm(density~treatment, data = gob.den.front)
hist(resid(tot.den.f))
#back reefs only
tot.den.b<-lm(density~treatment, data = gob.den.back)
hist(resid(tot.den.b))

#model selection with all data
#fixed factors only
mod1.egg<-lm(egg.count~treatment+week, data = gob.sub)
mod1a.egg<-lm(egg.count~treatment+week+deployment.day, data = gob.sub)
mod1b.egg.fb<-lm(egg.count~treatment+week+deployment.day+final.biomass, data = gob.sub)

#fixed factors with interactions
mod2.egg.int<-lm(egg.count~treatment*week, data = gob.sub)

mod13.egg<-lm(egg.count~treatment+week*deployment.day, data = gob.sub)
mod6.egg.fb<-lm(egg.count~treatment+week+final.biomass, data = gob.sub)


mod3.egg<-lmer(egg.count~treatment+(1|week), data = gob.sub)
mod4.egg<-lmer(egg.count~treatment+(1|week)+(1|deployment.day), data = gob.sub)
mod4a.egg<-lmer(egg.count~treatment+final.biomass+(1|week)+(1|deployment.day), data = gob.sub)#covariate of final biomass
mod4b.egg<-lmer(egg.count~treatment*week+(1|deployment.day), data = gob.sub)
mod4c.egg<-lmer(egg.count~treatment+week+(1|deployment.day), data = gob.sub)
mod4d.egg<-lmer(egg.count~week+(1|deployment.day), data = gob.sub)

mod5.egg.poiss<-glmer(egg.count~treatment+(1|week)+(1|deployment.day),family=poisson, data = gob.sub)
mod5a.egg.poiss<-glmer(egg.count~treatment+(1|week)+(1|deployment.day)+final.biomass,family=poisson, data = gob.sub)
#mod5a.egg.poiss<-glmer(egg.count~treatment+week+(1|deployment.day),family=gamma, data = gob.sub)

mod7.egg.int.fb<-lm(egg.count~treatment*week+final.biomass, data = gob.sub)
mod8.egg.3w.int.fb<-lm(egg.count~treatment*week*final.biomass, data = gob.sub)
mod9.egg.week.deploy.date<-lmer(egg.count~treatment+week+(1|deployment.day), data = gob.sub)
mod10.egg.3w.int.fb<-lmer(egg.count~treatment*week*final.biomass+(1|deployment.day), data = gob.sub)
mod11.egg.3w.fix.int.fb<-lmer(egg.count~treatment*week+final.biomass+(1|deployment.day), data = gob.sub)
mod12.egg.3w.fix.fb<-lmer(egg.count~treatment+week+final.biomass+(1|deployment.day), data = gob.sub)


anova(mod13.egg)
AIC(mod13.egg)
Anova(mod9.egg.week.deploy.date, type="II")

#really interested in the question of how reproduction might change over time
#so I want to look at week as a fixed factor

anova(mod1.egg)
anova(mod2.egg.int)

#model selection
#lm
anova(mod1.egg,mod2.egg.int,mod6.egg.fb)
AIC(mod1.egg)
AIC(mod2.egg.int)
AIC(mod6.egg.fb)
AIC(mod7.egg.int.fb)
AIC(mod8.egg.3w.int.fb)
hist(resid(mod8.egg.3w.int.fb))
AIC(mod10.egg.3w.int.fb)
hist(resid(mod10.egg.3w.int.fb))#seems like this one is best now that I changed week to a factor
AIC(mod11.egg.3w.fix.int.fb)
AIC(mod12.egg.3w.fix.fb)

#lmer
anova(mod3.egg,mod4.egg, mod4a.egg, mod10.egg.3w.int.fb, mod11.egg.3w.fix.int.fb,mod12.egg.3w.fix.fb)

hist(resid(mod4.egg))
qqnorm(resid(mod4.egg))
qqline(resid(mod4.egg))
boxplot(egg.count~treatment, data=gob.sub)
plot(mod4.egg)

hist(resid(mod4a.egg))

AIC(mod4a.egg)

Anova(mod4.egg, type = "II")

hist(resid(mod12.egg.3w.fix.fb))

#glmer
AIC(mod5.egg.poiss)
hist(resid(mod5.egg.poiss))
#not good..

##
#OK! I think my model should be model 11 
#it includes week as a fixed factor, and will show me if there is an interaction between egg counts and time
mod10.egg.3w.int.fb<-lmer(egg.count~treatment*week*final.biomass+(1|deployment.day), data = gob.sub)
mod11.egg.3w.fix.int.fb<-lmer(egg.count~treatment*week+final.biomass+(1|deployment.day), data = gob.sub)
mod12.egg.3w.fix.fb<-lmer(egg.count~treatment+week+final.biomass+(1|deployment.day), data = gob.sub)

anova(mod10.egg.3w.int.fb,mod11.egg.3w.fix.int.fb,mod12.egg.3w.fix.fb)
summary(mod10.egg.3w.int.fb)
anova(mod10.egg.3w.int.fb)
Anova(mod10.egg.3w.int.fb, type="II")

anova(mod11.egg.3w.fix.int.fb)
Anova(mod11.egg.3w.fix.int.fb, type="II")

Anova(mod12.egg.3w.fix.fb, type="II")

################################these are the results that I used in my presentation############

##mod 4b is the model that I want
#final biomass was not a good predictor becasue it was not applicable to early egg counts
hist(resid(mod4b.egg))
qqnorm(resid(mod4b.egg))
qqline(resid(mod4b.egg))

Anova(mod4b.egg, type="II")

anova(mod4b.egg,mod4c.egg)

hist(resid(mod4c.egg))
qqnorm(resid(mod4c.egg))
qqline(resid(mod4c.egg))
Anova(mod4c.egg,type="II")
plot(mod4c.egg)
boxplot(egg.count~treatment*week, data=gob.sub)

Anova(mod4d.egg, type="II")

############################# Number of nests per treatment per week#######

#all reefs combined
bargraph.CI(x.factor = treatment, response = nest.count, group = week, legend = TRUE, main="nest production per week", xlab="week",ylab="average number of eggs produced", data = gob.sub)
#just front eggs
bargraph.CI(x.factor = treatment, response = nest.count, group = week, legend = TRUE, main="front nest production per week", xlab="week",ylab="average number of eggs produced", data = gob.egg.front)
#just back eggs
bargraph.CI(x.factor = treatment, response = nest.count, group = week, legend = TRUE, main="back nest production per week", xlab="week",ylab="average number of eggs produced", data = gob.egg.back)


################### now want to look at average number of eggs per nest ##############

#grouped by week
bargraph.CI(x.factor = treatment, response = eggs.per.nest, group = week, legend = TRUE, main="egg production per nest", xlab="week",ylab="average number of eggs per reef", data = gob.sub)
#grouped by treatment
bargraph.CI(x.factor = week, response = eggs.per.nest, group = treatment, legend = TRUE, main="egg production per nest", xlab="week",ylab="average number of eggs per reef", data = gob.sub)



######################################
#                                    #      
#   Exploring log-transformed data   #
# WSN 11.7.18                        #
#                                    #
######################################

rm(list=ls())

library(sciplot)
library(lme4)
library(car)
library(lmerTest)

gob.sub<-read.csv("Data/10.29.18.eggs.csv")
gob.sub<-read.csv("Data/10.29.18.eggs.outliers.a.csv") #first cut out outliers removed
gob.sub<-read.csv("Data/10.29.18.eggs.outlier.b.csv") #second cut of outliers removed
gob.sub<-read.csv("Data/10.29.18.eggs.raw.outliers.csv")
gob.sub<-read.csv("Data/10.29.18.eggs.blah.csv")
gob.sub<-read.csv("Data/10.29.18.eggs.meh.csv")
#gob.sub$week<-as.factor(gob.sub$week)#need to make week a factor???
View(gob.sub)
head(gob.sub)
gob.sub$week<-as.factor(gob.sub$week)
gob.sub$treatment<-ordered(gob.sub$treatment,levels=c("Low","Medium","High"))

#IMPORTANT:
#dealing with average reproduction based on the average number of fish that I saw over the course
#of each week

#setting up log-transformed column
gob.sub$R.log<-log(gob.sub$egg.count+1)
gob.sub$Rage.log<-log(gob.sub$avg.egg.density.week+1)
View(gob.sub)

#model selection with R.log
#comparison to other model that I use in previous analysis
#uses raw egg counts, and includes deployment day as a factor
#this is what I used before
mod4b.egg<-lmer(eggs.prop.population.remaining~treatment*week+(1|deployment.day), data = gob.sub)
hist(resid(mod4b.egg))
qqnorm(resid(mod4b.egg))
qqline(resid(mod4b.egg))
anova(mod4b.egg)
summary(mod4b.egg)

mod.4b.egg.poiss<-glmer(egg.count~treatment*week+(1|deployment.day),family=poisson,data = gob.sub)
hist(resid(mod.4b.egg.poiss))
qqnorm(resid(mod.4b.egg.poiss))
qqline(resid(mod.4b.egg.poiss))
summary(mod.4b.egg.poiss)
Anova(mod.4b.egg.poiss,type = "II")


mod.log.egg<-lmer(R.log~treatment*week+(1|deployment.day), data = gob.sub)
hist(resid(mod.log.egg))


#egg counts based on densities taken on the sampling say seem to be more normal than
#egg counts averaged over teh whole week, not sure how I feel about that
#I think weekly densties are more realistic and more reliable
mod1.egg<-lmer(avg.egg.clutch.den~treatment*week+(1|deployment.day),data=gob.sub)#looks like the best model
hist(resid(mod1.egg))
qqnorm(resid(mod1.egg))
qqline(resid(mod1.egg))
anova(mod1.egg)
summary(mod1.egg)
plot(mod1.egg)
AIC(mod1.egg)

bargraph.CI(x.factor = week, response = avg.egg.clutch.den, group = treatment, Legend=TRUE, main="eggs per capita weekly density", data = gob.sub)
bargraph.CI(x.factor = treatment, response = egg.count, main="eggs per capita weekly density", data = gob.sub)

gob.sub$log.egg<-log(gob.sub$egg.count+1)

mod1.egg<-lm(egg.count~Treatment*week*avg.density.week,data=gob.sub)
hist(resid(mod1.egg))
qqnorm(resid(mod1.egg))
qqline(resid(mod1.egg))
anova(mod1.egg)
summary(mod1.egg)
plot(mod1.egg)
(mod1.egg)
AIC(mod1.egg)

gob.sub$log.density<-log(gob.sub$avg.egg.density.week+1)
mod2.egg<-aov(log.density~treatment,data=gob.sub)
hist(resid(mod2.egg))
qqnorm(resid(mod2.egg))
qqline(resid(mod2.egg))
anova(mod2.egg)
plot(mod2.egg)
(mod1.egg)
AIC(mod1.egg)


mod1.egg<-lm(avg.egg.density.week~treatment*week,data=gob.sub)
hist(resid(mod1.egg))
anova(mod1.egg)
plot(mod1.egg)
View(mod1.egg)

#number of eggs per clutch
mod1.egg<-lm(avg.egg.clutch.den~treatment*week,data=gob.sub)
hist(resid(mod1.egg))
anova(mod1.egg)
plot(mod1.egg)
View(mod1.egg)

#now lets look at egg count per week with log-transformed raw counts

mod5.egg<-lmer(R.log~treatment*week+(1|deployment.day), data = gob.sub)
hist(resid(mod4b.egg))
qqnorm(resid(mod4b.egg))
qqline(resid(mod4b.egg))
#looks fairly similar

#let's look at the egg counts averaged by weekly density, keeping dep.day in the model

mod6.egg<-lmer(avg.egg.density.week~treatment*week+(1|deployment.day), data = gob.sub)
hist(resid(mod6.egg))
qqnorm(resid(mod6.egg))
qqline(resid(mod6.egg))
#not great hist

#log of egg counts by average weekly density, no effect of deployment day
mod6a.egg<-lm(avg.egg.density.week~treatment*week, data = gob.sub)
hist(resid(mod6a.egg))
qqnorm(resid(mod6a.egg))
qqline(resid(mod6a.egg))


#now let's look at the log+1 number of eggs produced averaged by weekly density

mod7.egg<-lmer(Rage.log~treatment*week+(1|deployment.day), data = gob.sub)
hist(resid(mod7.egg))
qqnorm(resid(mod7.egg))
qqline(resid(mod7.egg))
#not so great

#poisson distributions
log.mod.6.poiss<-glmer(R.log~treatment*week+(1|deployment.day),family=poisson, data=gob.sub)
hist(resid(log.mod.6.poiss))

mod5.egg.poiss<-glmer(egg.count~treatment*week+(1|deployment.day),family=poisson, data = gob.sub)
hist(resid(mod5.egg.poiss))

#mod5.egg.gamma<-glm(Rage.log~treatment*week,family=Gamma(link="inverse"), data = gob.sub)

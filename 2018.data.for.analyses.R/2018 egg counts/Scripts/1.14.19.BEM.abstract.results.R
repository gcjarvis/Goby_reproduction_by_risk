######################################
#                                    #      
#   Egg counts with Trials 1-6 data  #
# 1/14/19                            #
# Putting together abstract for BEM  #
######################################

#want to incorporate density into the equation, but also might just go with raw counts?
#will get a df with reproduction per week, and avg. density per week, then run three separate 
#analyses:
#A) raw output: just look at results per treatment over time
#B) accounting for density: 1) divide output by average density
#                        2) use average density as a covariate with raw data

rm(list=ls())

library(sciplot)
library(lme4)
library(car)
library(lmerTest)
library(dplyr)
library(ggplot2)

gob.sub<-read.csv("Data/List of nests to count.1.10.19.csv")
gob.sub$Week<-as.factor(gob.sub$Week)
gob.sub$Treatment<-ordered(gob.sub$Treatment,levels=c("Low","Medium","High","Control"))
#density<-read.csv("C:\\Users\\George\\Desktop\\2018 summer\\2018 Goby\\2018 data for analyses, R\\2018 densities\\Data//density.1.14.19.csv")
#../goes up a level, instead of having to force R to go to a separate folder
density$Day<-as.factor(density$Day)

#total reproduction per week by reef per treatment
egg.per.week<-gob.sub %>%
  group_by(Trial,Reef,Week,Treatment) %>%
  summarize(Egg.count = sum(Egg.count))

#want to get density per reef per week 
den.per.week<-density %>%
  group_by(Trial,Reef,Week,Treatment) %>%
  summarize(Density = ceiling(mean(Density)))
#'ceiling' rounds the average density value, so you don't end up inflating reproduction averages
#rounds the density value up to the nearest whole number

#now want to join the two df's
den.and.egg<-left_join(egg.per.week,den.per.week,by=c("Trial","Week","Reef","Treatment"))
den.and.egg

#now going to calculate per capita output by dividing total egg count by density
den.avg.egg<-mutate(den.and.egg,per.capita.repro=Egg.count/Density)
den.avg.egg
View(den.avg.egg)

#adding column to df for log-transformed egg count data
#have to do log+1 because of times when egg counts were 0
den.avg.egg$ln.ec<-log(den.avg.egg$Egg.count+1)

#adding column for log+1 per capita output
den.avg.egg$ln.pc<-log(den.avg.egg$per.capita.repro+1)

#manipulating df variables

#order treatments
den.avg.egg$Treatment<-ordered(den.avg.egg$Treatment,levels=c("Low","Medium","High","Control"))
#store week as a factor for analyses
den.avg.egg$Week<-as.factor(den.avg.egg$Week)
#store trial as a factor for analyses
den.avg.egg$Trial<-as.factor(den.avg.egg$Trial)
den.avg.egg
df<-den.avg.egg

##########A) Unadjusted egg counts by treatment by week#############
#now that we have our df set up, we can go through model selection
mod1<-lm(Egg.count~Treatment, data=df)#not normal
hist(resid(mod1))
anova(mod1)
#log egg counts
mod1log<-lm(ln.ec~Treatment, data=df)#worse than non-logged counts!
hist(resid(mod1log))
anova(mod1log)

mod2<-lm(Egg.count~Treatment*Week, data=df)#more normal, but still not great
hist(resid(mod2))
qqnorm(resid(mod2))
qqline(resid(mod2))
anova(mod2)
summary(mod2)
#log
mod2log<-lm(ln.ec~Treatment*Week, data=df)#more normal, but still not great
hist(resid(mod2log))
qqnorm(resid(mod2log))
qqline(resid(mod2log))
anova(mod2log)
summary(mod2log)
#not going to waste any more time with log counts...looks terrible

mod3<-lm(Egg.count~Treatment*Week*Density, data=df)#better fit, 
hist(resid(mod3))
qqnorm(resid(mod3))
qqline(resid(mod3))
anova(mod3)
summary(mod3)
plot(mod3)#Ok fit and normailty, but ddefinitely a shape to the variances
#don't really want to start pulling data points as outliers

mod4<-lmer(Egg.count~Treatment*Week*Density+(1|Trial), data=df)
hist(resid(mod4))
qqnorm(resid(mod4))
qqline(resid(mod4))       
plot(mod4) #best fit so far, variances still look a bit funky because of some outlying points
anova(mod4)#density and week*density were sig.

mod5<-glmer(Egg.count~Treatment*Week*Density+(1|Trial),family=poisson,data=df)
hist(resid(mod5))
qqnorm(resid(mod5))
qqline(resid(mod5))       
plot(mod5)
Anova(mod5, type="II")
summary(mod5)
fixef(mod5)
ranef(mod5)
Anova(mod5)
#best fit so far, but it's tough to determine whether the variances are equal
#doing the plot function shows a clutster of points
#not sure how to analyze this, because when I try to run Anova (type=etc.) it doesn't seem to work

#########B: Per capita output by week###############
#now going to look at corrected metrics
mod1<-lm(per.capita.repro~Treatment, data=df)#not normal
hist(resid(mod1))
anova(mod1)
#log egg counts
mod1log<-lm(ln.pc~Treatment, data=df)#worse than non-logged counts!
hist(resid(mod1log))
anova(mod1log)

mod2<-lm(per.capita.repro~Treatment*Week, data=df)#more normal, but still not great
hist(resid(mod2))
qqnorm(resid(mod2))
qqline(resid(mod2))
anova(mod2) #indicates that week had an effect, which makes sense, lost fish over time, fewer eggs
summary(mod2)
#log
mod2log<-lm(ln.pc~Treatment*Week, data=df)#bad fit, likely becuase of all of the zeros
hist(resid(mod2log))
qqnorm(resid(mod2log))
qqline(resid(mod2log))
anova(mod2log)
summary(mod2log)
#not going to waste any more time with log counts...looks terrible

mod3<-lm(per.capita.repro~Treatment*Week*Density, data=df)#better fit, 
hist(resid(mod3))
qqnorm(resid(mod3))
qqline(resid(mod3))
anova(mod3) #treatment not significant here, only week, density barely nonsignificant
summary(mod3)
plot(mod3)#Ok fit and normailty, but ddefinitely a shape to the variances
#don't really want to start pulling data points as outliers

mod4<-lmer(per.capita.repro~Treatment*Week*Density+(1|Trial), data=df)
hist(resid(mod4))
qqnorm(resid(mod4))
qqline(resid(mod4))       
plot(mod4) #best fit so far, variances still look a bit funky because of some outlying points
anova(mod4)#nothing significant

mod5<-glmer(per.capita.repro~Treatment*Week*Density+(1|Trial),family=poisson,data=df)
hist(resid(mod5))
qqnorm(resid(mod5))
qqline(resid(mod5))       
plot(mod5)
Anova(mod5, type="II")
summary(mod5)
fixef(mod5)
ranef(mod5)
Anova(mod5)
#best fit so far, but it's tough to determine whether the variances are equal
#doing the plot function shows a clutster of points
#not sure how to analyze this, because when I try to run Anova (type=etc.) it doesn't seem to work
#not sure what to do for the Lmertest

####takeaway: it seems like the raw counts are a better fit, with the best fit using the raw counts
#and including density in the model as a covariate, and running it with a Poisson distribution
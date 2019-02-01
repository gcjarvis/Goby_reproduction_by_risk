#########################################
# 10/29/2018                            #  
# Density Analyses for WSN              #
# George Jarvis                         #
#########################################

#install.packages("tidyverse")

#clear workspace
rm(list=ls())

library(sciplot)
library(lme4)
library(car)
library(ggplot2)
library(lmerTest)
library(dplyr)

#quick density analyses, with just density data from 6/16,6/17, and 6/20
getwd()
gob.den<-read.csv("Data/June.density.clean.10.29.18.csv")
View(gob.den)
head(gob.den)

#There are 2 rows of data where no density was taken
gob.den<- gob.den[complete.cases(gob.den), ]#271 rows of data
#attach(gob.den)

#visualization of reefs pooled all together (front and back)
lineplot.CI(day,density,group=treatment, xlab="Days after deployment", ylab="Fish density per reef", data=gob.den)

#want to see if I can subset the data into front (1-9) and back (10-18) reefs
#good code here, taken from "discovering statistics using R" book
#back 9 reefs
gob.den.back<-gob.den[gob.den$reef>9, c("trial","date", "week", "day", "reef", "treatment", "density","release.day")]#I just want to get a subset of
# rows that only include data for reefs 1-9 and then 10-18.
View(gob.den.back)

#front 9 reefs
gob.den.front<-gob.den[gob.den$reef<10, c("trial","date", "week", "day", "reef", "treatment", "density","release.day")]#I just want to get a subset of
View(gob.den.front)

#Now want to see how the density data look separated by the two deployment dates (front, deployed first)
# (back, deployed second)

#visualization of front reefs only
lineplot.CI(day,density,group=treatment, data=gob.den.front, main="Densities Reefs 1-9", xlab="Day after deployment", ylab="Density (fish/reef)")
#really low after day 12, with a little bump back up on day 14

#visualization of back reefs only
lineplot.CI(day,density,group=treatment, data=gob.den.back, main="Densities Reefs 10-18", xlab="Day after deployment", ylab="Density (fish/reef)")

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
mod1.tot<-lm(density~treatment, data = gob.den)
mod2.tot<-lmer(density~treatment+(1|day), data = gob.den)
mod3.tot<-lmer(density~treatment+(1|day)+(1|release.day), data = gob.den)
mod4.tot.poiss<-glmer(density~treatment+(1|day)+(1|release.day),family=poisson, data = gob.den)

AIC(mod1.tot)
anova(mod2.tot,mod3.tot)
AIC(mod2.tot)
AIC(mod3.tot)
AIC(mod4.tot.poiss)
#looks like the model 3 with day and release day is the best model (AIC)

anova(mod3.tot)
summary(mod3.tot)

fixef(mod3.tot)
ranef(mod3.tot)

plot(mod3.tot)
hist(resid(mod3.tot))
Anova(mod3.tot, type="II")
Anova(mod3.tot, type="III")
qqnorm(resid(mod3.tot))
qqline(resid(mod3.tot))

#now going to look at the same models, but broken down by day
#front
#broke down by 
mod1.tot.f<-lm(density~treatment, data = gob.den.front)
mod2.tot.f<-lmer(density~treatment+(1|day), data = gob.den.front)
#mod3.tot.f<-lmer(density~treatment+(1|day), data = gob.den.front)
mod4.tot.poiss.f<-glmer(density~treatment+(1|day),family=poisson, data = gob.den.front)


###########figure for WSN talk

gob.sub<-read.csv("Data/June.density.10.29.18.csv")
gob.sub$Week<-as.factor(gob.sub$Week)
gob.sub$Treatment<-ordered(gob.sub$Treatment,levels=c("Low","Medium","High"))
gob.sub<- gob.sub[complete.cases(gob.sub), ]

gob.survey<-lmer(Density~Treatment*Week+(1|deployment.day),data=gob.sub)
hist(resid(gob.survey))
qqnorm(resid(gob.survey))
qqline(resid(gob.survey))
anova(gob.survey)
summary(gob.survey)

Anova(gob.survey,type="II")

viz.sur<-with(gob.sub, aggregate((Density), list(Week=Week,Treatment=Treatment), mean))
#now apply the se function to the 4th column [,3]
viz.sur$se<-with(gob.sub, aggregate((Density), list(Week=Week,Treatment=Treatment), function(x) sd(x)/sqrt(length(x))))[,3]
viz.sur

View(gob.sub)

p <- ggplot(viz.sur, aes(x=Week, y=x, group=Treatment, color=Treatment))+ 
  geom_linerange(aes(ymin=x-se, ymax=x+se), 
                 position=position_dodge(0)) +
  geom_line(size=0.75)+
  geom_line(aes(linetype=Treatment)) +
  geom_point(aes(shape=Treatment))+
  labs(x="Week", y = "Number of Gobies per Reef")+
  theme_classic() + 
  theme(axis.text.x=element_text(size=18, colour="black"),axis.text.y=element_text(size=18, colour="black"), axis.title=element_text(size=20,face="bold")) +
  theme(axis.title.y = element_text(size= 22, margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)), axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(breaks = seq(0,15,2)) +
  theme(text = element_text(family="Arial"))
p + scale_color_manual(values=c("#0072B2","#009E73","#D55E00")) + theme(legend.text=element_text(size=15)) + theme(legend.title =element_text(size=15, face="bold"))

p <- ggplot(gob.survey, aes(x=Week, y=x, group=Treatment, color=Treatment))+ 
  geom_linerange(aes(ymin=x-se, ymax=x+se), 
                 position=position_dodge(0)) +
  geom_line(size=0.75)+
  geom_line(aes(linetype=Treatment)) +
  geom_point(aes(shape=Treatment))+
  labs(x="Week", y = "Density per Reef")+
  theme_classic() + 
  theme(axis.text.x=element_text(size=20, colour="black"),axis.text.y=element_text(size=16, colour="black"), axis.title=element_text(size=22,face="bold")) +
  theme(axis.title.y = element_text(size= 22, margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)), axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0, 0)) +
  theme(text = element_text(family="Arial"))
p + scale_color_manual(values=c("#D55E00","#0072B2","#009E73")) + theme(legend.text=element_text(size=15)) + theme(legend.title =element_text(size=15, face="bold"))




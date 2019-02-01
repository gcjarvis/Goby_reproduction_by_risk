#########################################
# 12/3/2018                             #  
# Density Analyses                      #
# George Jarvis                         #
#########################################

#clear workspace
rm(list=ls())

library(sciplot)
library(lme4)
library(car)
library(ggplot2)
library(lmerTest)
library(dplyr)

#######subsetting data and looking at trial 5 densities##########
# as of 12/3/18 data entry for densities (3/4 of data from trial 5)
#filter function in dplyr is a good way to subset data for analyses

#importing data that includes trial 5 densities
gob.den<-read.csv("Data/density.11.28.18.csv")
gob.den$Trial<-as.factor(gob.den$Trial)
gob.den$Week<-as.factor(gob.den$Week)
gob.den.4<-filter(gob.den, Trial==4) # good way to subset data quickly
gob.den.5<-filter(gob.den, Trial==5) 
head(gob.den.5)

#trial 4 only
lineplot.CI(Week,Density,group=Treatment,legend = TRUE, xlab="Week", ylab="Fish density per reef", data=gob.den.4)
#no difference in the number of gobies seen per week, regardless of treatment

#trial 4 densitities only
#fewer fish seen over time, but equally among treatments
gob.survey<-lmer(Density~Treatment*Week+(1|deployment.day),data=gob.den.4)
hist(resid(gob.survey))
qqnorm(resid(gob.survey))
qqline(resid(gob.survey))
anova(gob.survey)
summary(gob.survey)

Anova(gob.survey,type="II")

#trial 5 only
lineplot.CI(Week,Density,group=Treatment,legend = TRUE, xlab="Week", ylab="Fish density per reef", data=gob.den.5)
#seems like much lower number of gobies seen in the high risk treatments compared to low and medium
#stats show that there are differences by treatment, and differences over time:
#fewer fish seen in high risk trt, and fewer fish seen over time

gob.survey<-lmer(Density~Treatment*Week+(1|deployment.day),data=gob.den.5)
hist(resid(gob.survey))
qqnorm(resid(gob.survey))
qqline(resid(gob.survey))
anova(gob.survey)
summary(gob.survey)

Anova(gob.survey,type="II")

#combined trials
lineplot.CI(Week,Density,group=Treatment,legend = TRUE, xlab="Week", ylab="Fish density per reef", data=gob.den)
#looks like the low densities from trial 5 are driving differences overall
#doesn't seem like trial acocunts for much of the variation in the data (when trial is random factor)
#seems like treatment*trial is significant though (when trial is included as a crossed, fixed effect)

#random effect, doesn't account for much of the total variation
gob.survey<-lmer(Density~Treatment*Week+(1|deployment.day)+(1|Trial),data=gob.den)
hist(resid(gob.survey))
qqnorm(resid(gob.survey))
qqline(resid(gob.survey))
anova(gob.survey)
summary(gob.survey)

#fixed, crossed effect, different effect of treatment depending on the trial
Anova(gob.survey,type="II")

gob.survey<-lmer(Density~Treatment*Week*Trial+(1|deployment.day),data=gob.den)
hist(resid(gob.survey))
qqnorm(resid(gob.survey))
qqline(resid(gob.survey))
anova(gob.survey)
summary(gob.survey)

Anova(gob.survey,type="II")

#######vizualization#####

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
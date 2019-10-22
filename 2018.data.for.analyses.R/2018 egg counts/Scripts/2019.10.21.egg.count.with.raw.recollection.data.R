# Description: reanalyzing, using raw number of gobies recollected as covariate
# Author: George C Jarvis
# Date: Mon Oct 21 16:34:41 2019
# Notes: I think this is more representative than using an average (20+ # reco/2),
#       so I want to go with this over a calculated metric
#       A) I went back through and redid the analyses, and it seems to be 
#           better this way, so I'm going to use that raw data as the covariate
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
#includes a cloumn ("Treatment") where uncaged and HR are coded as "High"
#also includes a column ("T6.comparison") where uncaged and high are separated
repro<-read.csv("Data/new.data.2019.9.30.csv", na.strings = "")

#mixed model####
mod.2<-lmer(Egg.count~ Treatment * Recollection + (1|Trial), data=repro)
hist(resid(mod.2))
qqnorm(resid(mod.2))
qqline(resid(mod.2))
anova(mod.2)
Anova(mod.2) 
fixef(mod.2)

#plotting relationship between recolelctions and egg count
plot(repro$Recollection, repro$Egg.count)
abline(lm(repro$Egg.count~repro$Recollection))
#positive relationship, but much more data for lower numbers of fish recollected

mod.ER<-lm(Egg.count~Recollection, data=repro)
hist(resid(mod.ER))
qqnorm(resid(mod.ER))
qqline(resid(mod.ER))
anova(mod.ER)
Anova(mod.ER)
summary(mod.ER)
#shows clear relationship, but only has an R-squared of 0.385
#this makes sense becasue of the large spread of repro for each recollection value

#making a B&W figure for this analysis####
repro$Treatment<-ordered(repro$Treatment,levels=c("Low","Medium","High"))
egg.anc<-with(repro, aggregate((Egg.count), list(Treatment=Treatment), mean))
egg.anc$se<-with(repro, aggregate((Egg.count), list(Treatment=Treatment), 
                                  function(x) sd(x)/sqrt(length(x))))[,2]

png(filename = "Output/19.10.22.gray.ancova.legend.4.png", width = 900, height = 600)

#plot seems a little more realistic when I don't round to the nearest whole number (ceiling)

#building from the bottom, not including shape by treatment
anc1<-ggplot(repro, aes(Recollection, Egg.count, shape=Treatment, linetype=Treatment, col=Treatment)) +
  geom_smooth(method="lm", se=FALSE, show.legend = TRUE)  +
  geom_point(size=3)+
  theme_classic()+
  labs(x="Number of Gobies Recollected", y="Total Reproduction Per Reef")+
  expand_limits(y=0)+
  scale_color_manual(values=c("black", "#666666", "grey"))+
  scale_linetype_manual(values=c("solid", "dashed", "twodash"))+
  theme(axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"), 
        axis.title=element_text(size=20))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), 
        axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0)), 
        axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(legend.text=element_text(size=18)) +
  theme(legend.title =element_text(size=20))+
  labs(color  = "Perceived Risk", linetype = "Perceived Risk", shape = "Perceived Risk")
anc1
#sort of lame that I had to code it like this to get the legend correc
dev.off()
# Description: grayscale plot for manuscript, density
# Author: George C Jarvis
# Date: Wed Oct 23 17:34:33 2019
# Notes: Will have to figure out how to signify that the last points for each
#       treatment are the number of fish that I recollected from each reef. That way,
#       I can conserve the number of figures that I have to include in the MS
# --------------

rm(list=ls())

#loading libraries
library(sciplot)
library(lme4)
library(car)
library(lmerTest)
library(dplyr)
library(ggplot2)
library(MASS)
library(nlme)
library(pwr)

#loading data
viz.surv<-read.csv("Data/density.2019.10.1.csv", na.strings = "")
viz.surv$Density<-as.numeric(viz.surv$Density)
viz.surv$den.max<-as.numeric(viz.surv$den.max)

#models, comparing density vs. den.max
#density
mod.1a<-glmer(Density ~ Treatment + (1|Trial),family = poisson, data=viz.surv)
hist(resid(mod.1a))
qqnorm(resid(mod.1a))
qqline(resid(mod.1a))
anova(mod.1a)
Anova(mod.1a)

#den.max
mod.1<-glmer(den.max ~ Treatment + (1|Trial),family=poisson, data=viz.surv)
hist(resid(mod.1))
qqnorm(resid(mod.1)) #not really normal, maybe there's a better transformation
qqline(resid(mod.1))
anova(mod.1)
Anova(mod.1) 
summary(mod.1)
#chi-squared shows an effect of treatment, showing more fish seen in HR/uncaged treament
fixed.effects(mod.1)
ranef(mod.1)

#plotting####
png(filename = "Output/19.10.23.densities.png", width = 900, height = 600)

#plot seems a little more realistic when I don't round to the nearest whole number (ceiling)

#building from the bottom, not including shape by treatment
anc1<-ggplot(repro, aes(Recollection, Egg.count, shape=Treatment, linetype=Treatment, col=Treatment)) +
  geom_smooth(method="lm", se=FALSE, show.legend = TRUE)  +
  geom_point(size=3)+
  theme_classic()+
  labs(x="Number of Gobies Recollected", y="Total Egg Production Per Reef")+
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
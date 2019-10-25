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
viz.surv<-read.csv("Data/density.2019.10.1.csv")
viz.surv<-na.omit(viz.surv)
viz.surv$Density<-as.numeric(viz.surv$Density)
viz.surv$den.max<-as.numeric(viz.surv$den.max)

#models, comparing density vs. den.max
#trying linear model, no poisson distribution
mod.a<-lmer(Density~Treatment*Day+(1|Trial), data=viz.surv)
hist(resid(mod.a))
qqnorm(resid(mod.a))
qqline(resid(mod.a))
anova(mod.a)
Anova(mod.a)

mod.b<-lmer(den.max~Treatment*Day+(1|Trial),data=viz.surv)#probably go with this model
hist(resid(mod.b))
qqnorm(resid(mod.b))
qqline(resid(mod.b))
anova(mod.b) 
Anova(mod.b)

#density
mod.1a<-glmer(Density ~ Treatment * Day + (1|Trial),family = poisson, data=viz.surv)
hist(resid(mod.1a))
qqnorm(resid(mod.1a))
qqline(resid(mod.1a))
anova(mod.1a)
Anova(mod.1a)

#den.max
mod.1<-glmer(den.max ~ Treatment * Day + (1|Trial),family=poisson, data=viz.surv)
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

#sort of a neat plot showing the raw data for visual surveys
png(filename = "Output/19.10.23.densities.all.data.not.for.pub.png", width = 900, height = 600)

#building from the bottom, not including shape by treatment
den1<-ggplot(viz.surv, aes(Day, Density, shape=Treatment, linetype=Treatment, col=Treatment)) +
  geom_smooth(se=TRUE, show.legend = TRUE)  +
  #geom_smooth(method="lm", se=FALSE, show.legend = TRUE)
  geom_point(size=3)+
  theme_classic()+
  labs(x="Day", y="Number of Gobies Seen")+
  expand_limits(y=0)+
  scale_color_manual(values=c("black", "#666666", "grey"))
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
den1
#sort of lame that I had to code it like this to get the legend correc
dev.off()

# 1) with density, not den.max

#average density per treatment, making df for ggplot
viz.surv$Treatment<-ordered(viz.surv$Treatment,levels=c("Low","Medium","High"))
den<-with(viz.surv, aggregate((Density), list(Day=Day,Treatment=Treatment), mean))
den$se<-with(viz.surv, aggregate((Density), list(Day=Day,Treatment=Treatment), 
                                          function(x) sd(x)/sqrt(length(x))))[,3]

png(filename = "Output/2019.10.24.den.manuscript.means.png", width = 1000, height = 500)

den.plot <- ggplot(den, aes(x=Day, y=x, shape=Treatment, color=Treatment, linetype=Treatment))+ 
  geom_linerange(aes(ymin=x-se, ymax=x+se), 
                 position=position_dodge(0)) +
  geom_point(size=3)+
  labs(x="Day", y = "Number of Fish Seen")+
  theme_classic() + 
  scale_color_manual(values=c("black", "#666666", "grey"))+
  scale_linetype_manual(values=c("solid", "dashed", "twodash"))+
  geom_line(aes(linetype=Treatment))+
  geom_line(size=1.0)+
  theme(axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"), 
        axis.title=element_text(size=20))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), 
        axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0)), 
        axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(legend.text=element_text(size=18)) +
  theme(legend.title =element_text(size=20))+
  scale_y_continuous(expand = c(0, 0),limits = c(0,21),breaks = c(5,10,15,20)) + scale_x_continuous(breaks=c(5, 10, 15, 20,25,30))+
  labs(color  = "Perceived Risk", linetype = "Perceived Risk", shape = "Perceived Risk")
den.plot
dev.off()

# 2) with den.max, ###change the df, and then replot####

#average density per treatment
viz.surv$Treatment<-ordered(viz.surv$Treatment,levels=c("Low","Medium","High"))
den<-with(viz.surv, aggregate((den.max), list(Day=Day,Treatment=Treatment), mean))
den$se<-with(viz.surv, aggregate((den.max), list(Day=Day,Treatment=Treatment), 
                                 function(x) sd(x)/sqrt(length(x))))[,3]

png(filename = "Output/2019.10.24.den.max.MS.means.png", width = 1000, height = 500)

den.plot <- ggplot(den, aes(x=Day, y=x, shape=Treatment, color=Treatment, linetype=Treatment))+ 
  geom_linerange(aes(ymin=x-se, ymax=x+se), 
                 position=position_dodge(0)) +
  geom_point(size=3)+
  labs(x="Day", y = "Number of Fish Seen")+
  theme_classic() + 
  scale_color_manual(values=c("black", "#666666", "grey"))+
  scale_linetype_manual(values=c("solid", "dashed", "twodash"))+
  geom_line(aes(linetype=Treatment))+
  geom_line(size=1.0)+
  theme(axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"), 
        axis.title=element_text(size=20))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), 
        axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0)), 
        axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(legend.text=element_text(size=18)) +
  theme(legend.title =element_text(size=20))+
  scale_y_continuous(expand = c(0, 0),limits = c(0,21),breaks = c(5,10,15,20)) + scale_x_continuous(breaks=c(5, 10, 15, 20,25,30))+
  labs(color  = "Perceived Risk", linetype = "Perceived Risk", shape = "Perceived Risk")
den.plot
dev.off()

#next thing to do####
#rerun the old analyses with the same data, but remove T6
#need to figure out what is going on with these data...and why
#there are inconsistencies with the plots

#2019.10.25 -- going to try and rerun the analyses with the data that I used
# when I did my defense/thesis, to find out what the problem is

#with old data, using density
gob.den<-read.csv("Data/density.2019.4.26.csv", na.strings = "")
gob.den$Treatment<-ordered(gob.den$Treatment,levels=c("Low","Medium","High","Control"))
gob.den$Density<-as.numeric(gob.den$Density)
retry<-with(gob.den, aggregate((Density), list(Day=Day,Treatment=Treatment), mean))
retry$se<-with(gob.den, aggregate((Density), list(Day=Day,Treatment=Treatment), 
                                 function(x) sd(x)/sqrt(length(x))))[,3]

#png(filename = "Output/2019.10.24.den.manuscript.means.png", width = 1000, height = 500)

den.plot.2 <- ggplot(retry, aes(x=Day, y=x, shape=Treatment, color=Treatment, linetype=Treatment))+ 
  geom_linerange(aes(ymin=x-se, ymax=x+se), 
                 position=position_dodge(0)) +
  geom_point(size=3)+
  labs(x="Day", y = "Number of Fish Seen")+
  theme_classic() + 
  scale_color_manual(values=c("black", "#666666", "grey","grey"))+
  scale_linetype_manual(values=c("solid", "dashed", "twodash", "dashed"))+
  geom_line(aes(linetype=Treatment))+
  geom_line(size=1.0)+
  theme(axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"), 
        axis.title=element_text(size=20))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), 
        axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0)), 
        axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(legend.text=element_text(size=18)) +
  theme(legend.title =element_text(size=20))+
  scale_y_continuous(expand = c(0, 0),limits = c(0,21),breaks = c(5,10,15,20)) + scale_x_continuous(breaks=c(5, 10, 15, 20,25,30))+
  labs(color  = "Perceived Risk", linetype = "Perceived Risk", shape = "Perceived Risk")
den.plot.2
dev.off()

#with old data using den.max

gob.den$den.max<-as.numeric(gob.den$den.max)
retry<-with(gob.den, aggregate((den.max), list(Day=Day,Treatment=Treatment), mean))
retry$se<-with(gob.den, aggregate((den.max), list(Day=Day,Treatment=Treatment), 
                                  function(x) sd(x)/sqrt(length(x))))[,3]

#png(filename = "Output/2019.10.24.den.manuscript.means.png", width = 1000, height = 500)

den.plot.2 <- ggplot(retry, aes(x=Day, y=x, shape=Treatment, color=Treatment, linetype=Treatment))+ 
  geom_linerange(aes(ymin=x-se, ymax=x+se), 
                 position=position_dodge(0)) +
  geom_point(size=3)+
  labs(x="Day", y = "Number of Fish Seen")+
  theme_classic() + 
  scale_color_manual(values=c("black", "#666666", "grey","grey"))+
  scale_linetype_manual(values=c("solid", "dashed", "twodash", "dashed"))+
  geom_line(aes(linetype=Treatment))+
  geom_line(size=1.0)+
  theme(axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"), 
        axis.title=element_text(size=20))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), 
        axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0)), 
        axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(legend.text=element_text(size=18)) +
  theme(legend.title =element_text(size=20))+
  scale_y_continuous(expand = c(0, 0),limits = c(0,21),breaks = c(5,10,15,20)) + scale_x_continuous(breaks=c(5, 10, 15, 20,25,30))+
  labs(color  = "Perceived Risk", linetype = "Perceived Risk", shape = "Perceived Risk")
den.plot.2
dev.off()

#breaking down by trial to see what's up, using density first, then den max second
#just trials 4 and 5

den.t4.5<-gob.den[(gob.den$Trial>3) & (gob.den$Trial<6), ]
retry.t4.5.den<-with(den.t4.5, aggregate((Density), list(Day=Day,Treatment=Treatment), mean))
retry.t4.5.den$se<-with(new.t4.5, aggregate((Density), list(Day=Day,Treatment=Treatment), 
                                        function(x) sd(x)/sqrt(length(x))))[,3]

#png(filename = "Output/2019.10.24.den.manuscript.means.png", width = 1000, height = 500)

den.plot.4 <- ggplot(retry.t4.5.den, aes(x=Day, y=x, shape=Treatment, color=Treatment, linetype=Treatment))+ 
  geom_linerange(aes(ymin=x-se, ymax=x+se), 
                 position=position_dodge(0)) +
  geom_point(size=3)+
  labs(x="Day", y = "Number of Fish Seen")+
  theme_classic() + 
  scale_color_manual(values=c("black", "#666666", "grey"))+
  scale_linetype_manual(values=c("solid", "dashed", "twodash"))+
  geom_line(aes(linetype=Treatment))+
  geom_line(size=1.0)+
  theme(axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"), 
        axis.title=element_text(size=20))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), 
        axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0)), 
        axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(legend.text=element_text(size=18)) +
  theme(legend.title =element_text(size=20))+
  scale_y_continuous(expand = c(0, 0),limits = c(0,21),breaks = c(5,10,15,20)) + scale_x_continuous(breaks=c(5, 10, 15, 20,25,30))+
  labs(color  = "Perceived Risk", linetype = "Perceived Risk", shape = "Perceived Risk")
den.plot.4
dev.off()

new.t4.5<-gob.den[(gob.den$Trial>3) & (gob.den$Trial<6), ]
retry.t4.5<-with(new.t4.5, aggregate((den.max), list(Day=Day,Treatment=Treatment), mean))
retry.t4.5$se<-with(new.t4.5, aggregate((den.max), list(Day=Day,Treatment=Treatment), 
                                  function(x) sd(x)/sqrt(length(x))))[,3]

#png(filename = "Output/2019.10.24.den.manuscript.means.png", width = 1000, height = 500)

den.plot.3 <- ggplot(retry.t4.5, aes(x=Day, y=x, shape=Treatment, color=Treatment, linetype=Treatment))+ 
  geom_linerange(aes(ymin=x-se, ymax=x+se), 
                 position=position_dodge(0)) +
  geom_point(size=3)+
  labs(x="Day", y = "Number of Fish Seen")+
  theme_classic() + 
  scale_color_manual(values=c("black", "#666666", "grey"))+
  scale_linetype_manual(values=c("solid", "dashed", "twodash"))+
  geom_line(aes(linetype=Treatment))+
  geom_line(size=1.0)+
  theme(axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"), 
        axis.title=element_text(size=20))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), 
        axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0)), 
        axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(legend.text=element_text(size=18)) +
  theme(legend.title =element_text(size=20))+
  scale_y_continuous(expand = c(0, 0),limits = c(0,21),breaks = c(5,10,15,20)) + scale_x_continuous(breaks=c(5, 10, 15, 20,25,30))+
  labs(color  = "Perceived Risk", linetype = "Perceived Risk", shape = "Perceived Risk")
den.plot.3
dev.off()
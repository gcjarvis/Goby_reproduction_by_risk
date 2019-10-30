# Description: final analyses for density for MS
# Author: George C Jarvis
# Date: Fri Oct 25 11:20:42 2019
# Notes: Code includes models, along with plots for manuscript. I used the 
#       metric for maximum number of fish that were on the reef, not just the 
#       number of fish that I saw on that day. In other words, if I saw more fish
#       on a later date for a single reef, I went back and made any observations
#       the maximum number that it could have been.
#       B)M. Steele wanted me to include the number of fish recollected from each reef
#         in this plot, but I can't do that because the number of fish seen is confounded
#         by treatment. In this case, the number seen was less than the number recollected,
#         at least for the high-risk treatment, given that the number of fish recollected 
#         was the same for all treatments
# --------------

rm(list=ls())

#loading packages
library(lme4)
library(lmerTest)
library(ggplot2)
library(dplyr)

#loading data
viz.surv<-read.csv("Data/density.2019.10.1.csv")
#omitting rows with NA values for density, where no survey was done
viz.surv<-na.omit(viz.surv)
viz.surv$den.max<-as.numeric(viz.surv$den.max)
#viz.surv$Day<-as.factor(viz.surv$Day)

#mixed model, with den.max as the response
mod.1<-lmer(den.max ~ Treatment*Day + (1|Trial), data=viz.surv)
hist(resid(mod.1))
qqnorm(resid(mod.1))
qqline(resid(mod.1))
anova(mod.1) #sig. effect of trt., and also of Day, no interactive effects
#L+M>H, and fewer fish seen over time, saw the same relative number of fish,
# regardless of treatment (no interactive effects, need to make that more clear)
summary(mod.1)
#chi-squared shows an effect of treatment, showing more fish seen in HR/uncaged treament
fixef(mod.1)
ranef(mod.1)

#repeated measures version of the same model
#mixed model, with den.max as the response
mod.1.nest<-lmer(den.max ~ Treatment*Day + (1|Treatment:Reef) + (1|Trial), data=viz.surv)
hist(resid(mod.1.nest))
qqnorm(resid(mod.1.nest))
qqline(resid(mod.1.nest))
anova(mod.1.nest) #sig. effect of trt., and also of Day, no interactive effects
#L+M>H, and fewer fish seen over time, saw the same relative number of fish,
# regardless of treatment (no interactive effects, need to make that more clear)
summary(mod.1.nest)
#chi-squared shows an effect of treatment, showing more fish seen in HR/uncaged treament
fixef(mod.1.nest)
ranef(mod.1.nest)

#mixed model, with den.max as the response, no repeated-measure, but including reef as a random effect
mod.reef<-lmer(den.max ~ Treatment*Day + (1|Reef) + (1|Trial), data=viz.surv)
hist(resid(mod.reef))
qqnorm(resid(mod.reef))
qqline(resid(mod.reef))
anova(mod.reef) #sig. effect of trt., and also of Day, no interactive effects
#L+M>H, and fewer fish seen over time, saw the same relative number of fish,
# regardless of treatment (no interactive effects, need to make that more clear)
summary(mod.reef)
#chi-squared shows an effect of treatment, showing more fish seen in HR/uncaged treament
fixef(mod.reef)
ranef(mod.reef)

#plotting
viz.surv$Treatment<-ordered(viz.surv$Treatment,levels=c("Low","Medium","High"))
den<-with(viz.surv, aggregate((den.max), list(Day=Day,Treatment=Treatment), mean))
den$se<-with(viz.surv, aggregate((den.max), list(Day=Day,Treatment=Treatment), 
                                 function(x) sd(x)/sqrt(length(x))))[,3]

png(filename = "Output/2019.10.24.den.max.MS.means.0.75.thickness.png", width = 1000, height = 500)

den.plot <- ggplot(den, aes(x=Day, y=x, shape=Treatment, color=Treatment, linetype=Treatment))+ 
  geom_linerange(aes(ymin=x-se, ymax=x+se), 
                 position=position_dodge(0)) +
  geom_point(size=3)+
  labs(x="Day", y = "Number of Fish Seen")+
  theme_classic() + 
  scale_color_manual(values=c("black", "#666666", "grey"))+
  scale_linetype_manual(values=c("solid", "dashed", "twodash"))+
  geom_line(aes(linetype=Treatment))+
  geom_line(size=0.75)+
  theme(axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"), 
        axis.title=element_text(size=20))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), 
        axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0)), 
        axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(legend.text=element_text(size=18)) +
  theme(legend.title =element_text(size=20))+
  scale_y_continuous(expand = c(0, 0),limits = c(0,18),breaks = c(4,8,12,16)) + scale_x_continuous(expand=c(0,0), limits= c(0,29),breaks=c(2,4,6,8,10,12,14,16,18,20,22,24,26,28))+
  labs(color  = "Perceived Risk", linetype = "Perceived Risk", shape = "Perceived Risk")
den.plot
dev.off()

#7, 14, 21,28

#with default line thickness, including space at the end of the plot to add in recollection data

png(filename = "Output/2019.10.24.den.max.MS.means.default.thickness.png", width = 1000, height = 500)
png(filename = "Output/2019.10.24.den.max.space.to.plot.recollections.png", width = 1200, height = 600, res=300)
#I like the 12x4 size for this one
png("Output/2019.10.28.test.300dpi.12x4.best.png", width = 12, height = 4, units = 'in', res = 300)
#taller figure to see the distinction in recollections
png("Output/2019.10.28.density.15x7.300dpi.combo.shorter.x.axis.png", width = 15, height = 7, units = 'in', res = 300)

den.plot <- ggplot(den, aes(x=Day, y=x, shape=Treatment, color=Treatment, linetype=Treatment))+ 
  geom_linerange(aes(ymin=x-se, ymax=x+se), 
                 position=position_dodge(0)) +
  geom_point(size=3)+
  labs(x="Day", y = "Number of Fish Seen")+
  theme_classic() + 
  scale_color_manual(values=c("black", "#666666", "grey"))+
  scale_linetype_manual(values=c("solid", "dashed", "twodash"))+
  geom_line(aes(linetype=Treatment))+
  geom_line(size=0.75)+
  theme(axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"), 
        axis.title=element_text(size=20))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), 
        axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0)), 
        axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(legend.text=element_text(size=18)) +
  theme(legend.title =element_text(size=20))+
  scale_y_continuous(expand = c(0, 0),limits = c(0,18),breaks = c(4,8,12,16)) + scale_x_continuous(expand=c(0,0), limits= c(0,30),breaks=c(2,4,6,8,10,12,14,16,18,20,22,24,26,28))+
  labs(color  = "Perceived Risk", linetype = "Perceived Risk", shape = "Perceived Risk")
den.plot
dev.off()

#going to figure out how to bring in recollection data, then add the single points to the figure

#adding the letter R to the label for the figure, just made the axis label different,
# -and then made the size of the test for the label the same as the other plot

reco<-read.csv("C:/Users/George/Desktop/2018 summer/2018 Goby/Goby_reproduction_by_risk/2018.data.for.analyses.R/2018 recollections/Data/2019.10.8.recollection.data.csv")

#setting up df for ggplot
reco$Treatment<-ordered(reco$Treatment,levels=c("Low","Medium","High"))
reco.fig<-with(reco, aggregate((Count), list(Treatment=Treatment), mean))
reco.fig$se<-with(reco, aggregate((Count), list(Treatment=Treatment), 
                                 function(x) sd(x)/sqrt(length(x))))[,2]

png("Output/2019.10.28.recollection.labels.png", width = 15, height = 7, units = 'in', res = 300)
reco.plot <- ggplot(reco.fig, aes(x="", y=x, shape=Treatment, color=Treatment, linetype=Treatment))+ 
  geom_linerange(aes(ymin=x-se, ymax=x+se), 
                 position=position_dodge(0)) +
  geom_point(size=3)+
  labs(x="", y = "Number of Fish Seen")+
  theme_classic() + 
  scale_color_manual(values=c("black", "#666666", "grey"))+
  scale_linetype_manual(values=c("solid", "dashed", "twodash"))+
  geom_line(aes(linetype=Treatment))+
  geom_line(size=0.75)+
  theme(axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"), 
        axis.title=element_text(size=20))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), 
        axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0)), 
        axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(legend.text=element_text(size=18)) +
  theme(legend.title =element_text(size=20))+
  scale_y_continuous(expand = c(0, 0),limits = c(0,18),breaks = c(4,8,12,16))+ #scale_x_continuous(expand=c(0,0), limits= c(0,32),breaks=c(2,4,6,8,10,12,14,16,18,20,22,24,26,28))+
  scale_x_discrete(labels=c("R"))+
  labs(color  = "Perceived Risk", linetype = "Perceived Risk", shape = "Perceived Risk")
reco.plot
dev.off()
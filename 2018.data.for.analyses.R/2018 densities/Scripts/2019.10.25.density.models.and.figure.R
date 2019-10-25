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

#loading data
viz.surv<-read.csv("Data/density.2019.10.1.csv")
#omitting rows with NA values for density, where no survey was done
viz.surv<-na.omit(viz.surv)
viz.surv$den.max<-as.numeric(viz.surv$den.max)
viz.surv$Day<-as.factor(viz.surv$Day)

#mixed model, with den.max as the 
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

7, 14, 21,28

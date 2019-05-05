# Description: Density analyses for thesis, decided to go with day, not week, so let's see if this has a big effect
# Author: George C Jarvis
# Date: Fri Apr 26 21:39:20 2019
# --------------

#I think if I go with day for all analyses, it's better

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
library(extrafont)
library(HH)

gob.den<-read.csv("Data/density.2019.4.26.csv")
gob.den<- na.omit(gob.den) #remove NA's

#subsetting by trial
den.2017<-gob.den[(gob.den$Trial<4), ]
den.2018.t4.5<-gob.den[(gob.den$Trial>3) & (gob.den$Trial<6), ]
den.2018.t6<-gob.den[(gob.den$Trial>5), ]

#linear models by treatment and day####

#2017######
mod1<-aov(Density~Treatment*Day, data=den.2017)
hist(resid(mod1))#more normal
qqnorm(resid(mod1))
qqline(resid(mod1))
anova(mod1)
plot(mod1)#again, looks like equal variance
TukeyHSD(mod1) #low and med are not diff from each other, but high is lower than both
tapply(den.2017$Density,den.2017$Treatment,mean)
#   High       Low    Medium 
# 8.293103 10.120690 10.275862 = 21% fewer fish seen on high-risk treatments

#2017, adding trial as a fixed factor, and nesting reef within treatment as a random factor to account for repeated measures
mod1<-lmer(Density~Treatment*Day + (1|Trial) + (1|Treatment:Reef), data=den.2017)
hist(resid(mod1))#more normal
qqnorm(resid(mod1))
qqline(resid(mod1))
anova(mod1)

mod1a<-lmer(Density~Treatment*Day*Trial + (1|Treatment:Reef), data=den.2017)
hist(resid(mod1a))#more normal
qqnorm(resid(mod1a))
qqline(resid(mod1a))
anova(mod1a)
Anova(mod1a)

#best model, the number of gobies seen on each day depended on the trial
#there was a margin. insig. interaction of treatment and trial, 
#there was a sig. diff in number of gobies seen over time

mod1b<-lmer(Density~Treatment*Day + Trial + (1|Treatment:Reef), data=den.2017)
hist(resid(mod1b))#more normal
qqnorm(resid(mod1b))
qqline(resid(mod1b))
anova(mod1b)

AIC(mod1,mod1a, mod1b)
ancova(Density~Treatment*Day, data=den.2017)
plot(mod1)#again, looks like equal variance

#seems like the best model is the one with trial included as a fixed factor
#interesting that when I run as a repeated measures, no sig. effect of treament

#2017 plot####
den.2017$Treatment<-ordered(den.2017$Treatment,levels=c("Low","Medium","High"))
viz.sur.2017<-with(den.2017, aggregate((Density), list(Day=Day,Treatment=Treatment), mean))
viz.sur.2017$se<-with(den.2017, aggregate((Density), list(Day=Day,Treatment=Treatment), 
                                          function(x) sd(x)/sqrt(length(x))))[,3]

png(filename = "Output/2017.den.thesis.png", width = 1000, height = 500)

p <- ggplot(viz.sur.2017, aes(x=Day, y=x, group=Treatment, color=Treatment))+ 
  geom_linerange(aes(ymin=x-se, ymax=x+se), 
                 position=position_dodge(0)) +
  geom_line(size=1.5)+
  geom_line(aes(linetype=Treatment)) +
  geom_point(aes(shape=Treatment))+
  labs(x="Day", y = "Number of Fish Seen")+
  theme_classic() + 
  theme(axis.text.x=element_text(size=28, colour="black"),axis.text.y=element_text(size=28, colour="black"), axis.title=element_text(size=35,face="bold")) +
  theme(axis.title.y = element_text(size= 34, margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)), axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0, 0),limits = c(0,18),breaks = c(4,8,12,16)) + scale_x_continuous(breaks=c(1,2,3,4,5,6))+
  theme(text = element_text(family="Arial"))
p + scale_color_manual(values=c("#0072B2","#009E73","#D55E00")) + theme(legend.text=element_text(size=20)) + theme(legend.title =element_text(size=28, face="bold"))


dev.off()

#2018.t4.5#####

#2018, adding trial as a fixed factor, and nesting reef within treatment as a random factor to account for repeated measures

#removed outlying data for day 18
mod2<-lmer(Density~Treatment*Day + (1|Trial) + (1|Treatment:Reef), data=gob.den)
hist(resid(mod2))#more normal
qqnorm(resid(mod2))
qqline(resid(mod2))
anova(mod2)
plot(mod2)

#best model based on AIC, going to go with this one
#sig. decrease in number of fish seen over time, but margnially insignigifcant 
# effect of treatment on number of fish seen

mod2a<-lmer(Density~Treatment*Day*Trial + (1|Treatment:Reef), data=gob.den)
hist(resid(mod2a))#more normal
qqnorm(resid(mod2a))
qqline(resid(mod2a))
anova(mod2a)
Anova(mod2a)
plot(mod2a)

mod2b<-lmer(Density~Treatment*Day + Trial + (1|Treatment:Reef), data=gob.den)
hist(resid(mod2b))#more normal
qqnorm(resid(mod2b))
qqline(resid(mod2b))
anova(mod2b)
Anova(mod2b)
plot(mod2b)

AIC(mod2,mod2a, mod2b)

mod2<-lmer(Density~Treatment*Day + (1|Trial) + (1|Treatment:Reef), data=den.2018.t4.5)
hist(resid(mod2))#more normal
qqnorm(resid(mod2))
qqline(resid(mod2))
anova(mod2)

mod2a<-lmer(Density~Treatment*Day*Trial + (1|Treatment:Reef), data=den.2018.t4.5)
hist(resid(mod2a))#more normal
qqnorm(resid(mod2a))
qqline(resid(mod2a))
anova(mod2a)
Anova(mod2a)
#says that this one is the best model, shows that depending on the trial, 
# I saw a different number of gobies on the reef, and that the effects of treatment
#on the number of gobies seen also vaired by trial

mod2b<-lmer(Density~Treatment*Day + Trial + (1|Treatment:Reef), data=den.2018.t4.5)
hist(resid(mod2b))#more normal
qqnorm(resid(mod2b))
qqline(resid(mod2b))
anova(mod2b)
Anova(mod2b)

AIC(mod2,mod2a, mod2b)
ancova(Density~Treatment*Day, data=den.2017)
plot(mod1)#again, looks like equal variance

#seems like the best model is the one with trial included as a fixed factor
#interesting that when I run as a repeated measures, no sig. effect of treament

#removed day 18, becuase it had some huge values, might have been a day where Hunter was recording??
gob.den<-read.csv("Data/density.2019.4.26.no.d18.t4.5.csv")
gob.den<- na.omit(gob.den)
den.2018.t4.5<-gob.den[(gob.den$Trial>3) & (gob.den$Trial<6), ]

mod2<-aov(Density~Treatment*Day, data=den.2018.t4.5)
hist(resid(mod2))#more normal
qqnorm(resid(mod2))
qqline(resid(mod2))
anova(mod2)
plot(mod2)#again, looks like equal variance
#TukeyHSD(mod2) #not working for this model for some reason
tapply(den.2018.t4.5$Density,den.2018.t4.5$Treatment,mean)
#    Low   Medium     High 
#8.364130 8.643243 6.693548 = 23% fewer fish seen on high-risk reefs


#2018 t4.5 plot#####
den.2018.t4.5$Treatment<-ordered(den.2018.t4.5$Treatment,levels=c("Low","Medium","High"))
viz.sur.2018.t4.5<-with(den.2018.t4.5, aggregate((Density), list(Day=Day,Treatment=Treatment), mean))
viz.sur.2018.t4.5$se<-with(den.2018.t4.5, aggregate((Density), list(Day=Day,Treatment=Treatment), 
                                          function(x) sd(x)/sqrt(length(x))))[,3]

png(filename = "Output/t4.5.2018.den.thesis.png", width = 1000, height = 500)#for presentations
png(filename = "Output/t4.5.2018.den.thesis..wide.png", width = 2000, height = 500)#for thesis (longer for bottom panel)

p <- ggplot(viz.sur.2018.t4.5, aes(x=Day, y=x, group=Treatment, color=Treatment))+ 
  geom_linerange(aes(ymin=x-se, ymax=x+se), 
                 position=position_dodge(0)) +
  geom_line(size=1.5)+
  geom_line(aes(linetype=Treatment)) +
  geom_point(aes(shape=Treatment))+
  labs(x="Day", y = "Number of Fish Seen")+
  theme_classic() + 
  theme(axis.text.x=element_text(size=28, colour="black"),axis.text.y=element_text(size=28, colour="black"), axis.title=element_text(size=35,face="bold")) +
  theme(axis.title.y = element_text(size= 34, margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)), axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0, 0),limits = c(0,18),breaks = c(4,8,12,16)) + scale_x_continuous(breaks=c(5, 10, 15, 20,25,30))+
  theme(text = element_text(family="Arial"))
p + scale_color_manual(values=c("#0072B2","#009E73","#D55E00")) + theme(legend.text=element_text(size=20)) + theme(legend.title =element_text(size=28, face="bold"))


dev.off()

#2018.t6#####

#2018.t6, adding trial as a fixed factor, and nesting reef within treatment as a random factor to account for repeated measures
#not surprising result, considering the fact that there's a big diff on day 4

mod1b<-lmer(Density~Treatment*Day + (1|Treatment:Reef), data=den.2018.t6)
hist(resid(mod1b))#more normal
qqnorm(resid(mod1b))
qqline(resid(mod1b))
anova(mod1b)

AIC(mod1,mod1a, mod1b)

#plotting
gob.den<-read.csv("Data/density.2019.4.26.no.d18.t4.5.csv")
gob.den<- na.omit(gob.den)
den.2018.t4.5<-gob.den[(gob.den$Trial>3) & (gob.den$Trial<6), ]

mod3<-aov(Density~Treatment*Day, data=den.2018.t6)
hist(resid(mod3))#more normal
qqnorm(resid(mod3))
qqline(resid(mod3))
anova(mod3)
plot(mod3)#again, looks like equal variance
TukeyHSD(mod3) #not working for this model for some reason
tapply(den.2018.t6$Density,den.2018.t6$Treatment,mean)
#   Low  Medium    High Control 
# 4.60    4.35    3.90    3.45 = significant interaction between day and treatment, driven by differences
# in high/control and low/med treatments, not sure what happened, because all cages were covered for first 24 hours


#2018 t6 plot#####
den.2018.t6$Treatment<-ordered(den.2018.t6$Treatment,levels=c("Low","Medium","High","Control"))
viz.sur.2018.t6<-with(den.2018.t6, aggregate((Density), list(Day=Day,Treatment=Treatment), mean))
viz.sur.2018.t6$se<-with(den.2018.t6, aggregate((Density), list(Day=Day,Treatment=Treatment), 
                                                    function(x) sd(x)/sqrt(length(x))))[,3]
viz.sur.2018.t6

png(filename = "Output/t6.2018.den.thesis.png", width = 1000, height = 500)

p <- ggplot(viz.sur.2018.t6, aes(x=Day, y=x, group=Treatment, color=Treatment))+ 
  geom_linerange(aes(ymin=x-se, ymax=x+se), 
                 position=position_dodge(0)) +
  geom_line(size=1.5)+
  geom_line(aes(linetype=Treatment)) +
  geom_point(aes(shape=Treatment))+
  labs(x="Day", y = "Number of Fish Seen")+
  theme_classic() + 
  theme(axis.text.x=element_text(size=28, colour="black"),axis.text.y=element_text(size=28, colour="black"), axis.title=element_text(size=35,face="bold")) +
  theme(axis.title.y = element_text(size= 34, margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)), axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0, 0),limits = c(0,18),breaks = c(4,8,12,16)) + scale_x_continuous(breaks=c(4, 7, 12, 14))+
  theme(text = element_text(family="Arial"))
p + scale_color_manual(values=c("#0072B2","#009E73","#D55E00","#899DA4")) + theme(legend.text=element_text(size=20)) + theme(legend.title =element_text(size=28, face="bold"))


dev.off()


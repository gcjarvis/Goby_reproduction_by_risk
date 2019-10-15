# Description: egg counts from ANCOVA for thesis
# Author: George C Jarvis
# Date: Sun Apr 28 11:49:51 2019
# --------------

rm(list=ls())

library(sciplot)
library(lme4)
library(car)
library(lmerTest)
library(dplyr)
library(plyr)
library(ggplot2)
library(MASS)
library(nlme)
library(pwr)
library(extrafont)

#importing datasets, needed to eliminate control manually
#for t1.2.3
egg.2017<-read.csv("Data/new.data.2019.4.23.t1.3.csv", na.strings = "") #uses adjusted counts for density
#for t4.5
egg.2018.t4.5<-read.csv("Data/new.data.2019.4.23.no.t6.csv", na.strings = "") #uses adjusted counts for density
View(egg.2018.t4.5)
#for t6
egg.2018.t6<-read.csv("Data/new.data.2019.4.23.csv", na.strings = "") #uses adjusted counts for density
egg.2018.t6<-egg.2018.t6[(egg.2018.t6$Trial>5), ]#subsetting t6 df, b/c lazy
mapvalues(egg.2018.t6, from = "Control", to = "Uncaged")

#plotting, not including anything about recollections here
#2017
egg.2017$Treatment<-ordered(egg.2017$Treatment,levels=c("Low","Medium","High"))
egg.anc.2017<-with(egg.2017, aggregate((Egg.count), list(Treatment=Treatment), mean))
egg.anc.2017$se<-with(egg.2017, aggregate((Egg.count), list(Treatment=Treatment), 
                                          function(x) sd(x)/sqrt(length(x))))[,2]

png(filename = "Output/2017.thesis.png", width = 700, height = 800)#thesis
png(filename = "Output/2017.bigger.thesis.png", width = 800, height = 700)#talk

egg.plot.2017<- ggplot(egg.anc.2017, aes(x=Treatment, y=x, fill=Treatment)) +
  geom_bar(stat="identity", colour= "black", width = 0.7, position="dodge")+ 
  scale_x_discrete(limits=c("Low","Medium","High"))+
  theme_classic() + theme(legend.position="none")  + #scale_fill_discrete(name="Sex",labels=c(" Male"," Female", " Transitional")) +
  theme(legend.key.size = unit(1.3,'line')) + 
  scale_fill_manual(values=c("#0072B2","#009E73","#D55E00")) + 
  theme(axis.text.x=element_text(size=32, colour="black"),axis.text.y=element_text(size=32, colour="black"), axis.title=element_text(size=37,face="bold")) +
  theme(axis.title.y = element_text(size= 37, margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 25, r = 0, b = 0, l = 0)), axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0, 0),limits = c(0, 4800))
egg.plot.2017 + geom_linerange(aes(ymin=x-se, ymax=x+se), size=0.5,   
                         position=position_dodge(.7)) + theme(text = element_text(family="Arial")) +
  labs(x="Risk Treatment", y="Number of Eggs Produced")

dev.off()

#2018.t4.5 
#note the change in y-axis scaling
egg.2018.t4.5$Treatment<-ordered(egg.2018.t4.5$Treatment,levels=c("Low","Medium","High"))
egg.anc.2018.t4.5<-with(egg.2018.t4.5, aggregate((Egg.count), list(Treatment=Treatment), mean))
egg.anc.2018.t4.5$se<-with(egg.2018.t4.5, aggregate((Egg.count), list(Treatment=Treatment), 
                                          function(x) sd(x)/sqrt(length(x))))[,2]

png(filename = "Output/2018.t4.5.thesis.png", width = 700, height = 800)#thesis
png(filename = "Output/2018.t4.5.bigger.thesis.png", width = 800, height = 700)#defense talk

egg.plot.2018<- ggplot(egg.anc.2018.t4.5, aes(x=Treatment, y=x, fill=Treatment)) +
  geom_bar(stat="identity", colour= "black", width = 0.7, position="dodge")+ 
  scale_x_discrete(limits=c("Low","Medium","High"))+
  theme_classic() + theme(legend.position="none")  + #scale_fill_discrete(name="Sex",labels=c(" Male"," Female", " Transitional")) +
  theme(legend.key.size = unit(1.3,'line')) + 
  scale_fill_manual(values=c("#0072B2","#009E73","#D55E00")) + 
  theme(axis.text.x=element_text(size=32, colour="black"),axis.text.y=element_text(size=32, colour="black"), axis.title=element_text(size=37,face="bold")) +
  theme(axis.title.y = element_text(size= 37, margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 25, r = 0, b = 0, l = 0)), axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0, 0),limits = c(0, 22000),breaks = c(5000,10000,15000,20000))
egg.plot.2018 + geom_linerange(aes(ymin=x-se, ymax=x+se), size=0.5,   
                               position=position_dodge(.7)) + theme(text = element_text(family="Arial")) +
  labs(x="Risk Treatment", y="Number of Eggs Produced")

dev.off()

#2018.t6
egg.2018.t6$Treatment<-ordered(egg.2018.t6$Treatment,levels=c("Low","Medium","High","Control"))
egg.anc.2018.t6<-with(egg.2018.t6, aggregate((Egg.count), list(Treatment=Treatment), mean))
egg.anc.2018.t6$se<-with(egg.2018.t6, aggregate((Egg.count), list(Treatment=Treatment), 
                                                    function(x) sd(x)/sqrt(length(x))))[,2]

png(filename = "Output/2018.t6.thesis.png", width = 700, height = 800)
png(filename = "Output/2018.t6.bigger.thesis.png", width = 800, height = 700)#talk

egg.plot.2018.t6<- ggplot(egg.anc.2018.t6, aes(x=Treatment, y=x, fill=Treatment)) +
  geom_bar(stat="identity", colour= "black", width = 0.7, position="dodge")+ 
  scale_x_discrete(limits=c("Low","Medium","High","Control"))+
  theme_classic() + theme(legend.position="none")  + #scale_fill_discrete(name="Sex",labels=c(" Male"," Female", " Transitional")) +
  theme(legend.key.size = unit(1.3,'line')) + 
  scale_fill_manual(values=c("#0072B2","#009E73","#D55E00", "#899DA4")) + 
  theme(axis.text.x=element_text(size=32, colour="black"),axis.text.y=element_text(size=32, colour="black"), axis.title=element_text(size=37,face="bold")) +
  theme(axis.title.y = element_text(size= 37, margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 25, r = 0, b = 0, l = 0)), axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0, 0),limits = c(0, 4800)) #,breaks = c(5000,10000,15000,20000))
egg.plot.2018.t6 + geom_linerange(aes(ymin=x-se, ymax=x+se), size=0.5,   
                               position=position_dodge(.7)) + theme(text = element_text(family="Arial")) +
  labs(x="Risk Treatment", y="Number of Eggs Produced")

dev.off()


####plots for thesis submission revisions

#control trt color code: , "#899DA4"  
  
#2017 separated by trial  
egg.2017$Treatment<-ordered(egg.2017$Treatment,levels=c("Low","Medium","High"))
egg.anc.2017<-with(egg.2017, aggregate((Egg.count), list(Treatment=Treatment), mean))
egg.anc.2017$se<-with(egg.2017, aggregate((Egg.count), list(Treatment=Treatment), 
                                          function(x) sd(x)/sqrt(length(x))))[,2]

egg.t1.3$Treatment<-ordered(egg.t1.3$Treatment, levels=c("Low","Medium","High"))

labels<-c("1" = "Trial 1", "2" = "Trial 2", "3" = "Trial 3")
labels

png(filename = "Output/2017.ancova.png", width = 800, height = 700)#thesis
#png(filename = "Output/2017.ancova.thesis.png", width = 800, height = 700)#talk
  
anc<-ggplot(egg.t1.3, aes(avg.inhab, Egg.count, color=Treatment, fill=Treatment)) +
  geom_smooth(method="lm", se=FALSE, show.legend = FALSE) +
  facet_grid(. ~ Trial,labeller=labeller(Trial = labels))  +
  geom_point(size=3) + 
  theme_classic() + 
  theme(strip.text = element_text(face="bold", size = 30))+
  xlab("Average Number of Gobies") +
  ylab("Reproduction per Reef") +
  expand_limits(y=0) +
  theme(axis.text.x=element_text(size=30, colour="black"),axis.text.y=element_text(size=30, colour="black"), axis.title=element_text(size=35,face="bold")) +
  theme(axis.title.y = element_text(size= 35, margin = margin(t = 0, r = 20, b = 0, l = 0)), 
        axis.title.x = element_text(margin = margin(t = 25, r = 0, b = 0, l = 0)), 
        axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))
anc + scale_color_manual(values=c("#0072B2","#009E73","#D55E00")) +
  theme(legend.text=element_text(size=25)) + 
  theme(legend.title =element_text(size=30, face="bold"))

dev.off()

#2018.t4.5 need to make a new labels object, to reflect only 2 trials? maybe can use same one

egg.t4.5$Treatment<-ordered(egg.t4.5$Treatment, levels=c("Low","Medium","High"))

labels.t4.5<-c("4" = "Trial 1", "5" = "Trial 2")
labels.t4.5

png(filename = "Output/2018.t.4.5.ancova.png", width = 900, height = 700)#thesis
#png(filename = "Output/2017.ancova.thesis.png", width = 800, height = 700)#talk

anc<-ggplot(egg.t4.5, aes(avg.inhab, Egg.count, color=Treatment, fill=Treatment)) +
  geom_smooth(method="lm", se=FALSE, show.legend = FALSE) +
  facet_grid(. ~ Trial,labeller=labeller(Trial = labels.t4.5))  +
  geom_point(size=3) + 
  theme_classic() + 
  theme(strip.text = element_text(face="bold", size = 30))+
  xlab("Average Number of Gobies") +
  ylab("Reproduction per Reef") +
  expand_limits(y=0) +
  theme(axis.text.x=element_text(size=30, colour="black"),axis.text.y=element_text(size=30, colour="black"), axis.title=element_text(size=35,face="bold")) +
  theme(axis.title.y = element_text(size= 35, margin = margin(t = 0, r = 20, b = 0, l = 0)), 
        axis.title.x = element_text(margin = margin(t = 25, r = 0, b = 0, l = 0)), 
        axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))
anc + scale_color_manual(values=c("#0072B2","#009E73","#D55E00")) +
  theme(legend.text=element_text(size=25)) + 
  theme(legend.title =element_text(size=30, face="bold"))

dev.off()

#not including trial for my written thesis, it would require me to go back and change
#all of my analyses and tables, etc. 

#2017 no trial, no faceting by trial
png(filename = "Output/2017.ancova.png", width = 800, height = 700)#thesis
#png(filename = "Output/2017.ancova.thesis.png", width = 800, height = 700)#talk

anc<-ggplot(egg.t1.3, aes(avg.inhab, Egg.count, color=Treatment, fill=Treatment)) +
  geom_smooth(method="lm", se=FALSE, show.legend = FALSE) +
  #facet_grid(. ~ Trial,labeller=labeller(Trial = labels))  +
  geom_point(size=3) + 
  theme_classic() + 
  #theme(strip.text = element_text(face="bold", size = 30))+
  xlab("Average Number of Gobies") +
  ylab("Reproduction per Reef") +
  expand_limits(y=0) +
  theme(axis.text.x=element_text(size=30, colour="black"),axis.text.y=element_text(size=30, colour="black"), axis.title=element_text(size=35,face="bold")) +
  theme(axis.title.y = element_text(size= 35, margin = margin(t = 0, r = 20, b = 0, l = 0)), 
        axis.title.x = element_text(margin = margin(t = 25, r = 0, b = 0, l = 0)), 
        axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))
anc + scale_color_manual(values=c("#0072B2","#009E73","#D55E00")) +
  theme(legend.text=element_text(size=25)) + 
  theme(legend.title =element_text(size=30, face="bold"))

dev.off()

#2018.t4.5 no facet, no trial

png(filename = "Output/2018.t.4.5.ancova.png", width = 800, height = 700)#thesis
#png(filename = "Output/2017.ancova.thesis.png", width = 800, height = 700)#talk

anc<-ggplot(egg.t4.5, aes(avg.inhab, Egg.count, color=Treatment, fill=Treatment)) +
  geom_smooth(method="lm", se=FALSE, show.legend = FALSE) +
  #facet_grid(. ~ Trial,labeller=labeller(Trial = labels.t4.5))  +
  geom_point(size=3) + 
  theme_classic() + 
  #theme(strip.text = element_text(face="bold", size = 30))+
  xlab("Average Number of Gobies") +
  ylab("Reproduction per Reef") +
  expand_limits(y=0) +
  theme(axis.text.x=element_text(size=30, colour="black"),axis.text.y=element_text(size=30, colour="black"), axis.title=element_text(size=35,face="bold")) +
  theme(axis.title.y = element_text(size= 35, margin = margin(t = 0, r = 20, b = 0, l = 0)), 
        axis.title.x = element_text(margin = margin(t = 25, r = 0, b = 0, l = 0)), 
        axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))
anc + scale_color_manual(values=c("#0072B2","#009E73","#D55E00")) +
  theme(legend.text=element_text(size=25)) + 
  theme(legend.title =element_text(size=30, face="bold"))

dev.off()

#2018.t6

egg.t6$Treatment<-ordered(egg.t6$Treatment, levels=c("Low","Medium","High","Control"))
library(plyr)
revalue(egg.t6$Treatment, c("Control"="Uncaged"))

mapvalues(x, from = c("beta", "gamma"), to = c("two", "three"))

png(filename = "Output/2018.t.6.ancova.png", width = 800, height = 700)#thesis
#png(filename = "Output/2017.ancova.thesis.png", width = 800, height = 700)#talk



anc<-ggplot(egg.t6, aes(avg.inhab, Egg.count, color=Treatment, fill=Treatment)) +
  geom_smooth(method="lm", se=FALSE, show.legend = FALSE) +
  #facet_grid(. ~ Trial,labeller=labeller(Trial = labels.t4.5))  +
  geom_point(size=3) + 
  theme_classic() + 
  #theme(strip.text = element_text(face="bold", size = 30))+
  xlab("Average Number of Gobies") +
  ylab("Reproduction per Reef") +
  expand_limits(y=0) +
  theme(axis.text.x=element_text(size=30, colour="black"),axis.text.y=element_text(size=30, colour="black"), axis.title=element_text(size=35,face="bold")) +
  theme(axis.title.y = element_text(size= 35, margin = margin(t = 0, r = 20, b = 0, l = 0)), 
        axis.title.x = element_text(margin = margin(t = 25, r = 0, b = 0, l = 0)), 
        axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))
anc + scale_color_manual(values=c("#0072B2","#009E73","#D55E00","#899DA4")) +
  theme(legend.text=element_text(size=25)) + 
  theme(legend.title =element_text(size=30, face="bold"))

dev.off()
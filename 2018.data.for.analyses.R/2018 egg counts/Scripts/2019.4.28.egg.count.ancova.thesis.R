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


#plotting, not including anything baout recollections here
#2017
egg.2017$Treatment<-ordered(egg.2017$Treatment,levels=c("Low","Medium","High"))
egg.anc.2017<-with(egg.2017, aggregate((Egg.count), list(Treatment=Treatment), mean))
egg.anc.2017$se<-with(egg.2017, aggregate((Egg.count), list(Treatment=Treatment), 
                                          function(x) sd(x)/sqrt(length(x))))[,2]

png(filename = "Output/2017.thesis.png", width = 700, height = 800)

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

png(filename = "Output/2018.t4.5.thesis.png", width = 700, height = 800)

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

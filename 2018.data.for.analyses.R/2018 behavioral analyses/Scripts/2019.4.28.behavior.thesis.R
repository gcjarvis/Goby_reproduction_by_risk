# Description: have to make a bunch of barplots for behavior analyses
# Author: George C Jarvis
# Date: Sun Apr 28 15:01:32 2019
# --------------

#note: only analyzed behaviors for exp 1 and and 2

rm(list=ls())

library(extrafont)
library(ggplot2)
library(sciplot)
library(dplyr)
library(plyr)
library(vegan)
library(psych)
library(lme4)

#setting up df's
# trials 1-3
behave<-read.csv("Data/2019.4.25.behavior.ss.csv") #doesn't include unnecessary NA's
b.2017.t1.2.3<-behave[(behave$Trial<4),]
b.2017.t1.2.3$Treatment<-ordered(b.2017.t1.2.3$Treatment,levels=c("Low","Medium","High"))

#trials 4 and 5
b.2018.t4.5<-behave[(behave$Trial>3) & (behave$Trial<6), ]
b.2018.t4.5$Treatment<-ordered(b.2018.t4.5$Treatment,levels=c("Low","Medium","High"))

df<-b.2017.t1.2.3
df<-b.2018.t4.5

#exposure time, plot proportions on a single axis, going to just change the df. no time to code it correctly
exp<-with(df, aggregate((proportion.exposed), list(Treatment=Treatment), mean))
exp$se<-with(df, aggregate((proportion.exposed), list(Treatment=Treatment), 
                                          function(x) sd(x)/sqrt(length(x))))[,2]

png(filename = "Output/2018.t.4.5.exp.thesis.png", width = 700, height = 800)

exp.plot<- ggplot(exp, aes(x=Treatment, y=x, fill=Treatment)) +
  geom_bar(stat="identity", colour= "black", width = 0.7, position="dodge")+ 
  scale_x_discrete(limits=c("Low","Medium","High"))+
  theme_classic() + theme(legend.position="none")  + #scale_fill_discrete(name="Sex",labels=c(" Male"," Female", " Transitional")) +
  theme(legend.key.size = unit(1.3,'line')) + 
  scale_fill_manual(values=c("#0072B2","#009E73","#D55E00")) + 
  theme(axis.text.x=element_text(size=32, colour="black"),axis.text.y=element_text(size=32, colour="black"), axis.title=element_text(size=37,face="bold")) +
  theme(axis.title.y = element_text(size= 37, margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 25, r = 0, b = 0, l = 0)), axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0, 0),limits = c(0, 0.85), breaks = c(0,0.2,0.4,0.6,0.8))
exp.plot + geom_linerange(aes(ymin=x-se, ymax=x+se), size=0.5,   
                               position=position_dodge(.7)) + theme(text = element_text(family="Arial")) +
  labs(x="Risk Treatment", y="Proportion of Time Exposed")

dev.off()

# distance moved

df<-b.2017.t1.2.3
df<-b.2018.t4.5

#exposure time, plot proportions on a single axis, going to just change the df. no time to code it correctly
exp<-with(df, aggregate((total.dist.moved), list(Treatment=Treatment), mean))
exp$se<-with(df, aggregate((total.dist.moved), list(Treatment=Treatment), 
                           function(x) sd(x)/sqrt(length(x))))[,2]

png(filename = "Output/2017.dist.thesis.png", width = 700, height = 800)

exp.plot<- ggplot(exp, aes(x=Treatment, y=x, fill=Treatment)) +
  geom_bar(stat="identity", colour= "black", width = 0.7, position="dodge")+ 
  scale_x_discrete(limits=c("Low","Medium","High"))+
  theme_classic() + theme(legend.position="none")  + #scale_fill_discrete(name="Sex",labels=c(" Male"," Female", " Transitional")) +
  theme(legend.key.size = unit(1.3,'line')) + 
  scale_fill_manual(values=c("#0072B2","#009E73","#D55E00")) + 
  theme(axis.text.x=element_text(size=32, colour="black"),axis.text.y=element_text(size=32, colour="black"), axis.title=element_text(size=37,face="bold")) +
  theme(axis.title.y = element_text(size= 37, margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 25, r = 0, b = 0, l = 0)), axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0, 0),limits = c(0, 280), breaks = c(50,100,150,200,250))
exp.plot + geom_linerange(aes(ymin=x-se, ymax=x+se), size=0.5,   
                          position=position_dodge(.7)) + theme(text = element_text(family="Arial")) +
  labs(x="Risk Treatment", y="Total Distance Moved (cm)")

dev.off()


# foraging rates

df<-b.2017.t1.2.3
df<-b.2018.t4.5

exp<-with(df, aggregate((bites.min), list(Treatment=Treatment), mean))
exp$se<-with(df, aggregate((bites.min), list(Treatment=Treatment), 
                           function(x) sd(x)/sqrt(length(x))))[,2]

png(filename = "Output/2018.bites.thesis.png", width = 700, height = 800)

exp.plot<- ggplot(exp, aes(x=Treatment, y=x, fill=Treatment)) +
  geom_bar(stat="identity", colour= "black", width = 0.7, position="dodge")+ 
  scale_x_discrete(limits=c("Low","Medium","High"))+
  theme_classic() + theme(legend.position="none")  + #scale_fill_discrete(name="Sex",labels=c(" Male"," Female", " Transitional")) +
  theme(legend.key.size = unit(1.3,'line')) + 
  scale_fill_manual(values=c("#0072B2","#009E73","#D55E00")) + 
  theme(axis.text.x=element_text(size=32, colour="black"),axis.text.y=element_text(size=32, colour="black"), axis.title=element_text(size=37,face="bold")) +
  theme(axis.title.y = element_text(size= 37, margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 25, r = 0, b = 0, l = 0)), axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0, 0),limits = c(0, 1.2), breaks = c(0,0.5,1.0))
exp.plot + geom_linerange(aes(ymin=x-se, ymax=x+se), size=0.5,   
                          position=position_dodge(.7)) + theme(text = element_text(family="Arial")) +
  labs(x="Risk Treatment", y="Bites per Minute")

dev.off()

# courtship

df<-b.2017.t1.2.3
df<-b.2018.t4.5

exp<-with(df, aggregate((courtship.min), list(Treatment=Treatment), mean))
exp$se<-with(df, aggregate((courtship.min), list(Treatment=Treatment), 
                           function(x) sd(x)/sqrt(length(x))))[,2]

png(filename = "Output/2018.court.min.thesis.png", width = 700, height = 800)

exp.plot<- ggplot(exp, aes(x=Treatment, y=x, fill=Treatment)) +
  geom_bar(stat="identity", colour= "black", width = 0.7, position="dodge")+ 
  scale_x_discrete(limits=c("Low","Medium","High"))+
  theme_classic() + theme(legend.position="none")  + #scale_fill_discrete(name="Sex",labels=c(" Male"," Female", " Transitional")) +
  theme(legend.key.size = unit(1.3,'line')) + 
  scale_fill_manual(values=c("#0072B2","#009E73","#D55E00")) + 
  theme(axis.text.x=element_text(size=32, colour="black"),axis.text.y=element_text(size=32, colour="black"), axis.title=element_text(size=37,face="bold")) +
  theme(axis.title.y = element_text(size= 37, margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 25, r = 0, b = 0, l = 0)), axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0, 0),limits = c(0, 0.18), breaks = c(0,0.05,0.1,0.15))
exp.plot + geom_linerange(aes(ymin=x-se, ymax=x+se), size=0.5,   
                          position=position_dodge(.7)) + theme(text = element_text(family="Arial")) +
  labs(x="Risk Treatment", y="Courtship Behaviors per Minute")

dev.off()

# movements per minute

df<-b.2017.t1.2.3
df<-b.2018.t4.5

exp<-with(df, aggregate((movements.min), list(Treatment=Treatment), mean))
exp$se<-with(df, aggregate((movements.min), list(Treatment=Treatment), 
                           function(x) sd(x)/sqrt(length(x))))[,2]

png(filename = "Output/2018.movement.rate.thesis.png", width = 700, height = 800)

exp.plot<- ggplot(exp, aes(x=Treatment, y=x, fill=Treatment)) +
  geom_bar(stat="identity", colour= "black", width = 0.7, position="dodge")+ 
  scale_x_discrete(limits=c("Low","Medium","High"))+
  theme_classic() + theme(legend.position="none")  + #scale_fill_discrete(name="Sex",labels=c(" Male"," Female", " Transitional")) +
  theme(legend.key.size = unit(1.3,'line')) + 
  scale_fill_manual(values=c("#0072B2","#009E73","#D55E00")) + 
  theme(axis.text.x=element_text(size=32, colour="black"),axis.text.y=element_text(size=32, colour="black"), axis.title=element_text(size=37,face="bold")) +
  theme(axis.title.y = element_text(size= 37, margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 25, r = 0, b = 0, l = 0)), axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0, 0),limits = c(0, 1.8), breaks = c(0,0.5,1.0,1.5))
exp.plot + geom_linerange(aes(ymin=x-se, ymax=x+se), size=0.5,   
                          position=position_dodge(.7)) + theme(text = element_text(family="Arial")) +
  labs(x="Risk Treatment", y="Movements per Minute")

dev.off()
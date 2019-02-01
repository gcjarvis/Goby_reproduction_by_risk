#egg counts

rm(list=ls())

library(sciplot)
library(lme4)
library(lmerTest)
library(car)
library(dplyr)
library(ggplot2)
library(extrafont)
library(plyr)

gob.sub<-read.csv("Data/10.29.18.eggs.csv")
gob.sub<-read.csv("Data/10.29.18.eggs.outliers.a.csv") #first cut out outliers removed
gob.sub<-read.csv("Data/10.29.18.eggs.outlier.b.csv") #second cut of outliers removed
gob.sub<-read.csv("Data/10.29.18.eggs.3.weeks.csv")
#gob.sub$week<-as.factor(gob.sub$week)#need to make week a factor???
View(gob.sub)
head(gob.sub)
gob.sub$week<-as.factor(gob.sub$week)
gob.sub$Treatment<-ordered(gob.sub$Treatment,levels=c("Low","Medium","High"))

egg.means<-with(gob.sub, aggregate((egg.count), list(Treatment=Treatment,week=week), mean))
egg.means
#now apply the se function to the 4th column [,3]
egg.means$se<-with(gob.sub, aggregate((egg.count), list(Treatment=Treatment,week=week), function(x) sd(x)/sqrt(length(x))))[,3]
egg.means

reco.plot<- ggplot(reco.means, aes(x=treatment, y=x, fill=Sex)) +
  geom_bar(stat="identity", colour= "black", width = 0.85, position="dodge")+ 
  scale_x_discrete(limits=c("Low","Medium","High"))+
  theme_classic() + theme(legend.position=c(0.82,0.82)) + #scale_fill_discrete(name="Sex",labels=c(" Male"," Female", " Transitional")) +
  theme(legend.key.size = unit(1.3,'line')) + 
  theme(legend.title=element_text(size=20) , legend.text=element_text(size=18)) + scale_fill_manual(values=c("#0072B2","#D55E00","#009E73")) + 
  theme(axis.text.x=element_text(size=20, colour="black"),axis.text.y=element_text(size=20, colour="black"), axis.title=element_text(size=25,face="bold")) +
  theme(axis.title.y = element_text(size= 25, margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 25, r = 0, b = 0, l = 0)), axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0, 0))
reco.plot + geom_linerange(aes(ymin=x-se, ymax=x+se), size=0.5,   
                           position=position_dodge(.85)) + theme(text = element_text(family="Arial")) +
  labs(x="Risk Treatment", y="Number of Gobies Recollected")

egg.wsn<- ggplot(egg.means, aes(x=week, y=x, fill=Treatment)) +
  geom_bar(stat="identity", colour= "black", width = 0.7, position="dodge")+ 
  scale_x_discrete(limits=c("1","2","3","4"))+
  theme_classic() + theme(legend.position="right")  + #scale_fill_discrete(name="Sex",labels=c(" Male"," Female", " Transitional")) +
  theme(legend.key.size = unit(1.3,'line')) + 
  theme(legend.title=element_text(size=20) , legend.text=element_text(size=18)) + scale_fill_manual(values=c("#0072B2","#009E73","#D55E00")) + 
  theme(axis.text.x=element_text(size=20, colour="black"),axis.text.y=element_text(size=20, colour="black"), axis.title=element_text(size=25,face="bold")) +
  theme(axis.title.y = element_text(size= 25, margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 25, r = 0, b = 0, l = 0)), axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0, 0))
egg.wsn + geom_linerange(aes(ymin=x-se, ymax=x+se), size=0.5,   
                           position=position_dodge(.7)) + theme(text = element_text(family="Arial")) +
  labs(x="Week", y="Total Reproduction per Reef")

#with just 3 weeks of data
egg.wsn<- ggplot(egg.means, aes(x=week, y=x, fill=treatment)) +
  geom_bar(stat="identity", colour= "black", width = 0.85, position="dodge")+ 
  scale_x_discrete(limits=c("1","2","3"))+
  theme_classic() + theme(legend.position=c(0.82,0.82)) + #scale_fill_discrete(name="Sex",labels=c(" Male"," Female", " Transitional")) +
  theme(legend.key.size = unit(1.3,'line')) + 
  theme(legend.title=element_text(size=20) , legend.text=element_text(size=18)) + scale_fill_manual(values=c("#0072B2","#D55E00","#009E73")) + 
  theme(axis.text.x=element_text(size=20, colour="black"),axis.text.y=element_text(size=20, colour="black"), axis.title=element_text(size=25,face="bold")) +
  theme(axis.title.y = element_text(size= 25, margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 25, r = 0, b = 0, l = 0)), axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0, 0))
egg.wsn + geom_linerange(aes(ymin=x-se, ymax=x+se), size=0.5,   
                         position=position_dodge(.85)) + theme(text = element_text(family="Arial")) +
  labs(x="Week", y="Total Reproduction per Reef")


########visual surveys

viz.sur<-with(goby.den.survey, aggregate((Proportion_pop), list(Days_after_deploy=Days_after_deploy,Treatment=Treatment), mean))
#now apply the se function to the 4th column [,3]
viz.sur$se<-with(goby.den.survey, aggregate((Proportion_pop), list(Days_after_deploy=Days_after_deploy,Treatment=Treatment), function(x) sd(x)/sqrt(length(x))))[,3]
viz.sur

p <- ggplot(viz.sur, aes(x=Days_after_deploy, y=x, group=Treatment, color=Treatment))+ 
  geom_linerange(aes(ymin=x-se, ymax=x+se), 
                 position=position_dodge(0)) +
  geom_line(size=0.75)+
  geom_line(aes(linetype=Treatment)) +
  geom_point(aes(shape=Treatment))+
  labs(x="Days After Deployment", y = "Porportion Observed")+
  theme_classic() + 
  theme(axis.text.x=element_text(size=20, colour="black"),axis.text.y=element_text(size=16, colour="black"), axis.title=element_text(size=22,face="bold")) +
  theme(axis.title.y = element_text(size= 22, margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)), axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0, 0)) +
  theme(text = element_text(family="Arial"))
p + scale_color_manual(values=c("#D55E00","#0072B2","#009E73")) + theme(legend.text=element_text(size=15)) + theme(legend.title =element_text(size=15, face="bold"))
# Description: plots for recollections, colored
# Author: George C Jarvis
# Date: Mon Apr 01 09:26:00 2019
# --------------


rm(list=ls())

#load packages
library(sciplot)
library(lme4)
library(lmerTest)
library(car)
library(plyr)
library(dplyr)
library(ggplot2)
library(extrafont)
library(tidyr)
library(wesanderson)

getwd()
#importing raw data
#tried to make the .csv file not contain any NA's 
reco.raw<-read.csv("Data/Recollections.all.trials.2.23.19.csv")
View(reco.raw)

#recollection by sex for each reef, includes immatures, maybe will want to
#--analyze that at some point
reco.wrang<-as.data.frame(reco.raw) %>%
  tidyr::gather(key = Sex.recollected, value = Count,-Year, -Trial, -Trial.duration,-Reef,-Treatment,-Initial.size,-Final.size,-Initial.sex,-Growth,-Final.sex,
                -Final.sex.female.transitional.equals.female,-final.Weight.g,-Deployment.day,-Immature)

#removing immatures from the dataset, because they weren't recollected
reco.wrang<-as.data.frame(reco.raw) %>%
  tidyr::gather(key = Sex.recollected, value = Count,-Year, -Trial, -Trial.duration,-Reef,-Treatment,-Initial.size,-Final.size,-Initial.sex,-Growth,-Final.sex,
                -Final.sex.female.transitional.equals.female,-final.Weight.g,-Deployment.day,-Immature)
View(reco.wrang)
#now just have to go back and count all of the fish per sex for each reef/trial

#tibl that includes sex.recollected + dep. day
reco.reef<-reco.wrang %>%
  group_by(Year,Trial,Reef,Treatment,Sex.recollected,Deployment.day) %>%
  summarize(Count = sum(Count))
View(reco.reef)

#simpler tibl that excludes sex and dep. day
reco.reef.simple<-reco.wrang %>%
  group_by(Year,Trial,Reef,Treatment) %>%
  summarize(Count = sum(Count))

#df exported to combine recollections with egg counts
write.csv(reco.reef.simple,"2019.4.23.recollection.data.csv")

View(reco.reef.simple)
#calculating proportion of initial population recollected
reco.reef.simple$Prop.reco<-(reco.reef.simple$Count/20)
View(reco.reef.simple)

#subsetting by trial
reco.t1.2.3<-reco.reef.simple[reco.reef.simple$Trial<4,]
reco.t4.5<-reco.reef.simple[(reco.reef.simple$Trial>3) & (reco.reef.simple$Trial<6), ]
reco.t6<-reco.reef.simple[reco.reef.simple$Trial==6,]

df<-reco.t1.2.3
df<-reco.t4.5
df<-reco.t6
View(df)

#models
mod1<-lm(Prop.reco~Treatment,data=df)
hist(resid(mod1))#not normal
qqnorm(resid(mod1))
qqline(resid(mod1))
anova(mod1)
boxplot(Prop.reco~Treatment,data=df)

#plotting#####

#trials 1-5 (only three treatment levels)
reco.means<-with(df, aggregate((Prop.reco), list(Treatment=Treatment), mean))
#reco.means
#now apply the se function to the 4th column [,3]
reco.means$se<-with(df, aggregate((Prop.reco), list(Treatment=Treatment), function(x) sd(x)/sqrt(length(x))))[,2]
#reco.means
reco.means$Treatment<-ordered(reco.means$Treatment,levels=c("Low","Medium","High"))

#t1.2.3
png(filename = "Output/reco.t.1.2.3.png", width = 700, height = 800)

reco.plot<- ggplot(reco.means, aes(x=Treatment, y=x, fill=Treatment)) +
  geom_bar(stat="identity", colour= "black", width = 0.85, position="dodge")+ 
  scale_x_discrete(limits=c("Low","Medium","High"))+
  theme_classic() + theme(legend.position="none") + #scale_fill_discrete(limits=c("Low","Medium", "High")) +
  #theme(legend.key.size = unit(1.3,'line')) + 
  #theme(legend.title=element_text(size=34) , legend.text=element_text(size=20)) 
  scale_fill_manual(values=c("#0072B2","#009E73","#D55E00")) + 
  theme(axis.text.x=element_text(size=32, colour="black"),axis.text.y=element_text(size=32, colour="black"), axis.title=element_text(size=37,face="bold")) +
  theme(axis.title.y = element_text(size= 37, margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 25, r = 0, b = 0, l = 0)), axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0, 0),limits = c(0,0.45))
reco.plot + geom_linerange(aes(ymin=x-se, ymax=x+se), size=0.5,   
                           position=position_dodge(.85)) + theme(text = element_text(family="Arial")) +
  labs(x="Risk Treatment", y="Proportion of Gobies Recollected")

dev.off()

#t4.5
png(filename = "Output/reco.t.4.5.png", width = 700, height = 800)

reco.plot<- ggplot(reco.means, aes(x=Treatment, y=x, fill=Treatment)) +
  geom_bar(stat="identity", colour= "black", width = 0.85, position="dodge")+ 
  scale_x_discrete(limits=c("Low","Medium","High"))+
  theme_classic() + theme(legend.position="none") + #scale_fill_discrete(limits=c("Low","Medium", "High")) +
  #theme(legend.key.size = unit(1.3,'line')) + 
  #theme(legend.title=element_text(size=34) , legend.text=element_text(size=20)) 
  scale_fill_manual(values=c("#0072B2","#009E73","#D55E00")) + 
  theme(axis.text.x=element_text(size=32, colour="black"),axis.text.y=element_text(size=32, colour="black"), axis.title=element_text(size=37,face="bold")) +
  theme(axis.title.y = element_text(size= 37, margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 25, r = 0, b = 0, l = 0)), axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0, 0),limits = c(0,0.45))
reco.plot + geom_linerange(aes(ymin=x-se, ymax=x+se), size=0.5,   
                           position=position_dodge(.85)) + theme(text = element_text(family="Arial")) +
  labs(x="Risk Treatment", y="Proportion of Gobies Recollected")

dev.off()

#trial 6
reco.means<-with(df, aggregate((Prop.reco), list(Treatment=Treatment), mean))
#reco.means
#now apply the se function to the 4th column [,3]
reco.means$se<-with(df, aggregate((Prop.reco), list(Treatment=Treatment), function(x) sd(x)/sqrt(length(x))))[,2]
#reco.means
reco.means$Treatment<-ordered(reco.means$Treatment,levels=c("Low","Medium","High","Control"))

#for presentation...maybe...for now I think I like the idea of putting them all on the same axis
png(filename = "Output/reco.t.6.pres.png", width = 700, height = 800)

reco.plot<- ggplot(reco.means, aes(x=Treatment, y=x, fill=Treatment)) +
  geom_bar(stat="identity", colour= "black", width = 0.85, position="dodge")+ 
  scale_x_discrete(limits=c("Low","Medium","High","Control"))+
  theme_classic() + theme(legend.position="none") + #scale_fill_discrete(limits=c("Low","Medium", "High")) +
  #theme(legend.key.size = unit(1.3,'line')) + 
  #theme(legend.title=element_text(size=34) , legend.text=element_text(size=20)) 
  scale_fill_manual(values=c("#0072B2","#009E73","#D55E00","#35274A")) + 
  theme(axis.text.x=element_text(size=32, colour="black"),axis.text.y=element_text(size=32, colour="black"), axis.title=element_text(size=37,face="bold")) +
  theme(axis.title.y = element_text(size= 37, margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 25, r = 0, b = 0, l = 0)), axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0, 0))
reco.plot + geom_linerange(aes(ymin=x-se, ymax=x+se), size=0.5,   
                           position=position_dodge(.85)) + theme(text = element_text(family="Arial")) +
  labs(x="Risk Treatment", y="Proportion of Gobies Recollected")

dev.off()

#for thesis (plotting all figures on same scale with one y-axis)
png(filename = "Output/reco.t.6.thesis.png", width = 700, height = 800)

reco.plot<- ggplot(reco.means, aes(x=Treatment, y=x, fill=Treatment)) +
  geom_bar(stat="identity", colour= "black", width = 0.85, position="dodge")+ 
  scale_x_discrete(limits=c("Low","Medium","High","Control"))+
  theme_classic() + theme(legend.position="none") + #scale_fill_discrete(limits=c("Low","Medium", "High")) +
  #theme(legend.key.size = unit(1.3,'line')) + 
  #theme(legend.title=element_text(size=34) , legend.text=element_text(size=20)) 
  scale_fill_manual(values=c("#0072B2","#009E73","#D55E00","#899DA4")) + 
  theme(axis.text.x=element_text(size=32, colour="black"),axis.text.y=element_text(size=32, colour="black"), axis.title=element_text(size=37,face="bold")) +
  theme(axis.title.y = element_text(size= 37, margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 25, r = 0, b = 0, l = 0)), axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0, 0),limits=c(0,0.45))
reco.plot + geom_linerange(aes(ymin=x-se, ymax=x+se), size=0.5,   
                           position=position_dodge(.85)) + theme(text = element_text(family="Arial")) +
  labs(x="Risk Treatment", y="Proportion of Gobies Recollected")

dev.off()
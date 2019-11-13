# Description: figures for MS, b&w
# Author: George C Jarvis
# Date: Wed Nov 13 21:19:41 2019
# Notes: will lay this out as a 5-paneled figure in PPT, aligning the x axes
#       but having different scales on the y axes
# --------------

#clear workspace
rm(list=ls())

#library
library(ggplot2)

#loading data
behave<-read.csv("Data/2019.10.25.behavior.includes.recollections.csv")
behave<-na.omit(behave)
behave$avg.inhab<-(ceiling((behave$Recollection+20)/2))
behave$Treatment<-ordered(behave$Treatment,levels=c("Low","Medium","High"))

#plotting
#have to set up different df for each behavior, shouldn't be too hard

#exposure time####
exp<-with(behave, aggregate((proportion.exposed), list(Treatment=Treatment), mean))
exp$se<-with(behave, aggregate((proportion.exposed), list(Treatment=Treatment), 
                           function(x) sd(x)/sqrt(length(x))))[,2]

png("Output/2019.11.13.exposure.9.5x5.5.300dpi.png", width = 6.5, height = 5.5, units = 'in', res = 300)

reco.plot<- ggplot(exp, aes(x=Treatment, y=x, fill=Treatment)) +
  geom_bar(stat="identity", colour= "black", width = 0.5, position="dodge")+ 
  scale_x_discrete(limits=c("Low","Medium","High"))+
  theme_classic() + 
  labs(x="Risk Treatment", y="Proportion of Time Exposed") +
  theme(legend.position="none") + 
  scale_fill_manual(values=c("grey", "grey", "grey")) +
  theme(axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"), 
        axis.title=element_text(size=20))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), 
        axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0)), 
        axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(legend.text=element_text(size=18)) +
  theme(legend.title =element_text(size=20))+
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0,0),limits = c(0,1.015))
reco.plot + geom_linerange(aes(ymin=x-se, ymax=x+se), size=0.5,   
                           position=position_dodge(.85)) + theme(text = element_text(family="Arial"))
#not sure why it's off the x axis now...has something to do with breaks
dev.off()

#total distance moved####
td<-with(behave, aggregate((total.dist.moved), list(Treatment=Treatment), mean))
td$se<-with(behave, aggregate((total.dist.moved), list(Treatment=Treatment), 
                               function(x) sd(x)/sqrt(length(x))))[,2]

png("Output/2019.11.13.td.moved.9.5x5.5.300dpi.png", width = 6.5, height = 5.5, units = 'in', res = 300)

reco.plot<- ggplot(td, aes(x=Treatment, y=x, fill=Treatment)) +
  geom_bar(stat="identity", colour= "black", width = 0.5, position="dodge")+ 
  scale_x_discrete(limits=c("Low","Medium","High"))+
  theme_classic() + 
  labs(x="Risk Treatment", y="Linear Distance Traveled (mm)") +
  theme(legend.position="none") + 
  scale_fill_manual(values=c("grey", "grey", "grey")) +
  theme(axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"), 
        axis.title=element_text(size=20))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), 
        axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0)), 
        axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(legend.text=element_text(size=18)) +
  theme(legend.title =element_text(size=20))+
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0,0),limits = c(0,205))
reco.plot + geom_linerange(aes(ymin=x-se, ymax=x+se), size=0.5,   
                           position=position_dodge(.85)) + theme(text = element_text(family="Arial"))
#not sure why it's off the x axis now...has something to do with breaks
dev.off()

#movements/min####
mm<-with(behave, aggregate((movements.min), list(Treatment=Treatment), mean))
mm$se<-with(behave, aggregate((movements.min), list(Treatment=Treatment), 
                              function(x) sd(x)/sqrt(length(x))))[,2]

png("Output/2019.11.13.movement.swims.9.5x5.5.300dpi.png", width = 6.5, height = 5.5, units = 'in', res = 300)

reco.plot<- ggplot(mm, aes(x=Treatment, y=x, fill=Treatment)) +
  geom_bar(stat="identity", colour= "black", width = 0.5, position="dodge")+ 
  scale_x_discrete(limits=c("Low","Medium","High"))+
  theme_classic() + 
  labs(x="Risk Treatment", y = bquote('Movement Rate'
(movements~min^-1)))+ # y="Movements per Minute") + #can't figure it out how to get parenthesis with superscript below axis label
  theme(legend.position="none") + 
  scale_fill_manual(values=c("grey", "grey", "grey")) +
  theme(axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"), 
        axis.title=element_text(size=20))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), 
        axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0)), 
        axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(legend.text=element_text(size=18)) +
  theme(legend.title =element_text(size=20))+
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0,0),limits = c(0,1.51))
reco.plot + geom_linerange(aes(ymin=x-se, ymax=x+se), size=0.5,   
                           position=position_dodge(.85)) + theme(text = element_text(family="Arial"))
#not sure why it's off the x axis now...has something to do with breaks
dev.off()

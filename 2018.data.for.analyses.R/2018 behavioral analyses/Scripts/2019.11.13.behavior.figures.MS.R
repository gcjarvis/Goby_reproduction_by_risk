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

png("Output/2019.2.1.exposure.9.5x5.5.300dpi.png", width = 6.5, height = 5.5, units = 'in', res = 300)

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
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0,0),limits = c(0,0.807))
reco.plot + geom_linerange(aes(ymin=x-se, ymax=x+se), size=0.5,   
                           position=position_dodge(.85)) + theme(text = element_text(family="Arial"))
#not sure why it's off the x axis now...has something to do with breaks
dev.off()

#this was the only behavior with a sig. trt. effect
#so it's likely the only one I will report percentages for in text


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

#foraging (bites per minute)####
fr<-with(behave, aggregate((bites.min), list(Treatment=Treatment), mean))
fr$se<-with(behave, aggregate((bites.min), list(Treatment=Treatment), 
                              function(x) sd(x)/sqrt(length(x))))[,2]

png("Output/2019.11.14.Bite.rate.no.parentheses.9.5x5.5.300dpi.png", width = 6.5, height = 5.5, units = 'in', res = 300)

reco.plot<- ggplot(fr, aes(x=Treatment, y=x, fill=Treatment)) +
  geom_bar(stat="identity", colour= "black", width = 0.5, position="dodge")+ 
  scale_x_discrete(limits=c("Low","Medium","High"))+
  theme_classic() + 
  labs(x="Risk Treatment", y = "Bites per Minute")+
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
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0,0),limits = c(0,1.1))
reco.plot + geom_linerange(aes(ymin=x-se, ymax=x+se), size=0.5,   
                           position=position_dodge(.85)) + theme(text = element_text(family="Arial"))
dev.off()

#movements per minute with rate in parentheses

png("Output/2019.11.14.bite.rate.with.parentheses.9.5x5.5.300dpi.png", width = 6.5, height = 5.5, units = 'in', res = 300)

reco.plot<- ggplot(fr, aes(x=Treatment, y=x, fill=Treatment)) +
  geom_bar(stat="identity", colour= "black", width = 0.5, position="dodge")+ 
  scale_x_discrete(limits=c("Low","Medium","High"))+
  theme_classic() + 
  labs(x="Risk Treatment",y=(expression(atop("Foraging", 
                                             paste((bites~min^-1))))))+
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
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0,0),limits = c(0,1.01))
reco.plot + geom_linerange(aes(ymin=x-se, ymax=x+se), size=0.5,   
                           position=position_dodge(.85)) + theme(text = element_text(family="Arial"))

dev.off()

#courtship rate (displays per minute)####
cr<-with(behave, aggregate((courtship.min), list(Treatment=Treatment), mean))
cr$se<-with(behave, aggregate((courtship.min), list(Treatment=Treatment), 
                              function(x) sd(x)/sqrt(length(x))))[,2]

png("Output/2019.11.14.courtship.rate.no.parentheses.9.5x5.5.300dpi.png", width = 6.5, height = 5.5, units = 'in', res = 300)

reco.plot<- ggplot(cr, aes(x=Treatment, y=x, fill=Treatment)) +
  geom_bar(stat="identity", colour= "black", width = 0.5, position="dodge")+ 
  scale_x_discrete(limits=c("Low","Medium","High"))+
  theme_classic() + 
  labs(x="Risk Treatment", y = "Courtship Displays per Minute")+
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
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0,0),limits = c(0,1.5))
reco.plot + geom_linerange(aes(ymin=x-se, ymax=x+se), size=0.5,   
                           position=position_dodge(.85)) + theme(text = element_text(family="Arial"))
dev.off()

#movements per minute with rate in parentheses

png("Output/2019.11.14.courtship.rate.with.parentheses.9.5x5.5.300dpi.png", width = 6.5, height = 5.5, units = 'in', res = 300)
png("Output/2019.2.1.interactions.with.parentheses.9.5x5.5.300dpi.png", width = 6.5, height = 5.5, units = 'in', res = 300)

reco.plot<- ggplot(cr, aes(x=Treatment, y=x, fill=Treatment)) +
  geom_bar(stat="identity", colour= "black", width = 0.5, position="dodge")+ 
  scale_x_discrete(limits=c("Low","Medium","High"))+
  theme_classic() + 
  labs(x="Risk Treatment",y=(expression(atop("Interactions with Conspecifics", 
                                             paste((displays~min^-1))))))+
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
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0,0),limits = c(0,0.151))
reco.plot + geom_linerange(aes(ymin=x-se, ymax=x+se), size=0.5,   
                           position=position_dodge(.85)) + theme(text = element_text(family="Arial"))

dev.off()

#movement rate (movements per minute), includes code for superscript####
mm<-with(behave, aggregate((movements.min), list(Treatment=Treatment), mean))
mm$se<-with(behave, aggregate((movements.min), list(Treatment=Treatment), 
                              function(x) sd(x)/sqrt(length(x))))[,2]

png("Output/2019.11.13.movement.swims.no.parentheses.9.5x5.5.300dpi.png", width = 6.5, height = 5.5, units = 'in', res = 300)

reco.plot<- ggplot(mm, aes(x=Treatment, y=x, fill=Treatment)) +
  geom_bar(stat="identity", colour= "black", width = 0.5, position="dodge")+ 
  scale_x_discrete(limits=c("Low","Medium","High"))+
  theme_classic() + 
  labs(x="Risk Treatment", y ="Movements per Minute")+
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
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0,0),limits = c(0,1.51),
                          labels = scales::number_format(accuracy = 0.01)) #changed to 2 decimal places for movement rate
reco.plot + geom_linerange(aes(ymin=x-se, ymax=x+se), size=0.5,   
                           position=position_dodge(.85)) + theme(text = element_text(family="Arial"))
dev.off()

#movements per minute with rate in parentheses

png("Output/2019.11.14.movement.swims.with.rate.parentheses.9.5x5.5.300dpi.png", width = 6.5, height = 5.5, units = 'in', res = 300)

reco.plot<- ggplot(mm, aes(x=Treatment, y=x, fill=Treatment)) +
  geom_bar(stat="identity", colour= "black", width = 0.5, position="dodge")+ 
  scale_x_discrete(limits=c("Low","Medium","High"))+
  theme_classic() +
  labs(x="Risk Treatment",y=(expression(atop("Movements min"^-1))))+ 
  #                                         paste((movements~min^-1))))))+
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
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0,0),limits = c(0,1.51),
                          labels = scales::number_format(accuracy = 0.01)) #changed to 2 decimal places for movement rate)
reco.plot + geom_linerange(aes(ymin=x-se, ymax=x+se), size=0.5,   
                           position=position_dodge(.85)) + theme(text = element_text(family="Arial"))

dev.off()


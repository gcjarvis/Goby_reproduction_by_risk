# Description: colored figures for talks
# Author: George C Jarvis
# Date: Sun Mar 31 21:05:06 2019
# --------------  

rm(list=ls())

library(extrafont)
library(ggplot2)
library(sciplot)
library(dplyr)
library(plyr)
library(vegan)
library(psych)
library(lme4)

#egg counts

egg.den.bio<-read.csv("Data/jarvis.egg.count.data.with.den.max.2019.3.6.csv") #uses adjusted counts for density
#updated data
egg.den.bio<-read.csv("Data/new.data.2019.4.16a.t1.5.csv", na.strings = "")

#trial 1-3
egg.2017.t1.2.3<-egg.den.bio[(egg.den.bio$Trial<4),]
egg.2017.t1.2.3$Treatment<-ordered(egg.2017.t1.2.3$Treatment,levels=c("Low","Medium","High"))

#trials 4 and 5
egg.2018.t4.5<-egg.den.bio[(egg.den.bio$Trial>3) & (egg.den.bio$Trial<6), ]
egg.2018.t4.5$Treatment<-ordered(egg.2018.t4.5$Treatment,levels=c("Low","Medium","High"))

#trial 6
egg.2017.t6<-egg.den.bio[(egg.den.bio$Trial>5),]
egg.2017.t6$Treatment<-ordered(egg.2017.t6$Treatment,levels=c("Low","Medium","High","Control"))


#setting up df's for plotting
df<-egg.2017.t1.2.3
df<-egg.2018.t4.5
df<-egg.2017.t6

#generic code using df
egg.means<-with(df, aggregate((Egg.count), list(Treatment=Treatment), mean))
egg.means
#now apply the se function to the 4th column [,3]
egg.means$se<-with(df, aggregate((Egg.count), list(Treatment=Treatment), function(x) sd(x)/sqrt(length(x))))[,2]
egg.means

#to get the figure to export as a png, make sure you change output file name
png(filename = "Output/barplot.t.1.2.3.png", width = 700, height = 800)

egg.wsn<- ggplot(egg.means, aes(x=Treatment, y=x, fill=Treatment)) +
  geom_bar(stat="identity", colour= "black", width = 0.7, position="dodge")+ 
  scale_x_discrete(limits=c("Low","Medium","High"))+
  theme_classic() + theme(legend.position="none")  + #scale_fill_discrete(name="Sex",labels=c(" Male"," Female", " Transitional")) +
  theme(legend.key.size = unit(1.3,'line')) + 
  scale_fill_manual(values=c("#0072B2","#009E73","#D55E00")) + 
  theme(axis.text.x=element_text(size=32, colour="black"),axis.text.y=element_text(size=32, colour="black"), axis.title=element_text(size=37,face="bold")) +
  theme(axis.title.y = element_text(size= 37, margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 25, r = 0, b = 0, l = 0)), axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0, 0),limits = c(0, 4800))
egg.wsn + geom_linerange(aes(ymin=x-se, ymax=x+se), size=0.5,   
                         position=position_dodge(.7)) + theme(text = element_text(family="Arial")) +
  labs(x="Risk Treatment", y="Total Reproduction per Reef")

dev.off()

#adding control treatment
png(filename = "Output/barplot.t.6.png", width = 700, height = 800)

egg.wsn<- ggplot(egg.means, aes(x=Treatment, y=x, fill=Treatment)) +
  geom_bar(stat="identity", colour= "black", width = 0.7, position="dodge")+ 
  scale_x_discrete(limits=c("Low","Medium","High","Control"))+
  theme_classic() + theme(legend.position="none")  + #scale_fill_discrete(name="Sex",labels=c(" Male"," Female", " Transitional")) +
  theme(legend.key.size = unit(1.3,'line')) + 
  scale_fill_manual(values=c("#0072B2","#009E73","#D55E00","red")) + 
  theme(axis.text.x=element_text(size=32, colour="black"),axis.text.y=element_text(size=32, colour="black"), axis.title=element_text(size=37,face="bold")) +
  theme(axis.title.y = element_text(size= 37, margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 25, r = 0, b = 0, l = 0)), axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0, 0),limits = c(0, 4800))
egg.wsn + geom_linerange(aes(ymin=x-se, ymax=x+se), size=0.5,   
                         position=position_dodge(.7)) + theme(text = element_text(family="Arial")) +
  labs(x="Risk Treatment", y="Total Reproduction per Reef")

dev.off()


egg.means<-with(egg.2018.t4.5, aggregate((Egg.count), list(Treatment=Treatment), mean))
egg.means
#now apply the se function to the 4th column [,3]
egg.means$se<-with(egg.2018.t4.5, aggregate((Egg.count), list(Treatment=Treatment), function(x) sd(x)/sqrt(length(x))))[,2]
egg.means

egg.wsn<- ggplot(egg.means, aes(x=Treatment, y=x, fill=Treatment)) +
  geom_bar(stat="identity", colour= "black", width = 0.7, position="dodge")+ 
  scale_x_discrete(limits=c("Low","Medium","High"))+
  theme_classic() + theme(legend.position="none")  + #scale_fill_discrete(name="Sex",labels=c(" Male"," Female", " Transitional")) +
  theme(legend.key.size = unit(1.3,'line')) + 
   scale_fill_manual(values=c("#0072B2","#009E73","#D55E00")) + 
  theme(axis.text.x=element_text(size=20, colour="black"),axis.text.y=element_text(size=20, colour="black"), axis.title=element_text(size=25,face="bold")) +
  theme(axis.title.y = element_text(size= 25, margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 25, r = 0, b = 0, l = 0)), axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0, 0))
egg.wsn + geom_linerange(aes(ymin=x-se, ymax=x+se), size=0.5,   
                         position=position_dodge(.7)) + theme(text = element_text(family="Arial")) +
  labs(x="Risk Treatment", y="test")#+
  #geom_point(data=egg.2018.t4.5,aes(Treatment,Egg.count,color=Treatment),position="dodge")+
  #scale_color_manual(values=c("black","black","black"))

########boxplot for BEM for egg counts (t4 and 5) with Base R############

#have to set up means df in order to overlay diamonds on plot
means <- tapply(egg.2018.t4.5$Egg.count,egg.2018.t4.5$Treatment,mean)
points(means,col="red",pch=18)

png(filename = "Output/baxplot.defense.t.6.png", width = 800, height = 700)
boxplot(Egg.count~Treatment,data=egg.2018.t4.5,
        col=c("#0072B2","#009E73","#D55E00"),
        ylab="Total Reproduction per Reef",xlab="Risk Treatment",
        cex.lab=2.0, cex.axis=1.65,
        par(family="sans"), par(font.lab=2), par(mar=c(5,6,4,1)+.1))
points(means,col="black",pch=23,bg="white",cex=2,lwd=1.5) 

dev.off()

#consider with notch?

boxplot(Egg.count~Treatment,data=egg.2018.t4.5,notch=TRUE,
        col=c("#0072B2","#009E73","#D55E00"),
        ylab="Total Reproduction per Reef",xlab="Risk Treatment",
        cex.lab=2.0, cex.axis=1.65,
        par(family="sans"), par(font.lab=2), par(mar=c(5,6,4,1)+.1))
        points(means,col="black",pch=23,bg="white",cex=2,lwd=1.5)
        

        
#lots of variation in output even within treatment

#thinking it might be easiest to try to recreate the plot from base R with my colors, axes labels, and Arial font

boxplot(len~supp*dose, data=ToothGrowth, notch=TRUE, 
        col=(c("gold","darkgreen")),
        main="Tooth Growth", xlab="Suppliment and Dose")


########trying boxplot for egg counts as well, tons of variation in the data
#not going to be able to do it for this talk, but want to eventually
# found good stuff here: http://www.sthda.com/english/wiki/ggplot2-box-plot-quick-start-guide-r-software-and-data-visualization

egg.wsn<- ggplot(egg.means, aes(x=Treatment, y=x, fill=Treatment)) +
  geom_boxplot(stat="identity", colour= "black", width = 0.7, position="dodge")+ 
  scale_x_discrete(limits=c("Low","Medium","High"))+
  theme_classic() + theme(legend.position="none")  + #scale_fill_discrete(name="Sex",labels=c(" Male"," Female", " Transitional")) +
  #theme(legend.key.size = unit(1.3,'line')) + 
  scale_fill_manual(values=c("#0072B2","#009E73","#D55E00")) + 
  theme(axis.text.x=element_text(size=20, colour="black"),axis.text.y=element_text(size=20, colour="black"), axis.title=element_text(size=25,face="bold")) +
  theme(axis.title.y = element_text(size= 25, margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 25, r = 0, b = 0, l = 0)), axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0, 0))


egg.wsn + geom_linerange(aes(ymin=x-se, ymax=x+se), size=0.5,   
                         position=position_dodge(.7)) + theme(text = element_text(family="Arial")) +
  labs(x="Risk Treatment", y="Number of Eggs per Reef")+
  geom_point(data=egg.2018.t4.5,aes(Treatment,Egg.count,color=Treatment),position="dodge")+
  scale_color_manual(values=c("black","black","black"))




egg.count<- ggplot(egg.means, aes(x=Treatment, y=x, fill=Treatment)) +
  geom_bar(stat="identity", colour= "black", width = 0.7, position="dodge")+ 
  scale_x_discrete(limits=c("Low","Medium","High"))+
  theme_classic() + theme(legend.position="none") + scale_fill_manual(values=c("#D55E00","#0072B2","#009E73")) + 
  theme(axis.text.x=element_text(size=25, colour="black"),axis.text.y=element_text(size=18, colour="black"), axis.title=element_text(size=25,face="bold")) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 25, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 22, r = 0, b = 0, l = 0)), axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0, 0))
egg.count + geom_linerange(aes(ymin=x-se, ymax=x+se), size=0.5,   
                           position=position_dodge(.9)) + theme(text = element_text(family="Arial")) +
  labs(x="Risk Treatment", y="Average # Eggs Produced per Reef")

########visual surveys#########

viz.sur<-with(egg.2018.t4.5, aggregate((Density), list(Week=Week,Treatment=Treatment), mean))
#now apply the se function to the 4th column [,3]
viz.sur$se<-with(egg.2018.t4.5, aggregate((Density), list(Week=Week,Treatment=Treatment), 
                                          function(x) sd(x)/sqrt(length(x))))[,3]
viz.sur

p <- ggplot(viz.sur, aes(x=Week, y=x, group=Treatment, color=Treatment))+ 
  geom_linerange(aes(ymin=x-se, ymax=x+se), 
                 position=position_dodge(0)) +
  geom_line(size=1.5)+
  geom_line(aes(linetype=Treatment)) +
  geom_point(aes(shape=Treatment))+
  labs(x="Week", y = "Number of Fish Per Reef")+
  theme_classic() + 
  theme(axis.text.x=element_text(size=28, colour="black"),axis.text.y=element_text(size=28, colour="black"), axis.title=element_text(size=35,face="bold")) +
  theme(axis.title.y = element_text(size= 34, margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)), axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0, 0)) +
  theme(text = element_text(family="Arial"))
p + scale_color_manual(values=c("#0072B2","#009E73","#D55E00")) + theme(legend.text=element_text(size=20)) + theme(legend.title =element_text(size=28, face="bold"))


#number of gobies recollected######
reco.means<-with(df, aggregate((Count), list(Treatment=Treatment), mean))
reco.means$se<-with(df, aggregate((Count), list(Treatment=Treatment), function(x) sd(x)/sqrt(length(x))))[,2]
reco.means

reco.means$Treatment<-ordered(reco.means$Treatment,levels=c("Low","Medium","High"))

png(filename = "Output/reco.bem.png", width = 900, height = 800)

reco.plot<- ggplot(reco.means, aes(x=Treatment, y=x, fill=Treatment)) +
  geom_bar(stat="identity", colour= "black", width = 0.85, position="dodge")+ 
  scale_x_discrete(limits=c("Low","Medium","High"))+
  theme_classic() + theme(legend.position="none") + #scale_fill_discrete(limits=c("Low","Medium", "High")) +
  #theme(legend.key.size = unit(1.3,'line')) + 
  #theme(legend.title=element_text(size=34) , legend.text=element_text(size=20)) 
  scale_fill_manual(values=c("#0072B2","#009E73","#D55E00")) + 
  theme(axis.text.x=element_text(size=30, colour="black"),axis.text.y=element_text(size=30, colour="black"), axis.title=element_text(size=35,face="bold")) +
  theme(axis.title.y = element_text(size= 35, margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 25, r = 0, b = 0, l = 0)), axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0, 0))
reco.plot + geom_linerange(aes(ymin=x-se, ymax=x+se), size=0.5,   
                           position=position_dodge(.85)) + theme(text = element_text(family="Arial")) +
  labs(x="Risk Treatment", y="Gobies Recollected per Reef")

dev.off()

#looking at proportion of zeros by treatment
plot(factor(Egg.count==0)~Treatment, data=egg.2018.t4.5)
plot(factor(Egg.count >1)~Treatment, data=egg.2018.t4.5,
     main="Proportion of reefs without eggs",
     xlab="Risk Treatment",
     cex.lab=2.0, cex.axis=1.65,
     par(family="sans"), par(font.lab=2))
points(means,col="black",pch=23,bg="white",cex=2,lwd=1.5))

###### behavior figures

exp.means<-with(behave, aggregate((proportion.exposed), list(Treatment=Treatment), mean))
exp.means
#now apply the se function to the 4th column [,3]
exp.means$se<-with(behave, aggregate((proportion.exposed), list(Treatment=Treatment), function(x) sd(x)/sqrt(length(x))))[,2]
exp.means

#to get the figure to export as a png
png(filename = "Output/exposure.wsn.png", width = 900, height = 800)

exp.wsn<- ggplot(exp.means, aes(x=Treatment, y=x, fill=Treatment)) +
  geom_bar(stat="identity", colour= "black", width = 0.7, position="dodge")+ 
  scale_x_discrete(limits=c("Low","Medium","High"))+
  theme_classic() + theme(legend.position="none")  + #scale_fill_discrete(name="Sex",labels=c(" Male"," Female", " Transitional")) +
  theme(legend.key.size = unit(1.3,'line')) + 
  scale_fill_manual(values=c("#0072B2","#009E73","#D55E00")) + 
  theme(axis.text.x=element_text(size=30, colour="black"),axis.text.y=element_text(size=30, colour="black"), axis.title=element_text(size=35,face="bold")) +
  theme(axis.title.y = element_text(size= 35, margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 25, r = 0, b = 0, l = 0)), axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0, 0))
exp.wsn + geom_linerange(aes(ymin=x-se, ymax=x+se), size=0.5,   
                         position=position_dodge(.7)) + theme(text = element_text(family="Arial")) +
  labs(x="Risk Treatment", y="Proportion of Time Exposed")
  
dev.off()

#courtship
plot(factor(courtship==0)~Treatment, data=behave)#60% of the data are 0, but equally across treatments
#quick model, not sure if I can really say much about this one statistically
mod1<-lm(courtship~Treatment+Density,family=binomial,data=behave.2018.t4.5)
qqnorm(resid(mod1))
hist(resid(mod1))

#subsetting for times when courtship was not eaul to zero
# based on variable values
behave.court <- behave[ which(behave$courtship>0), ]
View(behave.court)

behave$court.units<-(behave$courtship/5)

court.means<-with(behave, aggregate((court.units), list(Treatment=Treatment), mean))
#now apply the se function to the 4th column [,3]
court.means$se<-with(behave, aggregate((court.units), list(Treatment=Treatment), function(x) sd(x)/sqrt(length(x))))[,2]
court.means

#to get the figure to export as a png
png(filename = "Output/courtship.wsn.png", width = 900, height = 800)

court.wsn<- ggplot(court.means, aes(x=Treatment, y=x, fill=Treatment)) +
  geom_bar(stat="identity", colour= "black", width = 0.7, position="dodge")+ 
  scale_x_discrete(limits=c("Low","Medium","High"))+
  theme_classic() + theme(legend.position="none")  + #scale_fill_discrete(name="Sex",labels=c(" Male"," Female", " Transitional")) +
  theme(legend.key.size = unit(1.3,'line')) + 
  scale_fill_manual(values=c("#0072B2","#009E73","#D55E00")) + 
  theme(axis.text.x=element_text(size=30, colour="black"),axis.text.y=element_text(size=30, colour="black"), axis.title=element_text(size=35,face="bold")) +
  theme(axis.title.y = element_text(size= 35, margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 25, r = 0, b = 0, l = 0)), axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0, 0))
court.wsn + geom_linerange(aes(ymin=x-se, ymax=x+se), size=0.5,   
                         position=position_dodge(.7)) + theme(text = element_text(family="Arial")) +
  labs(x="Risk Treatment", y="Courtship Displays per Minute")

dev.off()

#foraging


behave$bites.min<-(behave$num.bites/5)
  
for.means<-with(behave, aggregate((bites.min), list(Treatment=Treatment), mean))
#now apply the se function to the 4th column [,3]
for.means$se<-with(behave, aggregate((bites.min), list(Treatment=Treatment), function(x) sd(x)/sqrt(length(x))))[,2]
for.means

#to get the figure to export as a png
png(filename = "Output/foraging.wsn.png", width = 900, height = 800)

forage.wsn<- ggplot(for.means, aes(x=Treatment, y=x, fill=Treatment)) +
  geom_bar(stat="identity", colour= "black", width = 0.7, position="dodge")+ 
  scale_x_discrete(limits=c("Low","Medium","High"))+
  theme_classic() + theme(legend.position="none")  + #scale_fill_discrete(name="Sex",labels=c(" Male"," Female", " Transitional")) +
  theme(legend.key.size = unit(1.3,'line')) + 
  scale_fill_manual(values=c("#0072B2","#009E73","#D55E00")) + 
  theme(axis.text.x=element_text(size=30, colour="black"),axis.text.y=element_text(size=30, colour="black"), axis.title=element_text(size=35,face="bold")) +
  theme(axis.title.y = element_text(size= 35, margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 25, r = 0, b = 0, l = 0)), axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0, 0))
forage.wsn + geom_linerange(aes(ymin=x-se, ymax=x+se), size=0.5,   
                           position=position_dodge(.7)) + theme(text = element_text(family="Arial")) +
  labs(x="Risk Treatment", y="Bites per Minute")

dev.off()


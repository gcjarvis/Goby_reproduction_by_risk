##############################################################################
# recollections data to get at growth rates                                  #
#   Hunter and George                                                        #
#    7/27/2018                                                               #
##############################################################################

rm(list=ls())

#load packages
library(sciplot)
library(lme4)
library(lmerTest)
library(car)
library(dplyr)
library(ggplot2)
library(extrafont)
#library(ggplot)

getwd()
sc<-read.csv("Data/Recollections_Trial_4.5.hunter.analyses.csv")
sc$Trial<-as.factor(sc$Trial)#made trial a factor
#only including adults that we tagged initially (immatures were listed as NA in datasheet)
sc<- sc[complete.cases(sc), ]#122 rows of data
sc$Growth<-as.numeric(sc$Growth)

sc.mod.g<-lm(Growth~Treatment, data=sc)
hist(resid(sc.mod.g))
qqnorm(resid(sc.mod.g))
boxplot(resid(sc.mod.g))
boxplot(Growth~Treatment, data=sc)

anova(sc.mod.g)
summary(sc.mod.g)

#plotting growth by risk
bargraph.CI(x.factor = Treatment, response = Growth, main="Growth by risk", data = sc)
#seems like no sig difference in growth among treatments

#adding in initial sex as a fixed factor
sc.mod.initial<-lm(Growth~Treatment+Initial.Sex, data=sc)
hist(resid(sc.mod.initial))
qqnorm(resid(sc.mod.initial))
boxplot(resid(sc.mod.initial))
boxplot(Growth~Treatment, data=sc)

anova(sc.mod.initial)
summary(sc.mod.initial)

#two-way ANOVA with initial sex
sc.mod.initial.cross<-lm(Growth~Treatment*Initial.Sex, data=sc)
hist(resid(sc.mod.initial.cross))
qqnorm(resid(sc.mod.initial.cross))
boxplot(resid(sc.mod.initial.cross))
#boxplot(Growth~Treatment, data=sc)
anova(sc.mod.initial.cross)
summary(sc.mod.initial.cross)

#two-way ANOVA with final sex
sc.mod.initial.cross.final.sex<-lm(Growth~Treatment*Sex, data=sc)
hist(resid(sc.mod.initial.cross.final.sex))
qqnorm(resid(sc.mod.initial.cross.final.sex))
boxplot(resid(sc.mod.initial.cross.final.sex))
#boxplot(Growth~Treatment, data=sc)
anova(sc.mod.initial.cross.final.sex)
summary(sc.mod.initial.cross.final.sex)

#model selection
anova(sc.mod.initial, sc.mod.g, sc.mod.initial.cross)

#factoring in initial sex
bargraph.CI(x.factor = Treatment, response = Growth, group=Initial.Sex, legend = TRUE, main="Growth by risk", data = sc)

#just running initial sex as the factor
bargraph.CI(x.factor = Initial.Sex, response = Growth, main="Growth by sex", data = sc)

#looking at the initial sizes of the fish among treatments to make sure one treatment didn't have a ton of huge fish (would have low growth rates if so)
bargraph.CI(x.factor = Treatment, response = Tag, main="Initial tag by treatment", data = sc)

bargraph.CI(x.factor = Initial.Sex, response = Treatment, main="Initial sex by treatment", data = sc)

#final size by risk
bargraph.CI(x.factor = Treatment, response = Final.Size, main="Final size by treatment", data = sc)

#growth by final sex
bargraph.CI(x.factor = Treatment, response = Growth, group=Sex, legend = TRUE, main="Growth by risk", data = sc)


######## ggplot figures ########
#now make ggplot tables
growth.fig<-(with(sc, aggregate((Growth), list(Treatment=Treatment,Initial.Sex), mean)))
#now apply the se function to the 3rd column [,3]

growth.fig$se<-with(sc, aggregate((Growth), list(Treatment=Treatment,Initial.Sex), function(x) sd(x)/sqrt(length(x))))[,3]
growth.fig

#egg.fig.1<-ddply(gobydata,c("Treatment"),summarise,mean=mean(Total_Eggs, na.rm=TRUE), SE=std.error(Total_Eggs, na.rm=TRUE))
#egg.fig.1

#no grid, no caps on error bars, NARROWER BARS, NO LEGEND
#egg.count<- ggplot(egg.fig, aes(x=Treatment, y=x, fill=Treatment)) +
  #geom_bar(stat="identity", colour="black",position="dodge")+ 
  #scale_x_discrete(limits=c("Low","Medium","High"))+
  #theme_classic()  
#egg.count + geom_linerange(aes(ymin=x-se, ymax=x+se), size=0.5,   
                           #position=position_dodge(.9)) + theme(text = element_text(size = 14)) +
  #labs(x="Risk Treatment", y="Proportion of fish observed Day 6")


#order of egg.count table is high, low, medium, so that's how you have to code the color order

##MASTER GGPLOT CODE USED FOR FIGURES --- this has all of the axes labels and colors as they should appear in the PPT
growth.plot<- ggplot(growth.fig, aes(x=Treatment, y=x, fill=Group.2)) +
  geom_bar(stat="identity", colour= "black", width = 0.7, position="dodge")+ 
  scale_x_discrete(limits=c("Low","Medium","High"))+
  theme_classic() + theme(legend.position="none") + scale_fill_manual(values=c("#D55E00","#009E73")) + 
  theme(axis.text.x=element_text(size=25, colour="black"),axis.text.y=element_text(size=18, colour="black"), axis.title=element_text(size=25,face="bold")) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 25, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 22, r = 0, b = 0, l = 0)), axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0, 0))
growth.plot + geom_linerange(aes(ymin=x-se, ymax=x+se)), size=0.5,   
                           position=position_dodge(.9)) + theme(text = element_text(family="Arial")) +
  labs(x="Risk Treatment", y="Average growth (mm)")

egg.count<- ggplot(egg.fig, aes(x=Treatment, y=x, fill=Treatment)) +
  geom_bar(stat="identity", colour= "black", width = 0.7, position="dodge")+ 
  scale_x_discrete(limits=c("Low","Medium","High"))+
  theme_classic() + theme(legend.position="none") + scale_fill_manual(values=c("#D55E00","#0072B2","#009E73")) + 
  theme(axis.text.x=element_text(size=25, colour="black"),axis.text.y=element_text(size=18, colour="black"), axis.title=element_text(size=25,face="bold")) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 25, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 22, r = 0, b = 0, l = 0)), axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0, 0))
egg.count + geom_linerange(aes(ymin=x-se, ymax=x+se), size=0.5,   
                           position=position_dodge(.9)) + theme(text = element_text(family="Arial")) +
  labs(x="Risk Treatment", y="Average # Eggs Produced per Reef")

head(sc)

sample(1:20, 16, replace=FALSE)

##data for hunter growth figure

#now make ggplot tables
growth.fig1<-(with(sc, aggregate((Growth), list(Treatment=Treatment), mean)))
#now apply the se function to the 3rd column [,3]

growth.fig1$se<-with(sc, aggregate((Growth), list(Treatment=Treatment), function(x) sd(x)/sqrt(length(x))))[,2]
growth.fig1


#means for initial sex
growth.fig.sex<-(with(sc, aggregate((Growth), list(Initial.Sex), mean)))
#now apply the se function to the 3rd column [,3]

growth.fig.sex$se<-with(sc, aggregate((Growth), list(Initial.Sex), function(x) sd(x)/sqrt(length(x))))[,2]
growth.fig.sex

##
sc.mod.g<-lm(Growth~Initial.Sex, data=sc)
hist(resid(sc.mod.g))
qqnorm(resid(sc.mod.g))
boxplot(resid(sc.mod.g))
boxplot(Growth~Treatment, data=sc)
anova(sc.mod.g)

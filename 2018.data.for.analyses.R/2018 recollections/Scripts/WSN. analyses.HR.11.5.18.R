#analyses for Hunter's poster
#11.5.18
#went back and only put down proportions of each sex present on each reef for each trial

rm(list=ls())

#load packages
library(sciplot)
library(lme4)
library(lmerTest)
library(car)
library(dplyr)
library(ggplot2)
library(extrafont)
library(ggplot2)
library(plyr)


getwd()
reco<-read.csv("Data/recollections.proportions.only.11.5.18.csv")
#only had one reef where I didn't recollect any fish
reco<- reco[complete.cases(reco), ]
reco$Sex<-ordered(reco$Sex,levels=c("Male","Female","Transitional"))
View(reco)

#subsetting the data for just proportion of females recollected
#reco.fem<-gob.sub[gob.sub$week==1, c("week", "deployment.day", "reef", "treatment", "egg.count", "initial.biomass","final.biomass","change.in.biomass","average.biomass","eggs.initial.biomass","eggs.avg.biomass","nest.count","eggs.per.nest")]
reco.fem<-reco[reco$sex=="female", c("trial", "sex", "reef", "treatment", "final.proportion")]
View(reco.fem)


bargraph.CI(x.factor = treatment, response = total.recollection, group=Sex, legend=TRUE, main="recollections by sex, trial 4",x.leg = 10, data = reco)

#model selection
reco1<-lm(final.proportion~treatment, data=reco)
hist(resid(reco1))
anova(reco1)
reco2<-lm(final.proportion~treatment*sex, data=reco)#looks like the best model
hist(resid(reco2))
qqnorm(resid(reco2))
qqline(resid(reco2))

anova(reco1,reco2)
AIC(reco1)
AIC(reco2)

#now including trial as a factor
reco3<-lmer(final.proportion~treatment*sex+(1|trial), data=reco)
hist(resid(reco3))
qqnorm(resid(reco3))
AIC(reco3)

#gamma distribution
#worried it might kickback because I have some zeros...yep, non-positive values...

#mod5.gamma<-glmer(final.proportion~treatment*sex+(1|trial), family=Gamma(link="inverse"), data=reco)
#anova(mod5.gamma)
#summary(mod5.gamma)

#######model that I'm going to use for Hunter's analyses#####
reco2<-lm(total.recollection~treatment*sex, data=reco)#looks like the best model
hist(resid(reco2))
qqnorm(resid(reco2))
qqline(resid(reco2))
boxplot(final.proportion~treatment, data=reco)

anova(reco2)

reco3<-lmer(final.proportion~treatment*sex+(1|trial), data=reco)#let's just go with this one
hist(resid(reco3))
qqnorm(resid(reco3))

anova(reco3)
summary(reco3)

#########now trying to see how the data look just with females
reco.fem1<-lm(final.proportion~treatment, data=reco.fem)#looks like the best model
hist(resid(reco.fem1))
qqnorm(resid(reco.fem1))
qqline(resid(reco.fem1))
#boxplot(final.proportion~treatment, data=reco.fem)

anova(reco.fem1)
bargraph.CI(x.factor = treatment, response = final.proportion, main="proportion female recollected", data = reco.fem)

reco.fem2<-lm(final.proportion~treatment*trial,data=reco.fem) #looks like the best model
hist(resid(reco.fem2))
qqnorm(resid(reco.fem2))
qqline(resid(reco.fem2))
anova(reco.fem2)
reco.fem3<-lmer(final.proportion~treatment+(1|trial),data=reco.fem)
hist(resid(reco.fem3))

##############Looking at number of fish recolleccted by treatment and sex33####

bargraph.CI(x.factor = treatment, response = total.recollection, group=sex, legend=TRUE, main="recollections by sex", data = reco)

#model selection
reco1.tot<-lm(total.recollection~treatment, data=reco)
hist(resid(reco1))
anova(reco1)

#model that I went with for my WSN talk:
reco2.tot<-lmer(total.recollection~treatment*Sex+(1|deployment.day), data=reco)#looks like the best model
hist(resid(reco2.tot))
qqnorm(resid(reco2.tot))
qqline(resid(reco2.tot))
anova(reco2.tot)
summary(reco2.tot)

reco.2.poiss<-glmer(total.recollection~treatment*Sex+(1|deployment.day),family=poisson, data=reco)
hist(resid(reco2.tot))
qqnorm(resid(reco2.tot))
qqline(resid(reco2.tot))
anova(reco2.tot)
summary(reco2.tot)



anova(reco1.tot,reco2.tot)
AIC(reco1.tot)
AIC(reco2.tot)

reco3.tot<-lmer(total.recollection~treatment*sex+(1|trial),data=reco)
hist(resid(reco3.tot))
AIC(reco3.tot)

reco3.poiss<-glmer(total.recollection~treatment*sex+(1|trial),family=poisson,data=reco)
hist(resid(reco3.poiss))
AIC(reco3.poiss)
qqnorm(resid(reco3.poiss))
qqline(resid(reco3.poiss))
boxplot(total.recollection~treatment*sex,data=reco)
boxplot(total.recollection~treatment,data=reco)
summary(reco3.poiss)

Anova(reco3.poiss, type="II")
Anova(reco3.poiss, type="III")

#getting means and standard error
reco.means<-with(reco, aggregate((total.recollection), list(Sex=Sex,treatment=treatment), mean))
reco.means

reco.means<-with(reco, aggregate((total.recollection), list(Sex=Sex,treatment=treatment), mean))
reco.means
#now apply the se function to the 4th column [,3]
reco.means$se<-with(reco, aggregate((total.recollection), list(Sex=Sex,treatment=treatment), function(x) sd(x)/sqrt(length(x))))[,3]
reco.means
#gob.sub$treatment<-ordered(gob.sub$treatment,levels=c("Low","Medium","High"))
#reco.means$sex<-ordered(reco.means$sex,levels=c("Male","Female","Transitional"))
reco.means

#viz.sur$se<-with(goby.den.survey, aggregate((Proportion_pop), list(Days_after_deploy=Days_after_deploy,Treatment=Treatment), function(x) sd(x)/sqrt(length(x))))[,3]
#viz.sur

##########plotting for WSN poster##########
##ggplot

#going to see if I can reorder the fish labels to get the two HASE together

#levels(spp.abund$Species)
#spp.abund$Species <- factor(spp.abund$Species, 
#levels = c("FHASE","MHASE", "FSEPU","PACL"))


spp.plot<- ggplot(spp.abund, aes(x=Risk, y=x, fill=Species)) +
  geom_bar(stat="identity", colour= "black", width = 0.7, position="dodge")+ 
  scale_x_discrete(limits=c("Low","Medium","High"))+
  theme_classic() + theme(legend.position="top") + #scale_fill_discrete(name="Species",labels=c("FHASE          ","MHASE          ", "SEPU          ","PACL          ")) +
  theme(legend.key.size = unit(1.5,'line')) + 
  theme(legend.title=element_text(size=15) , legend.text=element_text(size=10)) + scale_fill_manual(values=c("#D55E00","#0072B2","#CC79A7","#009E73")) + 
  theme(axis.text.x=element_text(size=18, colour="black"),axis.text.y=element_text(size=16, colour="black"), axis.title=element_text(size=20,face="bold")) +
  theme(axis.title.y = element_text(size= 20, margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)), axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0, 0))
spp.plot + geom_linerange(aes(ymin=x-se, ymax=x+se), size=0.5,   
                          position=position_dodge(.7)) + theme(text = element_text(family="Arial")) +
  labs(x="Access Level", y="Predator Abundance by Species")

#don't need to specify FSEPU, becuase we didn't see many MSEPU

reco.plot<- ggplot(reco.means, aes(x=treatment, y=x, fill=Sex)) +
  geom_bar(stat="identity", colour= "black", width = 0.85, position="dodge")+ 
  scale_x_discrete(limits=c("Low","Medium","High"))+
  theme_classic() + theme(legend.position="right") + #scale_fill_discrete(name="Sex",labels=c(" Male"," Female")) +
  theme(legend.key.size = unit(1.3,'line')) + 
  theme(legend.title=element_text(size=20) , legend.text=element_text(size=18)) + scale_fill_manual(values=c("#0072B2","#D55E00","#009E73")) + 
  theme(axis.text.x=element_text(size=20, colour="black"),axis.text.y=element_text(size=20, colour="black"), axis.title=element_text(size=25,face="bold")) +
  theme(axis.title.y = element_text(size= 25, margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 25, r = 0, b = 0, l = 0)), axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0, 0))
reco.plot + geom_linerange(aes(ymin=x-se, ymax=x+se), size=0.5,   
                          position=position_dodge(.85)) + theme(text = element_text(family="Arial")) +
  labs(x="Risk Treatment", y="Number of Gobies Recollected")

#making a histogram to show the distribution of males to females based on size

gob.hist<-read.csv("Data/natural.ratio.11.6.18.csv")
gob.hist<-gob.hist[complete.cases(gob.hist),]
#revisit this when I make the plot
#gob.hist$Sex<-ordered(gob.hist$Sex,levels=c("Male","Female"))
head(gob.hist)

#figuring out the group mean with plyr package, so I can plot lines on my hist

#mu <- ddply(df, "sex", summarise, grp.mean=mean(weight))
#head(mu)

gob.hist.sum <- ddply(gob.hist, "Sex", summarise, grp.mean=mean(Size))
gob.hist.sum

#plotting
ggplot(df, aes(x=weight, color=sex, fill=sex)) +
  geom_histogram(position="identity", alpha=0.5)+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=sex),
             linetype="dashed")+
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  labs(title="Weight histogram plot",x="Weight(kg)", y = "Count")+
  theme_classic()


#FUCK THIS!!!!!
histogram<-ggplot(gob.hist, aes(x=Size, color=Sex, fill=Sex)) +
  geom_histogram(position="identity", alpha=0.8,binwidth = 1)+
  scale_color_manual(values=c("#D55E00","#0072B2"))+
  scale_fill_manual(values=c("#D55E00","#0072B2"))+
  theme(axis.text.x=element_text(size=16, colour="black"),axis.text.y=element_text(size=16, colour="black"), axis.title=element_text(size=20,face="bold")) +
  theme(axis.title.y = element_text(size= 20, margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)), axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(legend.position=c(0.80,0.80)) +
  theme(legend.key.size = unit(1.5,'line')) + 
  theme(legend.title=element_text(size=20) , legend.text=element_text(size=15)) + 
  theme(scale_y_continuous(expand = c(0, 0)))
histogram+theme_classic()+theme(text = element_text(family="Arial")) +
  labs(x="Size (mm)", y="Frequency") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  theme(legend.position=c(0.80,0.80))+
  theme(legend.title=element_text(size=25) , legend.text=element_text(size=22))+
  theme(axis.text.x=element_text(size=25, colour="black"),axis.text.y=element_text(size=25, colour="black"), axis.title=element_text(size=30,face="bold")) +
  theme(axis.title.y = element_text(size= 30, margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 25, r = 0, b = 0, l = 0)), axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))+
  theme(legend.key.size = unit(1.5,'line'))
  
  


histogram+scale_x_continuous(name="Size",limits=c(0,44))

##recollection by treatment for WSN
#I went back and changed all of the transitional fish to female, b/c it seems like male-female sex change was rare
#didn't want to remove them, b/c they were still recollected

getwd()
reco<-read.csv("Data/recollections.proportions.only.11.5.18.george.transitional.changed.to.female.csv")
#only had one reef where I didn't recollect any fish
reco<- reco[complete.cases(reco), ]
reco$Sex<-ordered(reco$Sex,levels=c("Male","Female","Transitional"))
View(reco)

#model that I went with for my WSN talk:
reco2.tot<-lmer(total.recollection~treatment*Sex+(1|deployment.day), data=reco)#looks like the best model
hist(resid(reco2.tot))
qqnorm(resid(reco2.tot))
qqline(resid(reco2.tot))
anova(reco2.tot)
summary(reco2.tot)
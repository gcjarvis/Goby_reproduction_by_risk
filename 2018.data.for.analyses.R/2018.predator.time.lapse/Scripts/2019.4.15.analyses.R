# Description: have final dataset for ptl, condensed 
#   down by average counts and scores per reef
# Author: George C Jarvis
# Date: Mon Apr 15 09:34:47 2019
# --------------

#Notes: for these data, I got rid of rare spp.

rm(list=ls())

library(sciplot)
library(lme4)
library(car)
library(lmerTest)
library(dplyr)
library(ggplot2)
library(tidyr)
library(plyr)
library(vegan)
library(psych)
library(MASS)
library(boot)
library(effects)
library(simpleboot)

#these data do not contain any info about predator species
ptl<-read.csv("Data/2019.4.15.ptl.no.spp.csv") #short format (column for count and score)
#this is the formaat I have to use for stats
#ptl<-read.csv("Data/2019.4.15.ptl.wrangling.csv") #long-format, used for plotting
ptl<-read.csv("Data/2019.4.16.ptl.wrangling.csv") #long-format, used for plotting
ptl$Treatment<-ordered(ptl$Treatment, c("Low", "Medium","High","Control"))
ptl$measure<-ordered(ptl$measure, c("contain.pred", "sublethal","lethal"))

#subsetting data by trials, might have to analyze separately
ptl.1.3<- ptl[(ptl$Trial <4), ]
ptl.4.5<-ptl[(ptl$Trial>3) & (ptl$Trial<6), ]
ptl.6<-ptl[(ptl$Trial >5), ]

bargraph.CI(x.factor = Treatment, response = count, 
            legend=TRUE, main="counts per treatment, all trials pooled", 
            xlab="treatment", ylab="proportion of photos with predators", 
            x.leg=13, yleg=4500, data = ptl)

bargraph.CI(x.factor = Treatment, response = value, beside=TRUE,
            legend=levels(unique(ptl$measure)), main="scores, all trials pooled", 
            xlab="treatment", ylab="proportion of photos with threatening
            predators", 
            x.leg=13, yleg=4500, data = ptl)

#use this code if you want to plot count and score on the same plot
#this uses the long-format data with a combined factor of measure (count or score)

bargraph.CI(response=value, x.factor=Treatment, group = measure, data = ptl,
            xlab = "Risk Treatment", ylab = "Prop of photos with predators", cex.lab = 1.5, x.leg = 1,
            col = "black", angle = 45, cex.names = 1.25,
            density = c(0,20), legend = FALSE)

#not working the way I want it to in base R or sciplot

#trying in ggplot
ptl.means<-with(ptl, aggregate((value), list(Treatment=Treatment,measure=measure), mean))
ptl.means
#now apply the se function to the 4th column [,3]
ptl.means$se<-with(ptl, aggregate((value), list(Treatment=Treatment,measure=measure), function(x) sd(x)/sqrt(length(x))))[,3]
ptl.means

# Grouped
g<-ggplot(ptl.means, aes(fill=measure, y=x, x=Treatment)) + 
  geom_bar(position="dodge", stat="identity")
g + geom_linerange(aes(ymin=x-se, ymax=x+se), size=0.5,   
                         position=position_dodge(0.90)) + theme(text = element_text(family="Arial")) +
  labs(x="Risk Treatment", y="Proportion of photos")
  

#testing out some models to see if data are normal

#prop of photos with predators seen (focal species)

mod1<-lm(count~Treatment*Trial, data=ptl)
hist(resid(mod1))
qqnorm(resid(mod1))
boxplot(count~Treatment, data=ptl) #control has much less variance (low ss)
anova(mod1) #no diff in the prop of photos with predators seen per treatment
#also, no difference in the prop of photos with preds by trial (can pool)

#prop of photos with predators perceived as threat

mod2<-aov(score~Treatment, data=ptl)
hist(resid(mod2))
qqnorm(resid(mod2))
boxplot(score~Treatment, data=ptl)
anova(mod2) #difference overall, but likely driven by 0 in low treatment
#TukeyHSD(mod2,ptl$score~ptl$Treatment,ordered=TRUE)

tapply(ptl$score,ptl$Treatment,mean)
#    Control       High        Low     Medium 
#   0.05104167 0.06627667 0.00000000 0.03551525

#going to drop low from the model and see if there are still differences
ptl.sub.low<-subset(ptl,Treatment!="Low")
mod3<-lm(score~Treatment,data=ptl.sub.low)
hist(resid(mod3))
qqnorm(resid(mod3))
boxplot(score~Treatment, data=ptl.sub.low)
anova(mod3) #no diff between control, high, and med, as suspected

#now dropping control from the data to see if there's a diff between high and med
ptl.sub.cont<-subset(ptl.sub.low,Treatment!="Control")
mod4<-lm(score~Treatment,data=ptl.sub.cont)
hist(resid(mod4))
qqnorm(resid(mod4))
boxplot(score~Treatment, data=ptl.sub.cont)
anova(mod4) #no diff between high and medium, interesting

#takeaway: it seems that medium, high, and control were equally as risky
# all were more risky than low-risk treatments

#next step: go back and see if you need to make a distinction between lethal and 
# sublethal predation risk perceived, to show that risk was in fact higher
# in the high-risk treatments

#there has to be a way to distinguish between high and medium-risk, to 
# show that preds were in fact able to get closer, and they were

#it's not enough to show average scores by treatment, because then we can't tell
# if preds were hanging out in the same places in all treatments (likely the case)
# or if the averages are a result of some preds on the reef and some far away

#essentially, what are those averages showing? 5 + 1, average is 3, but that means 
# something different biologically that 3 + 3


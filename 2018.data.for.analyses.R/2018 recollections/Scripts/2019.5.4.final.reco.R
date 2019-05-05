# Description: incorporating trial, and using proportion of gobies recollected as the response variable, not number of gobies recollected
# Author: George C Jarvis
# Date: Sat May 04 14:44:30 2019
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

#models including trial as fixed and random effect

#2017
# trial as fixed factor
mod1<-lm(Prop.reco~Treatment*Trial,data=reco.t1.2.3)
hist(resid(mod1))
qqnorm(resid(mod1))
qqline(resid(mod1))
anova(mod1) #higher recollections in trial 3, but equal recollections
# within each trial, regardless of treatment
boxplot(Prop.reco~Treatment*Trial,data=reco.t1.2.3)

#trial as random factor

mod1a<-lmer(Prop.reco~Treatment+(1|Trial),data=reco.t1.2.3)
hist(resid(mod1a))#not normal
qqnorm(resid(mod1a))
qqline(resid(mod1a))
anova(mod1a)
boxplot(Prop.reco~Treatment, data=reco.t1.2.3)

#experiment 1 looks like there were higher recollections in Trial 3
# but those recollections did not vary by treatment (no interactive effects)
# also no effect of treatment on the proportion of gobies recollected
# including trial as a fixed factor seems to have a better fit for normality, 
# but there seems to be equal variance even without trial included

#2018.t.4.5
# trial as fixed factor
mod2<-lm(Prop.reco~Treatment*Trial, data=reco.t4.5)
hist(resid(mod2))
qqnorm(resid(mod2))
qqline(resid(mod2))
plot(mod2)
anova(mod2) 
boxplot(Prop.reco~Treatment*Trial,data=reco.t4.5) #as shown here

#trial as random factor

mod2a<-lmer(Prop.reco~Treatment+(1|Trial),data=reco.t4.5)
hist(resid(mod2a))#not normal
qqnorm(resid(mod2a))
qqline(resid(mod2a))
anova(mod2a)
boxplot(Prop.reco~Treatment, data=reco.t4.5)

#no difference in recollections by trial, but recollections were 
# more variable in trial 4 than in trial 5

#overall, doesn't seem to be too much of an effect of trial on recollections for
#Experiment 2, other than higher variation in proportion of rish recollected in trial 4

#2018.t6, only did one trial, so nothing to compare it to, just rerunning to 
# confirm old stats

reco.t6$log<-log(reco.t6$Prop.reco+1)

mod3<-lm(Prop.reco~Treatment,data=reco.t6)
hist(resid(mod3))#not normal
qqnorm(resid(mod3))
qqline(resid(mod3))
anova(mod3)
Anova(mod3)
boxplot(Prop.reco~Treatment, data=reco.t6)

#same as before, no diff


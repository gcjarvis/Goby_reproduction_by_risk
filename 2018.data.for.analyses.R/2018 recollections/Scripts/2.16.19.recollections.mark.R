######################################
# results for recollections for Mark #                     
#   George Jarvis                    #                                          
#    2.17.19                         #                                                               
######################################

#NOTE: the biomass data for each reef, including the biomass estimates for 2017 recollections, exists in
# the recollection .csv used in these analyses

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
library(ggplot)
library(tidyr)

getwd()
#importing raw data
#tried to make the .csv file not contain any NA's 
reco.raw<-read.csv("Data/Recollections.all.trials.2.23.19.csv")
#View(reco.raw)
#reco.raw<- reco.raw[complete.cases(reco.raw), ]
#reco.raw$growth<-as.integer(reco.raw$growth)
#reco.raw$weight<-as.integer(reco.raw$weight)
#manipulating raw data with dplyr
#trying to figure out why growth values are not the same as they are in the .csv file (10,11,12,etc.)
#growth is going to take some doing..

#let's try out total biomass for now

#total reproduction per week by reef per treatment

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
 
reco.reef<-reco.wrang %>%
  group_by(Year,Trial,Reef,Treatment,Sex.recollected,Deployment.day) %>%
  summarize(Count = sum(Count))
View(reco.reef)

reco.reef$Treatment<-as.character(reco.reef$Treatment)

#subsetting trials####
#2017 trials 1-3
reco.2017<-reco.reef[reco.reef$Trial<4,]
#2018 trials 4-5
reco.2018.t4.5<-reco.reef[c(reco.reef$Trial>3)&(reco.reef$Trial<6),]
#2018 trial 6
reco.2018.t6<-reco.reef[reco.reef$Trial==6,]

#transforming data
#2017 data aren't normally distributed
reco.2017$scount<-sqrt(reco.2017$Count)
reco.2017$logcount<-log(reco.2017$Count + 1)

#2018 t4 and t5 data aren't really either
reco.2018.t4.5$scount<-sqrt(reco.2018.t4.5$Count)
reco.2018.t4.5$logcount<-log(reco.2018.t4.5$Count + 1)

#2018 t6 likely needs transforming as well
reco.2018.t6$scount<-sqrt(reco.2018.t6$Count)
reco.2018.t6$logcount<-log(reco.2018.t6$Count+1)

#df's for analyses
df<-reco.reef
df<-reco.2017
df<-reco.2018.t4.5
df<-reco.2018.t6
View(df)

mod1<-lmer(Count~Treatment+(1|Sex.recollected),data=df)
hist(resid(mod1))#not normal
qqnorm(resid(mod1))
qqline(resid(mod1))
anova(mod1)
boxplot(Count~Treatment,data=df)

mod.log<-lm(logcount~Treatment,data = df)
hist(resid(mod.log))#not normal
qqnorm(resid(mod.log))
qqline(resid(mod.log))
anova(mod.log)

mod.log1<-lm(logcount~Treatment*Sex.recollected*Trial,data = df)
hist(resid(mod.log1))#not normal
qqnorm(resid(mod.log1))
qqline(resid(mod.log1))
anova(mod.log1)

mod.sqrt1<-lm(scount~Treatment*Sex.recollected, data=df)
hist(resid(mod.sqrt1))#not normal
qqnorm(resid(mod.sqrt1))
qqline(resid(mod.sqrt1))
anova(mod.sqrt1)
boxplot(scount~Treatment,data=df)

mod2<-lm(Count~Treatment*Sex.recollected,data=df)
hist(resid(mod2))#more normal than before
qqnorm(resid(mod2)) #still not really that normal
qqline(resid(mod2))
anova(mod2) #difference in the number of fish recollected based on sex

mod2a<-lm(Count~Treatment*Trial,data=df)
hist(resid(mod2a))
qqnorm(resid(mod2a)) 
qqline(resid(mod2a))
anova(mod2a)#no effect of trial on number of fish recollected

mod2ai<-lm(Count~Treatment*Trial*Sex.recollected,data=df)
hist(resid(mod2ai))
qqnorm(resid(mod2ai)) 
qqline(resid(mod2ai))
anova(mod2ai)#interactive effect of trial and sex of fish

mod2aii<-lm(Count~Treatment*Trial*Sex.recollected*Deployment.day,data=df)
hist(resid(mod2aii))
qqnorm(resid(mod2aii)) 
qqline(resid(mod2aii))
anova(mod2aii)#pretty big model,likely info in here that's not necessary

#most realistic model, might try a different dist, becuase current model 
#--assumptions looks a bit wonky, esp, normality, too many zeros
mod2aiii<-lmer(Count~Treatment*Sex.recollected+(1|Trial)+(1|Trial:Deployment.day),data=df)
hist(resid(mod2aiii))
qqnorm(resid(mod2aiii)) 
qqline(resid(mod2aiii))
anova(mod2aiii)
ranef(mod2aiii)
plot(mod2aiii)
boxplot(Count~Treatment,data=df)

#looks ok for 2018 t4 and t5 counts, so going to go ahead, 
#--drop interaction from model and rerun

mod2aa<-lmer(Count~Treatment+Sex.recollected+(1|Trial)+(1|Trial:Deployment.day),data=df)
hist(resid(mod2aa))
qqnorm(resid(mod2aa)) 
qqline(resid(mod2aa))
anova(mod2aa,type="II")
ranef(mod2aa)
plot(mod2aa)
boxplot(Count~Treatment,data=df)

mod2aiv<-glmer(Count~Treatment*Sex.recollected+(1|Trial)+(1|Trial:Deployment.day),family=poisson,data=df)
hist(resid(mod2aiv))
qqnorm(resid(mod2aiv)) 
qqline(resid(mod2aiv))
Anova(mod2aiv)
ranef(mod2aiv)
plot(mod2aiv)

#sqrt count
mod2aiii<-lmer(scount~Treatment*Sex.recollected+(1|Trial)+(1|Trial:Deployment.day),data=df)
hist(resid(mod2aiii))
qqnorm(resid(mod2aiii)) 
qqline(resid(mod2aiii))
anova(mod2aiii,type="I")
ranef(mod2aiii)
plot(mod2aiii)#slightly better for 2018 t4 and t5 data

#dropping int. term from model
mod2aai<-lmer(scount~Treatment+Sex.recollected+(1|Trial)+(1|Trial:Deployment.day),data=df)
hist(resid(mod2aai))
qqnorm(resid(mod2aai)) 
qqline(resid(mod2aai))
anova(mod2aai,type="II")
ranef(mod2aai)
plot(mod2aai)

mod2aiv<-glmer(scount~Treatment*Sex.recollected+(1|Trial)+(1|Trial:Deployment.day),family=poisson,data=df)
hist(resid(mod2aiv))
qqnorm(resid(mod2aiv)) 
qqline(resid(mod2aiv))
Anova(mod2aiv)
ranef(mod2aiv)
plot(mod2aiv)

#log + 1 counts
mod2aiii<-lmer(logcount~Treatment*Sex.recollected+(1|Trial)+(1|Trial:Deployment.day),data=df)
hist(resid(mod2aiii))
qqnorm(resid(mod2aiii)) 
qqline(resid(mod2aiii))
anova(mod2aiii)
ranef(mod2aiii)
plot(mod2aiii)

#dropping interaction term from model
mod2ax<-lmer(logcount~Treatment+Sex.recollected+(1|Trial)+(1|Trial:Deployment.day),data=df)
hist(resid(mod2ax))
qqnorm(resid(mod2ax)) 
qqline(resid(mod2ax))
anova(mod2ax)
ranef(mod2ax)
plot(mod2ax)

mod2aiv<-glmer(logcount~Treatment*Sex.recollected+(1|Trial)+(1|Trial:Deployment.day),family=poisson,data=df)
hist(resid(mod2aiv))
qqnorm(resid(mod2aiv)) 
qqline(resid(mod2aiv))
Anova(mod2aiv)
ranef(mod2aiv)
plot(mod2aiv)

#not sure what the best model will be for 2017 recollections, too many zeros

mod2b<-lmer(Count~Treatment*Sex.recollected+(1|Trial),data=df)
hist(resid(mod2b))
qqnorm(resid(mod2b)) 
qqline(resid(mod2b))
anova(mod2b)
ranef(mod2b)

mod2c<-lmer(Count~Treatment*Sex.recollected+(1|Trial)+,data=df)
hist(resid(mod2b))
qqnorm(resid(mod2b)) 
qqline(resid(mod2b))
anova(mod2b)
ranef(mod2b)

mod3<-lm(Count~Treatment*Sex.recollected+(1|Reef),data=df)#same results as mod2?
hist(resid(mod3))#more normal than before
qqnorm(resid(mod3)) #still not really that normal
qqline(resid(mod3))
anova(mod3)

mod4<-lm(Count~Treatment*Sex.recollected*Reef,data=df)
hist(resid(mod4))
qqnorm(resid(mod4)) 
qqline(resid(mod4))
anova(mod4)

#2018.t6 only analyses
mod2aiii<-lm(Count~Treatment*Sex.recollected,data=df)
hist(resid(mod2aiii))
qqnorm(resid(mod2aiii)) 
qqline(resid(mod2aiii))
anova(mod2aiii)
plot(mod2aiii)
boxplot(Count~Treatment,data=df)

#looks ok for 2018 t4 and t5 counts, so going to go ahead, 
#--drop interaction from model and rerun

mod2aa<-lm(Count~Treatment+Sex.recollected,data=df)
hist(resid(mod2aa))
qqnorm(resid(mod2aa)) 
qqline(resid(mod2aa))
anova(mod2aa)
plot(mod2aa)

#sqrt count
mod2aiii<-lm(scount~Treatment*Sex.recollected,data=df)
hist(resid(mod2aiii))
qqnorm(resid(mod2aiii)) 
qqline(resid(mod2aiii))
anova(mod2aiii)
plot(mod2aiii)#slightly better for 2018 t4 and t5 data

#dropping int. term from model
mod2aai<-lm(scount~Treatment+Sex.recollected,data=df)
hist(resid(mod2aai))
qqnorm(resid(mod2aai)) 
qqline(resid(mod2aai))
anova(mod2aai)
plot(mod2aai)
boxplot(scount~Treatment, data=df)

#log + 1 counts
mod2aiii<-lm(logcount~Treatment*Sex.recollected,data=df)#might be the best for testing of assumptions
hist(resid(mod2aiii))
qqnorm(resid(mod2aiii)) 
qqline(resid(mod2aiii))
anova(mod2aiii)
plot(mod2aiii)
boxplot(logcount~Treatment, data=df)

mod2aiv<-lm(logcount~Treatment+Sex.recollected,data=df)
hist(resid(mod2aiv))
qqnorm(resid(mod2aiv)) 
qqline(resid(mod2aiv))
anova(mod2aiv)
plot(mod2aiv)


# pick up here #####
#1. add in the deployment day as a random effect (can use dplyr to specify trials reefs for deployments days)
#2. try different distributions
#3. remove immatures from the data (not technically recollected)

bargraph.CI(x.factor = treatment, response = growth, group = final.sex, legend=TRUE, main="recollections by sex", data = reco.all)
#not surprising that the transitionals seemed to grow the most, bsed on fact that they change sex with increasing size

#plotting
#recollections by sex for each trial, combined dataset
bargraph.CI(x.factor = Trial, response = Count, group = Sex.recollected, legend=TRUE, main="All trials, recollection by sex",x.leg = 18, data = df)

#2017 t1-3
bargraph.CI(x.factor = Treatment, response = Count, group = Sex.recollected, legend=TRUE, main="2017 recollections by sex",x.leg = 10, data = df)

#2018 t4 and 5
#normal counts
bargraph.CI(x.factor = Treatment, response = Count, group = Sex.recollected, legend=TRUE, main="2018 trials 4 and 5 recollections by sex",x.leg = 10, data = df)
#found square root transformed model to be the best fit, assumptions-wise
bargraph.CI(x.factor = Treatment, response = scount, group = Sex.recollected, legend=TRUE, main="2018 trials 4 and 5 recollections by sex, sqrt trans",x.leg = 10, data = df)

#2018 t4 and 5
#normal counts
bargraph.CI(x.factor = Treatment, response = Count, group = Sex.recollected, legend=TRUE, main="2018 trial 6 recollections by sex",x.leg = 10, data = df)
#square root transformed count data
bargraph.CI(x.factor = Treatment, response = scount, group = Sex.recollected, legend=TRUE, main="2018 trial 6 recollections by sex, sqrt trans",x.leg = 10, data = df)
#found log + 1 transformed model to be the best fit, assumptions-wise
bargraph.CI(x.factor = Treatment, response = logcount, group = Sex.recollected, legend=TRUE, main="2018 trial 6 recollections by sex, log+1 trans",x.leg = 10, data = df)

#just looking at recollections without sex recollected as a factor
#2017
bargraph.CI(x.factor = Treatment, response = logcount, legend=TRUE, main="2017 trials 1-3 recollections, log+1 trans",x.leg = 10, data = df)
bargraph.CI(x.factor = Treatment, response = Count, legend=TRUE, main="2017 trials 1-3 recollections",x.leg = 10, data = df)

#2018.t4.5
bargraph.CI(x.factor = Treatment, response = logcount, legend=TRUE, main="2018 trial 4-5 recollections, log+1 trans",x.leg = 10, data = df)
bargraph.CI(x.factor = Treatment, response = Count, legend=TRUE, main="2018 trials 4-5 recollections",x.leg = 10, data = df)

#2018.t6
bargraph.CI(x.factor = Treatment, response = logcount, legend=TRUE, main="2018 trial 6 recollections, log+1 trans",x.leg = 10, data = df)
bargraph.CI(x.factor = Treatment, response = Count, legend=TRUE, main="2018 trial 6 recollections",x.leg = 10, data = df)



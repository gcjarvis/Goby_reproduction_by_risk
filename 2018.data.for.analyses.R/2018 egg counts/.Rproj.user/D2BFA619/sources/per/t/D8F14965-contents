################################################################
# sciplot code to get rough estimate of eggs by risk           #
#                                                              #
#    10/28/2018                                                #
################################################################

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
nests.risk<-read.csv("Data/10.29.18.eggs.csv")
head(nests.risk)

#adding a column to the data to log-transform +1 (for nests that had 0 eggs) the egg counts
nests.risk$log.egg<-(log(nests.risk$egg.count + 1))
View(nests.risk)

#model building with log-transformed data

log.mod1<-lm(log.egg~treatment, data=nests.risk)
hist(resid(log.mod1))#not so good
log.mod2<-lm(log.egg~treatment+week, data=nests.risk)
hist(resid(log.mod2))
log.mod3<-lm(log.egg~treatment*week, data=nests.risk)
hist(resid(log.mod3))
log.mod4<-lmer(log.egg~treatment+week+(1|deployment.day),data=nests.risk)
hist(resid(log.mod4))
log.mod5<-lmer(log.egg~treatment*week+(1|deployment.day), data=nests.risk)
hist(resid(log.mod5))
#log.mod6<-glmer(log.egg~treatment+week+(1|deployment.day), family=gamma, data=nests.risk)

anova(log.mod1,log.mod2,log.mod3)

#it doesn't seem like the log + 1 trans. data are normal



#want to see how number of eggs varies by treatment
#pooled data so I could get a total number of eggs laid per reef per treatment
#for example, if there were multiple TOLs with eggs I added them together
#I'm thinking that's the way to do it, instead of counting each TOL as a separate replicate, although I could do that

#want to see what production per week looks like
bargraph.CI(x.factor = week, group = treatment, legend = TRUE,xlab="Week", ylab="Average reproduction per nest", response = egg.count, data = nests.risk)
#seems like there's a lot of reproduction in the first two weeks, low and high look similar, with medium ramping up at the end

bargraph.CI(x.factor = treatment, response = egg.count,xlab="Risk treatment", ylab="Average reproduction per reef", data = nests.risk)
#doesn't look like there's a whole lot of difference in the number of eggs laid here. Will be interested to look at the average biomass per treatment

#now have standardized by average biomass (initial biomass/final biomass)

bargraph.CI(x.factor = week, group = treatment, legend = TRUE, response = eggs.avg.biomass, main="eggs per reef per week with biomass", data = nests.risk)
#similar trends as before, change in y axis

bargraph.CI(x.factor = treatment, response = eggs.avg.biomass, main="eggs per reef per avg. biomass", data = nests.risk)
#doesn't look like there's a whole lot of difference in the number of eggs laid here. Will be interested to look at the average biomass per treatment

bargraph.CI(x.factor = treatment, response = final.biomass, main="final biomass per reef", data = nests.risk)
#doesn't look like there's a whole lot of difference in the number of eggs laid here. Will be interested to look at the average biomass per treatment

#change in biomass by treatment
#positive values = increase in biomass
#negative values = decrese in biomass

bargraph.CI(x.factor = treatment, response = change.in.biomass, main="change in biomass per reef", data = nests.risk)
#not super informative here, and also probably not a grat way to represent the data visually

## now want to see if there are differences when I subset the data by 
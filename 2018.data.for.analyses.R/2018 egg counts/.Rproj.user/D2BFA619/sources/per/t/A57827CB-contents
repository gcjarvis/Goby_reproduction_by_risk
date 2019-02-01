######################################
#                                    #      
#   Egg counts with Trials 1-6 data  #
# 1/10/19                            #
#                                    #
######################################

rm(list=ls())

library(sciplot)
library(lme4)
library(car)
library(lmerTest)
library(dplyr)
library(ggplot2)

gob.sub<-read.csv("Data/List of nests to count.1.10.19.csv")
gob.sub$Week<-as.factor(gob.sub$Week)
gob.sub$Treatment<-ordered(gob.sub$Treatment,levels=c("Low","Medium","High","Control"))

bargraph.CI(x.factor = Week, response = Egg.count, group= Treatment, legend=TRUE, main="all trials", data = gob.sub)
bargraph.CI(x.factor = Treatment, response = Egg.count, group= Week, legend=TRUE, main="eggs by treatment", data = gob.sub)
bargraph.CI(x.factor = Treatment, response = Egg.count, legend=TRUE, main="eggs by treatment, no time", data = gob.sub)



egg.per.week<-gob.sub %>%
  group_by(Trial,Reef,Week,Treatment) %>%
  summarize(Egg.count = sum(Egg.count))
#not sure why this isn't working. Not going to spend more time on it

egg.per.week

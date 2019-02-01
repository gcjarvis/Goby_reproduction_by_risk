#########################################
# 1/3/2019                              #  
# Density Analyses                      #
# George Jarvis                         #
#########################################

#clear workspace
rm(list=ls())

library(sciplot)
library(lme4)
library(car)
library(ggplot2)
library(lmerTest)
library(dplyr)

gob.den<-read.csv("Data/density.1.3.19.csv")
gob.den$Trial<-as.factor(gob.den$Trial)
gob.den$Week<-as.factor(gob.den$Week)
gob.den.4<-filter(gob.den, Trial==4) # good way to subset data quickly
gob.den.5<-filter(gob.den, Trial==5) 
gob.den.6<-filter(gob.den, Trial==6)
gob.den.6$Trial<-as.factor(gob.den.6$Trial)
gob.den.6$Week<-as.factor(gob.den.6$Week)

#all together
lineplot.CI(Week,Density,group=Treatment,legend = TRUE,main="all trials", xlab="Week", ylab="Fish density per reef", data=gob.den)
#only trial 4
lineplot.CI(Week,Density,group=Treatment,legend = TRUE,main= "trial 4", xlab="Week", ylab="Fish density per reef", data=gob.den.4)
#only trial 5
lineplot.CI(Week,Density,group=Treatment,legend = TRUE,main= "trial 5", xlab="Week", ylab="Fish density per reef", data=gob.den.5)
#only 6
lineplot.CI(Week,Density,group=Treatment,legend = TRUE,main= "trial 6", xlab="Week", ylab="Fish density per reef", data=gob.den.6)

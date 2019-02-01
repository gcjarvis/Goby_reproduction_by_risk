#########################################
# 1/19/2019                             #  
# Temperature data                      #
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

temp<-read.csv("Data/2018.temp.data.working.1.19.20.csv")
temp$Trial<-as.factor(temp$Trial)
temp$Trial

#gob.den$Week<-as.factor(gob.den$Week)
#gob.den.4<-filter(gob.den, Trial==4) # good way to subset data quickly
#gob.den.5<-filter(gob.den, Trial==5) 
#gob.den.6<-filter(gob.den, Trial==6)
#gob.den.6$Trial<-as.factor(gob.den.6$Trial)
#gob.den.6$Week<-as.factor(gob.den.6$Week)

#as of 1.19.19, still need to go back and 1) remove the temperature data points from
#times after I removed the loggers from the reefs and 2) remove data after week two
#from trial 6

#it will be kind of hard to know the exact time that I removed them from the water,
#so to avoid any bias, I might just use data up until the morning of the last day
#that I shot photos of the reefs

#I think the resolution should be pretty good, even if I miss about 200 data points
#there will still be a ton of data to pull from

#I'm pretty sure I got the start times right, but it's interesting that the same didn't happen
#from end of trial 5 to start of trial 6

#all together
lineplot.CI(Week,Temp.F,group=Trial,legend = TRUE,main="Temp per trial over time", xlab="Week", ylab="temp (F)", data=temp)
#really cool - definitely got warmer from June -> August
lineplot.CI(Trial,Temp.F,group=Treatment,legend = TRUE,main="Temp per trial and trt", xlab="Trial", ylab="temp", data=temp)
#fairly consistent temps among treatments within each trial
lineplot.CI(Week,Temp.F,group=Treatment,legend = TRUE,main="Temp per trial and trt", xlab="week", ylab="temp", data=temp)
#kind of a useless figure

# Description: figuring out why density calculations are off
# Author: George C Jarvis
# Date: Fri Oct 25 10:23:39 2019
# Notes: I made a pivot table in excel and found the correct values
#     I think this has something to do with how I'm calculating means in my df
#     that I use for ggplot lineplot
# --------------

rm(list=ls())

#loading libraries
library(sciplot)
library(lme4)
library(car)
library(lmerTest)
library(dplyr)
library(ggplot2)
library(MASS)
library(nlme)
library(pwr)

#loading data, make sure not to do na.rm = ""...MY GOD!!!
#that was the issue
viz.surv<-read.csv("Data/density.2019.10.25.csv")
viz.surv$Density<-as.integer(viz.surv$Density)
viz.surv$den.max<-as.integer(viz.surv$den.max)
viz.surv$Day<-as.factor(viz.surv$Day)
viz.surv$Treatment<-ordered(viz.surv$Treatment,levels=c("Low","Medium","High"))


#quick plot
lineplot.CI(Day,den.max,group=Treatment,legend = TRUE,main="all trials", xlab="Day", ylab="Fish density per reef", data=viz.surv)
#gives me the same result as before, not sure why the values for day 1 aren't higher
#should be 14 for low, 13.8 for med, and 13.9 for high, but not less than 10

View(viz.surv)



den<-with(viz.surv, aggregate((den.max), list(Day=Day,Treatment=Treatment), mean))
den$se<-with(viz.surv, aggregate((den.max), list(Treatment=Treatment,Day=Day), 
                                 function(x) sd(x)/sqrt(length(x))))[,3]


lineplot.CI(Day,den.max,group=Treatment,legend = TRUE,main="all trials", xlab="Day", ylab="Fish density per reef", data=viz.surv)

#wanting to see how it looks as a bargraph by day
bargraph.CI(x.factor = Day, response = den.max, 
            group= Treatment, legend=TRUE, main="all trials, HR combined grouped by trial", 
            data = viz.surv)

#no grouping factor, seems like H>L=M, strange...with den.max
#goes against what Mark found...
bargraph.CI(x.factor = Treatment, response = den.max, 
            legend=TRUE, main="den.max, HR combined grouped by trial", 
            data = viz.surv)


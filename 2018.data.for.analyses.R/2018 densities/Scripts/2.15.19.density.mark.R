#########################################
# 2/15/2019                             #  
# Initial Density Analyses for mark     #
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

gob.den<-read.csv("Data/density.1.15.19.csv")
#gob.den$Trial<-as.factor(gob.den$Trial) #not able to subset data if you do this ahead of time???
gob.den$Week<-as.factor(gob.den$Week)
#den.avg.egg$Treatment<-ordered(den.avg.egg$Treatment,levels=c("Low","Medium","High","Control"))
gob.den$Treatment<-ordered(gob.den$Treatment,levels=c("Low","Medium","High","Control"))
#reordering the treatments doesn't change the order of the legend text, but it does control the order of the legend markers 
#it seems like line plot defaults to having the label with the highest value listed first (in this case "Medium")

#setting up df for different years (will also include all of the figures for 2018 trials separately)####
den.2017<-gob.den[gob.den$Trial<4, c("Trial","Day", "Reef", "Week", "Treatment", "Density")]
den.2017$Trial<-as.factor(den.2017$Trial) #need to do this for analyses
den.2017$Treatment<-ordered(den.2017$Treatment,levels=c("Low","Medium","High","Control"))
den.2018<-gob.den[gob.den$Trial>3, c("Trial","Day", "Reef", "Week", "Treatment", "Density")]
den.2018$Trial<-as.factor(den.2018$Trial)
den.2018$Treatment<-ordered(den.2018$Treatment,levels=c("Low","Medium","High","Control"))

#df names
#2a. 2017 + 2018 combo = gob.den
#2b. 2017 only = den.2017
#2c. 2018 only = den.2018

#2a. 2017 and 2018 combined
lineplot.CI(Week,Density,group=Treatment,legend = TRUE,main="2017 + 2018 densities by week and treatment", xlab="Week", ylab="Number of fish per reef (mean +/- se)", data=gob.den)
#2a.i consider doing barplots too
bargraph.CI(x.factor = Week, response = Density, group= Treatment, legend=TRUE, main="2017 + 2018 densities by week and treatment", xlab="Week", ylab="Number of fish per reef (mean +/- se)", x.leg=13, yleg=4500, data = gob.den)
#2a.ii pooled by treatment, densities over time
bargraph.CI(x.factor = Week, response = Density, legend=TRUE, main="2017 + 2018 density over time", xlab="Week", ylab="Number of fish per reef (mean +/- se)", x.leg=13, yleg=4500, data = gob.den)
#2a.iii pooled by week, densities by treatment
bargraph.CI(x.factor = Treatment, response = Density, legend=TRUE, main="2017 + 2018 densities by treatment", xlab="Week", ylab="Number of fish per reef (mean +/- se)", x.leg=13, yleg=4500, data = gob.den)
#2a.ix let's do by day too
lineplot.CI(Day,Density,group=Treatment,legend = TRUE,main="2017 + 2018 densities by day", xlab="Day", ylab="Number of fish per reef (mean +/- se)", data=gob.den)


#2b. 2017 only 
#not really useful to view by week for this trial...
#lineplot.CI(Week,Density,group=Treatment,legend = TRUE,main="2017 densities by week", xlab="Week", ylab="Number of fish per reef (mean +/- se)", data=den.2017)
#2b.i better to look at daily densities
lineplot.CI(Day,Density,group=Treatment,legend = TRUE,main="2017 densities by day and trial", xlab="Day", ylab="Number of fish per reef (mean +/- se)", data=den.2017)
#2b.ii barplot of densities by day, pooling by treatments
bargraph.CI(x.factor = Day, response = Density, legend=TRUE, main="2017 densities over time", xlab="Day", ylab="Number of fish per reef (mean +/- se)", x.leg=13, yleg=4500, data = den.2017)
#2b.iii barplot of densities by treatment, pooling by day
bargraph.CI(x.factor = Treatment, response = Density, legend=TRUE, main="2017 densities by treatment", xlab="Treatment", ylab="Number of fish per reef (mean +/- se)", x.leg=13, yleg=4500, data = den.2017)


#2c.line 2018 only, looks very similar to combined data I'd imagine
lineplot.CI(Week,Density,group=Treatment,legend = TRUE,main="2018 densities by week and treatment", xlab="Week", ylab="Number of fish per reef (mean +/- se)", data=den.2018)
#2c.i bargraph by week
bargraph.CI(x.factor = Week, response = Density, group= Treatment, legend=TRUE, main="2018 densities by week and treatment", xlab="Week", ylab="Number of fish per reef (mean +/- se)", x.leg=13, yleg=4500, data = den.2018)
#2c.ii bargraph of density by week, pooling by treatment
bargraph.CI(x.factor = Week, response = Density, legend=TRUE, main="2018 densities by week", xlab="Week", ylab="Number of fish per reef (mean +/- se)", x.leg=13, yleg=4500, data = den.2018)
#2c.iii bargraph of density by treatment, pooling by week
bargraph.CI(x.factor = Treatment, response = Density, legend=TRUE, main="2018 densities by treatment", xlab="Treatment", ylab="Number of fish per reef (mean +/- se)", x.leg=13, yleg=4500, data = den.2018)
#2c.ix.line lineplot by day, want to put this next to 2017 only to show the difference in density over time with diff. tagging methods
lineplot.CI(Day,Density,group=Treatment,legend = TRUE,main="2018 densities by day and treatment", xlab="Day", ylab="Number of fish per reef (mean +/- se)", data=den.2018)
#NOTE: I did not include any other figures after this in the update for Mark
#2c.x.bar barplot by day, want to put this next to 2017 only to show the difference in density over time with diff. tagging methods
#bargraph.CI(x.factor = Day, response = Density, group= Treatment, legend=TRUE, main="2018 densities by day and treatment", xlab="Week", ylab="Number of fish per reef (mean +/- se)", x.leg=115, yleg=4500, data = den.2018)
#2c.xi.bar.day barplot by day, pooled by treatment
#bargraph.CI(x.factor = Day, response = Density, legend=TRUE, main="2018 densities by day", xlab="Week", ylab="Number of fish per reef (mean +/- se)", x.leg=115, yleg=4500, data = den.2018)
#2c.xii.bar.trt barplot by treatment, pooled over time
#bargraph.CI(x.factor = Treatment, response = Density, legend=TRUE, main="2018 densities by treatment", xlab="Treatment", ylab="Number of fish per reef (mean +/- se)", x.leg=115, yleg=4500, data = den.2018)
#2c.x bargraph of density by day, not a useful plot, might also be a dup
#bargraph.CI(x.factor = Day, response = Density, group=Treatment, legend=TRUE, main="2018 densities by day", xlab="Day", ylab="Number of fish per reef (mean +/- se)", x.leg=115, yleg=4500, data = den.2018)

#only trial 4
lineplot.CI(Week,Density,group=Treatment,legend = TRUE,main= "trial 4", xlab="Week", ylab="Fish density per reef", data=gob.den.4)
#only trial 5
lineplot.CI(Week,Density,group=Treatment,legend = TRUE,main= "trial 5", xlab="Week", ylab="Fish density per reef", data=gob.den.5)
#only 6
lineplot.CI(Week,Density,group=Treatment,legend = TRUE,main= "trial 6", xlab="Week", ylab="Fish density per reef", data=gob.den.6)

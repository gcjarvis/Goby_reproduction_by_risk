#######################
# 2.20.19             #  
# Behavioral analyses #
# George Jarvis       #
#######################

#clear workspace
rm(list=ls())

library(sciplot)
library(lme4)
library(car)
library(lmerTest)
library(dplyr)
library(ggplot2)

#reading in entire raw datasheet (2017 + 2018 trials)
behave<-read.csv("Data/behavior.data.2.20.19.csv")
#ordering treatments for better graphical visualization
behave$Treatment<-ordered(behave$Treatment,levels=c("Low","Medium","High","Control"))

#subsetting data, only 2017
behave.2017<-behave[behave$Year==2017, c("Year","Trial", "Reef", "Week","Days.after.deployment","Observer","Density","Goby.ID", "Treatment", "Density","Gravidity","Nearest.TOL","Starting.location",
                                      "Starting.position","Total.time.exposed","proportion.exposed","Total.time.hidden","proportion.hidden","num.movement.swims","num.bites","num.times.chased",
                                      "chased.by.size.tag","num.times.chaser","chaser.of.size.tag","area.used.cm")]
#including data from 2018
behave.2018<-behave[behave$Year==2018, c("Year","Trial", "Reef", "Week","Days.after.deployment","Observer","Density","Goby.ID", "Treatment", "Density","Gravidity","Nearest.TOL","Starting.location",
                                      "Starting.position","Total.time.exposed","proportion.exposed","Total.time.hidden","proportion.hidden","num.movement.swims","num.bites","num.times.chased",
                                      "chased.by.size.tag","num.times.chaser","chaser.of.size.tag","area.used.cm")]
#preliminary results for Mark
#df names:
#a. behave: 2017 + 2018 data
#b. behave.2017: 2017 data
#c. behave.2018: 2018 data

#a. 2017 + 2018 data combined
#exposure time
#a.i proportion of time exposed by treatment, might be issues with ss, esp for control treatment
bargraph.CI(x.factor = Treatment, response = proportion.exposed, main="Exposure time by risk, 2018 + 2017 data", xlab="Treatment", ylab="Proportion of time exposed (mean +/- se)", data = behave)
#a.ii proportion of time exposed by treatment over time
bargraph.CI(x.factor = Week, response = proportion.exposed, group= Treatment, legend=TRUE, main="Exposure time by risk over time, 2018 + 2017 data", xlab="Week", ylab="Proportion of time exposed (mean +/- se)", x.leg=14, yleg=4500, data = behave)
### doesn't seem like exposure time changed over time, might not be the best way to look at it, also just might be unnecessary
#a.iii number of bites
bargraph.CI(x.factor = Treatment, response = num.bites, main="Number of bites by risk, 2018 + 2017 data", xlab="Treatment", ylab="Number of bites (mean +/- se)", data = behave)
#a.iv number of movement swims
bargraph.CI(x.factor = Treatment, response = num.movement.swims, main="Number of movement swims, 2018 + 2017 data", xlab="Treatment", ylab="Number of movemements (mean +/- se)", data = behave)

#b. 2017 data only
#b.i proportion of time exposed by treatment, might be issues with ss, esp for control treatment
bargraph.CI(x.factor = Treatment, response = proportion.exposed, main="Exposure time by risk, 2017 data only", xlab="Treatment", ylab="Proportion of time exposed (mean +/- se)", data = behave.2017)
#b.ii proportion of time exposed by treatment over time
bargraph.CI(x.factor = Week, response = proportion.exposed, group= Treatment, legend=TRUE, main="Exposure time by risk over time, 2017 data only", xlab="Week", ylab="Proportion of time exposed (mean +/- se)", x.leg=4, yleg=4500, data = behave.2017)
#b.iii number of bites
bargraph.CI(x.factor = Treatment, response = num.bites, main="Number of bites by risk, 2017 data only", xlab="Treatment", ylab="Number of bites (mean +/- se)", data = behave.2017)
#b.iv number of movement swims
bargraph.CI(x.factor = Treatment, response = num.movement.swims, main="Number of movement swims, 2017 data only", xlab="Treatment", ylab="Number of movemements (mean +/- se)", data = behave.2017)

#c. 2018 data only
#exposure time
#c.i proportion of time exposed by treatment, might be issues with ss, esp for control treatment
bargraph.CI(x.factor = Treatment, response = proportion.exposed, main="Exposure time by risk, 2018 data only", xlab="Treatment", ylab="Proportion of time exposed (mean +/- se)", data = behave.2018)
#c.ii proportion of time exposed by treatment over time
bargraph.CI(x.factor = Week, response = proportion.exposed, group= Treatment, legend=TRUE, main="Exposure time by risk over time, 2018 data only", xlab="Week", ylab="Proportion of time exposed (mean +/- se)", x.leg=14, yleg=4500, data = behave.2018)
### doesn't seem like exposure time changed over time, might not be the best way to look at it, also just might be unnecessary
#c.iii number of bites
bargraph.CI(x.factor = Treatment, response = num.bites, main="Number of bites by risk, 2018 data only", xlab="Treatment", ylab="Number of bites (mean +/- se)", data = behave.2018)
#c.iv number of movement swims
bargraph.CI(x.factor = Treatment, response = num.movement.swims, main="Number of movement swims, 2018 data only", xlab="Treatment", ylab="Number of movemements (mean +/- se)", data = behave.2018)

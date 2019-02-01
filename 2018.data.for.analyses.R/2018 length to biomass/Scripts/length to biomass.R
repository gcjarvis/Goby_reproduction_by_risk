#################################################
# Length to biomass curve for bluebanded gobies #                                                              #
#    7/26/2018                                  #
#   George Jarvis                               #
#################################################


#7/26/2018, took raw data in Excel file, removed the notes, and gravidity from data
#resaved as .csv file to work with

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
length.biomass<-read.csv("Data/natural ratio.csv")#603 rows of data
#only including cases where there is both length and biomass
length.biomass<- length.biomass[complete.cases(length.biomass), ]#426 rows of data

head(length.biomass)
plot(length.biomass$Size, length.biomass$Weight) #looks like apretty tight relationship
#definitely lacking data points for big gobies, but they are much less common, so that's understandable
#also might want to break down to gravid/non-gravid, and see what that looks like
#would be interesting to break down the data by sex (male/female) and also by gravids 
#I'm imagining three colored line curves on the same graph

lb.mod<-lm(Weight~Size, data=length.biomass)
anova(lb.mod)
summary(lb.mod)# 0.87 R-squared value for the data frame that contains all of the individuals
#would be interested to see how breaking down by sex would alter the curve
#will have to go back to the raw data and add in data for gravids

nests.risk$Trial<-as.factor(nests.risk$Trial)#made trial a factor
head(nests.risk)
# Description: combining data on reproduction and biomass/density recollected to send to Dustin
# Author: George C Jarvis
# Date: Tue Jul 14 11:46:00 2020
# Notes: need to add a few columns: duration of trial (days, 1= 14, 2 and 3= 7,4 = 28, 5= 27, 6=14),
#       density of females recollected, density total recollcted (mature males + females), total biomass female,
#       total biomass for the whole reef, reproduction per capita per day (per female density), and reproduction per female biomass per day
# --------------

rm(list=ls())

library(lme4)
library(lmerTest)
library(dplyr)
library(ggplot2)
library(emmeans) #for generating least-squares adjusted means from models (will likely help when I have to
## back-transform means if data transformation is needed)

#importing dataset####
repro<-read.csv("Data/goby_reproduction.2020.4.9.csv")

#intial data viz
pairs(repro)

#data manipulation####

#adding column for average number of inhabitants throughout the trial, rounded to nearest whole number of fish
repro$avg.inhab<-(ceiling((repro$Recollection+20)/2))

#adding egg/week variable to the dataset for comparisons among trials of differing length
repro<-repro %>%
  mutate(egg.week = ifelse(Trial<4, Egg.count/1,
                           ifelse(Trial == 4| Trial == 5, (ceiling(Egg.count/4)),
                                  ifelse(Trial == 6, (ceiling(Egg.count/2)), NA))))

#adding a column for year (as a proxy for tagging procedure), where trials 1-3 = 2017, and 4-6 = 2018
repro$Year <- ifelse(repro$Trial <=3, 2017, 2018)
#making Year and Trial factors
repro$Year<- as.factor(repro$Year)
repro$Trial<-as.factor(repro$Trial)

#sqrt-transforming egg counts to better satisfy assumptions
repro$sqrt.egg.week<-sqrt(repro$egg.week)

#subsetting data from Trial 6 for comparison of caged vs. uncaged treatments
#subset trial
repro.t6<-repro[repro$Trial==6,]
#subset treatments, caged and uncaged only to test for cage effects
## Note that Treatment factor for trial 6 models should be replaced with "T6.comparison"
repro.t6<-repro.t6[repro.t6$T6.comparison == "High" | repro.t6$T6.comparison == "Uncaged", ]

#exporting wrangled data for those that are not working with these data in R
#Data from Trials 1-5
#write.csv(repro,"Data\\egg_counts_after_data_wrangling.2020.4.7.csv", row.names = FALSE)
#Data from Trial 6 only
#write.csv(repro.t6,"Data\\Cage_effects_egg_counts_after_data_wrangling.2020.4.7.csv", row.names = FALSE)

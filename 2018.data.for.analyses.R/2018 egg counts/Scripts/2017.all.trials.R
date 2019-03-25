######################################
#                                    #      
#   Egg counts with 2017 data        #
#   all trials from 2017             #
# 2.23.19                            #
######################################

rm(list=ls())

library(sciplot)
library(lme4)
library(car)
library(lmerTest)
library(dplyr)
library(ggplot2)

#A. load data (egg counts and densities combined into one for raw counts and per capita)
#raw data for egg counts
#gob.sub<-read.csv("Data/jarvis.egg.count.data.19.2.22.csv")
#subset for just trials 1-3 from 2017
#egg.t123<-gob.sub[gob.sub$Trial<4,]
#thinking I might want to include daily counts (day 4 and day 7)

#Egg counts from 2017
gob.2017.raw<-read.csv("Data/2017.egg_counts_raw_data_19.2.23.csv")
#need to figure out the best way to account for density in my models
View(gob.2017.raw)

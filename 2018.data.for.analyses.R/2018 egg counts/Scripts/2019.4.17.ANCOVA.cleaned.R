# Description: ANCOVA results with visualization
# Author: George C Jarvis
# Date: Wed Apr 17 19:46:59 2019
# --------------

rm(list=ls())

library(sciplot)
library(lme4)
library(car)
library(lmerTest)
library(dplyr)
library(ggplot2)
library(MASS)
library(nlme)
library(pwr)
library(HH)#for ancova and plots

#data import

#2017.t.1.2.3
egg.den.bio<-read.csv("Data/new.data.2019.4.16a.no.t6.csv", na.strings = "") #uses adjusted counts for density
egg.den.bio<-na.omit(egg.den.bio)
egg.2017.t1.2.3<-egg.den.bio

#just t4.5
egg.den.bio<-read.csv("Data/new.data.2019.4.16a.t1.5.csv", na.strings = "") #uses adjusted counts for density
egg.den.bio<-na.omit(egg.den.bio)
egg.2018.t4.5<-egg.den.bio[(egg.den.bio$Trial>3) & (egg.den.bio$Trial<6), ]

#2018.t6
egg.den.bio<-read.csv("Data/new.data.2019.4.16.csv", na.strings = "") 
egg.2018.t6<-egg.den.bio[(egg.den.bio$Trial>5),]

#visualization

#2017
mod <- ancova(Egg.count ~ Treatment*Density, data=egg.2017.t1.2.3)
pred <- predict(mod)
#View(pred)

#plot
ggplot(data = cbind(egg.2017.t1.2.3, pred),
       aes(Density, Egg.count, color=Treatment)) + geom_point()  + stat_smooth(method="lm") +
  facet_grid(. ~ Treatment) + geom_line(aes(y=Density))

#2018.t.4.5
mod <- ancova(Egg.count ~ Treatment*Density, data=egg.2018.t4.5)
pred <- predict(mod)

#plot
ggplot(data = cbind(egg.2018.t4.5, pred),
       aes(Density, Egg.count, color=Treatment)) + geom_point()  + stat_smooth(method="lm") +
  facet_grid(. ~ Treatment) + geom_line(aes(y=Density))

#2018.t.6
mod <- ancova(Egg.count ~ Treatment*Density, data=egg.2018.t6)
pred <- predict(mod)

#plot
ggplot(data = cbind(egg.2018.t6, pred),
       aes(Density, Egg.count, color=Treatment)) + geom_point()  + stat_smooth(method="lm") +
  facet_grid(. ~ Treatment) + geom_line(aes(y=Density))
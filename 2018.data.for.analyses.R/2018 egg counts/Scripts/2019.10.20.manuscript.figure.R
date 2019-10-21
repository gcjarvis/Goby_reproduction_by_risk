# Description: plotting ANCOVA with all data combined into single plot
# Author: George C Jarvis
# Date: Sun Oct 20 07:56:37 2019
# Notes: will be the primary figure in my MS, no need to facet by exp. anymore
# --------------

rm(list=ls())

#loading packages
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

#loading data
repro<-read.csv("Data/new.data.2019.9.30.csv", na.strings = "")

#adding column for average density, rounded to nearest whole number of fish
repro$avg.inhab<-(ceiling((repro$Recollection+20)/2)) #looks a little weird in the plot
repro$avg.inhab.no.round<-((repro$Recollection+20)/2) #looking at it this way too


#plotting as ancova with all trials combined

#2018.t4.5 no facet, no trial

#setting up df
repro$Treatment<-ordered(repro$Treatment,levels=c("Low","Medium","High"))
egg.anc<-with(repro, aggregate((Egg.count), list(Treatment=Treatment), mean))
egg.anc$se<-with(repro, aggregate((Egg.count), list(Treatment=Treatment), 
                                        function(x) sd(x)/sqrt(length(x))))[,2]

png(filename = "Output/19.10.20.combined.ancova.color.not.rounded.reco.png", width = 800, height = 700)

#plot seems a little more realistic when I don't round to the nearest whole number (ceiling)

anc<-ggplot(repro, aes(avg.inhab.no.round, Egg.count, color=Treatment, fill=Treatment)) +
  geom_smooth(method="lm", se=FALSE, show.legend = FALSE) +
  #facet_grid(. ~ Trial,labeller=labeller(Trial = labels.t4.5))  +
  geom_point(size=3) + 
  theme_classic() + 
  #theme(strip.text = element_text(face="bold", size = 30))+
  xlab("Average Number of Gobies") +
  ylab("Reproduction per Reef") +
  expand_limits(y=0) +
  theme(axis.text.x=element_text(size=30, colour="black"),axis.text.y=element_text(size=30, colour="black"), axis.title=element_text(size=35,face="bold")) +
  theme(axis.title.y = element_text(size= 35, margin = margin(t = 0, r = 20, b = 0, l = 0)), 
        axis.title.x = element_text(margin = margin(t = 25, r = 0, b = 0, l = 0)), 
        axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))
anc + scale_color_manual(values=c("#0072B2","#009E73","#D55E00")) +
  theme(legend.text=element_text(size=25)) + 
  theme(legend.title =element_text(size=30, face="bold"))
  
dev.off()
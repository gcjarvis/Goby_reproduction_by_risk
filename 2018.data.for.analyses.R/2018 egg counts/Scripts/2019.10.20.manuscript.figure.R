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


#plotting as ancova with all trials combined COLOR...not for MS

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

#replotting with new dimensions and quality of image. NO COLOR. For MS
#I like the 12x4 size for this one
png("Output/2019.11.10.test.300dpi.12x4.egg.no.round.png", width = 12, height = 4, units = 'in', res = 300)
#taller figure to see the distinction in recollections
png("Output/2019.11.10.egg.10x6.300dpi.egg.rounded.png", width = 9, height = 5, units = 'in', res = 300)
png("Output/2019.11.12.treatments.ordered.egg.9.5x5.5.300dpi.egg.rounded.png", width = 9.5, height = 5.5, units = 'in', res = 300)

#building from the bottom, not including shape by treatment
anc1<-ggplot(repro, aes(avg.inhab, Egg.count, shape=Treatment, linetype=Treatment, col=Treatment)) +
  geom_smooth(method="lm", se=FALSE, show.legend = TRUE)  +
  geom_point(size=3)+
  theme_classic()+
  labs(x="Gobies per Reef", y="Total Egg Production Per Reef")+
  expand_limits(y=0)+
  scale_color_manual(values=c("black", "#666666", "grey"))+
  scale_linetype_manual(values=c("solid", "dashed", "twodash"))+
  theme(axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"), 
        axis.title=element_text(size=20))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), 
        axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0)), 
        axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(legend.text=element_text(size=18)) +
  theme(legend.title =element_text(size=20))+
  scale_x_continuous(breaks=c(10,12,14,16,18,20))+
  labs(color  = "Perceived Risk", linetype = "Perceived Risk", shape = "Perceived Risk")
anc1
#sort of lame that I had to code it like this to get the legend correc
dev.off()

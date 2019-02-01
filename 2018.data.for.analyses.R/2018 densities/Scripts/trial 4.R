#######################
# 6/20/2018           #  
# Density Analyses    #
# George Jarvis       #
#######################

#clear workspace
rm(list=ls())

#quick density analyses, with just density data from 6/16,6/17, and 6/20
getwd()
gob.den<-read.csv("Data/2018.density.surveys.csv")
head(gob.den)

library(sciplot)
bargraph.CI(data = gob.den, #specify data frame
            x.factor = Risk, #factor on x axis
            response = Density)
#no effect of risk on density
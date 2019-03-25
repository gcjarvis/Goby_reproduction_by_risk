# Description: Making quick tables to export to excel, then word
# Author: George C Jarvis
# Date: Mon Mar 18 21:45:24 2019
# --------------

#Description: Neat, quick way to export tables from R with a more professional look
#Source: https://www.youtube.com/watch?v=0A_MsElAJYo

#Exporting Anova table
#example data from PTL
ptl<-read.csv("Data/2019.3.18.ptl.csv")
#make a df of the output you want, model coefficients, for example
mod1<-lm(Count~Treatment, data=ptl)
mod.anova<-summary(mod1)
write.table(mod.anova,"clipboard",sep="\t") #copies this

#2019.3.18: can't get write.table to work with lm or aov output...

#stargazer instead, work pretty well for now
library(stargazer)
stargazer(mod.anova,type = "text")


stargazer(mydata, type = "text", title="Descriptive statistics", digits=1, out="table1.txt")
stargazer(mod1, type = "text", title="Descriptive statistics", digits=1)

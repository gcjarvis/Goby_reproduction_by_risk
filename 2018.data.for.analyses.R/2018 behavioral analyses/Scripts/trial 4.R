#######################
# 6/10/18             #  
# Behavioral analyses #
# George Jarvis       #
#######################

#clear workspace
rm(list=ls())
go.b<-read.csv("Data/behavior.data.6.10.18.csv")
head(go.b)

#quick plotting with sciplot
library(sciplot)
#png(file = "graphs/combined 2018.png", height = 300, width = 700) #this opens your graphics device
#par(family='Century Gothic', mar=c(5.1,6.1,1.1,1.1))
bargraph.CI(data = go.b, #specify data frame
            x.factor = risk, #factor on x axis
            response = num.movement.swims), #response variable, just run the code to this point if you want to see quick,dirty results

            #group = Crab, #indexing factor--for legend
            #col = c("turquoise", "turquoise4"), #color in order of legend
            space = c(0,0.3), #spaces b/w bars--0within levels, 0.3 b/w levels of x factor
            legend = TRUE,x.leg = 0.1 , #null is FALSE, x.leg specifies legend placement on x axis
            #leg.lab = c("Low", "Medium","High"), #specifies labels in legend for grouping variable in alphabetical order
            ylim = c(0,12), 
            #names.arg = c("High", "Low", "Medium")) #names of x axis levels
abline(h=0) 
dev.off()
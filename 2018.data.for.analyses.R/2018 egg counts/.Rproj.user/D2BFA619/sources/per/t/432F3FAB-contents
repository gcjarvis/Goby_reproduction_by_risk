################################################################
# sciplot code to get rough estimate of egg laying by risk     #
#                                                              #
#    6/25/2018                                                 #
################################################################

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
nests.risk<-read.csv("Data/nest_counts_prelim_6.25.18.csv")
head(nests.risk)

#subsetting the data for only position "back" nests b/c/ densities held for longer
nests.risk<-nests.risk[nests.risk$Block=="back",]
head(nests.risk)

#linear models
###############average number of nests by treatment##############
nest.trt<-lm(Nest.count~Treatment, data=nests.risk)
hist(resid(nest.trt)) #looks good
qqnorm(resid(nest.trt))
qqline(resid(nest.trt))
plot(nest.trt)
boxplot(Nest.count~Treatment, data=nests.risk)

summary(nest.trt) # significant difference in number of nests by treatment (p=0.042)

#plotting nests by risk
bargraph.CI(x.factor = Treatment, response = Nest.count, main="Nests per reef", data = nests.risk)
#looks like the sig. difference is driven by decreased nests in the medium risk treatment

############average number of nests by treatment*week######
#not sure if week would be a fixed or random effect, will try both
#week as a fixed factor, two-way crossed factor anova
nest.trt.week.fixed<-lm(Nest.count~Treatment*Week, data=nests.risk)
hist(resid(nest.trt.week.fixed))
qqnorm(resid(nest.trt.week.fixed))
qqline(resid(nest.trt.week.fixed))
summary(nest.trt.week.fixed) #no interaction, would run as separate ANOVAs

bargraph.CI(x.factor = Treatment, response = Nest.count, group= Week, legend= TRUE, main="Nests per reef per week", data = nests.risk)

#week as a random factor
nest.trt.week.ran<-lmer(Nest.count~Treatment+(1|Week), data=nests.risk)
hist(resid(nest.trt.week.ran))
qqnorm(resid(nest.trt.week.ran))
qqline(resid(nest.trt.week.ran))
summary(nest.trt.week.ran)
fixef(nest.trt.week.ran)

Anova(nest.trt.week.ran, type="II")
Anova(nest.trt.week.ran, Type="III")

nests.week<-lm(Nest.count~Week, data=nests.risk)
summary(nests.week)
bargraph.CI(x.factor = Week, response = Nest.count, legend= TRUE, main="Nests per week", data = nests.risk)
#no real difference in nest production over time, going to try and group by risk (re-visualizing the previous sciplot)

bargraph.CI(x.factor = Week, response = Nest.count, group= Treatment, legend= TRUE, main="Nests per week by treatment", data = nests.risk)
#not seeing any of the weeks have the trend that I saw last summer...
#medium risk seems to be producing the fewest eggs, but that is based off of nest counts
#this DOESN'T ACCOUNT FOR NESTS THAT WERE THERE MULTIPLE WEEKS IN A ROW
#Have to go back and figure out which nests were actually new

bargraph.CI(x.factor = Block, response = Nest.count, group= Treatment, legend= TRUE, main="Nests by treatment and position", data = nests.risk)

bargraph.CI(x.factor = Block, response = Nest.count, legend= TRUE, main="Nests by treatment and position", data = nests.risk)

nest.week<-lmer(Nest.count~Treatment+(1|Block), data=nests.risk)
summary(nest.week)
ranef(nest.week)

Anova(nest.week, type="II")
Anova(nest.week, type="III") #more reproduction in back reefs than front reefs

#########running with the same analyses I ran with my previous exp###########

#model building with nonlinear distributions
#nest.poiss<-glmer(Egg.count~Treatment+(1|Position)+(1|Trial), family=poisson, data=egg.day1)
#summary(mod.poiss)
#fixef(mod.poiss)

nest.poiss<-glmer(Nest.count~Treatment+(1|Block)+(1|Week), family=poisson, data=nests.risk)
summary(nest.poiss)
fixef(nest.poiss)
ranef(nest.poiss)

plot(nest.poiss)
hist(resid(nest.poiss))

qqnorm(resid(nest.poiss))
qqline(resid(nest.poiss))

Anova(nest.poiss, type="III") #for right-skewed data
Anova(nest.poiss, type="II") #test that I used for my previous analyses, the chi-squared value is much lower
# than when my level of replication was each TOL (went from 4383.9 (with TOL rep) to 691.55 (with reef rep))
# Surprisingly, it seems like the treatment is still significant, but this will be much harder to explain to an audience
# it seems like there is a lot of variation that is explained by the position and trial



#plotting
bargraph.CI(data = Data2018, #specify data frame
            x.factor = Measurement, #factor on x axis
            response = Length, #response variable
            group = Crab, #indexing factor--for legend
            legend = TRUE)

bargraph.CI(data = Data2018, #specify data frame
            x.factor = Measurement, #factor on x axis
            response = Length, #response variable
            group = Crab, #indexing factor--for legend
            legend = TRUE)

#pumps out a plot with means and SE, but also shows the customization aspects for the figure
#this code is almost the same as the barplot funciton in R, but it takes out the second step of
#...calculating SE and applying it to the barplot function
png(file = "graphs/combined 2018.png", height = 300, width = 700) #this opens your graphics device
par(family='Century Gothic', mar=c(5.1,6.1,1.1,1.1))
bargraph.CI(data = Data2018, #specify data frame
            x.factor = Measurement, #factor on x axis
            response = Length, #response variable
            group = Crab, #indexing factor--for legend
            col = c("turquoise", "turquoise4"), #color in order of legend
            space = c(0,0.3), #spaces b/w bars--0within levels, 0.3 b/w levels of x factor
            legend = TRUE,x.leg = 0.1 , #null is FALSE, x.leg specifies legend placement on x axis
            leg.lab = c("No Crab", "Crab"), #specifies labels in legend for grouping variable in alphabetical order
            ylim = c(0,12), 
            names.arg = c("Aperture Width", "Height", "Width")) #names of x axis levels
abline(h=0) 
dev.off()


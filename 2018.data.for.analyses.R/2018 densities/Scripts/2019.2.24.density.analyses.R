#########################################
# 2/24/2019                             #  
# Trials 4 and 5 density analyses       #
# George Jarvis                         #
#########################################

#clear workspace
rm(list=ls())

library(sciplot)
library(lme4)
library(car)
library(ggplot2)
library(lmerTest)
library(dplyr)
library(tidyr)
library(pwr)

gob.den<-read.csv("Data/density.2019.3.6.csv") #contains max.density counts now as well. 33% of all observations were
gob.den<-mutate(gob.den, Caging=ifelse(Treatment=="High","Uncaged",
                                                    ifelse(Treatment=="Control","Uncaged","Caged")))
gob.den$Caging<-as.factor(gob.den$Caging)
#--changed to reflect maximum densities on the reefs (den.max)
#gob.den$Trial<-as.factor(gob.den$Trial) #not able to subset data if you do this ahead of time???
#gob.den$Week<-as.factor(gob.den$Week)
#den.avg.egg$Treatment<-ordered(den.avg.egg$Treatment,levels=c("Low","Medium","High","Control"))
gob.den$Treatment<-ordered(gob.den$Treatment,levels=c("Low","Medium","High","Control"))
#reordering the treatments doesn't change the order of the legend text, but it does control the order of the legend markers 
#it seems like line plot defaults to having the label with the highest value listed first (in this case "Medium")

#setting up df for different years (will also include all of the figures for 2018 trials separately)####
den.t1.2.3<-gob.den[gob.den$Trial<4, ]
den.t1.2.3$Trial<-as.factor(den.t1.2.3$Trial) #need to do this for analyses
den.t1.2.3$Treatment<-ordered(den.t1.2.3$Treatment,levels=c("Low","Medium","High","Control"))
den.t4.5<-gob.den[(gob.den$Trial>3)&(gob.den$Trial<6), ]
den.t4.5$Trial<-as.factor(den.t4.5$Trial)
den.t6<-gob.den[gob.den$Trial>5,]
Treatment<-ordered(den.2018$Treatment,levels=c("Low","Medium","High","Control"))

#calculating means for abstracts
#price.cut <- diamonds %>%
  #group_by(cut) %>%
  #summarize(Mean = mean(price, na.rm=TRUE))

den.means.trt<-df %>%
  group_by(Treatment) %>%
  summarize(den.mean=mean(Density, na.rm=TRUE))
den.means.trt

#2017 - high risk densities were 19% lower overall in comp. to low and med
#2018 t4.5 - ""21% less in comp. to high
#2018 t6. - no diff in number of fish seen by treatment

#df assignments for analyses
df<-den.t1.2.3
df<-den.t4.5
df<-den.t6
#just looking at week 1 for trial 6
df<-den.t6[(den.t6$Week==1),]

#adding a column with cage vs uncaged distinstion for T6, I think I would be able to see
#differences there if it wasn't muddled up with HR and CC treatments being so similar

#notice how I did this here, was able to use dplyr
#-- to enter two contitional "ifelse" statements in the same pipe
#-- Had to make sure that the last condition (if trt not equal to "High" or "Control")
#-- equal to the value "Uncaged". Had to put it at the end of the string, 
#-- not after each "ifelse" statement #proud

#gob.den <- gob.den %>% 
gob.den.cage.uncaged<-mutate(gob.den, Caging=ifelse(Treatment=="High","Uncaged",
                       ifelse(Treatment=="Control","Uncaged","Caged")))
View(gob.den)



#model building
mod1<-aov(Density~Caging, data=df)
hist(resid(mod1))
qqnorm(resid(mod1))
qqline(resid(mod1))
anova(mod1)
summary(mod1)
TukeyHSD(mod1,"Caging", ordered=TRUE)
plot(TukeyHSD(mod1,"Caging"))

#this is what I will need to do for trial 1-3, 2017, week isn't great becuase it was only 1
mod2<-aov(Density~Treatment*Day, data=df)
hist(resid(mod2))
qqnorm(resid(mod2))
qqline(resid(mod2))
anova(mod2)
summary(mod2)
TukeyHSD(mod2,"Treatment",ordered=TRUE)
plot(TukeyHSD(mod2,"Treatment")) #no diff, must be a normality thing

mod2<-lm(Density~Treatment*Week, data=df)
hist(resid(mod2))
qqnorm(resid(mod2))
qqline(resid(mod2))
anova(mod2)
summary(mod2)#for trial 6, no sig. interaction term
plot(mod2)
boxplot(Density~Treatment, data=df) #not exactly equal variance for t6

#for 2018, trial 6, removing interaction term and reanalyzing
mod2i<-lm(Density~Treatment+Week, data=df)
hist(resid(mod2i))
qqnorm(resid(mod2i))
qqline(resid(mod2i))
anova(mod2i)
summary(mod2i)
plot(mod2i)

#consider looking at the outlying data to see if that makes a difference

#mod2ii<-lm(Density~Treatment*Week+(1|Week:Day), data=df)#didn't work b/c not every day exists
#--within each week
#hist(resid(mod2ii))
#qqnorm(resid(mod2ii))
#qqline(resid(mod2ii))
#anova(mod2ii)
#summary(mod2ii)

#building off of 2017 data
mod2a<-lm(Density~Treatment*Day*deployment.day, data=df)
hist(resid(mod2a))
qqnorm(resid(mod2a))
qqline(resid(mod2a))
anova(mod2a)
summary(mod2a)#not the best model based on AIC value (compared to mod1 and mod2)
plot(mod2a)#maybe unequal variance in the medium treatment, but not too bad, maybe some outliers too

mod2a<-lm(Density~Treatment*Week*deployment.day, data=df)
hist(resid(mod2a))
qqnorm(resid(mod2a))
qqline(resid(mod2a))
anova(mod2a)
summary(mod2a)

mod2b<-lm(Density~Treatment*Day+(1|Trial), data=df)
hist(resid(mod2b))
qqnorm(resid(mod2b))
qqline(resid(mod2b))
anova(mod2b)
summary(mod2b)

mod2b<-lm(Density~Treatment*Week+(1|Trial), data=df)
hist(resid(mod2b))
qqnorm(resid(mod2b))
qqline(resid(mod2b))
anova(mod2b)
summary(mod2b)

mod2c<-lm(Density~Treatment*Day*Trial, data=df)
hist(resid(mod2c))
qqnorm(resid(mod2c))
qqline(resid(mod2c))
anova(mod2c)
summary(mod2c)#best model based on AIC, but I'm not sure I want to have
##  trial as a fixed effect, didn't really expect to see any differences by trial
##  but it appears like there were differences among trials

mod2c<-lm(Density~Treatment*Week*Trial, data=df)
hist(resid(mod2c))
qqnorm(resid(mod2c))
qqline(resid(mod2c))
anova(mod2c)
summary(mod2c)

mod2d<-lm(Density~Treatment*Day*deployment.day+(1|Trial), data=df)
hist(resid(mod2d))
qqnorm(resid(mod2d))
qqline(resid(mod2d))
anova(mod2d)
summary(mod2d)

#2017 data only, running mixed model with nested term
mod2di<-lmer(Density~Treatment*Day+(1|Trial)+(1|Trial:deployment.day), data=df)
hist(resid(mod2di))
qqnorm(resid(mod2di))
qqline(resid(mod2di))
anova(mod2di)
summary(mod2di)
ranef(mod2di)
#no interaction between treatment and day, so dropping int. from model and running type II anova

mod2dii<-lmer(Density~Treatment+Day+(1|Trial)+(1|Trial:deployment.day), data=df)
hist(resid(mod2dii))
qqnorm(resid(mod2dii))
qqline(resid(mod2dii))
anova(mod2dii,type="II")
summary(mod2dii)
ranef(mod2dii)#going to go with this model for analyses on 19.2.27

mod2d<-lm(Density~Treatment*Week*deployment.day+(1|Trial), data=df)
hist(resid(mod2d))
qqnorm(resid(mod2d))
qqline(resid(mod2d))
anova(mod2d)
summary(mod2d)

mod2e<-lm(Density~Treatment*Day*deployment.day+Trial, data=df)
hist(resid(mod2e))
qqnorm(resid(mod2e))
qqline(resid(mod2e))
anova(mod2e)
summary(mod2e)

mod2e<-lm(Density~Treatment*Week*deployment.day+Trial, data=df)
hist(resid(mod2e))
qqnorm(resid(mod2e))
qqline(resid(mod2e))
anova(mod2e)
summary(mod2e)

mod2f<-lm(Density~Treatment*Day*Trial+deployment.day, data=df)
hist(resid(mod2f))
qqnorm(resid(mod2f))
qqline(resid(mod2f))
anova(mod2f)
summary(mod2f)

mod2f<-lm(Density~Treatment*Week*Trial+deployment.day, data=df)
hist(resid(mod2f))
qqnorm(resid(mod2f))
qqline(resid(mod2f))
anova(mod2f)
summary(mod2f)

mod2fi<-lmer(Density~Treatment*Week+(1|Trial:deployment.day), data=df)
hist(resid(mod2fi))
qqnorm(resid(mod2fi))
qqline(resid(mod2fi))
anova(mod2fi)
summary(mod2fi)
ranef(mod2fi)

#nesting deployment day within trial
#--might consider making dep. days into unique modifiers (or just run them as a nested factor within trial)

#Rationale for seting up mixed model with trial and dep. day nested within trial as random effects
#I don't think I should include trial as a fixed factor (no reason to explian why trial would have an effect on visual counts)
#--by this logic, I should analyze the effect of treatment and week on density
#--instead, should run deployment day nested within trial, and include trial as a random effect

#--did that, and found that there was no interaction when run as a Type I ANOVA
mod2fii<-lmer(Density~Treatment*Week*Trial+(1|Trial:deployment.day), data=df)
hist(resid(mod2fii))
qqnorm(resid(mod2fii))
qqline(resid(mod2fii))
anova(mod2fii,type="I")#best model based on AIC values
summary(mod2fii)
ranef(mod2fii)#not such a difference in trial 5, but much more fish were seen during trial
#-- 4 if they were deployed on the second day (more than twice as many fish)

#testing out different types of ANOVA, turns out that order of factors in model matters
mod2fiii<-lmer(Density~Treatment*Week+(1|Trial)+(1|Trial:deployment.day), data=df)
hist(resid(mod2fiii))
qqnorm(resid(mod2fiii))
qqline(resid(mod2fiii))
anova(mod2fiii,type="I")#best model based on AIC values
summary(mod2fiii)
ranef(mod2fiii)

mod2fiv<-lmer(Density~Week*Treatment+(1|Trial)+(1|Trial:deployment.day), data=df)
hist(resid(mod2fiv))
qqnorm(resid(mod2fiv))
qqline(resid(mod2fiv))
anova(mod2fiv,type="I")#best model based on AIC values
summary(mod2fiv)
ranef(mod2fiv)

#dropping interaction term of treatment and week (not significant), and running as a 
#--type II ANOVA
#this is the model that I went with for results for trials 5 and 6
mod2fv<-lmer(Density~Week+Treatment+(1|Trial)+(1|Trial:deployment.day), data=df)
hist(resid(mod2fv))
qqnorm(resid(mod2fv))
qqline(resid(mod2fv))
anova(mod2fv,type="II")
summary(mod2fv)
ranef(mod2fv)

#trial 6 mixed model, all deployed on the same day! go back up top to simpler models
mod2fvi<-lmer(Density~Week*Treatment, data=df)
hist(resid(mod2fvi))
qqnorm(resid(mod2fivi))
qqline(resid(mod2fivi))
anova(mod2fvi,type="I")#best model based on AIC values
summary(mod2fvi)
ranef(mod2fivi)

mod2fv<-lmer(Density~Treatment*Week+(1|Trial)+(1|Trial:deployment.day), data=df)
hist(resid(mod2fv))
qqnorm(resid(mod2fv))
qqline(resid(mod2fv))
anova(mod2fv,type="I")
summary(mod2fv)
ranef(mod2fv)

#no interaction, so dropped the interaction from the model and ran it as a type II anova
mod2fvi<-lmer(Density~Treatment+Week+(1|Trial)+(1|Trial:deployment.day), data=df)
hist(resid(mod2fvi))
qqnorm(resid(mod2fvi))
qqline(resid(mod2fvi))
anova(mod2fvi,type="II")
summary(mod2fvi)
ranef(mod2fvi)

#reordering factors to get ranef
mod2fvi<-lmer(Density~Week+Treatment+(1|Trial)+(1|Trial:deployment.day), data=df)
hist(resid(mod2fvi))
qqnorm(resid(mod2fvi))
qqline(resid(mod2fvi))
anova(mod2fvi,type="II")
summary(mod2fvi)
ranef(mod2fvi)

mod2g<-lmer(Density~Treatment*Day*Trial+(1|deployment.day), data=df)
hist(resid(mod2g))
qqnorm(resid(mod2g))
qqline(resid(mod2g))
anova(mod2g)
summary(mod2g)

mod2g<-lmer(Density~Treatment*Week*Trial+(1|deployment.day), data=df)
hist(resid(mod2g))
qqnorm(resid(mod2g))
qqline(resid(mod2g))
anova(mod2g)
summary(mod2g)

mod2h<-lm(Density~Treatment*Day*Trial*deployment.day, data=df)
hist(resid(mod2h))
qqnorm(resid(mod2h))
qqline(resid(mod2h))
anova(mod2h)
summary(mod2h)#best model based on AIC value, but I'm really interested in 
#--number of fish seen based on treatment over time

mod2h<-lm(Density~Treatment*Week*Trial*deployment.day, data=df)
hist(resid(mod2h))
qqnorm(resid(mod2h))
qqline(resid(mod2h))
anova(mod2h)
summary(mod2h)

mod2i<-lm(Density~Treatment*Day+Trial+deployment.day, data=df)
hist(resid(mod2i))
qqnorm(resid(mod2i))
qqline(resid(mod2i))
anova(mod2i)
summary(mod2i)

mod2i<-lm(Density~Treatment*Week+Trial+deployment.day, data=df)
hist(resid(mod2i))
qqnorm(resid(mod2i))
qqline(resid(mod2i))
anova(mod2i)
summary(mod2i)

mod2j<-lmer(Density~Treatment*Day+(1|Trial)+(1|deployment.day), data=df)
hist(resid(mod2j))
qqnorm(resid(mod2j))
qqline(resid(mod2j))
anova(mod2j)
summary(mod2j)

mod2j<-lmer(Density~Treatment*Week+(1|Trial)+(1|deployment.day), data=df)
hist(resid(mod2j))
qqnorm(resid(mod2j))
qqline(resid(mod2j))
anova(mod2j)
summary(mod2j)
ranef(mod2j)

#checking AIC values per model for 2017
AIC(mod1,mod2,mod2a,mod2c,mod2d,mod2e,mod2f,mod2g,mod2h,mod2i,mod2j,mod2fi,mod2fii,mod2fiii,
    mod2fiv,mod2fv)

anova(mod2,mod2a,mod1)
AIC(mod1)
AIC(mod2)
AIC(mod2a)
AIC(mod2b)
AIC(mod2c)
AIC(mod2d)

mod2e<-lm(Density~Day*Trial, data=df)
hist(resid(mod2e))
qqnorm(resid(mod2e))
qqline(resid(mod2e))
anova(mod2e)
summary(mod2e)

mod3<-lm(Density~Treatment*Week, data=df)
hist(resid(mod3))
qqnorm(resid(mod3))
qqline(resid(mod3))
anova(mod3)
summary(mod3)

gob.survey<-lmer(Density~Treatment*Week+(1|deployment.day),data=gob.den.5)
hist(resid(gob.survey))
qqnorm(resid(gob.survey))
qqline(resid(gob.survey))
anova(gob.survey)
summary(gob.survey)

#2017 plot
df$Treatment<-ordered(df$Treatment, c("Low","Medium","High"))
df$Treatment<-ordered(df$Treatment, c("Low","Medium","High","Control"))
df$Trial<-ordered(df$Trial, c("1","2","3"))
lineplot.CI(Day,Density,group=Treatment,legend = TRUE,main="2017 densities all trials by treatment", xlab="Day", ylab="Number of fish per reef (mean +/- se)", data=df)
lineplot.CI(Day,Density,group=Trial,legend = TRUE,main="2017 densities by day and trial", xlab="Day", ylab="Number of fish per reef (mean +/- se)", data=df)
bargraph.CI(x.factor = Week, response = Egg.count, group= Treatment, legend=TRUE, main="2017 + 2018 total counts", xlab="Week", ylab="Eggs per reef (mean +/- se)", x.leg=13, yleg=4500, data = df)
bargraph.CI(x.factor = Trial, response = Density, group= Treatment, legend=TRUE, main="2017 densities by trial and treatment", xlab="Trial", ylab="Number of fish seen per reef (mean +/- se)", data = df)
bargraph.CI(x.factor = Trial, response = Density, group= deployment.day, legend=FALSE, main="2017 densities by deployment day and trial", xlab="Trial", ylab="Number of fish seen per reef (mean +/- se)", data = df)
bargraph.CI(x.factor = Day, response = Density, group= Treatment, legend=TRUE, main="2017, trials 1-3 density by treatment and day", xlab="Day", ylab="Number of fish per reef (mean +/- se)", data = df)
#bargraph.CI(x.factor = Trial, response = Density, legend=TRUE, main="2017, trials 1-3 density by treatment and day", xlab="Day", ylab="Number of fish per reef (mean +/- se)", data = df)

#2018 trial 4 and 5 plots
bargraph.CI(x.factor = Trial, response = Density, group= deployment.day, legend=TRUE, main="2018,trials 4 and 5, densities by deployment day", xlab="Trial", ylab="Number of fish seen per reef (mean +/- se)", data = df)
bargraph.CI(x.factor = Week, response = Density, group= Treatment, legend=TRUE, main="2018, trials 4 and 5 density by treatment and week", xlab="Week", ylab="Number of fish per reef (mean +/- se)", data = df)
lineplot.CI(Week,Density,group=Treatment,legend = TRUE,main="2018, trials 4 and 5, density by treatment and week", xlab="Week", ylab="Number of fish per reef (mean +/- se)", data=df)

#2018 trial 6 plots
bargraph.CI(x.factor = Week, response = Density, group=Treatment, legend=TRUE, main="2018,trial 6, densities by treatment and week", xlab="Week", ylab="Number of fish seen per reef (mean +/- se)",yleg=8, data = df)
lineplot.CI(Week,Density,group=Treatment,legend = TRUE,main="2018,trial 6, density by treatment and week", xlab="Week", ylab="Number of fish per reef (mean +/- se)", data=df)
bargraph.CI(x.factor = Day, response = Density, group= Treatment, legend=TRUE, main="2018, trials 4 and 5 density by treatment and day", xlab="Day", ylab="Number of fish per reef (mean +/- se)", data = df)

#2018 trial 6, just week 1
bargraph.CI(x.factor = Day, response = Density,group = Treatment, legend=TRUE, main="2018,trial 6, densities by treatment and day, week 1 data only", xlab="Day", ylab="Number of fish seen per reef (mean +/- se)",yleg=11, data = df)

#looking at caged vs. uncaged within week 1
bargraph.CI(x.factor = Day, response = Density,group = Caging, legend=TRUE, main="2018,trial 6, densities by treatment and day, week 1 data only", xlab="Day", ylab="Number of fish seen per reef (mean +/- se)",yleg=11, data = df)


#using den.max ("adjusted density counts") instead of raw den counts####
#df assignments for analyses
df<-den.t1.2.3
df<-den.t4.5
df<-den.t6
#just looking at week 1 for trial 6
df<-den.t6[(den.t6$Week==1),]

#model building
mod1<-lm(den.max~Treatment, data=df)
hist(resid(mod1))
qqnorm(resid(mod1))
qqline(resid(mod1))
anova(mod1)
summary(mod1)

#this is what I will need to do for trial 1-3, 2017, week isn't great becuase it was only 1
mod2<-lm(den.max~Treatment*Day, data=df)
hist(resid(mod2))
qqnorm(resid(mod2))
qqline(resid(mod2))
anova(mod2)
summary(mod2)

#dropping interaction term from model and rerunning
mod2a<-lm(den.max~Treatment+Day, data=df)
hist(resid(mod2))
qqnorm(resid(mod2))
qqline(resid(mod2))
anova(mod2)
summary(mod2)

mod2<-lm(den.max~Treatment*Week, data=df)
hist(resid(mod2))
qqnorm(resid(mod2))
qqline(resid(mod2))
anova(mod2)
summary(mod2)#for trial 6, no sig. interaction term
plot(mod2)

#for 2018, trial 6, removing interaction term and reanalyzing
mod2i<-lm(den.max~Treatment+Week, data=df)
hist(resid(mod2i))
qqnorm(resid(mod2i))
qqline(resid(mod2i))
anova(mod2i)
summary(mod2i)
plot(mod2i)

#trial 6 power analysis for linear model
pwr.f2.test(u =3, v = 72, f2 = 0.05, sig.level = 0.05, power = )


#consider looking at the outlying data to see if that makes a difference

#mod2ii<-lm(den.max~Treatment*Week+(1|Week:Day), data=df)#didn't work b/c not every day exists
#--within each week
#hist(resid(mod2ii))
#qqnorm(resid(mod2ii))
#qqline(resid(mod2ii))
#anova(mod2ii)
#summary(mod2ii)

#building off of 2017 data
mod2a<-lm(den.max~Treatment*Day*deployment.day, data=df)
hist(resid(mod2a))
qqnorm(resid(mod2a))
qqline(resid(mod2a))
anova(mod2a)
summary(mod2a)#not the best model based on AIC value (compared to mod1 and mod2)
plot(mod2a)#maybe unequal variance in the medium treatment, but not too bad, maybe some outliers too

mod2a<-lm(den.max~Treatment*Week*deployment.day, data=df)
hist(resid(mod2a))
qqnorm(resid(mod2a))
qqline(resid(mod2a))
anova(mod2a)
summary(mod2a)

mod2b<-lm(den.max~Treatment*Day+(1|Trial), data=df)
hist(resid(mod2b))
qqnorm(resid(mod2b))
qqline(resid(mod2b))
anova(mod2b)
summary(mod2b)

mod2b<-lm(den.max~Treatment*Week+(1|Trial), data=df)
hist(resid(mod2b))
qqnorm(resid(mod2b))
qqline(resid(mod2b))
anova(mod2b)
summary(mod2b)

mod2c<-lm(den.max~Treatment*Day*Trial, data=df)
hist(resid(mod2c))
qqnorm(resid(mod2c))
qqline(resid(mod2c))
anova(mod2c)
summary(mod2c)#best model based on AIC, but I'm not sure I want to have
##  trial as a fixed effect, didn't really expect to see any differences by trial
##  but it appears like there were differences among trials

mod2c<-lm(den.max~Treatment*Week*Trial, data=df)
hist(resid(mod2c))
qqnorm(resid(mod2c))
qqline(resid(mod2c))
anova(mod2c)
summary(mod2c)

mod2d<-lm(den.max~Treatment*Day*deployment.day+(1|Trial), data=df)
hist(resid(mod2d))
qqnorm(resid(mod2d))
qqline(resid(mod2d))
anova(mod2d)
summary(mod2d)

#2017 data only, running mixed model with nested term
mod2di<-lmer(den.max~Treatment*Day+(1|Trial)+(1|Trial:deployment.day), data=df)
hist(resid(mod2di))
qqnorm(resid(mod2di))
qqline(resid(mod2di))
anova(mod2di)
summary(mod2di)
ranef(mod2di)
#no interaction between treatment and day, so dropping int. from model and running type II anova

mod2dii<-lmer(den.max~Treatment+Day+(1|Trial)+(1|Trial:deployment.day), data=df)
hist(resid(mod2dii))
qqnorm(resid(mod2dii))
qqline(resid(mod2dii))
anova(mod2dii,type="II")
summary(mod2dii)
ranef(mod2dii)#going to go with this model for analyses on 19.2.27

mod2d<-lm(den.max~Treatment*Week*deployment.day+(1|Trial), data=df)
hist(resid(mod2d))
qqnorm(resid(mod2d))
qqline(resid(mod2d))
anova(mod2d)
summary(mod2d)

mod2e<-lm(den.max~Treatment*Day*deployment.day+Trial, data=df)
hist(resid(mod2e))
qqnorm(resid(mod2e))
qqline(resid(mod2e))
anova(mod2e)
summary(mod2e)

mod2e<-lm(den.max~Treatment*Week*deployment.day+Trial, data=df)
hist(resid(mod2e))
qqnorm(resid(mod2e))
qqline(resid(mod2e))
anova(mod2e)
summary(mod2e)

mod2f<-lm(den.max~Treatment*Day*Trial+deployment.day, data=df)
hist(resid(mod2f))
qqnorm(resid(mod2f))
qqline(resid(mod2f))
anova(mod2f)
summary(mod2f)

mod2f<-lm(den.max~Treatment*Week*Trial+deployment.day, data=df)
hist(resid(mod2f))
qqnorm(resid(mod2f))
qqline(resid(mod2f))
anova(mod2f)
summary(mod2f)

mod2fi<-lmer(den.max~Treatment*Week+(1|Trial:deployment.day), data=df)
hist(resid(mod2fi))
qqnorm(resid(mod2fi))
qqline(resid(mod2fi))
anova(mod2fi)
summary(mod2fi)
ranef(mod2fi)

#nesting deployment day within trial
#--might consider making dep. days into unique modifiers (or just run them as a nested factor within trial)

#Rationale for seting up mixed model with trial and dep. day nested within trial as random effects
#I don't think I should include trial as a fixed factor (no reason to explian why trial would have an effect on visual counts)
#--by this logic, I should analyze the effect of treatment and week on den.max
#--instead, should run deployment day nested within trial, and include trial as a random effect

#--did that, and found that there was no interaction when run as a Type I ANOVA
mod2fii<-lmer(den.max~Treatment*Week*Trial+(1|Trial:deployment.day), data=df)
hist(resid(mod2fii))
qqnorm(resid(mod2fii))
qqline(resid(mod2fii))
anova(mod2fii,type="I")#best model based on AIC values
summary(mod2fii)
ranef(mod2fii)#not such a difference in trial 5, but much more fish were seen during trial
#-- 4 if they were deployed on the second day (more than twice as many fish)

#testing out different types of ANOVA, turns out that order of factors in model matters
mod2fiii<-lmer(den.max~Treatment*Week+(1|Trial)+(1|Trial:deployment.day), data=df)
hist(resid(mod2fiii))
qqnorm(resid(mod2fiii))
qqline(resid(mod2fiii))
anova(mod2fiii,type="I")#best model based on AIC values
summary(mod2fiii)
ranef(mod2fiii)

mod2fiv<-lmer(den.max~Week*Treatment+(1|Trial)+(1|Trial:deployment.day), data=df)
hist(resid(mod2fiv))
qqnorm(resid(mod2fiv))
qqline(resid(mod2fiv))
anova(mod2fiv,type="I")#best model based on AIC values
summary(mod2fiv)
ranef(mod2fiv)

#dropping interaction term of treatment and week (not significant), and running as a 
#--type II ANOVA
#this is the model that I went with for results for trials 5 and 6
mod2fv<-lmer(den.max~Week+Treatment+(1|Trial)+(1|Trial:deployment.day), data=df)
hist(resid(mod2fv))
qqnorm(resid(mod2fv))
qqline(resid(mod2fv))
anova(mod2fv,type="II")
summary(mod2fv)
ranef(mod2fv)

#trial 6 mixed model, all deployed on the same day! go back up top to simpler models
mod2fvi<-lmer(den.max~Week*Treatment, data=df)
hist(resid(mod2fvi))
qqnorm(resid(mod2fivi))
qqline(resid(mod2fivi))
anova(mod2fvi,type="I")#best model based on AIC values
summary(mod2fvi)
ranef(mod2fivi)

mod2fv<-lmer(den.max~Treatment*Week+(1|Trial)+(1|Trial:deployment.day), data=df)
hist(resid(mod2fv))
qqnorm(resid(mod2fv))
qqline(resid(mod2fv))
anova(mod2fv,type="I")
summary(mod2fv)
ranef(mod2fv)

#no interaction, so dropped the interaction from the model and ran it as a type II anova
mod2fvi<-lmer(den.max~Treatment+Week+(1|Trial)+(1|Trial:deployment.day), data=df)
hist(resid(mod2fvi))
qqnorm(resid(mod2fvi))
qqline(resid(mod2fvi))
anova(mod2fvi,type="II")
summary(mod2fvi)
ranef(mod2fvi)

#reordering factors to get ranef
mod2fvi<-lmer(den.max~Week+Treatment+(1|Trial)+(1|Trial:deployment.day), data=df)
hist(resid(mod2fvi))
qqnorm(resid(mod2fvi))
qqline(resid(mod2fvi))
anova(mod2fvi,type="II")
summary(mod2fvi)
ranef(mod2fvi)

mod2g<-lmer(den.max~Treatment*Day*Trial+(1|deployment.day), data=df)
hist(resid(mod2g))
qqnorm(resid(mod2g))
qqline(resid(mod2g))
anova(mod2g)
summary(mod2g)

mod2g<-lmer(den.max~Treatment*Week*Trial+(1|deployment.day), data=df)
hist(resid(mod2g))
qqnorm(resid(mod2g))
qqline(resid(mod2g))
anova(mod2g)
summary(mod2g)

mod2h<-lm(den.max~Treatment*Day*Trial*deployment.day, data=df)
hist(resid(mod2h))
qqnorm(resid(mod2h))
qqline(resid(mod2h))
anova(mod2h)
summary(mod2h)#best model based on AIC value, but I'm really interested in 
#--number of fish seen based on treatment over time

mod2h<-lm(den.max~Treatment*Week*Trial*deployment.day, data=df)
hist(resid(mod2h))
qqnorm(resid(mod2h))
qqline(resid(mod2h))
anova(mod2h)
summary(mod2h)

mod2i<-lm(den.max~Treatment*Day+Trial+deployment.day, data=df)
hist(resid(mod2i))
qqnorm(resid(mod2i))
qqline(resid(mod2i))
anova(mod2i)
summary(mod2i)

mod2i<-lm(den.max~Treatment*Week+Trial+deployment.day, data=df)
hist(resid(mod2i))
qqnorm(resid(mod2i))
qqline(resid(mod2i))
anova(mod2i)
summary(mod2i)

mod2j<-lmer(den.max~Treatment*Day+(1|Trial)+(1|deployment.day), data=df)
hist(resid(mod2j))
qqnorm(resid(mod2j))
qqline(resid(mod2j))
anova(mod2j)
summary(mod2j)

mod2j<-lmer(den.max~Treatment*Week+(1|Trial)+(1|deployment.day), data=df)
hist(resid(mod2j))
qqnorm(resid(mod2j))
qqline(resid(mod2j))
anova(mod2j)
summary(mod2j)
ranef(mod2j)

#checking AIC values per model for 2017
AIC(mod1,mod2,mod2a,mod2c,mod2d,mod2e,mod2f,mod2g,mod2h,mod2i,mod2j,mod2fi,mod2fii,mod2fiii,
    mod2fiv,mod2fv)

anova(mod2,mod2a,mod1)
AIC(mod1)
AIC(mod2)
AIC(mod2a)
AIC(mod2b)
AIC(mod2c)
AIC(mod2d)

mod2e<-lm(den.max~Day*Trial, data=df)
hist(resid(mod2e))
qqnorm(resid(mod2e))
qqline(resid(mod2e))
anova(mod2e)
summary(mod2e)

mod3<-lm(den.max~Treatment*Week, data=df)
hist(resid(mod3))
qqnorm(resid(mod3))
qqline(resid(mod3))
anova(mod3)
summary(mod3)

gob.survey<-lmer(den.max~Treatment*Week+(1|deployment.day),data=gob.den.5)
hist(resid(gob.survey))
qqnorm(resid(gob.survey))
qqline(resid(gob.survey))
anova(gob.survey)
summary(gob.survey)

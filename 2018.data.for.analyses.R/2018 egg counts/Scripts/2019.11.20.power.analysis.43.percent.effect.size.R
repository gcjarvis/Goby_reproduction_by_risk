# Description: calculation for power analysis with 43% reduction in output
# Author: George C Jarvis
# Date: Wed Nov 20 07:46:50 2019
# Notes: this effect size of 43% reduction under high risk environment was taken from (Mukherjee et al. 2014)
# --------------

library(lme4)
library(lmerTest)
library(nlme)
library(car)
library(pwr)
library(tidyverse)
library(simr) #testing effect sizes of mixed models
# Update: 2020.1.30:
#this package works with models built from lme4 models (e.g., lmer, not lme)
#will see if this works

options(contrasts = c("contr.sum","contr.poly")) #this is important, run before ANOVA, will set SS to type III

#loading data####
repro<-read.csv("Data/new.data.2019.9.30.csv")
repro<-na.omit(repro) # no NA's to omit

#data manipulation####
#adding column for average density, rounded to nearest whole number of fish
repro$avg.inhab<-(ceiling((repro$Recollection+20)/2))

#adding a column for year (as a proxy for tagging procedure), where trials 1-3 = 2017, and 4-6 = 2018
repro$Year <- ifelse(repro$Trial <=3, 2017, 2018)
#want it as a factor? Going to make a variable with it as a factor, run the model, and see if I get different results
repro$Year.fact<- as.factor(repro$Year)

#adding egg/week variable to the dataset
repro<-repro %>%
  mutate(egg.week = ifelse(Trial<4, Egg.count/1,
                           ifelse(Trial == 4| Trial == 5, (ceiling(Egg.count/4)),
                                  ifelse(Trial == 6, (ceiling(Egg.count/2)), NA))))

# original code for power analysis, think this is wrong####
pwr.f2.test(u = 2, v = 99.09, f2 = 0.43, sig.level = 0.05, power = NULL)


#2020.1.8.retesting

pwr.f2.test(u = 2, v = 99, f2 = NULL, sig.level = 0.05, power = .80)

#ran model as linear model, with trial included as a fixed effect

pwr.f2.test(u = 2, v = 86, f2 = NULL, sig.level = 0.05, power = .80)

#still very high power to detect differences based on our sample size (~100% chance)

#parameters:

#u degrees of freedom for numerator
#v degrees of freedom for denominator
#f2 effect size
#sig.level Significance level (Type I error probability)
#power Power of test (1 minus Type II error probability)


pwr.f2.test(u = NULL, v = NULL, f2 = NULL, sig.level = 0.05, power = NULL)

#2020.3.7 update after emailing CSUN R wizards#####

#tricky part is figuring out what the correct error term is so that I have
# - the right denom df

#going to mess around with it

#based off of 2, 99 df (num, denom)

pwr.f2.test(u = 2, v = 6, f2 = 0.20, sig.level = 0.05, power = NULL)

#2020.1.8 update, using fully-reduced model from "2019.12.12.testing.out.l.schuster.code" script####

#trying out new package to test for power of mixed models

#SIMR: Power Analysis for Generalised Linear Mixed Models by Simulation
#NOTE: I don't think any of this is right

# I'm not sure how to do a power analysis in R with a mixed model
# I think doing it as a linear model with no random effects is wrong

#1) now that I have the model, I can test the effects, might not work with nlme though...
#idea is that you can change the effect size of the fixed effect

#original model from other script

mod2.2.luk<-lme(egg.week~(Treatment*Year.fact)+ Treatment+
                  avg.inhab+Year.fact,random=~1|Trial,repro,method="REML")

#ordered treatments so I can see which is which in output
mod2.2.luk<-lme(egg.week~(treatment.ordered*Year.fact)+ treatment.ordered+
                  avg.inhab+Year.fact,random=~1|Trial,repro,method="REML")

#trying to simplify model, just using treatment by trial nested within year
#note: this was my final model, so this is the one that I will use in the simR package
mod2.2.simp<-lme(egg.week~Treatment*Year,random=~1|Trial,repro,method="REML")

levels(repro$treatment.ordered)
levels(repro$Treatment)

fixef(mod2.2.simp)


fixef(mod2.2.simp)["Treatment"]
ranef(mod2.2.luk)
## x
## -0.1148147
fixef(model1)["x"] <- -0.05

fixef(model1)["x"]
## x
## -0.1148147
fixef(model1)["x"] <- -0.05


#power analysis with simR package. tried out on 2020.1.30#####

#note: this was my final model for egg counts, 
# - so this is the one that I will use in the simR package

#example dataset "simdata"

#The data set is representative of environmental monitoring data, with a response variable z (e.g. bird abun- dance) 
# - measured at 10 levels of the continuous fixed effect variable x (e.g. study year) 
# -for three groups g (e.g. study site). There is also a continuous response variable 
# -y, which is not used in this tutorial.

View(simdata)

#poisson model for example with glmer
model1 <- glmer(z ~ x + (1|g), family="poisson", data=simdata)
hist(resid(model1))
qqnorm(resid(model1))
qqline(resid(model1))
summary(model1) #looking at "Estimate", it appears that the effect size of x is -0.11

#going to try with my model (nlme package)
mod2.2.simp<-lme(egg.week~Treatment*Year+avg.inhab,random=~1|Trial,repro,method="REML")
anova(mod2.2.simp, type='marginal')
summary(mod2.2.simp)
#"Value column from this output is the same thing as "Estimate" column when run as a model with lme4 package

#going to reun model with lme4 package
mod1.3<-lmer(egg.week~Treatment*Year+(1|Year:Trial), data=repro)
hist(resid(mod1.3))
qqnorm(resid(mod1.3))
qqline(resid(mod1.3))
anova(mod1.3)
Anova(mod1.3)
summary(mod1.3) #"Estimate" column is in this output. I can calculate effect sizes this way
# or just look at the effect sizes that I calculated from the LS means package

#need to reformat my data so it looks like the data in the template
# - for this 

#exporting wrangled data
write.csv(repro,"Data\\wegg.power.analysis.2020.1.30.csv", row.names = FALSE)
#changed categorical factor of treatment to numeric (L=1, M=2, H=3)
#renamed column "Treatment.power"

#importing new dataset
repower<-read.csv(file = "Data/wegg.power.analysis.2020.1.30.csv")
repower<-na.omit(repower)

#rerunning models with new dataset 

# nlme package
mod1.3<-lmer(egg.week~Treatment.power*Year+(1|Year:Trial), data=repower)
hist(resid(mod1.3))
qqnorm(resid(mod1.3))
qqline(resid(mod1.3))
anova(mod1.3)
Anova(mod1.3)
summary(mod1.3) #"Estimate" column is in this output. I can calculate effect sizes this way
# or just look at the effect sizes that I calculated from the LS means package

fixef(mod1.3)["Treatment.power"]
#Fixed effects:
#                       Estimate Std. Error t value
#(Intercept)           765258.64 2828646.53   0.271
#Treatment.power      -147910.7

#that represents a 0.19 effect size (intercept/trt.pwr)

#close to the ~20% effect size I saw between low and medium?? Not sure if that's correct

#checking power of model before I change any of the effect sizes
#template

powerSim(model1)#worked, but took a lot of iterations

powerSim(mod1.3)#didn't work, trouble finding fixed effect in model
#will remove year?
mod1.3a<-lmer(egg.week~Treatment.power+(1|Year:Trial), data=repower)
abline(fixef(mod1.3a))
#worked, but will likely take some time, and also, this is not the model
# - that I ultimately want to run the power analysis on

powerSim(mod1.3a)
#Power for predictor 'Treatment.power', 
#(95% confidence interval):=======|
#  15.70% (13.50, 18.11)

#15% power to detect power based on the effect size that we saw (19% decrease)

#now going to change effect size for a 50% decrease in output for treatment
#-389238.82 = 50% effect ((-147910.75/0.19)=(x/0.5)); solve for x = -389238.82

fixef(mod1.3a)["Treatment.power"] <- -389238.82

#rerunning power analysis to see what my power would be with that effect size
powerSim(mod1.3a)
#Power for predictor 'Treatment.power', 
#(95% confidence interval):=======|
# 100.0% (99.63, 100.0)


powerSi

fixef(mod1.3)

fixef(mod1.3)

#setting the effect size for treatment
#template (from green and macleod, 2016 paper; in mendeley library)
fixef(model1)["x"]
x
## -0.1148147 fixef(model1)["x"] <- -0.05

fixef(model1)["x"] <- -0.05

#thid package, to test for power when there is more than one fixed effect
library(powerMediation)
powerInteract(nTotal, a, b, effsize, alpha = 0.05, nTests = 1)
powerInteract(110, 2, 3, 0.8, alpha = 0.05, nTests = 1)
#86% chance to detect a difference with an effect size of 0.8

powerInteract(110, 2, 3, 0.20, alpha = 0.05, nTests = 1)
#with what we actually saw, we had an 11% chance to see an effect

# retrying with code from simR package####
# code found at: https://humburg.github.io/Power-Analysis/simr_power_analysis.html

#rationale: I was having trouble with simR before because I was constructing the model 
# - incorrectly. I need to rsonctruct the model with more explicit parameters

#reduced model that I used for my analyses
#NOTE: this was constructed with nlme package, need to change to lmer
mod2.2.luk<-lme(egg.week~(Treatment*Year.fact)+ Treatment+
                  avg.inhab+Year.fact,random=~1|Trial,repro,method="REML")
anova(mod2.2.luk)
summary(mod2.2.luk)


#rationale: even though the model construction is different, the results should be fairly similar
#- Also, if it allows me to be able to run the simR package, then this is the best I can do

#making lmer model within the simR package
#template
model <- makeLmer(y ~ treat*time + (1|class/id), fixef=fixed, VarCorr=rand, sigma=res, data=covars)
model

#my model
mod.sim<- makeLmer(egg.week~(Treatment*Year.fact)+Treatment+avg.inhab+
                     Year.fact+(1|Year.fact/Trial),data=repro)
#doesn't like it like that, but I can run the model and add the parameters that
# - the code is lacking (efect sizes and variance)

#running the lmer model
mod2.2<-lmer(egg.week~Treatment*Year.fact+avg.inhab+(1|Year.fact:Trial),
             data=repro)
hist(resid(mod2.2))
qqnorm(resid(mod2.2))
qqline(resid(mod2.2))
anova(mod2.2)
Anova(mod2.2)
ranef(mod2.2)
summary(mod2.2)

#want to get the variance of the slopes (mainly Year.fact)
#temlplate
(summary(m)$coefficients[2,2])**2

(summary(mod2.2)$coefficients[4,2])**2# [4,2] gets me to std.error for year
#trying another way
#template
vcov(m)[2,2]

#mine
vcov(mod2.2)[4,2] #not the same values...but I feel more confident in this value
# - just based off of the egg counts (other value is ~100k, no way that could be it)

#trying to calculate the variance for year
#going to calculate it by hand, based on se and n=2 (2 years)

321.8912*(sqrt(2)) # = 455.2229, doesn't seem to be right


mod.view<-(summary(mod2.2)$coefficients) #not what I want, but great way to export a table!
View(mod.view)

#setting the effect sizes based on summary of lmer model

#template (done ahead of time typically)

## Intercept and slopes for intervention, time1, time2, intervention:time1, intervention:time2
fixed <- c(5, 0, 0.1, 0.2, 1, 0.9)

## Random intercepts for participants clustered by class
rand <- list(0.5, 0.1)

## residual variance
res <- 2

#want to see which treatments correspond to the summary table

#my version (found difference between each estimate, relative to the intercept)

fixed<-c(-6274.33,-6414.55,-5888.02,-6058.39,-5537.64,-6089.48,
         -6127.15)

#random intercepts for trials, clustered by year (have thesse values)
rand<-list(546.68241,-94.38653,-452.29588,290.42689,325.57428,-616.00117)
# Error in fn(theta) : theta size mismatch (6!=2)

# I don't think it likes that I have all these random intercepts, only wants 2 values
# - for the vector

#random intercept variance for trial nested within year
#intercept for nested factor is first, then year.fact (used same value as the one for fixedef)
#need to get the variance for Year.fact
var(Year.fact)

rand1<-list(413133,103614) #decided to go with the value for the variance of the slope for year
(summary(mod2.2)$coefficients[4,2])**2 #103614

#residual variance
res<- 3769288

#rerunning model with specified parameters given from model
#mod.sim<- makeLmer(egg.week~(Treatment*Year.fact)+Treatment+avg.inhab+
#                     Year.fact+(1|Year.fact/Trial),fixef = fixed, VarCorr = rand1,
#                     sigma = res, data=repro)

#using rand1 for varcorr
mod.sim<- makeLmer(egg.week~(Treatment*Year.fact)+Treatment+avg.inhab+
                    Year.fact+(1|Year.fact/Trial),fixef = fixed, VarCorr = rand1,
                    sigma = res, data=repro)
mod.sim

#power analysis
#interested in the effect of treatment here, so have to compare to alternative
# - models by specifying the fixed factors that we are not interested in 
# - (i.e. Year.Fact, and avg.inhab) with the "fcompare" fuction

#template
sim_treat <- powerSim(model, nsim=100, test = fcompare(y~time))
sim_treat

#mine
sim_treat <- powerSim(mod.sim, nsim=100, test = fcompare(egg.week~Year.fact))
sim_treat

# Description: testing lukas's code for nlme nested design vs. my previous models
# Author: George C Jarvis
# Date: Thu Dec 12 23:12:13 2019
# Notes: Lukas seems to think that this will be better for my nested ANCOVA models
#       for egg counts. We'll see.  I want to compare results for both models
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
library(vegan)
library(multcomp)

options(contrasts = c("contr.sum","contr.poly")) #this is important, run before ANOVA, will set SS to type III

#importing dataset, adding number of gobies on each reef, ordering treatments####
repro<-read.csv("Data/new.data.2019.9.30.csv")
repro<-na.omit(repro) # no NA's to omit

pairs(repro) #pretty cool to see correlations

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

#exporting wrangled data
write.csv(repro,"Data\\egg.counts.2019.12.23.csv", row.names = FALSE)

#subsetting data to do t6 comparison between HR caged and uncaged treatments
#subset trial
repro.t6<-repro[repro$Trial==6,]
#subset treatments, caged and uncaged only (T6.comparison)

#subsetting by multiple character factors within a variable. Awesome code!
repro.t6.comp<-repro.t6[repro.t6$T6.comparison == "High" | repro.t6$T6.comparison == "Uncaged", ]
View(repro.t6.comp)

#modeling####
#original model with lme4, with weekly output per reef as response variable, shown is the reduced model 
# - (from "2019.7.12.adding.year.factor" script)
mod2.2<-lmer(egg.week~Treatment*Year.fact+avg.inhab+(1|Year.fact:Trial),
             data=repro)
hist(resid(mod2.2))
qqnorm(resid(mod2.2))
qqline(resid(mod2.2))
anova(mod2.2)
Anova(mod2.2)
summary(mod2.2) # same results, will check out a model comparison?

#Lukas's model with nlme, let's do the full model and then reduce if needed
mod2.luk<-lme(egg.week~Treatment*Year.fact*avg.inhab,random=~1|Trial,repro,method="REML")
summary(mod2.luk)
anova(mod2.luk, type='marginal')

#running reduced model
mod2.1.luk<-lme(egg.week~(Treatment*Year.fact)+(avg.inhab*Year.fact)+(Treatment*avg.inhab)+ Treatment+
               avg.inhab+Year.fact,random=~1|Trial,repro,method="REML")
summary(mod2.1.luk)
anova(mod2.1.luk, type='marginal')

#running further reduced model
mod2.2.luk<-lme(egg.week~(Treatment*Year.fact)+ Treatment+
                  avg.inhab+Year.fact,random=~1|Trial,repro,method="REML")
summary(mod2.2.luk)
anova(mod2.2.luk, type='marginal')

anova(mod2.luk,mod2.1.luk,mod2.2.luk) #aic value is highest for the least-reduced model, but 
# also might not be meaningful because I'm comparing models with different fixed effects

#NOTE: this doesn't seem to work...
#I'm skeptical that my nested factor isn't being considered correctly in R
#I'm going to rerun the reduced model (mod2.2.luk) in nlme, and this time specify the nested term,
# - as opposed to "random=~1|Trial", where it's assumed that trial is nested within year
#best case scenario, I get the same output, in which case I think I'd feel better about nesting explicitly

mod2.3.luk<-lme(egg.week~(Treatment*Year.fact)+ Treatment+
                  avg.inhab+Year.fact,random=~1|Year.fact/Trial,repro,method="REML")
summary(mod2.3.luk)
anova(mod2.3.luk, type='marginal') #NaNs produced, R doesn't like it

# lme(Thickness ~ 1, random= ~1 | Lot/Wafer, data=Oxide) 

#will try another way

mod2.4.luk<-lme(egg.week~(Treatment*Year.fact)+ Treatment+
                  avg.inhab+Year.fact,random= 
                  list(~1|Year.fact,~1|Trial),repro,method="REML")
summary(mod2.4.luk)
anova(mod2.4.luk, type='marginal') #NaNs as well, not sure if I'll be able to specify

# lme(Thickness ~ 1, random=list(~1|Lot, ~1|Wafer), data=Oxide)

#2020.3.8 after watching video on mixed effects models####

#trial as a factor
repro$Trial.factor<-as.factor(repro$Trial)

mod2.4.luk<-lme(egg.week~(Treatment*Year.fact)+avg.inhab,random= 
                  ~1|Trial/Year.fact,repro,method="REML")
summary(mod2.4.luk)
anova(mod2.4.luk, type='marginal') #NaNs as well, not sure if I'll be able to specify

#as all fixed factors?
mod.fix.all<-aov(egg.week~Treatment*Year*Trial+avg.inhab, data=repro)
summary(mod.fix.all)

mod.fix.all.fact<-aov(egg.week~Treatment*Year*Trial.factor+avg.inhab, data=repro)
summary(mod.fix.all.fact)

#output is the same, which shows that R was at least nesting correctly

#Arithmetic vs. LS means (takes model means into account)####
if(!require(FSA)){install.packages("FSA")}
if(!require(psych)){install.packages("psych")}
if(!require(lsmeans)){install.packages("lsmeans")}
if(!require(car)){install.packages("car")}

library(FSA) #arithmetic means
library(psych)
library(lsmeans)
library(car)

#arithmetic means
library(FSA)

Summarize(egg.week~Treatment,
          data=repro,
          digits=3)
#nice, because it gives sample sizes for each treatment

#Summarize(Height ~ Classroom,
#          data=Data,
#          digits=3)

#LS means, taking model into consideration
mod2.2.luk<-lme(egg.week~(Treatment*Year.fact)+ Treatment+
                  avg.inhab+Year.fact,random=~1|Trial,repro,method="REML")

#making sure that the factor is nesting correctly, found code from here
# - https://stat.ethz.ch/pipermail/r-sig-mixed-models/2013q3/020645.html
# I got the same results when I ran it this way, so I'm assuming that R
# - is nesting Trial within year.fact properly
mod2.2.1.luk<-lme(fixed= egg.week~(Treatment*Year.fact)+ Treatment+
                avg.inhab+Year.fact,random=~1|Trial,
                weights = varIdent((form=~1|Year.fact)) ,repro,method="REML")


summary(mod2.2.1.luk)
anova(mod2.2.1.luk, type='marginal')

library(lsmeans)

lsmeans(mod2.2.luk,
        pairwise~Treatment*avg.inhab,
        adjust="tukey") #gives warning that results may be misleading 
#                         because of interactions

#the lsmeans are not too different from the arithmetic means

#model = lm(Height ~ Classroom + Sex + Classroom:Sex,
#           data = Data)

#library(lsmeans)

#lsmeans(model,
#        pairwise ~ Classroom,
#        adjust="tukey")
#

#testing for differences in output between HR caged and uncaged treatments####

#NOTES: 1) trial 6 only, using "repro.t6.comp" df, and "T6.comparison" for treatments
# 2) it's no longer a nested ANCOVA, it's just an ANCOVA with trt and avg.inhab,
# - so it's just a linear model
# 3) will add the reduced model results to the table for egg counts,
# - and the full model results to the supplementary table for egg counts

#full model
modt6.luk<-lm(egg.week~T6.comparison*avg.inhab,data=repro.t6.comp)
summary(modt6.luk)
anova(modt6.luk)

#reducing because no effect of covariate
modt6.luk<-lm(egg.week~T6.comparison,data=repro.t6.comp)
summary(modt6.luk)
anova(modt6.luk)

#added these stats to tables

#testing trial as a fixed effect####
#not sure it makes sense to run it this way
mod3.t<-lm(egg.week~Treatment*Year.fact*Trial*avg.inhab,repro)
anova(mod3.t)

#trial had an effect, want to see graphically
repro$treatment.ordered<-ordered(repro$Treatment, levels=c("Low","Medium","High"))
#reproduction by treatment
bargraph.CI(x.factor = avg.inhab, response = egg.week, group = Trial,
            legend=TRUE, main="reproduction per avg.inhab", 
            data = repro)

#citations for nlme package and maybe multcomp, 
# - though I'm not sure I want to do multiple comparisons with my models
# - but I doubt I'll get much kickback from reviewers for using that

#let's try multcomp package with my new models (nlme)###
citation(package = "nlme")
citation(package="multcomp")
citation(package="lsmeans")

#running multcomp with lme model for egg counts

#summary(glht(YOUR MODEL, linfct=mcp(YOUR FIXED FACTOR="Tukey")))
#not sure how I feel about that
summary(glht(mod2.2.luk, linfct=mcp(Treatment="Tukey")))

# I don't think I have to do this, because I'm not doing post-hoc tests
# I'm just going to say that there were overall differences when there were,
# - and which factors in the model contributed to those differences

####2020.3.14 log liklihood testing for significance for random nested factor####

#rationale: going to comapare the model results from the full model
# - with that of a reduced model (i.e. one that doesn't include the random factor)
# - of trial nested within year

#got the code from L. Schuster

#mixed model (includes random factor, use nlme package, "full model")
# NOTE: use method= 'ML'

#year as numeric, trial as integer (different result from others)
egg.mm<-lme(egg.week~Treatment*Year+avg.inhab,random=~1|Trial,repro,method='ML')
summary(egg.mm)

#year as factor, trial as integer (same result as egg.mmff)
egg.mmf<-lme(egg.week~Treatment*Year.fact+avg.inhab,random=~1|Trial,repro,method='ML')
summary(egg.mmf)

#year as factor, trial as factor (same results as egg.mmf)
egg.mmff<-lme(egg.week~Treatment*Year.fact+avg.inhab,random=~1|Trial.fact,repro,method='ML')
summary(egg.mmff)

#not sure about the best way include year in the model - as factor or numeric?

#for now, will move ahead with log-liklihood testing

#trying out egg.mm first
hist(resid(egg.mm)) #looks good
logLik(egg.mm) #log-lik: -988.4227 (df=9)

#egg.mmf
hist(resid(egg.mmf)) # same
logLik(egg.mmf) #log-lik: -988.4227 (df=9)

#egg.mmff
hist(resid(egg.mmff)) # same
logLik(egg.mmff) #log-lik: -988.4227 (df=9)

#now going to run the reduced models without the random effect

#year as numeric, trial as integer (different result from others)
regg<-lm(egg.week~Treatment*Year+avg.inhab,repro)
summary(regg)

#year as factor, trial as integer (same result as egg.mmff)
regg.f<-lm(egg.week~Treatment*Year.fact+avg.inhab,repro)
summary(regg.f)

#no need to run another version because no trial vs. trial.fact term

#testing logLik of each of the reduced models

#trying out egg.mm first

hist(resid(regg)) #looks good, slightly worse than mixed model with 
# - year as numeric
logLik(regg) #log-lik: -989.345 (df=8) #looks very similar to other model

#egg.mmf
hist(resid(regg.f)) # same
logLik(regg.f) #log-lik: -989.345 (df=8) #same

#compare logLik of full model to that of reduced model
# NOTE: year and trial as numeric and integer, respectively

2 * (logLik(egg.mm) - logLik(regg)) #chi2 = 1.844638
pchisq(2 * (logLik(egg.mm) - logLik(regg)), df=1, lower.tail = F)
# chi2 = 1.844638, p = 0.1744083

#same analysis, but with year as trial as factors

2 * (logLik(egg.mmff) - logLik(regg.f)) #chi2 = 1.844638
pchisq(2 * (logLik(egg.mmff) - logLik(regg.f)), df=1, lower.tail = F)
# chi2 = 1.844638, p = 0.1744083

#same result = thank god

#testing out code from L. Shuster to see if there are interactions with
# - the random and fixed effects

options(contrasts= c("contr.sum","contr.poly"))

repro$Random <- paste0(repro$Year, repro$Trial) # this is a little trick to specify a nested factor, 
                                             # it will combine year and trial into one, then just
                                             # use it as the random effect below (Random)

m <- lmer(egg.week ~ Treatment*Year*avg.inhab + (1|Random) + (1|Random:Treatment) +
            (1|Random:avg.inhab) + (1|Random:Treatment:avg.inhab), REML=F, repro)
summary(m)
anova(m, type=3) # check the df!! The output should be the same as when using the nlme package
logLik(m)

m.1 <- update(m, .~. -(1|Random:Treatment:avg.inhab))
summary(m.1)
anova(m.1, type=3)
logLik(m.1)

2*(logLik(m) - logLik(m.1)) # Chi2 =  -2.185038e-08
pchisq(2*(logLik(m) - logLik(m.1)), df = 1, lower.tail=F) # P = 1

#trying agin with year.fact

repro$Random <- paste0(repro$Year.fact, repro$Trial) # this is a little trick to specify a nested factor, 
# it will combine year and trial into one, then just
# use it as the random effect below (Random)

m.new.with.int <- lmer(egg.week ~ Treatment*Year.fact*avg.inhab + (1|Random) + (1|Random:Treatment), REML=F, repro)
hist(resid(m.new.with.int))
summary(m.new.with.int)
anova(m, type=3) # check the df!! The output should be the same as when using the nlme package
logLik(m)

m.1 <- update(m, .~. -(1|Random:Treatment:avg.inhab))
summary(m.1)
anova(m.1, type=3)
logLik(m.1)

2*(logLik(m) - logLik(m.1)) # Chi2 =  -2.185038e-08
pchisq(2*(logLik(m) - logLik(m.1)), df = 1, lower.tail=F) # P = 1

#visualizing data from each trial individually####

library(patchwork) #allows me to easily arrange multiple plots
# - for visual comparison

#subsetting into different df's based on trial

repro$treat.ordered<-ordered(repro$Treatment,levels=c("Low","Medium","High"))

# based on variable values
t1 <- repro[ which(repro$Trial==1),]

p1<-ggplot(t1, aes(avg.inhab, egg.week, shape=treat.ordered, linetype=treat.ordered, col=treat.ordered)) +
  geom_smooth(method="lm", se=FALSE, show.legend = TRUE)  +
  geom_point(size=3)+
  theme_classic()+
  labs(x="Gobies per Reef",y=(expression(atop("Reproductive Output", 
                                              paste((eggs~laid~reef^-1))))))+
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
#  scale_x_continuous(breaks=c(10,12,14,16,18,20)) + scale_y_continuous(limits = c(0,40000))+
  labs(color  = "Perceived Risk", linetype = "Perceived Risk", shape = "Perceived Risk")
p1

t2<- repro[ which(repro$Trial==2),]

p2<-ggplot(t2, aes(avg.inhab, egg.week, shape=treat.ordered, linetype=treat.ordered, col=treat.ordered)) +
  geom_smooth(method="lm", se=FALSE, show.legend = TRUE)  +
  geom_point(size=3)+
  theme_classic()+
  labs(x="Gobies per Reef",y=(expression(atop("Reproductive Output", 
                                              paste((eggs~laid~reef^-1))))))+
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
  #  scale_x_continuous(breaks=c(10,12,14,16,18,20)) + scale_y_continuous(limits = c(0,40000))+
  labs(color  = "Perceived Risk", linetype = "Perceived Risk", shape = "Perceived Risk")
p2

t3<- repro[ which(repro$Trial==3),]

p3<-ggplot(t3, aes(avg.inhab, egg.week, shape=treat.ordered, linetype=treat.ordered, col=treat.ordered)) +
  geom_smooth(method="lm", se=FALSE, show.legend = TRUE)  +
  geom_point(size=3)+
  theme_classic()+
  labs(x="Gobies per Reef",y=(expression(atop("Reproductive Output", 
                                              paste((eggs~laid~reef^-1))))))+
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
  #  scale_x_continuous(breaks=c(10,12,14,16,18,20)) + scale_y_continuous(limits = c(0,40000))+
  labs(color  = "Perceived Risk", linetype = "Perceived Risk", shape = "Perceived Risk")
p3

t4<- repro[ which(repro$Trial==4),]

p4<-ggplot(t4, aes(avg.inhab, egg.week, shape=treat.ordered, linetype=treat.ordered, col=treat.ordered)) +
  geom_smooth(method="lm", se=FALSE, show.legend = TRUE)  +
  geom_point(size=3)+
  theme_classic()+
  labs(x="Gobies per Reef",y=(expression(atop("Reproductive Output", 
                                              paste((eggs~laid~reef^-1))))))+
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
  #  scale_x_continuous(breaks=c(10,12,14,16,18,20)) + scale_y_continuous(limits = c(0,40000))+
  labs(color  = "Perceived Risk", linetype = "Perceived Risk", shape = "Perceived Risk")
p4

t5<- repro[ which(repro$Trial==5),]

p5<-ggplot(t5, aes(avg.inhab, egg.week, shape=treat.ordered, linetype=treat.ordered, col=treat.ordered)) +
  geom_smooth(method="lm", se=FALSE, show.legend = TRUE)  +
  geom_point(size=3)+
  theme_classic()+
  labs(x="Gobies per Reef",y=(expression(atop("Reproductive Output", 
                                              paste((eggs~laid~reef^-1))))))+
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
  #  scale_x_continuous(breaks=c(10,12,14,16,18,20)) + scale_y_continuous(limits = c(0,40000))+
  labs(color  = "Perceived Risk", linetype = "Perceived Risk", shape = "Perceived Risk")
p5

t6<- repro[ which(repro$Trial==6),]

p6<-ggplot(t6, aes(avg.inhab, egg.week, shape=treat.ordered, linetype=treat.ordered, col=treat.ordered)) +
  geom_smooth(method="lm", se=FALSE, show.legend = TRUE)  +
  geom_point(size=3)+
  theme_classic()+
  labs(x="Gobies per Reef",y=(expression(atop("Reproductive Output", 
                                              paste((eggs~laid~reef^-1))))))+
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
  #  scale_x_continuous(breaks=c(10,12,14,16,18,20)) + scale_y_continuous(limits = c(0,40000))+
  labs(color  = "Perceived Risk", linetype = "Perceived Risk", shape = "Perceived Risk")
p6

#patching them together in one frame

p1 + p2 + p3 + p4 + p5 + p6 # too big for frame

png("Output/repro.trials.stacked.egg.week.png", width = 15, height = 10, units = 'in', res = 300)
(p1 / p2 / p3) | (p4 / p5 / p6)
dev.off()


# Description: Script for reproductive output from Jarvis and Steele 
# Author: George C Jarvis
# Date: Thu Apr 09 15:14:09 2020
# Notes:
# --------------

rm(list=ls())

library(lme4)
library(lmerTest)
library(dplyr)
library(ggplot2)
library(visreg) #visualizing linear models

options(contrasts = c("contr.sum","contr.poly")) #this is important, run before ANOVA, will set SS to type III

#importing dataset, adding number of gobies on each reef, ordering treatments####
repro<-read.csv("Data/goby_reproduction.2020.4.9.csv")

#basic data viz
pairs(repro) #pretty cool to see correlations

#data manipulation####

#sqrt-transforming egg counts

#adding column for average density, rounded to nearest whole number of fish
repro$avg.inhab<-(ceiling((repro$Recollection+20)/2))

#adding a column for year (as a proxy for tagging procedure), where trials 1-3 = 2017, and 4-6 = 2018
repro$Year <- ifelse(repro$Trial <=3, 2017, 2018)
#want it as a factor? Going to make a variable with it as a factor, run the model, and see if I get different results
repro$Year.fact<- as.factor(repro$Year)
repro$Trial.fact<-as.factor(repro$Trial)

#adding egg/week variable to the dataset
repro<-repro %>%
  mutate(egg.week = ifelse(Trial<4, Egg.count/1,
                           ifelse(Trial == 4| Trial == 5, (ceiling(Egg.count/4)),
                                  ifelse(Trial == 6, (ceiling(Egg.count/2)), NA))))

#exporting wrangled data
write.csv(repro,"Data\\egg.counts.2019.12.23.csv", row.names = FALSE)
write.csv(repro,"Data\\egg.counts.2020.4.7.csv", row.names = FALSE)

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

# trying to see if I can reproduce the SYSTAT results in R 2020.3.19####

#using lme4 package
moda<-lmer(egg.week~Treatment*Year.fact*avg.inhab+(1|Year.fact/Trial),data=repro)
hist(resid(moda))
qqnorm(resid(moda))
qqline(resid(moda))
anova(moda)
Anova(moda)
summary(moda) # not doing it right, denDF are incorrrect

#going to try out some different models with nlme package
# - in Zurr et al. 2009 book on mixed models in R

mod1<-lme(egg.week~Treatment*Year, random = ~1 + Year|Trial, data = repro)
summary(mod1)
anova(mod1)

#testing out other models with the info from B Bolker's mixed model markdown####
#rationale: want to see if I can recreate SYSTAT output in R, meaning I want to 
# - see if I can get the correct error df for the fixed effect of treatment

# - https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html, specifically in the
# - "model definition --> "model speciofication" section

# intercept varying among sites and among blocks within sites (nested random effects)

emod<-lmer(egg.week~Treatment*Year+avg.inhab+(Treatment|Trial/Year),repro,method='ML')
summary(emod) #got some error messages
anova(emod)# not the same, still a lot of error df in treatment

emod<-lmer(egg.week~Treatment*Year*avg.inhab+(1|Treatment:Trial),repro,method='ML')


emod1<-lmer(egg.week~Treatment*Year+avg.inhab+(1|Treatment:Year:Trial),repro,REML=FALSE)
summary(emod1)
anova(emod1)

emod2<-lme(egg.week~Treatment*Year+avg.inhab, random=~1+Treatment|Trial/Year,repro,method='ML')
summary(emod2)
anova(emod2)

#attempting to remodel and compare model output between SYSTAT and R (full models first)

emod<-lmer(egg.week~Treatment*Year*avg.inhab+(1|Trial/Year),repro,REML=FALSE)
summary(emod)
anova(emod)

emoda<-lmer(egg.week~Treatment*Year*avg.inhab+(Treatment|Trial:Year),repro,REML = FALSE)
summary(emoda)

emodb<-lmer(egg.week~Treatment*Year*avg.inhab+(Treatment|Trial/Year),repro,REML = FALSE)
summary(emodb)
anova(emodb)

emodc<-lmer(egg.week~Treatment*Year*avg.inhab+(Treatment*Trial+(Treatment|Trial:Year)),repro,REML = FALSE)
summary(emodc)
anova(emodc)

#2020.4.5 testing out new nesting code in lmer ####

#start with full model, and also want to include nested interactions with covariate

emodx<-lmer(egg.week~Treatment*Year*avg.inhab+(Treatment|Year/Trial),repro,REML=FALSE)
summary(emodx)
anova(emodx)

#same thing but with different syntax, because sometimes the nesting operator does not
## work correctly

emodxi<-lmer(egg.week~Treatment*Year*avg.inhab+(Treatment|Year) + (Treatment|Year:Trial),
             repro,REML=FALSE)
summary(emodxi)
anova(emodxi)
Anova(emodxi)

#specifying the interactions of treatment x trial(year)
## and adding in the effects of the interaction between avg.ihab*Treatment*Trial

#---

#email from Steve on 2020.3.31

#So, because you listed every trial with a unique identifier in the structure of the data file, it is already explicitly crossed with risk and nested  in Year. So, the question is to how to specify the interaction between Treatment x Trial(Year). Either of these should work:
#  (1 | Trial:Treatment)   # or if that doesn't work in case it doesn't like an interaction term as a grouping factor (sometimes that works for me, sometimes it doesn't), then do the following:

#new_variable<- Trial:Treatment   # make a new variable that is a compound of the other two
#(1 | new_variable)     # now you have a single grouping factor with all the levels of the interaction.

#Looking at the plots you sent weeks ago for how the risk treatments show different slopes with goby number in different trials, then arguably you should be estimating separate slopes [vs goby number] for risk treatments in different trials. You can estimate random slopes that are either correlated with the intercepts of each trial, or independently of each trial intercept in the following way:
#  (Goby | new_variable)   # This gives you correlated slopes and intercepts for each risk treatment in each trial

#(Goby || new_variable)   # This gives you uncorrelated slopes and intercepts and is therefore equal to
#(1 | new_variable) + (0 + Goby | new_variable)  # it can also be written as
#(1 | new_variable) + (Goby - 1 | new_variable)

#---

#converting trial to a factor to use in first suggestion from Steve
repro$Trial.factor<-as.factor(repro$Trial)

new_variable<- repro$Trial.factor:repro$Treatment #not as simple as this
View(new_variable)

emodx1<-lmer(egg.week~Treatment*Year*avg.inhab+(1|new_variable),repro,REML=TRUE)
summary(emodx1)
anova(emodx1)

# adding in avg.inhab, testing correlated vs. uncorrelated slopes and intercepts

#correlated slopes
emodx1<-lmer(egg.week~Treatment*Year*avg.inhab+(1|new_variable)
             + (avg.inhab|new_variable),repro,REML=TRUE)
summary(emodx1)
anova(emodx1)

#uncorrelated slopes
emodx2<-lmer(egg.week~Treatment*Year*avg.inhab+(1|new_variable)
             + (avg.inhab||new_variable),repro,REML=TRUE)
summary(emodx2)
anova(emodx2)

#will try with the interaction as the grouping factor: (1 | Trial:Treatment) 

emodxii<-lmer(egg.week~Treatment*Year*avg.inhab+(1|Trial:Treatment),repro,REML=TRUE)
summary(emodxii)
anova(emodxii)

#curious, now year is significant? Might not stay that way if I reduce the model

#adding in the effect of avg.inhab, first with correlated, then with uncorrelated intercepts

emodxiii<-lmer(egg.week~Treatment*Year*avg.inhab+(1|Trial:Treatment) + (avg.inhab|Trial:Treatment),
               repro,REML=TRUE)
summary(emodxiii)
anova(emodxiii)

#including avg.inhab with uncorrelated slopes and intercepts

emodxiv<-lmer(egg.week~Treatment*Year*avg.inhab+(1|Trial:Treatment) + (avg.inhab||Trial:Treatment),
              repro,REML=TRUE)
summary(emodxiv)
anova(emodxiv)

# 2020.4.6 chat with lukas#####
#new variable that combines trial and year
repro$Random <- paste0(repro$Trial.fact, repro$Year.fact)
View(repro$Random)
repro$Random.fact<-as.factor(repro$Random)

#full model
#raw egg counts
lukas.mod<-lmer(egg.week~Treatment*Year*avg.inhab+(1|Random.fact)+(1|Treatment:Random.fact)
                + (1|avg.inhab:Treatment:Random.fact)+(1|avg.inhab:Random.fact),repro,REML=F)

#sqrt-transformed egg counts
lukas.mod<-lmer(sqrt(egg.week)~Treatment*Year*avg.inhab+(1|Random.fact)+(1|Treatment:Random.fact)
                + (1|avg.inhab:Treatment:Random.fact)+(1|avg.inhab:Random.fact),repro,REML=F)

hist(resid(lukas.mod))
qqnorm(resid(lukas.mod))
qqline(resid(lukas.mod))
plot(lukas.mod)

anova(lukas.mod, type = 3)
summary(lukas.mod)

#dropping 3-way interaction with random effects
lukas.mod1 <- update(lukas.mod, .~. -(1|avg.inhab:Treatment:Random.fact))
hist(resid(lukas.mod1))
qqnorm(resid(lukas.mod1))
qqline(resid(lukas.mod1))
plot(lukas.mod1)
anova(lukas.mod1)
summary(lukas.mod1)

#chi-square test comparing models, log-likelihood
logLik(lukas.mod) #'log Lik.' -986.9818 (df=17)

logLik(lukas.mod1) #'log Lik.' -986.9818 (df=16)

#compare logLik of full model to that of reduced model

2 * (logLik(lukas.mod) - logLik(lukas.mod1)) #chi2 = -2.185584e-08 (df=17)
pchisq(2 * (logLik(lukas.mod) - logLik(lukas.mod1)), df=1, lower.tail = F)
# chi2 = -2.185584e-08, p = 1
anova(lukas.mod)

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



#next removing the treatment*Random (trial(year) interaction)

#--------2020.4.6 clean slate------------####
#can't do anything until I figure out if my model is being analyzed correctly in R

#using information from https://stat.ethz.ch/pipermail/r-sig-mixed-models/2014q4/022840.html
## on how to code an interaction between fixed and random effects

#using lmer, given that Trial(Year) is already specified, given the way that I coded Trial

#need to start simple, but should still be using the "anova" function with my models to see if I can get similar
## results to what I'm seeing in SYSTAT, at least in terms of denDF

#comparing results with data for raw egg counts first, denDF shouldn't change based on transformed vs. raw data

me<-lmer(egg.week ~ Treatment*Year.fact*avg.inhab + (1|Year.fact:Trial.fact) + (1|Treatment:Year.fact:Trial.fact) +
           (1|avg.inhab:Year.fact:Trial.fact) + (1|Treatment:avg.inhab:Year.fact:Trial.fact), REML=F, repro)
summary(me)
anova(me) #Okay, this seems to be estimating denDF for year, avg.inhab, and avg.inhab*year sort of correctly
#but there are tons of error messages
coef(me)

me1.1<-lmer(egg.week ~ Treatment*Year*avg.inhab + (1|Trial) + (1|Treatment:Trial) +
              (avg.inhab||Trial) + (avg.inhab||Treatment:Trial), REML=TRUE, repro)
summary(me1.1)
anova(me1.1)
coef(me1.1)

#going to try some other notation
me1<-lmer(egg.week ~ Treatment*Year*avg.inhab + (1|Random) + (0+Treatment|Random) +
            (1|avg.inhab:Random) + (0+Treatment|avg.inhab:Random), REML=F, repro)
summary(me1)
anova(me1)

#more notation changes
me2<-lmer(egg.week ~ Treatment*Year*avg.inhab + (1|Trial) + (1|Treatment:Trial) +
            (1|Trial:avg.inhab), REML=TRUE, repro)
summary(me2)
anova(me2)

#more
coef(me2)

#I sort of want to try and build my model from the bottom up to see if I can affect the df for treatment
m1<-lmer(egg.week~Treatment+(1|Trial),repro,REML=F)
summary(m1)
anova(m1)
coef(m1)

m2<-lmer(egg.week~Treatment*Year+(1|Trial)+(0+Treatment|Trial),repro,REML=F)
summary(m2)
anova(m2)
coef(m2)

#maybe things get resolved with a fully-reduced model? no...

#I found a source that explains that when the Treatment*Trial(Year) term doesn't have that large of an effect, then the error for that
## term has less of a proportional effect over teh residual effect

#Because the variance associated with the interaction is essentially zero (in the presence of the subnum random main-effect) 
##the interaction term has no effect on the calculation of denominator degrees of freedom, F-values and p-values:

#anova(model3, type=1)
#Type I Analysis of Variance Table with Satterthwaite's method
#Sum Sq Mean Sq NumDF DenDF F value Pr(>F)
#group           12065.3 12065.3     1    18  2.4334 0.1362
#direction        1951.8  1951.8     1  5169  0.3936 0.5304
#group:direction 11552.2 11552.2     1  5169  2.3299 0.1270

#However, subnum:direction is the enclosing error stratum for subnum 
## so if we remove subnum all the associated SSQ falls back into subnum:direction



mred<-lmer(egg.week~Treatment+(1|Trial),repro)
hist(resid(mred))
qqnorm(resid(mred))
qqline(resid(mred))
summary(mred)
anova(mred, type = 1)
coef(mred)

# reading through the Bates et al. 2015 paper on lme4, this could be where I break through ####

#going through the different modules (n=4) in lme4, not sure what this will show, but might help me understand which parameters
## are important in my model?

#not sure what these do

#parsedFormula <- lFormula(formula = Reaction ~ Days + (Days | Subject), + data = sleepstudy)

#module 1 - formula module
parsedFormula <- lFormula(formula = egg.week~Treatment*Year.fact*avg.inhab + (1|Trial.fact) + (1|Treatment:Trial.fact) +
                            (1|Trial:avg.inhab)+(1|avg.inhab:Treatment:Trial.fact), REML=F, repro)
summary(parsedFormula)

# 2 - objective function module
#R> devianceFunction <- do.call(mkLmerDevfun, parsedFormula)

devianceFunction <- do.call(mkLmerDevfun, parsedFormula)

# reanalyzing data with log-likelihood estimates and chi-square tests####
#found out that I have to include Year as "Year.fact", in order to make sure that model is scaled correctly

#model structure 1: all slopes are the same (effcets are the same among trials), but intercepts are random (magnitude differs)
## Note: have to use REML=FALSE if you want to compare models with log-likelihood ratio tests (Pinheiro & Bates, 2000; Bolker et al., 2009)

me<-lmer(egg.week ~ Treatment*Year.fact*avg.inhab + (1|Trial.fact) + (1|Treatment:Trial.fact) +
           (1|Trial.fact:avg.inhab)+(1|avg.inhab:Treatment:Trial.fact), REML=F, repro)
hist(resid(me))
qqnorm(resid(me))
qqline(resid(me))
plot(me)

#trying out performance package for data visualization
library(performance)
png("Output/2020.4.8.model.testing.png", width = 9.5, height = 5.5, units = 'in', res = 300)
check_model(me)
dev.off()

#eh...not so great...

me2<-update(me, .~. -(1|avg.inhab:Treatment:Trial.fact))
summary(me2)
anova(me2) #Okay, this seems to be estimating denDF for year, avg.inhab, and avg.inhab*year sort of correctly
#but there are tons of error messages
coef(me2)

anova(me,me2) #no difference in models when three-waty interaction with random intercept is removed

#next logical removal is the 1|treatment:trial term

me3<-update(me2, .~. -(1|Treatment:Trial.fact))
summary(me3)
anova(me3)

anova(me2,me3) #no difference

#next logical removal is the 1|Trial:avg.inhab term, but I'm skeptical to take that out b/c it accounts for 12% of the
## residual variance in random effects (in fact, there's more variance attributed to this term than trial alone)
## fortunately, can test whwether it makes sense to remove that term with log-likelihood test comapring m3 with m4

#the key question will be whether to evaluate the effect of avg.inhab with random slope (different effect of covariate among trials),
## and random intercept (same effect of covariate, but different magnitude among trials)

#should compare model with avg.inhab|Trial.fact vs. 1|Trial:avg.inhab and see what log-likelihood shows

me4<-update(me3, .~. -(1|Trial:avg.inhab)) #all interactions with covariate that were N.S. were removed, including random effects
## in fact, it may put more variance in the overall model, including variance for fixed effects
summary(me4) # trial accounts for 6% of the residual variance
anova(me4)

anova(me3,me4) # P = 0.3939, doesn't suggest any differnce due to the dropping of the random effect of Trial:avg.inhab

#so, decided to leave the random factor of 1|Trial in the model, after removing all other nonsignificant randm effects at P>0.25
## as tested for with log-likelihood tests
#NOTE: this model structure assumes all random effects had random intercepts (1|Random...), but not random slopes among trials

#okay, now doing the same thing with fixed factors, going to start with highest-order interactions with the covariate

me5<-update(me4, .~. -(Treatment:Year.fact:avg.inhab))
summary(me5) # trial accounts for 6% of the residual variance
anova(me5)

anova(me4,me5)# P = 0.5013

#now removing Year.fact:avg.inhab (the order of removal seems arbitrary, but that's the next one)
me6<-update(me5, .~. -(Year.fact:avg.inhab))
summary(me6) # trial variance is still 6% of residual variance for random effects
anova(me6)

anova(me5,me6)# P = 0.779

#now want to take out Treatment:avg.inhab

me7<-update(me6, .~. -(Treatment:avg.inhab))
summary(me7) # trial variance is still 6% of residual variance for random effects
#NOTE: this summary shows me the estimate for the effect of avg inhab (741.827),
## indicating that for every goby added, there's an increase of about 740 eggs on average

#the trouble that I'm anticipating is that there might be a case where one of the years or the interactions
## is significant, not just the term itself (e.g. avg.inhab)

# you'll see that Treatment1 and Treatment2 are in the model, so I guess I can show what treatment si driving the
## trend? We'll cross that bridge when we come to it (likely behavior analyses)

anova(me7)

anova(me6,me7)# P = 0.7064

#dropping the fixed effects that I have reduced the model to do the log-likelihood estimates
#NOTE: m7 is the fully-reduced model, so just have to jeep iterating that model +/- individual fixed factors

# - avg.inhab, but keeping in other fixed factors

## should see a sig. log-likelihood result
me8<-update(me7, .~. -(avg.inhab))
summary(me8) # residual variance for random effects just shot up a ton (soaked up all variance from avg.inhab)
anova(me8)

anova(me7,me8)#chisq = 38.021      df = 1      p = 7e-10, super significant, so don't want to remove that term

# - Treatment*Year.fact

me9<-update(me7, .~. -(Treatment:Year.fact))
summary(me9)
anova(me9)

anova(me7,me9) #chisq = 0.5296     df= 2     p = 0.7674

# - Treatment

me10<-update(me7,.~. -(Treatment))
summary(me10)
anova(me10)

anova(me7,me10) #chisq = 0     df = 0         p = 1 #this is acceptable, amanda has this in her paper for a factor

#same analysis with log-likelihood estimates, but with sqrt.eggs as response####
## all I have to do is change the orignial model

me<-lmer(sqrt(egg.week) ~ Treatment*Year.fact*avg.inhab + (1|Trial.fact) + (1|Treatment:Trial.fact) +
           (1|Trial.fact:avg.inhab)+(1|avg.inhab:Treatment:Trial.fact), REML=F, repro)
hist(resid(me)) #I honestly think that 
qqnorm(resid(me))
qqline(resid(me))
plot(me)

summary(me)
anova(me) #Okay, this seems to be estimating denDF for year, avg.inhab, and avg.inhab*year sort of correctly
#but there are tons of error messages
coef(me)

me2<-update(me, .~. -(1|avg.inhab:Treatment:Trial.fact))
summary(me2)
anova(me2) #Okay, this seems to be estimating denDF for year, avg.inhab, and avg.inhab*year sort of correctly
#but there are tons of error messages
coef(me2)

anova(me,me2) #no difference in models when three-waty interaction with random intercept is removed

#next logical removal is the 1|treatment:trial term

me3<-update(me2, .~. -(1|Treatment:Trial.fact))
summary(me3)
anova(me3)

anova(me2,me3) #no difference

#next logical removal is the 1|Trial:avg.inhab term, but I'm skeptical to take that out b/c it accounts for 12% of the
## residual variance in random effects (in fact, there's more variance attributed to this term than trial alone)
## fortunately, can test whwether it makes sense to remove that term with log-likelihood test comapring m3 with m4

#the key question will be whether to evaluate the effect of avg.inhab with random slope (different effect of covariate among trials),
## and random intercept (same effect of covariate, but different magnitude among trials)

#should compare model with avg.inhab|Trial.fact vs. 1|Trial:avg.inhab and see what log-likelihood shows

me4<-update(me3, .~. -(1|Trial:avg.inhab)) #all interactions with covariate that were N.S. were removed, including random effects
## in fact, it may put more variance in the overall model, including variance for fixed effects
summary(me4) # trial accounts for 6% of the residual variance
anova(me4)

anova(me3,me4) # P = 0.3939, doesn't suggest any differnce due to the dropping of the random effect of Trial:avg.inhab

#so, decided to leave the random factor of 1|Trial in the model, after removing all other nonsignificant randm effects at P>0.25
## as tested for with log-likelihood tests
#NOTE: this model structure assumes all random effects had random intercepts (1|Random...), but not random slopes among trials

#okay, now doing the same thing with fixed factors, going to start with highest-order interactions with the covariate

me5<-update(me4, .~. -(Treatment:Year.fact:avg.inhab))
summary(me5) # trial accounts for 6% of the residual variance
anova(me5)

anova(me4,me5)# P = 0.5013

#now removing Year.fact:avg.inhab (the order of removal seems arbitrary, but that's the next one)
me6<-update(me5, .~. -(Year.fact:avg.inhab))
summary(me6) # trial variance is still 6% of residual variance for random effects
anova(me6)

anova(me5,me6)# P = 0.779

#now want to take out Treatment:avg.inhab

me7<-update(me6, .~. -(Treatment:avg.inhab))
summary(me7) # trial variance is still 6% of residual variance for random effects
#NOTE: this summary shows me the estimate for the effect of avg inhab (741.827),
## indicating that for every goby added, there's an increase of about 740 eggs on average

#the trouble that I'm anticipating is that there might be a case where one of the years or the interactions
## is significant, not just the term itself (e.g. avg.inhab)

# you'll see that Treatment1 and Treatment2 are in the model, so I guess I can show what treatment si driving the
## trend? We'll cross that bridge when we come to it (likely behavior analyses)

anova(me7)

anova(me6,me7)# P = 0.7064

#dropping the fixed effects that I have reduced the model to do the log-likelihood estimates
#NOTE: m7 is the fully-reduced model, so just have to jeep iterating that model +/- individual fixed factors

# - avg.inhab, but keeping in other fixed factors

## should see a sig. log-likelihood result
me8<-update(me7, .~. -(avg.inhab))
summary(me8) # residual variance for random effects just shot up a ton (soaked up all variance from avg.inhab)
anova(me8)

anova(me7,me8)#chisq = 38.021      df = 1      p = 7e-10, super significant, so don't want to remove that term

# - Treatment*Year.fact

me9<-update(me7, .~. -(Treatment:Year.fact))
summary(me9)
anova(me9)

anova(me7,me9) #chisq = 0.5296     df= 2     p = 0.7674

# - Treatment

me10<-update(me7,.~. -(Treatment))
summary(me10)
anova(me10)

anova(me7,me10)

#not done yet...

#next step will be to go back and rewrite the model where the slope for avg.inhab AND the intercept for avg.inhab is random

#the point that I'm not sure about is whether to code the interactions as avg.inhab|Trial, or as avg.inhab|avg.inhab:Trial

# it seems to me that it is redundant, and that I might want to include avg.inhab|Trial in addition to 1|avg.inhab:Trial

#stock code to be applied to my data
lme(y ~ time * tx, 
    random = list(therapist = ~time * tx, 
                  subjects = ~time),
    data=df)

newmod<-lme(egg.week~Treatment*Year.fact*avg.inhab, 
            random=list(Trial.fact=~1,Trial.fact:avg.inhab=~1),data=repro, method="ML")

#comparing lme to lmer models ####
# have to put all of the interactions into single terms so that I can list them in 
repro$trial_inhab_treatment<-paste0(repro$Trial.fact,repro$avg.inhab,repro$Treatment)
repro$trial_treatment<-paste0(repro$Trial.fact,repro$Treatment)
repro$trial_inhab<-paste0(repro$Trial.fact,repro$avg.inhab)
#View(repro$trial_inhab)

#maybe if avg.inhab were a factor?
repro$trial_inhab_factor_treatment<-paste0(repro$Trial.fact,as.factor(repro$avg.inhab),repro$Treatment)
repro$trial_treatment<-paste0(repro$Trial.fact,repro$Treatment)
repro$trial_inhab_factor<-paste0(repro$Trial.fact,as.factor(repro$avg.inhab))

#full model
newmod<-lme(egg.week~Treatment*Year.fact*avg.inhab, 
            random=list(Trial.fact=~1,trial_treatment=~1,trial_inhab=~1,trial_inhab_treatment=~1),data=repro, method="ML")
summary(newmod)
anova(newmod)
#seems like there's something going on with the trial_inhab term, and that's why the models aren't the same, and also why all terms
## that have avg.inhab included in them do not have lower dendf

newmod1<-lme(egg.week~Treatment*Year.fact*avg.inhab, 
             random=list(Trial.fact=~1,trial_treatment=~1,trial_inhab_factor=~1,trial_inhab_factor_treatment=~1),data=repro, method="ML")
summary(newmod1)
anova(newmod1)

#now want to go back and chack the full model with all of the random effects to see if denom df are correct!

#lmer full model
me<-lmer(egg.week ~ Treatment*Year.fact*avg.inhab + (1|Trial.fact) + (1|Treatment:Trial.fact) +
           (1|Trial.fact:avg.inhab)+(1|avg.inhab:Treatment:Trial.fact), REML=F, repro)
summary(me)
anova(me) #Okay, this seems to be estimating denDF for year, avg.inhab, and avg.inhab*year sort of correctly
#but there are tons of error messages
#coef(me)

#SUMMARY: not entirely the same, although it does have a closer df 

boxplot(repro$egg.week~repro$avg.inhab)

n<-1883/(sqrt(379033+135203+3245737))
n

View(repro)

#testing out visreg package with my model to try and figure out why the interept is acting funny

#anyway, not going to put too much thought into this, but it seems weird
library(visreg)

model <- lmer(sqrt(egg.week) ~ Treatment + avg.inhab + Year.fact + Treatment:Year.fact +(1|Trial.fact), data=repro)
visreg(model, "Treatment", by="Year.fact")
visreg(model, "Treatment")

#not sure why, when I include avg.inhab in the model,
## that it shows values for reproduction that are less than 0?
## It might have something to do with the structure of the model?
## It's only noticeable in the plot for 2017, which shows negative values

or (if no interaction):
  visreg(model, ‘Treatment’)



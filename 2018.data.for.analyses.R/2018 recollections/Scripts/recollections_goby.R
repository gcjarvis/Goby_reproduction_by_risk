# Description: Script for recollections (i.e. survival) from Jarvis and Steele
# Author: George C Jarvis
# Date: Sun May 24 13:48:41 2020
# Notes:
# --------------

rm(list=ls())

library(lme4)
library(lmerTest)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(emmeans) #for generating least-squares adjusted means from models

options(contrasts = c("contr.sum","contr.poly")) #this is important, run before ANOVA, will set SS to type III

# import dataset ####
reco<-read.csv("Data/2019.10.8.recollection.data.csv")

# data manipulation ####

#adding column for survivorship, dividing recollections by 20 (initial number of fish on reefs)
reco$Survivorship<-reco$Count/20

#changing trial and year to factors
reco$Trial<- as.factor(reco$Trial)
reco$Year<- as.factor(reco$Year)


#subsetting data from trial 6 to compare recollections between HR caged and uncaged treatments
reco.t6<-reco[reco$Trial==6,]

#subset treatments, caged and uncaged only (T6.comparison)
#subsetting by multiple character factors within a variable. Awesome code!
reco.t6.comp<-reco.t6[reco.t6$T6.comparison == "High" | reco.t6$T6.comparison == "Control", ]
head(reco.t6.comp)
view(reco.t6.comp)

#only need to do a t-test for trial 6
# are there differences in survivorship between caged and uncaged reefs?

#need to reframe the data for trial 6, putting it into long format

reco.t6.comp<-as_tibble(reco.t6.comp)
reco.t6.comp

t6.comp<-reco.t6.comp %>% pivot_wider(names_from = T6.comparison, values_from = Survivorship)
View(t6.comp)

#exporting data, then reimporting back into R
write.csv(t6.comp,"Data\\2020_24_5_recollections_t6_comp.csv", row.names = FALSE)

#removed NA's and empty rows, also removed all unnecessary columns

t6.comp.wrangled<- read.csv("Data\\2020_24_5_recollections_t6_comp.csv")
View(t6.comp.wrangled)

# analyses ####

#Trials 1-6:
# analyzing mixed models with log-likelihood estimates and chi-square tests.
# The full model will include treatment, year, the interaction (all fixed), 
## plus trial nested within year, and treatment x trial nested within year as random effects.
# I am planning to keep Trial in the model, but treatment x trial will likely be removed if it is not significant
#NOTE: Trial 6 data were included in this analysis, but all high-risk reefs (caged and uncaged)
## were considered "High risk" in the analysis
# Start with full model, then reduce model first by non-significant random effects
# Will keep Treatement, Year, and T x Y as fixed effects, because those are the factors I'm interested in

#Trial 6 only:
# testing for differences in survival between high-risk caged and uncaged treatments
# t-test

# 1a. All trials ####
# full model
re<-lmer(Survivorship ~ Treatment*Year + (1|Trial) + (1|Treatment:Trial), REML=F, reco)
hist(resid(re))
qqnorm(resid(re))
qqline(resid(re))
plot(re)
summary(re) #no extra variance explained by random effect of treatment x trial
anova(re)

#removing two-way interaction of random effect
re2<-update(re, .~. -(1|Treatment:Trial))
summary(re2)
anova(re2)

anova(re,re2) #no difference in models when two-way interaction with random intercept is removed, taking it out

#final model is re2 (Survivorship ~ Treatment*Year + (1|Trial))
# means for survival from final model
#me8 is the final model that I will likely end up with (those are all of the fixed and random factors that I'm interested in)
hist(resid(re2))#not bad, slightly right-skewed
#LS-adjusted means from model
emmeans(re2, pairwise~Treatment) #warning message re: interactions, but I think it's okay
boxplot(Survivorship~Treatment,data=reco)# variances don't look too bad, and medians look pretty good to me (i.e. match up with LS-means for the most part)

# now can test for effects individually

# - trial
re3<-update(re, .~. -(1|Trial))
summary(re3)
anova(re3)

anova(re3, re) #there were differences by trial, but not the important factor here, 
# so I will keep it in the model because it explains some of the variance

# - treatment:year
re4<- update(re2, .~. -(Treatment:Year))
summary(re4)
anova(re4)

anova(re4, re2)

# - year, can't also have the interaction term of treatment:year, so working from re4
re5<- update(re4, .~. -(Year))
summary(re5)
anova(re5)

anova(re5, re4)

# - treatment, working off of re4 for same reason as before
re6<- update(re4, .~. -(Treatment))
summary(re6)
anova(re6)

anova(re6, re4)

# 1b. T-test for trial 6 (survival for high-risk caged vs. uncaged) ####

t.test(t6.comp.wrangled$Control,t6.comp.wrangled$High)
# no sig. difference in survival based on caging 
# proportion surviving: high-risk caged = 0.19; high-risk uncaged = 0.14)
# reported these means for caging effects in the paper

# plotting (all trials, but not caging effects in trial 6) ####

survival<-with(reco, aggregate((Survivorship), list(Treatment=Treatment), mean))
survival$se<-with(reco, aggregate((Survivorship), list(Treatment=Treatment), function(x) sd(x)/sqrt(length(x))))[,2]

png("Output/2019.11.17.survivorship.6.5x5.5.300dpi.png", width = 6.5, height = 5.5, units = 'in', res = 300)

reco.plot<- ggplot(survival, aes(x=Treatment, y=x, fill=Treatment)) +
  geom_bar(stat="identity", colour= "black", width = 0.5, position="dodge")+ 
  scale_x_discrete(limits=c("Low","Medium","High"))+
  theme_classic() + 
  labs(x="Risk Treatment", y="Survivorship") +
  theme(legend.position="none") + 
  scale_fill_manual(values=c("grey", "grey", "grey")) +
  theme(axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"), 
        axis.title=element_text(size=20))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), 
        axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0)), 
        axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(legend.text=element_text(size=18)) +
  theme(legend.title =element_text(size=20))+
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0,0),limits = c(0,0.31),
                                                             labels = scales::number_format(accuracy = 0.01))
reco.plot + geom_linerange(aes(ymin=x-se, ymax=x+se), size=0.5,   
                           position=position_dodge(.85)) + theme(text = element_text(family="Arial"))
dev.off()


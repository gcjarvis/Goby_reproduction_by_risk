# Description: script for analyses and plot for survival from Jarvis and Steele
# Author: George C Jarvis
# Date: Tue Jun 15 17:30:31 2021
# Notes:
# --------------

library(lme4)
library(lmerTest)
library(tidyverse)
library(car)
library(emmeans) #for estimating marginal means from models
library(sciplot)

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
#head(reco.t6.comp)
#view(reco.t6.comp)

#only need to do a t-test for trial 6
# are there differences in survivorship between caged and uncaged reefs?

#need to reframe the data for trial 6, putting it into long format

reco.t6.comp<-as_tibble(reco.t6.comp)
#reco.t6.comp

t6.comp<-reco.t6.comp %>% pivot_wider(names_from = T6.comparison, values_from = Survivorship)
#View(t6.comp)

# simplifying df for t test in trial 6: did recollections differ between high-risk treatments with partial cages, and those with no cages?

# (control = high-risk uncaged, high = high-risk with partial cages)

t6.HR_uncaged <- t6.comp %>% select(Control) %>% filter(!is.na(Control))
t6.HR_partial_cage <- t6.comp %>% select(High) %>% filter(!is.na(High))

t6.comp.wrangled <- as.data.frame(c(t6.HR_uncaged,t6.HR_partial_cage))
View(t6.comp.wrangled)

# analyses ####

#Notes:

#Trials 1-6:

# analyzing mixed models with likelihood ratio tests ('lme4' for model construction, and 'car' package for analysis of deviance tests).
# The full model will include treatment, year, the interaction (all fixed), plus trial as a random effect

#NOTE: Trial 6 data were included in this analysis, but all high-risk reefs (with and without cages) were considered "High risk" in the analysis

# Start with full model, then reduce model first by non-significant random effects

#Trial 6 only:

# testing for differences in survival between high-risk caged and uncaged treatments
# t-test

# 1a. All trials ####

re<-lmer(Survivorship ~ Treatment + Year + (Treatment*Year) + (1|Trial), REML=F, reco)
hist(resid(re))
qqnorm(resid(re))
qqline(resid(re))
summary(re) #no extra variance explained by random effect of treatment x trial
Anova(re, type = "III")

# results:

# > Anova(re, type = "III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: Survivorship
# Chisq Df Pr(>Chisq)    
# (Intercept)    68.8702  1     <2e-16 ***
#   Treatment       0.4396  2     0.8027    
# Year            0.1135  1     0.7362    
# Treatment:Year  0.9213  2     0.6309    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# removing treatment*year

re1<-lmer(Survivorship ~ Treatment + Year+ (1|Trial), REML=F, reco)
Anova(re1, type = "III")

# results:

# > Anova(re1, type = "III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: Survivorship
# Chisq Df Pr(>Chisq)    
# (Intercept) 66.7706  1   3.05e-16 ***
#   Treatment    0.4710  2     0.7902    
# Year         0.0958  1     0.7569    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# testing the random effect of trial:

# making dummy variable with randomized values for trial

set.seed(595)
reco$trial_rand <- sample(reco$Trial)

# running model with trial and dummy variable included as random effects

re_dr<-lmer(Survivorship~ Treatment + Year +
                 (1|Trial) + (1|trial_rand), data = reco)
summary(re_dr) # trial explain 3% of the residual variance, low, but think it's worth including
Anova(re_dr, type = "III")

# removing trial effect, then running likelihood ratio test

re_nt <- lmer(Survivorship ~ Treatment + Year + (1|trial_rand), data = reco)
summary(re_nt)

2*(logLik(re_dr) - logLik(re_nt)) # Chi2 =  10.68
pchisq(2*(logLik(re_dr) - logLik(re_nt)), df = 1, lower.tail=F) # P < 0.01 df = 1

ranef(re1)

# results:

# > ranef(re1)
# $Trial
# (Intercept)
# 1 -0.02543350
# 2 -0.04811743
# 3  0.07355093
# 4  0.06175552
# 5  0.01432549
# 6 -0.07608101

# checking emmeans

emmeans(re1, pairwise ~ Treatment)

# $emmeans
# Treatment emmean     SE   df lower.CL upper.CL
# High       0.233 0.0422 18.6    0.144    0.321
# Low        0.248 0.0432 20.9    0.158    0.338
# Medium     0.258 0.0432 20.9    0.168    0.348
# 
# Results are averaged over the levels of: Year 
# Degrees-of-freedom method: kenward-roger 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast      estimate     SE  df t.ratio p.value
# High - Low     -0.0152 0.0375 107 -0.405  0.9138 
# High - Medium  -0.0252 0.0375 107 -0.671  0.7806 
# Low - Medium   -0.0100 0.0385 106 -0.260  0.9634 
# 
# Results are averaged over the levels of: Year 
# Degrees-of-freedom method: kenward-roger 
# P value adjustment: tukey method for comparing a family of 3 estimates 

# plotting

# estimated marginal means

re1.emm <- emmeans(re1, ~ Treatment)
plot(re1.emm)

# plots of raw data by treatment

bargraph.CI(x.factor = Treatment, response = Survivorship, main="raw data - survival vs. treatment", 
            xlab="Treatment", ylab="survival", data = reco)

# raw means seem to be congruent with marginal means, plot raw means by treatment

# 1b. T-test for trial 6 (survival for high-risk caged vs. uncaged) ####

t.test(t6.comp.wrangled$Control,t6.comp.wrangled$High)

# results:

# > t.test(t6.comp.wrangled$Control,t6.comp.wrangled$High)
# 
# Welch Two Sample t-test
# 
# data:  t6.comp.wrangled$Control and t6.comp.wrangled$High
# t = -0.96225, df = 6.1132, p-value = 0.3724
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.1765767  0.0765767
# sample estimates:
#   mean of x mean of y 
# 0.14      0.19 

# no sig. difference in survival based on caging 
# reported these means for caging effects in the paper

# plotting (all trials, but not caging effects in trial 6) ####

survival<-with(reco, aggregate((Survivorship), list(Treatment=Treatment), mean))
survival$se<-with(reco, aggregate((Survivorship), list(Treatment=Treatment), function(x) sd(x)/sqrt(length(x))))[,2]

#png("Output/2019.11.17.survivorship.6.5x5.5.300dpi.png", width = 6.5, height = 5.5, units = 'in', res = 300)

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
#dev.off()


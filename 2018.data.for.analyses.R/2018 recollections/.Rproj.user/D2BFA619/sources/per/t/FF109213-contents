##############################################################################
# recollections data to get at proportion of population that changed sex     #
#                                                                            #
#   Hunter and George                                                        #
#    7/27/2018                                                               #
##############################################################################

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
sc<-read.csv("Data/Recollections_Trial_4.5.hunter.analyses.csv")
sc$Trial<-as.factor(sc$Trial)#made trial a factor
#only including adults that we tagged initially (immatures were listed as NA in datasheet)
sc<- sc[complete.cases(sc), ]#122 rows of data
View(sc)

head(sc)
tail(sc)

sc.mod<-lm(Sex.change~Treatment, data=sc)
hist(resid(sc.mod))
qqnorm(resid(sc.mod))
boxplot(resid(sc.mod))

anova(sc.mod)
summary(sc.mod)


hist(resid())

#plotting prop sex change by risk
bargraph.CI(x.factor = Treatment, response = Sex.change, main="Proportion sex change by risk", data = sc)
#seems like no sig difference in prop sex change among treatments

#let's do the same plot but grouped by trial
bargraph.CI(x.factor = Treatment, response = Sex.change, group = Trial, legend=TRUE, main="Proportion sex change by risk", data = sc)

bargraph.CI(x.factor = Trial, response = Sex.change, main="Proportion sex change by trial", data = sc)
#running lm to see if there are sig.differences in sex between trials
sex.change.trial.mod<-lm(Sex.change~Trial, data=sc)
hist(resid(sc.mod))
qqnorm(resid(sc.mod))
boxplot(resid(sc.mod))

anova(sex.change.trial.mod)
summary(sex.change.trial.mod)

#tukey test
library(agricolae)
HSD.test(mod1.day.new,"Treatment", group=TRUE,console=TRUE)

#logit transformation
#logit(p, percents=range.p[2] > 1, adjust) #base code, will revisit
logit(sc$Sex.change, percents=range [2] > 1, adjust)

#glm model

sex.change.glm<-glmer(Sex.change~Treatment*Initial.Sex+(1|Trial),family=binomial(), data=sc)
hist(resid(sex.change.glm))
qqnorm(resid(sex.change.glm))
qqline(resid(sex.change.glm))

glmmformula <- update(modelformula, . ~ . + (1|TRTID10))
glmmformula <- update(Sex.change~Treatment*Initial.Sex+(1|Trial), data=sc)
sex.change.nb <- glmer.nb(sex.change.glm, data = scaled.mydata)

summary(sex.change.glm)
anova(sex.change.glm)
Anova(sex.change.glm, type="II")
Anova(sex.change.glm, type="III")

#nbglmm <- glmer.nb(glmmformula, data = scaled.mydata)




##getting means and SE for uhnter
sex.fig<-(with(sc, aggregate((Sex.change), list(Treatment=Treatment), mean)))
#now apply the se function to the 3rd column [,3]

sex.fig$se<-with(sc, aggregate((Sex.change), list(Treatment=Treatment), function(x) sd(x)/sqrt(length(x))))[,2]
sex.fig

sc.mod<-lm(Sex.change~Treatment, data=sc)
hist(resid(sc.mod))
qqnorm(resid(sc.mod))
boxplot(resid(sc.mod))
anova(sc.mod)
summary(sc.mod)

hist(resid())

######trying a logit transformation######
#logit=log( p / (1 - p))

###plotting####
sex.df<-with(sc, aggregate((Sex.change), list(Treatment=Treatment), mean))
#now apply the se function to the 4th column [,3]
sex.df$se<-with(sc, aggregate((Sex.change), list(Treatment=Treatment), function(x) sd(x)/sqrt(length(x))))[,2]
sex.df
#gob.sub$treatment<-ordered(gob.sub$treatment,levels=c("Low","Medium","High"))

sex.change.plot<- ggplot(sex.df, aes(x=Treatment, y=x, fill=Treatment)) +
  geom_bar(stat="identity", colour= "black", width = 0.65, position="dodge")+ 
  scale_x_discrete(limits=c("Low","Medium","High"))+
  theme_classic() + 
  scale_fill_manual(values=c("#0072B2","#D55E00","#009E73")) + 
  theme(axis.text.x=element_text(size=20, colour="black"),axis.text.y=element_text(size=20, colour="black"), axis.title=element_text(size=20,face="bold")) +
  theme(axis.title.y = element_text(size= 20, margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)), axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.ticks.x = element_blank()) + scale_y_continuous(expand = c(0, 0))
sex.change.plot + geom_linerange(aes(ymin=x-se, ymax=x+se), size=0.5,   
                           position=position_dodge(.7)) + theme(text = element_text(family="Arial")) +
  labs(x="Risk Treatment", y="Proportional Sex Change per Reef")

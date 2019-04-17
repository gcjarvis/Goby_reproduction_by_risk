######################################
#                                    #      
#   Egg counts with Trials 3 ways    #
#   a. Combined                      #
#   b. 2017 data only                #
#   c. 2018 data only                #
# 2/14/19                            #
######################################

#note: if you're having trouble with packages, go to this file path
# C:\Users\George\Documents\R\win-library\3.5
# then delete the folder for the package that you're struggling with
# then reinstall the packages
#   you may have to do this for the dependencies of packages as well

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

#A. load data (egg counts and densities combined into one for raw counts and per capita)
#raw data for egg counts
gob.sub<-read.csv("Data/List of nests to count.1.10.19.csv")
gob.sub$Week<-as.factor(gob.sub$Week)
gob.sub$Treatment<-ordered(gob.sub$Treatment,levels=c("Low","Medium","High","Control"))
gob.sub$Trial<-as.factor(gob.sub$Trial)
#bringing in density
density<-read.csv("Data/density.2019.3.6.csv")
density$Trial<-as.factor(density$Trial)
#density<-read.csv("C:\\Users\\George\\Desktop\\2018 summer\\2018 Goby\\2018 data for analyses, R\\2018 densities\\Data//density.1.14.19.csv")
#../goes up a level, instead of having to force R to go to a separate folder
density$Week<-as.factor(density$Week)
density$Day<-as.factor(density$Day)
#needed to select only the columns I needed for analyses to deal with NA's
density<-density[,c("Trial","Date","Week","Day","Reef","Treatment","Density","den.max","deployment.day")]
#complete cases only, which allows me to claculate average weekly densities
# I was having issues with Na's in my datasets if I had a single NA values for missed densities
density<-density[complete.cases(density),] #lost 3 rows
#also want to bring in biomass recollected per reef
reco.raw<-read.csv("Data/recollection.data.for.biomass.2.21.19.csv")
reco.raw$Trial<-as.factor(reco.raw$Trial)

#B.setting up df's####

#total reproduction per week by reef per treatment
egg.per.week<-gob.sub %>%
  group_by(Trial,Reef,Week,Treatment) %>%
  summarize(Egg.count = sum(Egg.count))
#egg.per.week

#density per reef per week with raw counts
den.per.week<-density %>%
  group_by(Trial,Reef,Week,Treatment,deployment.day) %>%
  summarize(Density = ceiling(mean(Density)))

#density per reef per week with den.max ("adjusted") counts
den.per.week<-density %>%
  group_by(Trial,Reef,Week,Treatment,deployment.day) %>%
  summarize(Density = ceiling(mean(den.max)))
#den.per.week
#'ceiling' rounds the average density value, so you don't end up inflating reproduction averages
#rounds the density value up to the nearest whole number

#recollections and biomass
#clean.reco<-reco.raw %>%
  #group_by(Trial,Reef,Treatment,Initial.biomass.total,Initial.female.biomass, Initial.male.biomass,
           #Final.biomass.total,Final.biomass.female,Final.biomass.male,Avg.biomass.total,
           #Avg.fem.biomass,Avg.male.biomass)

#trying a different way...has to be a better way to do this, but I found a workaround
clean.reco<-reco.raw[,c("Trial", "Reef", "Treatment", "Initial.biomass.total", 
                       "Initial.female.biomass","Initial.male.biomass","Final.biomass.total",
                       "Final.biomass.female","Final.biomass.male","Avg.biomass.total",
                       "Avg.fem.biomass","Avg.male.biomass")]
clean.reco<-clean.reco[complete.cases(clean.reco),]
#View(clean.reco)

#now want to join the two df's
den.and.egg<-left_join(egg.per.week,den.per.week,by=c("Trial","Week","Reef","Treatment"))
den.and.egg#error because I don't have "control" treatment in all trials
tail(den.and.egg)#wanted to see what the error message from line 44 was referring to
#I think the coercion is okay in this case

#now going to calculate per capita output by dividing total egg count by density
den.avg.egg<-mutate(den.and.egg,per.capita.repro= ceiling(Egg.count/Density))#ceiling for this as well (doesn't make sense to have less than a whole egg)
#might consider changing the ceiling portion of this line of code^

den.avg.egg
#View(den.avg.egg)
#have one case where there were no fish and no eggs = 0/0 or NaN in R's calculations
#want to correct that with the following function
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
# now want to apply that function to my df to make that Nan value = 0
den.avg.egg[is.nan(den.avg.egg)] <- 0
#View(den.avg.egg)
#NOTE: den.avg.egg is a df that doesn't contain any of the biomass data

#bring in biomass (total, avg. female, avg.male) recollected
egg.den.bio<-left_join(den.avg.egg,clean.reco,by=c("Trial","Reef","Treatment"))
egg.den.bio
#View(egg.den.bio)
#have to make trial an integer to subset it
egg.den.bio$Trial<-as.integer(egg.den.bio$Trial)
#NOTE: I have the same biomass values for all weeks that were sampled
#   this might be an issue later down the line, but I'll talk to Nyssa about it
#   I'm thinking there might be a way to summarize the data to get one row for 
#   biomass (e.g. average.total, average.male, avg.fem), but I'm not sure how
#   that would get analyzed in the linear model...
#   maybe make a separate object? Not sure, because it won't be applied to all levels
#   for week
#needed to add a column that assigned years to trial numbers
#reference code (this is how to select based on two separate columns, price and volume)
#df <- df %>% 
  #mutate(Portfolio=ifelse((Price<P1)&(Volume<V1),11,NA)) %>% 
  #mutate(Portfolio=ifelse((Price<P1)&(Volume<V2)& is.na(Portfolio),12,Portfolio))

egg.den.bio <- egg.den.bio %>% 
  mutate(Year=ifelse((Trial<4),2017,2018))
#View(egg.den.bio)
egg.den.bio$Treatment<-ordered(egg.den.bio$Treatment,levels=c("Low","Medium","High","Control"))

#want to get output per average total biomass 
egg.den.bio<-mutate(egg.den.bio,egg.per.tot.avg.biomass=ceiling(Egg.count/Avg.biomass.total))#ceiling for this as well (doesn't make sense to have less than a whole egg)
#View(egg.den.bio)

#exporting .csv file for Mark
write.csv(egg.den.bio,'jarvis.egg.count.data.with.den.max.2019.3.6.csv')

#C. subsetting df for a. b. and c. scenarios (see box above)####

#loading .csv file most updated as of 2019.3.6 #######
#egg.den.bio<-read.csv("Data/jarvis.egg.count.data.2019.2.26.csv") #uses raw counts for density

egg.den.bio<-read.csv("Data/jarvis.egg.count.data.with.den.max.2019.3.6.csv") #uses adjusted counts for density

#new dataset as of 2019.4.16, summed egg counts by reef####
#did this for one month, and also only have average densities
#over the whole trial
egg.den.bio<-read.csv("Data/new.data.2019.4.16.csv", na.strings = "") #uses adjusted counts for density

#now data with just t1-3, NEED TO DO NA.STRINGS = "" FOR ANCOVA!!
egg.den.bio<-read.csv("Data/new.data.2019.4.16a.no.t6.csv", na.strings = "") #uses adjusted counts for density
levels(egg.den.bio$Treatment)
View(egg.den.bio)

#just t4.5
egg.den.bio<-read.csv("Data/new.data.2019.4.16a.t1.5.csv", na.strings = "") #uses adjusted counts for density
egg.den.bio<-na.omit(egg.den.bio)
egg.2018.t4.5<-egg.den.bio[(egg.den.bio$Trial>3) & (egg.den.bio$Trial<6), ]

#a. combined df is for combined counts
#den.avg.egg
egg.den.bio<-na.omit(egg.den.bio)
#b. 2018 only
egg.2017.t1.2.3<-egg.den.bio[egg.den.bio$Trial<4,]
#egg.2017.t1.2.3 <- egg.2017.t1.2.3[egg.2017.t1.2.3$Treatment = "Low" |"Medium" |"High"),] # select only rows with any value < 8
egg.2017.t1.2.3<-subset(egg.2017.t1.2.3,Treatment!="Control")
levels(egg.2017.t1.2.3$Treatment)
View(egg.2017.t1.2.3)
egg.2017.t1.2.3<-na.omit(egg.2017.t1.2.3)

subset(dataframe, A==B & E!=0)
egg.2017.t1.2.3<-egg.2017.t1.2.3[!(egg.2017.t1.2.3$Treatment=="Control"),]
d<-d[!(d$A=="B" & d$E==0),]

egg.2017.t1.2.3$Treatment<-ordered(egg.2017.t1.2.3$Treatment,levels=c("Low","Medium","High"))
egg.2017
tail(egg.2018)


#c. 2018 t4 and t5 only only
egg.2018.t4.5<-egg.den.bio[(egg.den.bio$Trial>3) & (egg.den.bio$Trial<6), ]
egg.2018$Treatment<-ordered(egg.2018$Treatment,c("Low","Medium","High"))
#egg.2018<-den.avg.egg[den.avg.egg$Trial>3, c("Trial", "Reef", "Week", "Treatment", "Egg.count", "Density","per.capita.repro")]
#egg.2018<-egg.den.bio[egg.den.bio$Trial>3,]
#egg.2018$Treatment<-ordered(egg.2018$Treatment,levels=c("Low","Medium","High","Control"))
#egg.2018
#tail(egg.2018)#notice that I only ran 2-wk trials for trial 6
#I'm not sure if I mentioned that in the methods that I sent to Mark, maybe yes!

#seeing if there is a treatment effect within week 3 and 4
egg.2018.wk3.4<-egg.2018.t4.5[(egg.2018.t4.5$Week>1),]
df<-egg.2018.wk3.4

#need to make a df that just contains the last trial (trial 6)
egg.2018.t6<-egg.den.bio[egg.den.bio$Trial==6,]
egg.2018.t6

#2.14.19 --> log-transforming count data didn't seem to make it better in analyses
#adding column to df for log-transformed egg count data
#have to do log+1 because of times when egg counts were 0
#den.avg.egg$ln.ec<-log(den.avg.egg$Egg.count+1)

#adding column for log+1 per capita output
#den.avg.egg$ln.pc<-log(den.avg.egg$per.capita.repro+1)


#manipulating df variables######

#order treatments
den.avg.egg$Treatment<-ordered(den.avg.egg$Treatment,levels=c("Low","Medium","High","Control"))
#store week as a factor for analyses
den.avg.egg$Week<-as.factor(den.avg.egg$Week)
#store trial as a factor for analyses
den.avg.egg$Trial<-as.factor(den.avg.egg$Trial)
den.avg.egg
#df<-den.avg.egg
df<-egg.den.bio #combined
df<-egg.2017.t1.2.3#2017 only
#df$Trial<-as.factor(df$Trial)
df<-egg.2018.t4.5 #2018 only
df<-egg.2018.t6
#df$Treatment<-factor(df$Treatment, ordered=FALSE)#unordered factor so I can interpret results
df$Treatment<-ordered(df$Treatment,levels=c("Low","Medium","High"))
df$Treatment<-ordered(df$Treatment,levels=c("Low","Medium","High","Control"))
df<-egg.2018.wk3.4
#wanted to see if raw data showed anything different...not really
#I don't think it's worth it to go down that rabbit hole, I think it's best to
#--just keep it as total reproduction per reef per week, or just total reproduction over
#--the course of a month
df<-gob.sub

#log-transform counts for new data (2019.4.16)
#thought this would be better for t6 counts
df$log<-log(df$Egg.count+1)#...no

##########A) Unadjusted egg counts by treatment by week#############
#2017.t1.2.3 analyses
mod1<-lm(Egg.count~Treatment*Density,data=df)#not normal
#coef(summary(mod1))[,4]
hist(resid(mod1))
qqnorm(resid(mod1))
qqline(resid(mod1))
anova(mod1)
Anova(mod1, type="II")
summary(mod1)
plot(mod1)#looks like equal variance though
boxplot(Egg.count~Treatment*Density, data=df)
View(df)

library(HH)
ancova(Egg.count ~ Treatment + Density, data=df)

#with interactive effects of covariate
ancova(Egg.count ~ Treatment * Density, data=df)
       
mod2<-lm(Egg.count~Treatment*Deployment.day, data=df)
hist(resid(mod2))#more normal
qqnorm(resid(mod2))
qqline(resid(mod2))
anova(mod2)
plot(mod2)#again, looks like equal variance

mod3<-lm(Egg.count~Treatment*deployment.day*Trial, data=df)
hist(resid(mod3))#more normal
qqnorm(resid(mod3))
qqline(resid(mod3))
anova(mod3)
plot(mod3)

mod4<-lmer(Egg.count~Treatment*deployment.day+(1|Trial),data=df)
hist(resid(mod4))#more normal
qqnorm(resid(mod4))
qqline(resid(mod4))
anova(mod4,type = "II") #shows same thing here as before
Anova(mod4) #go with this output? Type II Wald chi squared test
plot(mod4)
ranef(mod4)

mod5<-lm(Egg.count~Treatment*deployment.day*Density*Trial,data=df)
hist(resid(mod5))
qqnorm(resid(mod5))
qqline(resid(mod5))
anova(mod5)
summary(mod5)
plot(mod5)#doesn't look so great...

mod6<-lm(Egg.count~Treatment*Density*Week,data=df)
hist(resid(mod6))
qqnorm(resid(mod6))
qqline(resid(mod6))
anova(mod6)
summary(mod6)
plot(mod6)

mod7<-lm(Egg.count~Treatment*Density*Trial,data=df)
hist(resid(mod7))
qqnorm(resid(mod7))
qqline(resid(mod7))
anova(mod7)
summary(mod7)
plot(mod7)

mod8<-lm(Egg.count~Treatment*Trial+Density,data=df)#bad
hist(resid(mod8))
qqnorm(resid(mod8))
qqline(resid(mod8))
anova(mod8)
summary(mod8)
plot(mod8)

mod9<-lm(Egg.count~Treatment*Trial*Deployment.day+Density,data=df)
hist(resid(mod9))
qqnorm(resid(mod9))
qqline(resid(mod9))
anova(mod9)
summary(mod9)
plot(mod9)

mod10<-lm(Egg.count~Treatment*Density*Deployment.day,data=df)
hist(resid(mod10))
qqnorm(resid(mod10))
qqline(resid(mod10))
anova(mod10)
summary(mod10)
plot(mod10)

#this is what I think the model should be
mod11<-lmer(Egg.count~Treatment*Density+(1|Trial)+(1|deployment.day),data=df)
hist(resid(mod11))
qqnorm(resid(mod11))
qqline(resid(mod11))
anova(mod11)
anova(mod11,type='II')
Anova(mod11)
summary(mod11)
plot(mod11)
ranef(mod11)

mod11a<-lmer(Egg.count~Treatment*Density+(1|Trial)+deployment.day,data=df)
hist(resid(mod11a))
qqnorm(resid(mod11a))
qqline(resid(mod11a))
anova(mod11a)
anova(mod11a,type='I')
Anova(mod11a)
summary(mod11a)
plot(mod11a)
ranef(mod11a)

#now adding trial as random effect and nesting dep. day
#--I think it has to be that way, because I don't want all day 1's being treated the same
#--among trials
mod11ai<-lmer(Egg.count~Treatment*Density+(1|Trial)+(1|Deployment.day),data=df)
hist(resid(mod11ai))
qqnorm(resid(mod11ai))
qqline(resid(mod11ai))
anova(mod11ai)
anova(mod11ai,type='I')
Anova(mod11ai)
summary(mod11ai)
plot(mod11ai)
ranef(mod11ai)
View(egg.2017.t1.2.3)

mod11aii<-lmer(Egg.count~Treatment+Density+(1|Trial)+(1|Trial:Deployment.day),data=df)
hist(resid(mod11aii))
qqnorm(resid(mod11aii))
qqline(resid(mod11aii))
anova(mod11aii)
anova(mod11aii,type='II')
Anova(mod11aii)
summary(mod11aii)
plot(mod11aii)
ranef(mod11aii)

#strange thing, but I can't get ranef output unless I switch order of factors in 
#--model
mod11aii<-lmer(Egg.count~Density+Treatment+(1|Trial)+(1|Trial:Deployment.day),data=df)
hist(resid(mod11aii))
qqnorm(resid(mod11aii))
qqline(resid(mod11aii))
anova(mod11aii)
anova(mod11aii,type='II')
Anova(mod11aii)
summary(mod11aii)
plot(mod11aii)
ranef(mod11aii)

AIC(mod11a,mod11ai,mod11aii)

#poisson dist? Doesn't look better than normal dist, although data don't look super normal
mod11aiii<-glmer(Egg.count~Treatment+Density+(1|Trial)+(1|Deployment.day),family=poisson,data=df)
hist(resid(mod11aiii))
qqnorm(resid(mod11aiii))
qqline(resid(mod11aiii))
Anova(mod11aiii)
summary(mod11aiii)
plot(mod11aiii)
ranef(mod11aiii)

#now going to work with 2018, t4 and t5 data nesting dep. day
mod11aiv<-lmer(Egg.count~Treatment*Density*Week+(1|Trial)+(1|Trial:Deployment.day),data=df)
hist(resid(mod11aiv))
qqnorm(resid(mod11aiv))
qqline(resid(mod11aiv))
Anova(mod11aiv)
anova(mod11aiv,type='I')
Anova(mod11aiv)
summary(mod11aiv)
plot(mod11aiv)
ranef(mod11aiv)

#reordering factors so I can get random effects
mod11aiv<-lmer(Egg.count~Treatment*Week*Density+(1|Trial)+(1|Trial:Deployment.day),data=df)
hist(resid(mod11aiv))
qqnorm(resid(mod11aiv))
qqline(resid(mod11aiv))
Anova(mod11aiv)
anova(mod11aiv,type='I')
Anova(mod11aiv)
summary(mod11aiv)
plot(mod11aiv)
ranef(mod11aiv)

#looking at 2018 t6 only, no factor of deployment day or trial (only 1 for both)
mod11av<-lm(Egg.count~Treatment*Week*Density,data=df)
hist(resid(mod11av))
qqnorm(resid(mod11av))
qqline(resid(mod11av))
#Anova(mod11av)
anova(mod11av)
summary(mod11av)
plot(mod11av)
#ranef(mod11av)

mod11avx<-lm(Egg.count~Treatment+Week*Density,data=df)
hist(resid(mod11avx))
qqnorm(resid(mod11avx))
qqline(resid(mod11avx))
anova(mod11avx)
summary(mod11avx)
plot(mod11avx)
#ranef(mod11av)

#taking out interactions, becasue none were significant
mod11avi<-lm(Egg.count~Treatment+Week+Density,data=df)
hist(resid(mod11avi))
qqnorm(resid(mod11avi))
qqline(resid(mod11avi))
#Anova(mod11av)
anova(mod11avi)
summary(mod11avi)
plot(mod11avi)
#ranef(mod11av)

mod11b<-lm(Egg.count~Treatment*Density+Trial+Deployment.day,data=df)
hist(resid(mod11b))
qqnorm(resid(mod11b))
qqline(resid(mod11b))
anova(mod11b)
anova(mod11b,type='II')
Anova(mod11b)
summary(mod11b)
plot(mod11b)
ranef(mod11b)

mod12<-lmer(Egg.count~Treatment*Density+(1|Trial),data=df)
hist(resid(mod12))
qqnorm(resid(mod12))
qqline(resid(mod12))
anova(mod12)
Anova(mod12)
summary(mod12)
plot(mod12)
ranef(mod12)

mod13<-lm(Egg.count~Treatment+Density+Trial,data=df)
hist(resid(mod13))
qqnorm(resid(mod13))
qqline(resid(mod13))
anova(mod13)
Anova(mod13)
summary(mod13)
plot(mod13)
ranef(mod13)

AIC(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod11a,mod11b,mod12,mod13)

#2018 trials 4 and 5 egg counts ###
#now that we have our df set up, we can go through model selection
mod1<-lm(Egg.count~Treatment, data=df)#not normal
hist(resid(mod1))
qqnorm(resid(mod1))
qqline(resid(mod1))
anova(mod1)
#log egg counts
mod1log<-lm(ln.ec~Treatment, data=df)#worse than non-logged counts!
hist(resid(mod1log))
anova(mod1log)

mod2<-lm(Egg.count~Treatment*Week, data=df)#more normal, but still not great
hist(resid(mod2))
qqnorm(resid(mod2))
qqline(resid(mod2))
anova(mod2)
summary(mod2)
#log
mod2log<-lm(ln.ec~Treatment*Week, data=df)#more normal, but still not great
hist(resid(mod2log))
qqnorm(resid(mod2log))
qqline(resid(mod2log))
anova(mod2log)
summary(mod2log)
#not going to waste any more time with log counts...looks terrible

#pulling density out of the other model, because it was sig.
mod2a<-lm(Egg.count~Treatment*Density, data=df) 
hist(resid(mod2a))
qqnorm(resid(mod2a))
qqline(resid(mod2a))
anova(mod2a)
summary(mod2a)
plot(mod2a)

mod3<-lm(Egg.count~Treatment*Week*Density, data=df)#better fit, 
hist(resid(mod3))
qqnorm(resid(mod3))
qqline(resid(mod3))
anova(mod3)#showing interaction between treatment and week
#-- means that the treatment effect on reproductive output varied by week
Anova(mod3,type="II")
summary(mod3)
plot(mod3)#Ok fit and normailty, but definitely a shape to the variances
#don't really want to start pulling data points as outliers
#this could be the best one for 2018, trials 4 and 5

mod4<-lmer(Egg.count~Treatment*Week*Density+(1|Trial), data=df)
hist(resid(mod4))
qqnorm(resid(mod4))
qqline(resid(mod4))       
plot(mod4) #best fit so far, variances still look a bit funky because of some outlying points
anova(mod4)#density and week*density were sig.
#looks a little better here, might consider bringing in deployment day
summary(mod4)

#figuring out if density is better as a interactive factor or covariate
AIC(mod4,mod4a)#mod4 is better based on AIC value

mod4a<-lmer(Egg.count~Treatment*Week*Density+(1|Trial)+(1|Deployment.day), data=df)
hist(resid(mod4a))
qqnorm(resid(mod4a))
qqline(resid(mod4a))
anova(mod4a)
summary(mod4a)
plot(mod4a)

mod4b<-lmer(Egg.count~Treatment*Week*Density+(1|Trial)+Deployment.day, data=df)
hist(resid(mod4b))
qqnorm(resid(mod4b))
qqline(resid(mod4b))
anova(mod4b)
summary(mod4b)
plot(mod4b)#slightly better than model than mod4a

(AIC(mod4a,mod4b))

mod4c<-lm(Egg.count~Treatment*Week*Density+Deployment.day, data=df)
hist(resid(mod4c))
qqnorm(resid(mod4c))
qqline(resid(mod4c))
anova(mod4c)
summary(mod4c)
plot(mod4c)#better than model that doesn't include deployment day, 
AIC(mod3,mod4c)

mod4d<-lm(Egg.count~Treatment*Week*Density+Deployment.day+Trial, data=df)
hist(resid(mod4d))
qqnorm(resid(mod4d))
qqline(resid(mod4d))
anova(mod4d)
summary(mod4d)
plot(mod4d)#better than model that doesn't include deployment day
AIC(mod3,mod4c,mod4d)

mod4e<-lmer(Egg.count~Treatment*Week*Density+Trial+(1|Deployment.day), data=df)
hist(resid(mod4e))
qqnorm(resid(mod4e))
qqline(resid(mod4e))
anova(mod4e)
summary(mod4e)
plot(mod4e)#looks like reefs deployed on d2 laid fewer eggs, based on all other factors in the model
AIC(mod3,mod4c,mod4d,mod4e)#best model so far in terms of AIC value
ranef(mod4e)

#dropping trial from the model
mod4f<-lmer(Egg.count~Treatment*Week*Density+(1|Deployment.day), data=df)
hist(resid(mod4f))
qqnorm(resid(mod4f))
qqline(resid(mod4f))
anova(mod4f)
summary(mod4f)
plot(mod4f)
ranef(mod4f)
AIC(mod3,mod4c,mod4d,mod4e,mod4f)#dropping trial didn't improve AIC value
#--wondering what to include in the model here, it seems like there weren't any differences just based
#--on trial (t4 and 5), so why include it in the model? AIC value is pretty close, and I think it increases power
#--to not include unnecessary factors

#just 2018 trial 6 data, deployed all fish on the same day, deployment day not a factor
mod1<-lm(Egg.count~Treatment,data=df)
hist(resid(mod1))
qqnorm(resid(mod1))
qqline(resid(mod1))
anova(mod1)
plot(mod1) #not super equal variance, much smaller in the high treatment

mod2<-lm(Egg.count~Treatment*Density,data=df)
hist(resid(mod2))
qqnorm(resid(mod2))
qqline(resid(mod2))
anova(mod2)
plot(mod2)

mod3<-lm(Egg.count~Treatment*Density*Week,data=df)
hist(resid(mod3))
qqnorm(resid(mod3))
qqline(resid(mod3))
anova(mod3)
plot(mod3)

mod3a<-lm(Egg.count~Treatment*Week+Density,data=df)
hist(resid(mod3a))
qqnorm(resid(mod3a))
qqline(resid(mod3a))
anova(mod3a)
plot(mod3a)#nest one based on AIC value, to include density as a covariate

mod4<-lm(Egg.count~Treatment*Density*Week+Position,data=df)
hist(resid(mod4))
qqnorm(resid(mod4))
qqline(resid(mod4))
anova(mod4)
plot(mod4)

mod5<-lmer(Egg.count~Treatment*Density*Week+(1|Position),data=df)
hist(resid(mod5))
qqnorm(resid(mod5))
qqline(resid(mod5))
anova(mod5)
plot(mod5)
ranef(mod5)#seems like lower reproduction in the front of the array

AIC(mod1,mod2,mod3,mod3a,mod4)#mod 5 was calculated differently, so not a great comparison to include

#including week as a random effect

mod4e<-lmer(Egg.count~Treatment*Week*Density+Trial+(1|Deployment.day), data=df)
hist(resid(mod4e))
qqnorm(resid(mod4e))
qqline(resid(mod4e))
anova(mod4e)
summary(mod4e)



#bringing year into the mix
mod2a<-lm(Egg.count~Treatment*Week*Density*Year, data=df)#better fit, 
hist(resid(mod2a))
qqnorm(resid(mod2a))
qqline(resid(mod2a))
anova(mod2a)
summary(mod2a)
plot(mod2a)

#just trial 6, no year
mod2a<-lm(Egg.count~Treatment*Week*Density*Final.biomass.total, data=df)#better fit, 
hist(resid(mod2a))
qqnorm(resid(mod2a))
qqline(resid(mod2a))
anova(mod2a)
summary(mod2a)
plot(mod2a)

#t6 with biomass data (didn't set up model)

mod2ai<-lm(Egg.count~Treatment*Week*Density+(1|Year)+(1|Trial), data=df)#better fit, 
hist(resid(mod2ai))
qqnorm(resid(mod2ai))
qqline(resid(mod2ai))
anova(mod2ai)
summary(mod2ai)
plot(mod2ai)
#fixef(mod2ai)
#ranef(mod2ai)

#year and GLMM, poisson dist, keep getting error message
mod2aii<-glmer(Egg.count~Treatment*Week*Density+(1|Year)+(1|Trial),family=poisson, data=df)#better fit, 
hist(resid(mod2aii))
qqnorm(resid(mod2aii))
qqline(resid(mod2aii))
anova(mod2aii)
summary(mod2aii)
plot(mod2aii)
fixef(mod2aii)
ranef(mod2aii)

mod4<-lmer(Egg.count~Treatment*Week*Density+(1|Trial), data=df)
hist(resid(mod4))
qqnorm(resid(mod4))
qqline(resid(mod4))       
plot(mod4) #best fit so far, variances still look a bit funky because of some outlying points
anova(mod4)#density and week*density were sig.

mod5<-glmer(Egg.count~Treatment*Week*Density+(1|Trial),family=poisson,data=df)
hist(resid(mod5))
qqnorm(resid(mod5))
qqline(resid(mod5))       
plot(mod5)
Anova(mod5, type="II")
summary(mod5)
fixef(mod5)
ranef(mod5)
Anova(mod5)
#best fit so far, but it's tough to determine whether the variances are equal
#doing the plot function shows a clutster of points
#not sure how to analyze this, because when I try to run Anova (type=etc.) it doesn't seem to work

#reproduction per avg.total.biomass#####
mod2<-lm(egg.per.tot.avg.biomass~Treatment*Week, data=df)#more normal, but still not great
hist(resid(mod2))
qqnorm(resid(mod2))
qqline(resid(mod2))
anova(mod2)
summary(mod2)

mod2b<-lm(egg.per.tot.avg.biomass~Treatment*Week*Density, data=df)#more normal, but still not great
hist(resid(mod2b))
qqnorm(resid(mod2b))
qqline(resid(mod2b))
anova(mod2b)
summary(mod2b)

mod2aii<-glmer(egg.per.tot.avg.biomass~Treatment*Week*Density+(1|Year)+(1|Trial),family=poisson, data=df)#better fit, 
hist(resid(mod2aii))
qqnorm(resid(mod2aii))
qqline(resid(mod2aii))
anova(mod2aii)
summary(mod2aii)
plot(mod2aii)
fixef(mod2aii)
ranef(mod2aii)

#########B: Per capita output by week###############
#now going to look at corrected metrics
mod1<-lm(per.capita.repro~Treatment, data=df)#not normal
hist(resid(mod1))
anova(mod1)
#log egg counts
mod1log<-lm(ln.pc~Treatment, data=df)#worse than non-logged counts!
hist(resid(mod1log))
anova(mod1log)

mod2<-lm(per.capita.repro~Treatment*Week, data=df)#more normal, but still not great
hist(resid(mod2))
qqnorm(resid(mod2))
qqline(resid(mod2))
anova(mod2) #indicates that week had an effect, which makes sense, lost fish over time, fewer eggs
summary(mod2)
#log
mod2log<-lm(ln.pc~Treatment*Week, data=df)#bad fit, likely becuase of all of the zeros
hist(resid(mod2log))
qqnorm(resid(mod2log))
qqline(resid(mod2log))
anova(mod2log)
summary(mod2log)
#not going to waste any more time with log counts...looks terrible

mod3<-lm(per.capita.repro~Treatment*Week*Density, data=df)#better fit, 
hist(resid(mod3))
qqnorm(resid(mod3))
qqline(resid(mod3))
anova(mod3) #treatment not significant here, only week, density barely nonsignificant
summary(mod3)
plot(mod3)#Ok fit and normailty, but ddefinitely a shape to the variances
#don't really want to start pulling data points as outliers

mod4<-lmer(per.capita.repro~Treatment*Week*Density+(1|Trial)+(1|Trial:deployment.day), data=df)
hist(resid(mod4))
qqnorm(resid(mod4))
qqline(resid(mod4))       
plot(mod4) #best fit so far, variances still look a bit funky because of some outlying points
anova(mod4)#nothing significant

mod5<-glmer(per.capita.repro~Treatment*Week*Density+(1|Trial),family=poisson,data=df)
hist(resid(mod5))
qqnorm(resid(mod5))
qqline(resid(mod5))       
plot(mod5)
Anova(mod5, type="II")
summary(mod5)
fixef(mod5)
ranef(mod5)
Anova(mod5)
#best fit so far, but it's tough to determine whether the variances are equal
#doing the plot function shows a clutster of points
#not sure how to analyze this, because when I try to run Anova (type=etc.) it doesn't seem to work
#not sure what to do for the Lmertest

####takeaway: it seems like the raw counts are a better fit, with the best fit using the raw counts
#and including density in the model as a covariate, and running it with a Poisson distribution

#X. Data visualization for results with Mark, raw counts and per capita ########
#df titles:
#a. den.avg.egg
#b. egg.2017
#c. egg.2018

df$Treatment<-ordered(df$Treatment,c("Low","Medium","High","Control"))

#x. looking at total output per year
bargraph.CI(x.factor = Year, response = Egg.count, group= Treatment, legend=TRUE, main="Egg counts by year and treatment", xlab="Year", ylab="Eggs per reef (mean +/- se)", x.leg=9, yleg=4500, data = egg.den.bio)
#y. total output by week and treatment with just t6 data
bargraph.CI(x.factor = Week, response = Egg.count, group= Treatment, legend=TRUE, main="Egg counts by week and treatment, trial 6 only", xlab="Week", ylab="Eggs per reef (mean +/- se)", x.leg=9, yleg=4500, data = egg.t6)
#z.total output by week and treatment using avg. output per total avg.biomass throughout the trial
bargraph.CI(x.factor = Week, response = egg.per.tot.avg.biomass, group= Treatment, legend=TRUE, main="Output per avg. total biomass by week and treatment, 2017 + 2018", xlab="Week", ylab="Eggs per avg. total biomass (mean +/- se)", x.leg=9, yleg=4500, data = df)


#a. combined by week total raw counts
bargraph.CI(x.factor = Week, response = Egg.count, group= Treatment, legend=TRUE, main="2017 + 2018 total counts", xlab="Week", ylab="Eggs per reef (mean +/- se)", x.leg=13, yleg=4500, data = den.avg.egg)
#a.i eggs produced over time, regardless of treatment
bargraph.CI(x.factor = Week, response = Egg.count, legend=TRUE, main="2017 + 2018 total counts over time", xlab="Week", ylab="Eggs per reef (mean +/- se)", x.leg=13, yleg=4500, data = den.avg.egg)
#a.ii eggs produced per treatment, regardless of time
bargraph.CI(x.factor = Treatment, response = Egg.count, main="2017 total counts per treatment", xlab="Treatment", ylab="Eggs per reef (mean +/- se)", x.leg=13, yleg=4500, data = df)

#b. 2017 only by week total raw counts
#egg counts by treatment
bargraph.CI(x.factor = Week, response = Egg.count,group=Treatment, legend=TRUE, main="Trial 6 egg counts", xlab="Week", ylab="Eggs per reef (mean +/- se)", data = df)
#egg counts by treatment and deployment day
bargraph.CI(x.factor = Deployment.day, response = Egg.count, group= Treatment, legend=TRUE, main="2017 total counts", xlab="Deployment day", ylab="Eggs per reef (mean +/- se)", data = egg.2017.t1.2.3)
#looking at eggs by deployment day
bargraph.CI(x.factor = Deployment.day, response = Egg.count, main="2017 total counts", xlab="Deployment day", ylab="Eggs per reef (mean +/- se)", x.leg=13, yleg=4500, data = egg.2017.t1.2.3)
#there's no need to do the other figures because there was only one week of sampling
#b.i eggs produced over time, regardless of treatment
#bargraph.CI(x.factor = Week, response = Egg.count, legend=TRUE, main="2017 + 2018 total counts over time", xlab="Week", ylab="Eggs per reef (mean +/- se)", x.leg=13, yleg=4500, data = den.avg.egg)
#b.ii eggs produced per treatment, regardless of time
#bargraph.CI(x.factor = Treatment, response = Egg.count, legend=TRUE, main="2017 + 2018 total counts per treatment", xlab="Treatment", ylab="Eggs per reef (mean +/- se)", x.leg=13, yleg=4500, data = den.avg.egg)

#density by treatment
bargraph.CI(x.factor = Treatment, response = Density, legend=TRUE, main="2017 density counts", xlab="Treatment", ylab="Fish per reef (mean +/- se)", data = egg.2017.t1.2.3)
#--maybe a difference between high and low, but not a big difference between medium and low
bargraph.CI(x.factor = Trial, response = Egg.count, legend=TRUE, main="2017 eggs by trial", xlab="Trial", ylab="Eggs per reef (mean +/- se)", data = egg.2017.t1.2.3)
#--doesn't look like there are any differences in output among trials in 2017

#c. 2018 only by week total raw counts
bargraph.CI(x.factor = Week, response = Egg.count, group= Treatment, legend=TRUE, main="2018 total counts per treatment over time", xlab="Week", ylab="Eggs per reef (mean +/- se)", x.leg=13, yleg=4500, data = df)
#a.i eggs produced over time, regardless of treatment
bargraph.CI(x.factor = Week, response = Egg.count, legend=TRUE, main="2018 total counts over time", xlab="Week", ylab="Eggs per reef (mean +/- se)", x.leg=13, yleg=4500, data = egg.2018)
#a.ii eggs produced per treatment, regardless of time
bargraph.CI(x.factor = Treatment, response = Egg.count, legend=TRUE, main="2018 t6 total counts per treatment", xlab="Treatment", ylab="Eggs per reef (mean +/- se)", x.leg=13, yleg=4500, data = df)

#now going to do per capita output for each of the previous figures####
#want to see how densities may have altered reproductive output
#divided total counts by average weekly densities (average of with an average = no no)
#really only have to change the response variable
#2a. combined by week total raw counts (only 1 week of sampling)
bargraph.CI(x.factor = Week, response = per.capita.repro, group= Treatment, legend=TRUE, main="2018, t4.5 per capita reproduction by risk over time", xlab="Week", ylab="Eggs per fish per reef (mean +/- se)", x.leg=10, yleg=4500, data = df)
#2a.i eggs produced over time, regardless of treatment
bargraph.CI(x.factor = Week, response = per.capita.repro, legend=TRUE, main="2017 + 2018 per capita counts over time", xlab="Week", ylab="Eggs per fish per reef (mean +/- se)", x.leg=13, yleg=4500, data = den.avg.egg)
#2a.ii eggs produced per treatment, regardless of time
bargraph.CI(x.factor = Treatment, response = per.capita.repro, legend=TRUE, main="2017 + 2018 per capita counts per treatment", xlab="Treatment", ylab="Eggs per fish per reef (mean +/- se)", x.leg=13, yleg=4500, data = den.avg.egg)

#means by treatment and time
per.capita.t4.5<-df %>%
  group_by(Treatment) %>%
  summarize(per.capita = ceiling(mean(per.capita.repro)))
per.capita.t4.5

#2b. 2017 only by week total raw counts (only 1 week of sampling)
bargraph.CI(x.factor = Week, response = per.capita.repro, group= Treatment, legend=TRUE, main="2017 per capita counts", xlab="Week", ylab="Eggs per fish per reef (mean +/- se)", x.leg=13, yleg=4500, data = egg.2017)
#there's no need to do the other figures because there was only one week of sampling

#2c. 2018 only by week total raw counts
bargraph.CI(x.factor = Week, response = per.capita.repro, group= Treatment, legend=TRUE, main="2018 per capita counts per treatment over time", xlab="Week", ylab="Eggs per fish per reef (mean +/- se)", x.leg=13, yleg=4500, data = egg.2018)
#2c.i eggs produced over time, regardless of treatment
bargraph.CI(x.factor = Week, response = per.capita.repro, legend=TRUE, main="2018 per capita counts over time", xlab="Week", ylab="Eggs per fish per reef (mean +/- se)", x.leg=13, yleg=4500, data = egg.2018)
#2c.ii eggs produced per treatment, regardless of time
bargraph.CI(x.factor = Treatment, response = per.capita.repro, legend=TRUE, main="2018 per capita counts per treatment", xlab="Treatment", ylab="Eggs per fish per reef (mean +/- se)", x.leg=13, yleg=4500, data = egg.2018)


##not sure why, but "Deployment.day" worked before, and now it has to be coded as "deployment.day"####

df<-egg.2017.t1.2.3#2017 only
#df$Trial<-as.factor(df$Trial)
df<-egg.2018.t4.5 #2018 only
df<-egg.2018.t6

#2017.t1.2.3 analyses
mod1<-lm(Egg.count~Treatment, data=df)#not normal
hist(resid(mod1))
qqnorm(resid(mod1))
qqline(resid(mod1))
anova(mod1)
plot(mod1)#looks like equal variance though

mod2<-lm(Egg.count~Treatment*deployment.day, data=df)
hist(resid(mod2))#more normal
qqnorm(resid(mod2))
qqline(resid(mod2))
anova(mod2)
plot(mod2)#again, looks like equal variance

mod3<-lm(Egg.count~Treatment*deployment.day*Trial, data=df)
hist(resid(mod3))#more normal
qqnorm(resid(mod3))
qqline(resid(mod3))
anova(mod3)
plot(mod3)

mod4<-lmer(Egg.count~Treatment*deployment.day+(1|Trial),data=df)
hist(resid(mod4))#more normal
qqnorm(resid(mod4))
qqline(resid(mod4))
anova(mod4,type = "II") #shows same thing here as before
Anova(mod4) #go with this output? Type II Wald chi squared test
plot(mod4)
ranef(mod4)

mod5<-lm(Egg.count~Treatment*deployment.day*Density*Trial,data=df)
hist(resid(mod5))
qqnorm(resid(mod5))
qqline(resid(mod5))
anova(mod5)
summary(mod5)
plot(mod5)#doesn't look so great...

mod6<-lm(Egg.count~Treatment*Density,data=df)
hist(resid(mod6))
qqnorm(resid(mod6))
qqline(resid(mod6))
anova(mod6)
summary(mod6)
plot(mod6)

mod7<-lm(Egg.count~Treatment*Density*Trial,data=df)
hist(resid(mod7))
qqnorm(resid(mod7))
qqline(resid(mod7))
anova(mod7)
summary(mod7)
plot(mod7)

mod8<-lm(Egg.count~Treatment*Trial+Density,data=df)#bad
hist(resid(mod8))
qqnorm(resid(mod8))
qqline(resid(mod8))
anova(mod8)
summary(mod8)
plot(mod8)

mod9<-lm(Egg.count~Treatment*Trial*deployment.day+Density,data=df)
hist(resid(mod9))
qqnorm(resid(mod9))
qqline(resid(mod9))
anova(mod9)
summary(mod9)
plot(mod9)

mod10<-lm(Egg.count~Treatment*Density,data=df)
hist(resid(mod10))
qqnorm(resid(mod10))
qqline(resid(mod10))
anova(mod10)
summary(mod10)
plot(mod10)
boxplot(Egg.count~Treatment,data=df)

mod10a<-lmer(Egg.count~Treatment*Density+(1|deployment.day),data=df)
hist(resid(mod10a))
qqnorm(resid(mod10a))
qqline(resid(mod10a))
anova(mod10a)
summary(mod10a)
plot(mod10a)

#this is what I think the model should be
mod11<-lmer(Egg.count~Treatment*Density+(1|Trial)+(1|deployment.day),data=df)
hist(resid(mod11))
qqnorm(resid(mod11))
qqline(resid(mod11))
anova(mod11)
anova(mod11,type='II')
Anova(mod11)
summary(mod11)
plot(mod11)
ranef(mod11)

mod11a<-lmer(Egg.count~Treatment*Density+(1|Trial)+deployment.day,data=df)
hist(resid(mod11a))
qqnorm(resid(mod11a))
qqline(resid(mod11a))
anova(mod11a)
anova(mod11a,type='I')
Anova(mod11a)
summary(mod11a)
plot(mod11a)
ranef(mod11a)

#now adding trial as random effect and nesting dep. day
#--I think it has to be that way, because I don't want all day 1's being treated the same
#--among trials
mod11ai<-lmer(Egg.count~Treatment*Density+(1|Trial)+(1|deployment.day),data=df)
hist(resid(mod11ai))
qqnorm(resid(mod11ai))
qqline(resid(mod11ai))
anova(mod11ai)
anova(mod11ai,type='I')
Anova(mod11ai)
summary(mod11ai)
plot(mod11ai)
ranef(mod11ai)
View(egg.2017.t1.2.3)

mod11aii<-lmer(Egg.count~Treatment+Density+(1|Trial)+(1|Trial:deployment.day),data=df)
hist(resid(mod11aii))
qqnorm(resid(mod11aii))
qqline(resid(mod11aii))
anova(mod11aii)
anova(mod11aii,type='II')
Anova(mod11aii)
summary(mod11aii)
plot(mod11aii)
ranef(mod11aii)

#strange thing, but I can't get ranef output unless I switch order of factors in 
#--model
mod11aii<-lmer(Egg.count~Density+Treatment+(1|Trial)+(1|Trial:deployment.day),data=df)
hist(resid(mod11aii))
qqnorm(resid(mod11aii))
qqline(resid(mod11aii))
anova(mod11aii)
anova(mod11aii,type='II')
Anova(mod11aii)
summary(mod11aii)
plot(mod11aii)
ranef(mod11aii)

AIC(mod11a,mod11ai,mod11aii)

#poisson dist? Doesn't look better than normal dist, although data don't look super normal
mod11aiii<-glmer(Egg.count~Treatment+Density+(1|Trial)+(1|deployment.day),family=poisson,data=df)
hist(resid(mod11aiii))
qqnorm(resid(mod11aiii))
qqline(resid(mod11aiii))
Anova(mod11aiii)
summary(mod11aiii)
plot(mod11aiii)
ranef(mod11aiii)

#now going to work with 2018, t4 and t5 data nesting dep. day
mod11aiv<-lmer(Egg.count~Treatment*Density*Week+(1|Trial)+(1|Trial:deployment.day),data=df)
hist(resid(mod11aiv))
qqnorm(resid(mod11aiv))
qqline(resid(mod11aiv))
Anova(mod11aiv)
anova(mod11aiv,type='I')
Anova(mod11aiv)
summary(mod11aiv)
plot(mod11aiv)
ranef(mod11aiv)

#reordering factors so I can get random effects
mod11aiv<-lmer(Egg.count~Treatment*Week*Density+(1|Trial)+(1|Trial:deployment.day),data=df)
hist(resid(mod11aiv))
qqnorm(resid(mod11aiv))
qqline(resid(mod11aiv))
Anova(mod11aiv)
anova(mod11aiv,type='I')
Anova(mod11aiv)
summary(mod11aiv)
plot(mod11aiv)
ranef(mod11aiv)

#looking at 2018 t6 only, no factor of deployment day or trial (only 1 for both)
mod11av<-lm(Egg.count~Treatment*Week*Density,data=df)
hist(resid(mod11av))
qqnorm(resid(mod11av))
qqline(resid(mod11av))
#Anova(mod11av)
anova(mod11av)
summary(mod11av)
plot(mod11av)
#ranef(mod11av)

mod11avx<-lm(Egg.count~Treatment+Week*Density,data=df)
hist(resid(mod11avx))
qqnorm(resid(mod11avx))
qqline(resid(mod11avx))
anova(mod11avx)
summary(mod11avx)
plot(mod11avx)
#ranef(mod11av)

#taking out interactions, becasue none were significant
mod11avi<-lm(Egg.count~Treatment+Week+Density,data=df)
hist(resid(mod11avi))
qqnorm(resid(mod11avi))
qqline(resid(mod11avi))
#Anova(mod11av)
anova(mod11avi)
summary(mod11avi)
plot(mod11avi)
#ranef(mod11av)

mod11b<-lm(Egg.count~Treatment*Density+Trial+deployment.day,data=df)
hist(resid(mod11b))
qqnorm(resid(mod11b))
qqline(resid(mod11b))
anova(mod11b)
anova(mod11b,type='II')
Anova(mod11b)
summary(mod11b)
plot(mod11b)
ranef(mod11b)

mod12<-lmer(Egg.count~Treatment*Density+(1|Trial),data=df)
hist(resid(mod12))
qqnorm(resid(mod12))
qqline(resid(mod12))
anova(mod12)
Anova(mod12)
summary(mod12)
plot(mod12)
ranef(mod12)

mod13<-lm(Egg.count~Treatment+Density+Trial,data=df)
hist(resid(mod13))
qqnorm(resid(mod13))
qqline(resid(mod13))
anova(mod13)
Anova(mod13)
summary(mod13)
plot(mod13)
ranef(mod13)

AIC(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod11a,mod11b,mod12,mod13)

#2018 trials 4 and 5 egg counts ###
#now that we have our df set up, we can go through model selection
mod1<-lm(Egg.count~Treatment, data=df)#not normal
hist(resid(mod1))
qqnorm(resid(mod1))
qqline(resid(mod1))
anova(mod1)
#log egg counts
mod1log<-lm(ln.ec~Treatment, data=df)#worse than non-logged counts!
hist(resid(mod1log))
anova(mod1log)

mod2<-lm(Egg.count~Treatment*Week, data=df)#more normal, but still not great
hist(resid(mod2))
qqnorm(resid(mod2))
qqline(resid(mod2))
anova(mod2)
summary(mod2)
#log
mod2log<-lm(ln.ec~Treatment*Week, data=df)#more normal, but still not great
hist(resid(mod2log))
qqnorm(resid(mod2log))
qqline(resid(mod2log))
anova(mod2log)
summary(mod2log)
#not going to waste any more time with log counts...looks terrible

#pulling density out of the other model, because it was sig.
mod2a<-lm(Egg.count~Treatment*Density, data=df) 
hist(resid(mod2a))
qqnorm(resid(mod2a))
qqline(resid(mod2a))
anova(mod2a)
summary(mod2a)
plot(mod2a)

mod3<-lm(Egg.count~Treatment*Week*Density, data=df)#better fit, 
hist(resid(mod3))
qqnorm(resid(mod3))
qqline(resid(mod3))
anova(mod3)#showing interaction between treatment and week
#-- means that the treatment effect on reproductive output varied by week
Anova(mod3,type="II")
summary(mod3)
plot(mod3)#Ok fit and normailty, but definitely a shape to the variances
#don't really want to start pulling data points as outliers
#this could be the best one for 2018, trials 4 and 5

mod4<-lmer(Egg.count~Treatment*Week*Density+(1|Trial), data=df)
hist(resid(mod4))
qqnorm(resid(mod4))
qqline(resid(mod4))       
plot(mod4) #best fit so far, variances still look a bit funky because of some outlying points
anova(mod4)#density and week*density were sig.
#looks a little better here, might consider bringing in deployment day
summary(mod4)

#figuring out if density is better as a interactive factor or covariate
AIC(mod4,mod4a)#mod4 is better based on AIC value

mod4a<-lmer(Egg.count~Treatment*Week*Density+(1|Trial)+(1|deployment.day), data=df)
hist(resid(mod4a))
qqnorm(resid(mod4a))
qqline(resid(mod4a))
anova(mod4a)
summary(mod4a)
plot(mod4a)

mod4b<-lmer(Egg.count~Treatment*Week*Density+(1|Trial)+deployment.day, data=df)
hist(resid(mod4b))
qqnorm(resid(mod4b))
qqline(resid(mod4b))
anova(mod4b)
summary(mod4b)
plot(mod4b)#slightly better than model than mod4a

(AIC(mod4a,mod4b))

mod4c<-lm(Egg.count~Treatment*Week*Density+deployment.day, data=df)
hist(resid(mod4c))
qqnorm(resid(mod4c))
qqline(resid(mod4c))
anova(mod4c)
summary(mod4c)
plot(mod4c)#better than model that doesn't include deployment day, 
AIC(mod3,mod4c)

mod4d<-lm(Egg.count~Treatment*Week*Density+deployment.day+Trial, data=df)
hist(resid(mod4d))
qqnorm(resid(mod4d))
qqline(resid(mod4d))
anova(mod4d)
summary(mod4d)
plot(mod4d)#better than model that doesn't include deployment day
AIC(mod3,mod4c,mod4d)

mod4e<-lmer(Egg.count~Treatment*Week*Density+Trial+(1|deployment.day), data=df)
hist(resid(mod4e))
qqnorm(resid(mod4e))
qqline(resid(mod4e))
anova(mod4e)
summary(mod4e)
plot(mod4e)#looks like reefs deployed on d2 laid fewer eggs, based on all other factors in the model
AIC(mod3,mod4c,mod4d,mod4e)#best model so far in terms of AIC value
ranef(mod4e)

#dropping trial from the model
mod4f<-lmer(Egg.count~Treatment*Week*Density+(1|deployment.day), data=df)
hist(resid(mod4f))
qqnorm(resid(mod4f))
qqline(resid(mod4f))
anova(mod4f)
summary(mod4f)
plot(mod4f)
ranef(mod4f)
AIC(mod3,mod4c,mod4d,mod4e,mod4f)#dropping trial didn't improve AIC value
#--wondering what to include in the model here, it seems like there weren't any differences just based
#--on trial (t4 and 5), so why include it in the model? AIC value is pretty close, and I think it increases power
#--to not include unnecessary factors

#just 2018 trial 6 data, deployed all fish on the same day, deployment day not a factor
mod1<-lm(Egg.count~Treatment,data=df)
hist(resid(mod1))
qqnorm(resid(mod1))
qqline(resid(mod1))
anova(mod1)
plot(mod1) #not super equal variance, much smaller in the high treatment

mod2<-lm(Egg.count~Treatment*Density,data=df)
hist(resid(mod2))
qqnorm(resid(mod2))
qqline(resid(mod2))
anova(mod2)
plot(mod2)

mod3<-lm(Egg.count~Treatment*Density*Week,data=df)
hist(resid(mod3))
qqnorm(resid(mod3))
qqline(resid(mod3))
anova(mod3)
plot(mod3)

mod3a<-lm(Egg.count~Treatment*Week+Density,data=df)
hist(resid(mod3a))
qqnorm(resid(mod3a))
qqline(resid(mod3a))
anova(mod3a)
plot(mod3a)#nest one based on AIC value, to include density as a covariate

mod4<-lm(Egg.count~Treatment*Density*Week+Position,data=df)
hist(resid(mod4))
qqnorm(resid(mod4))
qqline(resid(mod4))
anova(mod4)
plot(mod4)

mod5<-lmer(Egg.count~Treatment*Density*Week+(1|Position),data=df)
hist(resid(mod5))
qqnorm(resid(mod5))
qqline(resid(mod5))
anova(mod5)
plot(mod5)
ranef(mod5)#seems like lower reproduction in the front of the array

AIC(mod1,mod2,mod3,mod3a,mod4)#mod 5 was calculated differently, so not a great comparison to include

#including week as a random effect

mod4e<-lmer(Egg.count~Treatment*Week*Density+Trial+(1|deployment.day), data=df)
hist(resid(mod4e))
qqnorm(resid(mod4e))
qqline(resid(mod4e))
anova(mod4e)
summary(mod4e)



#bringing year into the mix
mod2a<-lm(Egg.count~Treatment*Week*Density*Year, data=df)#better fit, 
hist(resid(mod2a))
qqnorm(resid(mod2a))
qqline(resid(mod2a))
anova(mod2a)
summary(mod2a)
plot(mod2a)

#just trial 6, no year
mod2a<-lm(Egg.count~Treatment*Week*Density*Final.biomass.total, data=df)#better fit, 
hist(resid(mod2a))
qqnorm(resid(mod2a))
qqline(resid(mod2a))
anova(mod2a)
summary(mod2a)
plot(mod2a)

#t6 with biomass data (didn't set up model)

mod2ai<-lm(Egg.count~Treatment*Week*Density+(1|Year)+(1|Trial), data=df)#better fit, 
hist(resid(mod2ai))
qqnorm(resid(mod2ai))
qqline(resid(mod2ai))
anova(mod2ai)
summary(mod2ai)
plot(mod2ai)
#fixef(mod2ai)
#ranef(mod2ai)

#year and GLMM, poisson dist, keep getting error message
mod2aii<-glmer(Egg.count~Treatment*Week*Density+(1|Year)+(1|Trial),family=poisson, data=df)#better fit, 
hist(resid(mod2aii))
qqnorm(resid(mod2aii))
qqline(resid(mod2aii))
anova(mod2aii)
summary(mod2aii)
plot(mod2aii)
fixef(mod2aii)
ranef(mod2aii)

mod4<-lmer(Egg.count~Treatment*Week*Density+(1|Trial), data=df)
hist(resid(mod4))
qqnorm(resid(mod4))
qqline(resid(mod4))       
plot(mod4) #best fit so far, variances still look a bit funky because of some outlying points
anova(mod4)#density and week*density were sig.

mod5<-glmer(Egg.count~Treatment*Week*Density+(1|Trial),family=poisson,data=df)
hist(resid(mod5))
qqnorm(resid(mod5))
qqline(resid(mod5))       
plot(mod5)
Anova(mod5, type="II")
summary(mod5)
fixef(mod5)
ranef(mod5)
Anova(mod5)
#best fit so far, but it's tough to determine whether the variances are equal
#doing the plot function shows a clutster of points
#not sure how to analyze this, because when I try to run Anova (type=etc.) it doesn't seem to work

#reproduction per avg.total.biomass#####
mod2<-lm(egg.per.tot.avg.biomass~Treatment*Week, data=df)#more normal, but still not great
hist(resid(mod2))
qqnorm(resid(mod2))
qqline(resid(mod2))
anova(mod2)
summary(mod2)

mod2b<-lm(egg.per.tot.avg.biomass~Treatment*Week*Density, data=df)#more normal, but still not great
hist(resid(mod2b))
qqnorm(resid(mod2b))
qqline(resid(mod2b))
anova(mod2b)
summary(mod2b)

mod2aii<-glmer(egg.per.tot.avg.biomass~Treatment*Week*Density+(1|Year)+(1|Trial),family=poisson, data=df)#better fit, 
hist(resid(mod2aii))
qqnorm(resid(mod2aii))
qqline(resid(mod2aii))
anova(mod2aii)
summary(mod2aii)
plot(mod2aii)
fixef(mod2aii)
ranef(mod2aii)

#########B: Per capita output by week###############
#now going to look at corrected metrics
mod1<-lm(per.capita.repro~Treatment, data=df)#not normal
hist(resid(mod1))
anova(mod1)
#log egg counts
mod1log<-lm(ln.pc~Treatment, data=df)#worse than non-logged counts!
hist(resid(mod1log))
anova(mod1log)

mod2<-lm(per.capita.repro~Treatment*Week, data=df)#more normal, but still not great
hist(resid(mod2))
qqnorm(resid(mod2))
qqline(resid(mod2))
anova(mod2) #indicates that week had an effect, which makes sense, lost fish over time, fewer eggs
summary(mod2)
#log
mod2log<-lm(ln.pc~Treatment*Week, data=df)#bad fit, likely becuase of all of the zeros
hist(resid(mod2log))
qqnorm(resid(mod2log))
qqline(resid(mod2log))
anova(mod2log)
summary(mod2log)
#not going to waste any more time with log counts...looks terrible

mod3<-lm(per.capita.repro~Treatment*Week*Density, data=df)#better fit, 
hist(resid(mod3))
qqnorm(resid(mod3))
qqline(resid(mod3))
anova(mod3) #treatment not significant here, only week, density barely nonsignificant
summary(mod3)
plot(mod3)#Ok fit and normailty, but ddefinitely a shape to the variances
#don't really want to start pulling data points as outliers

mod4<-lmer(per.capita.repro~Treatment*Week*Density+(1|Trial), data=df)
hist(resid(mod4))
qqnorm(resid(mod4))
qqline(resid(mod4))       
plot(mod4) #best fit so far, variances still look a bit funky because of some outlying points
anova(mod4)#nothing significant

mod5<-glmer(per.capita.repro~Treatment*Week*Density+(1|Trial),family=poisson,data=df)
hist(resid(mod5))
qqnorm(resid(mod5))
qqline(resid(mod5))       
plot(mod5)
Anova(mod5, type="II")
summary(mod5)
fixef(mod5)
ranef(mod5)
Anova(mod5)
#best fit so far, but it's tough to determine whether the variances are equal
#doing the plot function shows a clutster of points
#not sure how to analyze this, because when I try to run Anova (type=etc.) it doesn't seem to work
#not sure what to do for the Lmertest

#tests after chatting with mark on 2019.3.7####
#run models three diferent ways: 
#1) ~Treatment*Density (no trial, dep.day, or week (2018 only) as factors)
#2) ~Treatment*Density*Trial*deployment.day (trial and dep.day as fixed)
#3) as is, with Trial and deployment.day as random factors

#decided to keep density in the model with treatment
#--because it seems to correlate positively with density
plot(Egg.count~Density,main="egg count vs. density, all data", data=df)
abline(lm(df$Egg.count~df$Density), col="red")# regression line (y~x)
lines(lowess(df$Density,df$Egg.count), col="blue") # lowess line (x,y)



#2017
#egg vs. density: no real relationship
plot(Egg.count~Density,main="egg count vs. density, 2017", data=df)
abline(lm(df$Egg.count~df$Density), col="red")# regression line (y~x)
lines(lowess(df$Density,df$Egg.count), col="blue") # lowess line (x,y)

#1) no difference by treatment, didn't include density
#GOING WITH THIS ONE
mod10<-lm(Egg.count~Treatment,data=df)
hist(resid(mod10))
qqnorm(resid(mod10))
qqline(resid(mod10))
anova(mod10) 
summary(mod10)
plot(mod10)
boxplot(Egg.count~Treatment,data=df) 
plot(Egg.count~Density,data=df) #total reproduction increases with density
abline(lm(df$Egg.count~df$Density))
plot(per.capita.repro~Density,data=df) #per.capita repro becomes less variable with density

#2) adding deployment day into the model as fixed factor
#don't think it improved the model, slightly better fit overall (R^2 value)
#--but still not great
mod10a<-lm(Egg.count~Treatment*Trial,data=df)
hist(resid(mod10a))
qqnorm(resid(mod10a))
qqline(resid(mod10a))
anova(mod10a)
summary(mod10a)
plot(mod10a)

#2a) now as random factor
#doesn't improve the model very much
mod10ai<-lm(Egg.count~Treatment+deployment.day,data=df)
hist(resid(mod10ai))
qqnorm(resid(mod10ai))
qqline(resid(mod10ai))
anova(mod10ai)
summary(mod10ai)
plot(mod10ai)
ranef(mod10ai)

#3) no effect of treatment or density on reproductive output
mod11aiv<-lmer(Egg.count~Treatment+(1|Trial:deployment.day),data=df)
hist(resid(mod11aiv))
qqnorm(resid(mod11aiv))
qqline(resid(mod11aiv))
Anova(mod11aiv)
anova(mod11aiv,type='I')
Anova(mod11aiv)
summary(mod11aiv)
plot(mod11aiv)
ranef(mod11aiv)

#2018.t4.5
#1)  no effect of treatment, lower total reproduction with lower density
mod10<-lm(Egg.count~Treatment*Density,data=df)
hist(resid(mod10))
qqnorm(resid(mod10))
qqline(resid(mod10))
anova(mod10) 
summary(mod10)
plot(mod10)
boxplot(Egg.count~Treatment,data=df) 
####kind of like this model for trial 4+5

#tukey test to look for differences, but only works if there are sig differences in the model
require(graphics)

summary(fm1 <- aov(breaks ~ wool + tension, data = warpbreaks))
TukeyHSD(fm1, "tension", ordered = TRUE)
plot(TukeyHSD(fm1, "tension"))
# }
summary(mod10<-aov(Egg.count~Week,data=df))
TukeyHSD(mod10, c("Week"), ordered = TRUE)
df$Week<-as.factor(df$Week)

#trying emmeans package to determine differences among treatments
#don't trust this, I don't think you can do tukey test with lmer, just mathematically
#library(emmeans)
#emmeans(model, list(pairwise ~ Group), adjust = "tukey")

#emmeans(mod10, list(pairwise ~ Treatment), adjust = "tukey")

#1a) sig ineraction for trt:week and trt:density, but not consistent among weeks
#both changes in direction and changes in relative position to each other
mod10a<-lm(Egg.count~Treatment*Density*Week,data=df)
hist(resid(mod10a))
qqnorm(resid(mod10a))
qqline(resid(mod10a))
anova(mod10a) 
summary(mod10a)
plot(mod10a)
boxplot(Egg.count~Treatment,data=df) #equal variance


#1b) no sig. effect of trial
mod10b<-lm(Egg.count~Treatment*Density*Week*Trial,data=df)
hist(resid(mod10b))
qqnorm(resid(mod10b))
qqline(resid(mod10b))
anova(mod10b) 
summary(mod10b)
plot(mod10b)
#boxplot(Egg.count~Treatment,data=df)

#1c) sig.effect of dep. day
mod10b<-lm(Egg.count~Treatment*Density*Week*deployment.day,data=df)
hist(resid(mod10b))
qqnorm(resid(mod10b))
qqline(resid(mod10b))
anova(mod10b) 
summary(mod10b)
plot(mod10b)
tapply(df$Egg.count,df$deployment.day,mean)
tapply(df$Egg.count,df$Treatment,mean)
tapply(df$Egg.count,df$Trial,mean)


#1d) removed week, dep day not sig. now, density driving trends
mod10bi<-lm(Egg.count~Treatment*Density*deployment.day,data=df)
hist(resid(mod10bi))
qqnorm(resid(mod10bi))
qqline(resid(mod10bi))
anova(mod10bi) 
summary(mod10bi)
plot(mod10bi)
tapply(df$Egg.count,df$deployment.day,mean)

#1e) dep.day as random factor
#WENT WITH THIS ONE
mod10bii<-lmer(Egg.count~Treatment*Density+(1|Trial:deployment.day),data=df)
hist(resid(mod10bii))
qqnorm(resid(mod10bii))
qqline(resid(mod10bii))
anova(mod10bii,type="I") 
summary(mod10bii)
ranef(mod10bii)
plot(mod10bii)

#2) 
mod10a<-lm(Egg.count~Treatment*Density*Trial*deployment.day,data=df)
hist(resid(mod10a))
qqnorm(resid(mod10a))
qqline(resid(mod10a))
anova(mod10a)
summary(mod10a)
plot(mod10a)

#3) no effect of treatment or density on reproductive output
mod11aiv<-lmer(Egg.count~Treatment*Density*Week+(1|Trial)+(1|Trial:deployment.day),data=df)
hist(resid(mod11aiv))
qqnorm(resid(mod11aiv))
qqline(resid(mod11aiv))
Anova(mod11aiv)
anova(mod11aiv,type='I')
Anova(mod11aiv)
summary(mod11aiv)
plot(mod11aiv)
ranef(mod11aiv)

#2018 t6
df<-egg.2018.t6
df$Treatment<-ordered(df$Treatment,levels=c("Low","Medium","High","Control"))

#egg vs. density: seems like repro increased with density
plot(Egg.count~Density,main="egg count vs. density, 2018 t6", data=df)
abline(lm(df$Egg.count~df$Density), col="red")# regression line (y~x)
lines(lowess(df$Density,df$Egg.count), col="blue") # lowess line (x,y)

#1) 
mod10<-lm(Egg.count~Treatment,data=df)
hist(resid(mod10))
qqnorm(resid(mod10))
qqline(resid(mod10))
anova(mod10) 
summary(mod10)
plot(mod10)
boxplot(Egg.count~Treatment,data=df) 
plot(Egg.count~Density,data=df) #total reproduction increases with density
abline(lm(df$Egg.count~df$Density))
plot(per.capita.repro~Density,data=df) #per.capita repro becomes less variable with density

#2) adding density to the model
#WENT WITH THIS ONE
mod10a<-lm(Egg.count~Treatment*Density,data=df)
hist(resid(mod10a))
qqnorm(resid(mod10a))
qqline(resid(mod10a))
anova(mod10a)
summary(mod10a)
plot(mod10a)

#2a) now adding week as a factor, because reproduction decreased over time
#not exactly sure how I feel about this, because density decreased by week
mod10ai<-lmer(Egg.count~Treatment*Density+(1|Week),data=df)
hist(resid(mod10ai))
qqnorm(resid(mod10ai))
qqline(resid(mod10ai))
anova(mod10ai)
summary(mod10ai)
plot(mod10ai)
ranef(mod10ai)

mod10aii<-lm(Egg.count~Treatment*Density*Week,data=df)
hist(resid(mod10aii))
qqnorm(resid(mod10aii))
qqline(resid(mod10aii))
anova(mod10aii)
summary(mod10aii)
plot(mod10aii)
ranef(mod10aii)

AIC(mod10a,mod10ai,mod10aii)

#### power analysis to determine how many reefs I would need to see
# sig. effect among treatments, given my current data

#for general linear models (which I am using)
pwr.f2.test(u =, v = , f2 = , sig.level = , power = )
#u=num df, v=den df, f2=effect size, sig.level = 0.05, power = ??)


#trial 1.2.3
pwr.f2.test(u =2, v = 51, f2 = 0.17, sig.level = 0.05, power = )


##looking at means with new dataset
tapply(df$Egg.count,df$Treatment,mean)

bargraph.CI(x.factor = Treatment, response = Egg.count, legend=FALSE, main="T1.2.3", xlab="Treatment", ylab="Eggs per reef (mean +/- se)", x.leg=9, yleg=4500, data = df)


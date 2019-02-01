######################################
#                                    #      
#   Egg counts with Trials 4+5 data  #
# 12/15/18                           #
#                                    #
######################################

#good notes for overally trends (as of 12/15/18):

#1.14.19 - want to come back and visit this, but I ned to figure out how to combine all trials into a single df

#low reproduction in weeks 1 and 3, but it ramps up in weeks 2 and 4
#highest per capita output in high risk trt
#highest reproduction in weeks 2 and 4, it's interesting to see this, but I think it's density dependent

rm(list=ls())

library(sciplot)
library(lme4)
library(car)
library(lmerTest)
library(dplyr)
library(ggplot2)

#gob.sub<-read.csv("Data/List of nests to count.12.15.18.csv") #just trials 4 and 5, no controal treatment
gob.sub<-read.csv("Data/List of nests to count.1.13.19.csv")
#subsetting these data to only reflect trials 4-6
gob.sub.46<-gob.sub[gob.sub$Trial>3, c("Trial","Week", "Reef","Treatment", "Egg.count")]
gob.sub.46$Week<-as.factor(gob.sub.46$Week)
gob.sub.46$Treatment<-ordered(gob.sub.46$Treatment,levels=c("Low","Medium","High","Control"))
View(gob.sub)
#density<-read.csv("C:\\Users\\George\\Desktop\\2018 summer\\2018 Goby\\2018 data for analyses, R\\2018 densities\\Data//density.12.4.18.csv")
density<-read.csv("C:\\Users\\George\\Desktop\\2018 summer\\2018 Goby\\2018 data for analyses, R\\2018 densities\\Data//density.1.14.19.csv")
density$Day<-as.factor(density$Day)
density$Week<-as.factor(density$Week)
density$Treatment<-ordered(density$Treatment,levels=c("Low","Medium","High","Control"))

View(density)
#this version includes information about BEG
#density<-read.csv("C:\\Users\\George\\Desktop\\2018 summer\\2018 Goby\\2018 data for analyses, R\\2018 densities\\Data//density.12.4.18.messing.withR.deleted.BEG.info.csv")
#this infomation doesn't
#had to enter double dashes in order to get it to read the path correctly
#gob.sub$week<-as.factor(gob.sub$week)#need to make week a factor???
View(gob.sub)
head(gob.sub)
gob.sub$Week<-as.factor(gob.sub$Week)
gob.sub$Treatment<-ordered(gob.sub$Treatment,levels=c("Low","Medium","High"))
#gob.sub$Treatment<-as.factor(gob.sub$Treatment)
#gob.egg.back$treatment

#analyses done on 12.15.18, before I had a chance to sum total reproduction by week
#need to figure out how to take raw data (with rows for individual TOL's/reef) and wrangle it, while keeping
#information for all other variables
#more importantly, how do I incorporate one single value for density for each reef?

bargraph.CI(x.factor = Week, response = Egg.count, group= Treatment, legend=TRUE, main="eggs", data = gob.sub)
#really interesting to see what is going on in week 3 with the high-risk trt. 
#it seems like there is 
bargraph.CI(x.factor = Treatment, response = Egg.count, main="eggs", data = gob.sub)# low, med, then high, althgugh likely
#not statistically different
bargraph.CI(x.factor = Week, response = Egg.count, main="eggs", data = gob.sub) #reproduction decreases over time
#bargraph.CI(x.factor = Treatment, response = Egg.count, group= Week, legend=TRUE, main="eggs", data = gob.sub)

##trying to figure out how to wrangle data in R
# SAMPLE: Create table...
xtabs(formula=count~name + sex, data=data)

xtabs(formula=Egg.count~Treatment, data=gob.sub)#this over-simplifies my dataframe, but it gives me what I want
#when I ony need to look at 1 or 2 variables.

#####mesing with dplyr to wrangle raw data

#this works, and it's kind of nice that it only brings in the variables that you want to look at
egg.per.week<-gob.sub %>%
  group_by(Trial,Reef,Week,Treatment) %>%
  summarize(eggs.per.reef = sum(Egg.count))#you can name the summarized variable something different
#in this case, "eggs.per.reef"

egg.per.week<-gob.sub %>%
  group_by(Trial,Reef,Week,Treatment) %>%
  summarize(Egg.count = sum(Egg.count))#or you can name it the same thing as before - they give the same result

#the only reason I could see for using the same variable name, is to simplify having to make quick figures where you
#already  the code

egg.per.week

bargraph.CI(x.factor = Week, response = eggs.per.reef, group= Treatment, legend=TRUE, main="eggs", data = egg.per.week)
#really interesting to see what is going on in week 3 with the high-risk trt. 
#it seems like there is 
bargraph.CI(x.factor = Treatment, response = Egg.count, main="eggs", data = gob.sub)# low, med, then high, althgugh likely
#not statistically different
bargraph.CI(x.factor = Week, response = Egg.count, main="eggs", data = gob.sub) #reproduction decreases over time
#bargraph.CI(x.factor = Treatment, response = Egg.count, group= Week, legend=TRUE, main="eggs", data = gob.sub)


#trying to see if I can work around density

#emailing with Nyssa on 12/15/18:
#re: how to get density in the df
#You can put density = mean(density) (or whatever the density column name is) 
#to get the mean density or sum (density) to get the total for your week by group in the summarize line. 

#Or, if you just want to add all the data that are numeric then you could use the summarize_if() function
#and tell it to do the same function to all the numeric columns.

#Or, if density is in a different dataframe then use the left_join() function to merge the 
#two dataframes and it will put the appropriate density value for all treatments, weeks, etc 
#as long as the columns and factors are named identically.

#trying on my own with option 3. I was running into issues because there are some dates where I surveyed densities, but not
#nests for eggs. Becasue of that, I had to remove "Date" from the joining 

left_join(gob.sub,density,by=c("Trial","Week","Day","Reef","Treatment"))
#left_join(density,gob.sub,by="Egg.count")

m.count.den<-inner_join(gob.sub,density,by=c("Trial","Week","Day","Reef","Treatment"))#this worked, and it put the same value
#for each density of there were multiple TOL's used for the reef

#going to combine trials 4 through 6 with density values for these same trials
m.count.den1<-inner_join(gob.sub.46,density,by=c("Trial","Week","Day","Reef","Treatment"))#this worked, and it put the same value
#for each density of there were multiple TOL's used for the reef

egg.per.week<-m.count.den %>%
  group_by(Trial,Reef,Week,Treatment,Density) %>%
  summarize(eggs.per.reef = sum(Egg.count))

View(egg.per.week)


m.count.den<-left_join(gob.sub,density,by=c("Trial","Week","Day","Reef","Treatment"))#this worked, and it put the same value


########going to work with density now
#want to get a df that has average density by week
#wanted to get the average density per reef per week (mean(x)), 
#but I also wanted to round up to the nearest whole number (ceiling(x))
den.per.week<-density %>%
  group_by(Trial,Reef,Week,Treatment) %>%
  summarize(Density = ceiling(mean(Density))) #including day kind of messes things up, b/c
#there were some days where I got densities but not egg counts
View(den.per.week)
#for some reason, I'm losing 27 rows of data when I run the dplyr function, and I'm not sure which ones it is

View(density)
View(den.per.week)

######cleaning up code for what worked############
#FYI, this is the code you want to use, not the other ones up above, things got wonky up there

#summing egg counts by week
egg.per.week<-gob.sub %>%
  group_by(Trial,Reef,Week,Treatment) %>%
  summarize(Egg.count = sum(Egg.count))

egg.per.week

#averaging densities by week
den.per.week<-density %>%
  group_by(Trial,Reef,Week,Treatment) %>%
  summarize(Density = ceiling(mean(Density)))

den.per.week

#going to join the two of them now:

den.and.egg<-left_join(egg.per.week,den.per.week,by=c("Trial","Week","Reef","Treatment"))
#left_join(density,gob.sub,by="Egg.count")

#now going to calculate per capita output by dividing total egg count by density
den.avg.egg<-mutate(den.and.egg,per.capita.repro=Egg.count/Density)
den.avg.egg
View(den.avg.egg)

#manipulating df variables

#order treatments
den.avg.egg$Treatment<-ordered(den.avg.egg$Treatment,levels=c("Low","Medium","High"))
#store week as a factor for analyses
den.avg.egg$Week<-as.factor(den.avg.egg$Week)
den.avg.egg

####quick visualization: plotting average per capita output per reef per treatment over time

bargraph.CI(x.factor = Week, response = per.capita.repro, group= Treatment, legend=TRUE, main="eggs", data = den.avg.egg)
#low reproduction in weeks 1 and 3, but it ramps up in weeks 2 and 4
bargraph.CI(x.factor = Treatment, response = per.capita.repro, main="eggs", data = den.avg.egg)# low, med, then high, althgugh likely
#highest per capita output in high risk trt
bargraph.CI(x.factor = Week, response = per.capita.repro, main="eggs", data = den.avg.egg) #reproduction decreases over time
#highest reproduction in weeks 2 and 4, it's interesting to see this, but I think it's density dependent

#next:

#1) model selection
#2) figure out if trials can be pooled, of if they need to be run separately
#3) run trial 4 and 5 separately, get plots, and look at them (use filter function to get just trial 4, and just trial 5)
#4) look at density data over time to see if there are differences between trials as well
#5) figure out what to do with temp data - same as density data? Just get average temp per week to see if it's affecting output?

#question: if we're not seeing differences in output, is that due to more nests with fewer eggs in each nest for high risk
#vs. low vs. med? What is the mechanism? Egg spreading vs. putting all eggs in one basket? could be an interesting approach
#plus it doesn't take that much time to look at, and it could be an interesting part of the story

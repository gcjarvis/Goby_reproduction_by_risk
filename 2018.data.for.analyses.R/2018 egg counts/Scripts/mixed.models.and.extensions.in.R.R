# Description: practice with mixed models from Zurr et al. 2009
# Author: George C Jarvis
# Date: Wed Mar 18 16:37:47 2020
# Notes: datasets are nested in the data for egg counts, but will have to map to them through
# - the 'read.table' function
# --------------

rm(list=ls())

# 2.1.1 Clevland dotplots - grouping points together and looking for heterogeneity

Nereis <- read.table(file = "Data/mixed.model.practice.zurr.et.al/Nereis.txt", header = TRUE, dec = ".")
names(Nereis)

# data exploration with dotplots
dotchart(Nereis$concentration,
         ylab= "Order of observations",
         xlab= "Concentration", main= "Cleveland dotplot")

#now adding a grouping factor
dotchart(Nereis$concentration,
         groups = factor(Nereis$nutrient),
         ylab= "Nutrient",
         xlab= "Concentration", main= "Cleveland dotplot",
         pch= Nereis$nutrient)

#2.1.2 Pairplots - can look for correllations among variables

pairs(Nereis)

#2.1.3 Boxplots - used to compare variances and compare median values

boxplot(concentration ~ factor(nutrient),
        varwidth = TRUE, xlab = "nutrient", #varwidth sets width of boxplot to sample size
        main = "Boxplot of concentration conditional on
        nutrient", ylab = "concentration", data = Nereis)

#2.1.4 xyplot from lattice package - used to compare trends among different samples
#more appropriate with more than 2 explanatory variables

# in this example, are relationships between niotrogen and age the same for all whale teeth sampled?

library(lattice)
TeethNitrogen <- read.table(file = "Data/mixed.model.practice.zurr.et.al/TeethNitrogen.txt", header = TRUE, dec = ".")

names(TeethNitrogen)

pairs(TeethNitrogen) #seems to be an increasing relationship between age and nitrogen

xyplot(X15N ~ Age|factor(Tooth), type="l",
       xlab = "Estimated age", col=1,
       ylab= expression(paste(delta^{15},"N")),
       strip = function(bg='white',...)
       strip.default(bg='white',...),
       data= TeethNitrogen)

# 2.3.6 testing assumptions of linear models

# wedge clam length and mass data throughout the year
# - the relationship between length and mass depends on the month measured
# - we want to include month 

Clams <- read.table(file = "Data/mixed.model.practice.zurr.et.al/Clams.txt", header = TRUE, dec = ".")

pairs(Clams)
#nonlinear relationship between length and mass, so log transformed both
# - BUT only because we want to analyze with linear regression

Clams$fMONTH<-factor(Clams$MONTH)
coplot(LNAFD~LNLENGTH|fMONTH, data=Clams) #relationship changes throughout the season

#coplot(egg.week~Treatment|Trial, data = repro)

#modeling clam regression
M1<- lm(LNAFD~LNLENGTH*fMONTH, data= Clams)
drop1(M1,test="F") #shows us that there is a sig. interaction, sp we move to next assumptions

#model validation: residuals vs. fitted values for homogeneity, QQ-plot or hist. for
# - normality, residuals vs. explanatory variable to check independence

op<-par(mfrow= c(2,2), mar=c(5,4,1,2))
plot(M1, add.smooth = FALSE, which=1)
E<-resid(M1)
hist(E,xlab="Residuals", main="")
plot(Clams$LNLENGTH,E,xlab="Log(Length)",
     ylab="Residuals")
plot(Clams$fMONTH,E,xlab="Month",ylab="Residuals")
par(op)

#teeth data again, same exploration
TN<-TeethNitrogen
M2<-lm(X15N~Age, subset = (TN$Tooth == "Moby"), data=TN) #that subset function is handy...
op<-par(mfrow = c(2,2))
plot(M2, add.smooth = FALSE) #also nice function to get all in one plane
par(op)

#it looks like there is a pattern in the left panels, 
# - suggesting non-independence of observations

#examine with linear regression

N.Moby<-TN$X15N[TN$Tooth=="Moby"]
Age.Moby<-TN$Age[TN$Tooth=="Moby"]

plot(y=N.Moby,x=Age.Moby,
     xlab="Estimated Moby Age",
     ylab=expression(paste(delta^{15},"N Moby")))
abline(M2)

summary(M2)
#output tells us what the equation for the regression relating nitrogen content and age
# - yi= y-intercept estimate + estimate for factor (slope)* xi

#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
#(Intercept) 11.748940   0.163559   71.83   <2e-16 ***
#Age          0.113794   0.006186   18.40   <2e-16 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 0.4859 on 40 degrees of freedom
#Multiple R-squared:  0.8943,	Adjusted R-squared:  0.8917 
#F-statistic: 338.4 on 1 and 40 DF,  p-value: < 2.2e-16


# so from this model summary we get:
# y(i)= 11.749 + 0.114 * age(i)

# estimated sloe and intercept are significantly different from 0 @5% level
# model explains 89% of the variation
# the estimator for standard deviation is equal to 0.486 (listed as "residual standard error" in summary)

#never knew that's what the output meant!



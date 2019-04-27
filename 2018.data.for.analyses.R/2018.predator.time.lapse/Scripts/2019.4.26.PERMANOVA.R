# Description: PERMANOVA for reporting of results
# Author: George C Jarvis
# Date: Fri Apr 26 13:22:39 2019
# --------------

#NOTE: these are the data for observations of predators by species, including the 

rm(list=ls())

library(sciplot)
library(lme4)
library(car)
library(lmerTest)
library(dplyr)
library(ggplot2)
library(tidyr)
library(plyr)
library(vegan)
library(psych)
library(vegan)


#importing raw data###
#old df ptl<-read.csv("Data/2019.3.19.ptl.csv")
#new df as of 2019.4.26
ptl<-read.csv("Data/2019.4.26a.ptl.csv")
#arranging by treatment for PERMANOVA + PermDisp
ptl<-arrange(ptl, Treatment)
View(ptl)
colSums(is.na(ptl))# no NA's

#1)permanova for number of number of photos that contained a pred#####
permanovamodel1<-adonis(contain.pred~Treatment, 
                        data = ptl, permutations = 100,
                        method="euclidean", by= "terms")
permanovamodel1 #suggests differences in location of means in multivariate space

#permdisp for counts
## euclidean distances between samples
dis <- vegdist(ptl$contain.pred, method = "euclidean")

#including each factor, along with the number of replicates for each (#factor, #replicates)
groups <- factor(c(rep(1,54),rep(2,207),rep(3,170),rep(4,154)), labels =
  c("Control","High","Low","Medium"))

## Calculate multivariate dispersions
mod <- betadisper(dis, groups)
mod

## Perform test
anova(mod)

## Permutation test for F
permutest(mod, pairwise = TRUE)

## Tukey's Honest Significant Differences
(mod.HSD <- TukeyHSD(mod))
plot(mod.HSD)

## Plot the groups and distances to centroids on the
## first two PCoA axes
plot(mod)

#2)sublethal predation, without altering dataset to exclude low risk####
#removing non-comparable treatments (i.e. low risk)
ptl<-read.csv("Data/2019.4.26a.ptl.no.low.trt.csv")
#permanova 
permanovamodel2<-adonis(sublethal~Treatment, 
                        data = ptl, permutations = 100,
                        method="euclidean", by= "terms")
permanovamodel2 #no differences in location of datapoints in multivariate space

#permdisp for counts
## euclidean distances between samples
dis2 <- vegdist(ptl$sublethal, method = "euclidean")

#with data for low-risk still included:
#including each factor, along with the number of replicates for each (#factor, #replicates)
groups <- factor(c(rep(1,54),rep(2,207),rep(3,170),rep(4,154)), labels =
                   c("Control","High","Low","Medium"))

#without data for low risk 
groups <- factor(c(rep(1,54),rep(2,207),rep(3,154)), labels =
                   c("Control","High","Medium"))

## Calculate multivariate dispersions
mod2 <- betadisper(dis2, groups) #suggests none here
mod2

## Perform test
anova(mod2)

## Permutation test for F
permutest(mod2, pairwise = TRUE)

## Tukey's Honest Significant Differences
(mod.HSD2 <- TukeyHSD(mod2))
plot(mod.HSD2)

## Plot the groups and distances to centroids on the
## first two PCoA axes
plot(mod2)

#3)lethal predation, without data for low and medium-risk treatments########
#removing non-comparable treatments (i.e. low and medium risk)
ptl<-read.csv("Data/2019.4.26a.ptl.no.low.no.med.trt.csv")
#permanova 
permanovamodel3<-adonis(lethal~Treatment, 
                        data = ptl, permutations = 100,
                        method="euclidean", by= "terms")
permanovamodel3 #no differences in location of datapoints in multivariate space

#permdisp for counts
## euclidean distances between samples
dis3 <- vegdist(ptl$lethal, method = "euclidean")

#without data for low and medium risk
#including each factor, along with the number of replicates for each (#factor, #replicates)
groups <- factor(c(rep(1,54),rep(2,207)), labels =
                   c("Control","High"))

## Calculate multivariate dispersions
mod3 <- betadisper(dis3, groups) #suggests none here
mod3

## Perform test
anova(mod3)

## Permutation test for F
permutest(mod3, pairwise = TRUE)

## Tukey's Honest Significant Differences
(mod.HSD3 <- TukeyHSD(mod3))
plot(mod.HSD3)

## Plot the groups and distances to centroids on the
## first two PCoA axes
plot(mod2)


#nmds plot, not important for now######

library(vegan)


PercentTotal<-read.csv('data/inverts_common_nMDS.csv')
View(PercentTotal)

#create	the	ordination	output	using	bray	curtis
ord<-metaMDS(PercentTotal[,-1],k=2,	distance='bray')	


#	if	it	does	not	converge	add	more iterations
ord<-metaMDS(PercentTotal[,-1],k=2,	distance='bray',	trymax	=	100)	#add	more	iterations
#let's	look	at	the	2D	stress.	Is	it	<	0.3?	
ord$stress
#	It	is	0.2	which	is	"good/ok".	So,	let's	continue
#	Let's	look	at	the	stress	plot

stressplot(ord)
#	looks	like	a	good	fit,	want	to	minimize	scatter
#	basic	plot
ordiplot(ord)	#	dots	represent	sites	(tide	pools	in	Nyssa's	case)	and	
#+	represents	species
#	add	text

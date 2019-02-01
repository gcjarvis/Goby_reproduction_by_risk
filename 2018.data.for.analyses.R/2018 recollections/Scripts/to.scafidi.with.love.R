#scafidi's analyses
#for use with WSN


library(ggplot2)
library(extrafont)
library(plyr)

abund.g<-read.csv("Data/Invert_abundance_forWSN.csv")

abund.g$trans<-(abund.g$Total.Abundance+0.1)^(1/3)

egg.means<-with(gob.sub, aggregate((egg.count), list(treatment=treatment,week=week), mean))
egg.means
#now apply the se function to the 4th column [,3]
reco.means$se<-with(reco, aggregate((total.recollection), list(Sex=Sex,treatment=treatment), function(x) sd(x)/sqrt(length(x))))[,3]

abund.means<-with(abund.g, aggregate((Total.Abundance), list(Algal.Species), mean))
abund.means
#now apply the se function to the 4th column [,3]
abund.means$se<-with(abund.g, aggregate((Total.Abundance), list(Algal.Species), function(x) sd(x)/sqrt(length(x))))[,2]

#################################################
# Length to biomass curve for bluebanded gobies #                                                              #
#    10/28/2018                                 #
#   George Jarvis                               #
#################################################

#want to get the equation of the line from a linear regression, 
#then figure out the initial biomass for reefs
#can take the average of the initial and final biomass values in order to
#account for differences in recollections among treatments

#example of how to get a model and to get equation of the line for the plot
model <- lm(y ~ x, data = mydf)
abline(model, col = "red")
summary(model)

length.biomass<-read.csv("Data/10.28.18.length to biomass.csv")
length.biomass<-length.biomass[,-3]#not sure why there was a random third column, so I deleted it
head(length.biomass)

lb.mod<-lm(weight~size, data=length.biomass)
summary(lb.mod)
plot(lb.mod)

LBplot<-plot(weight~size, data=length.biomass)
abline(lb.mod, col="red")

#trying to figure out the best way to plot a power function with nls (nonlinear least squares):
#taken from: http://www.css.cornell.edu/faculty/dgr2/teach/R/R_CurveFit.pdf

#plotting data to a generic power curve:
#stock:
> plot(y ~ x, main = "Known cubic, with noise")
> s <- seq(0, 1, length = 100)
> lines(s, s^3, lty = 2, col = "green")

#my data...doesn't seem to want to plot a fake surve to my data (lty=line type):
plot(weight ~ size, main = "Known cubic, with noise", data=length.biomass)
s <- seq(0, 40, length = 1000)
lines(s, s^3, lty = 1, col = "red")


#looking at a fitted curve using the nls function:
#stock:
m <- nls(y ~ I(x^power), data = ds, start = list(power = 1),
         + trace = T)
#my code:
m <- nls(weight ~ I(size^power), data = length.biomass, start = list(power = 1), trace = T)
m
summary(m)
#seems to work, but it's not as good as the fit that I got in excel
#after I checked the known values against the calculated values

#another shot with nls...not working:
#stock:
nls(y~b*x^z,start = list(b = 0, z = 1),data=your database)
#my code:
nls(weight~b*size^z,start = list(b = 0, z = 1),data=length.biomass)
#: Error in nlsModel(formula, mf, start, wts) : 
#     singular gradient matrix at initial parameter estimates


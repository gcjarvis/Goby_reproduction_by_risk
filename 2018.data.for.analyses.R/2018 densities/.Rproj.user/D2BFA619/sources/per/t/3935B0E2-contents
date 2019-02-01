##########################
#    Sciplot master code #
#    By: Lansing Perng   #
#    4/29/2018           #
##########################

rm(list=ls())

bargraph.CI(data = Data, #specify data frame
            x.factor = Teacher, #factor on x axis
            response = Score, #response variable
            group = Gender, #indexing factor--for legend
            col = c("pink", "lightblue1"), #color in order on legend, ordered alphabetically
            #border = NA, #removes border
            xlab = "Teacher", #label x axis
            ylab = "Scores", #label y axis
            main = "Scores by Teacher and Gender of Student", #title
            font.main = 1, #make title not bold
            space = c(0,0.3), #spaces b/w bars--0 within levels, 0.3 b/w levels of factor on x axis
            legend = TRUE, #default is FALSE
            x.leg = 7.5, #x.leg specifies legend placement on x axis, since there are no numbers, this takes some trial and error
            leg.lab = c("Female", "Male"), #specifies labels for indexing factor in alphabetical order
            cex.leg = 1.1, #expansion factor for legend--multiply default size by 1.1
            cex.lab = 1.2, #expansion factor of x and y axis labels
            cex.names = 1.1, #expansion factor for names of factor on x axis (4 levels of Teacher)
            cex.axis = 1, #expansion factor for numbered axis labels
            cex.main = 1.2, #expansion factor for title
            ylim = c(0,100), #sets y limits: min and max 
            names.arg = c("Jack", "John", "Jill", "Jane")) #names of x axis levels
abline(h=0) #adds horizontal line at y = 0-->adds x axis line
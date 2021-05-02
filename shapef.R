# Working with shapefile in R as df
library(rgdal)

mySHP <- file.choose()

myFILE <- readOGR(mySHP)

mode(myFILE)
length(myFILE)

library(ggplot2)

myDF <- fortify(myFILE)

class(myDF)

head(myDF)


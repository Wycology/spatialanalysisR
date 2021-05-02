# Working with shapefile in R as df
library(rgdal)

mySHP <- file.choose()

myFILE <- readOGR(mySHP)

mode(myFILE)

# Working with shapefile in R as df
library(rgdal) # library for working with shapefile in R

mySHP <- base::file.choose() # Choosing the shapefile from directory    

myFILE <- rgdal::readOGR(mySHP) # Reading the file into R

base::mode(myFILE) # Checking the mode of the file
base::length(myFILE) # Checking the number of features in the file

library(ggplot2) # Loading the package with fortify function

myDF <- ggplot2::fortify(myFILE) # Running the fortify function on the object to create dataframe

base::class(myDF) # Confirming the class of the dataframe

head(myDF) # the first few rows of the dataframe

nrow(myDF) # checking the number of rows

summary(myDF) # Unnecessary summary of the dataframe

library(tidyverse)
tibble::glimpse(myDF)

ggplot(myDF, aes(long, lat, group = group)) +
  coord_map('mercator') +
  geom_path()
# sinusoidal, cylequalarea, cylindrical, rectangular, gall, mollweide, gilbert,
# azequidistant, gnomonic, perspective, orthographic, stereographic, laue, fisheye,
# newyorker, conic, simpleconic, lambert, albers, bonne, polyconic, aitoff, lagrange,
# bicentric, elliptic, globular, vandergrinten, eisenlohr, guyou, square, tetra,
# harrison, trapezoidal, lune, mecca, homing, sp\_mercator, sp\_albers
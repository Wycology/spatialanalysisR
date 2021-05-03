# Working with shapefile in R as df
library(rgdal) # library for working with shapefile in R

mySHP <- file.choose() # Choosing the shapefile from directory    

myFILE <- readOGR(mySHP) # Reading the file into R

mode(myFILE) # Checking the mode of the file
length(myFILE) # Checking the number of features in the file

library(ggplot2) # Loading the package with fortify function

myDF <- fortify(myFILE) # Running the fortify function on the object to create dataframe

class(myDF) # Confirming the class of the dataframe

head(myDF) # the first few rows of the dataframe

nrow(myDF) # checking the number of rows

myDF

summary(myDF)

library(tidyverse)
glimpse(myDF)

ggplot(myDF, aes(long, lat, group = group)) +
  coord_map('mercator') +
  geom_path()
# sinusoidal, cylequalarea, cylindrical, rectangular, gall, mollweide, gilbert,
# azequidistant, gnomonic, perspective, orthographic, stereographic, laue, fisheye,
# newyorker, conic, simpleconic, lambert, albers, bonne, polyconic, aitoff, lagrange,
# bicentric, elliptic, globular, vandergrinten, eisenlohr, guyou, square, tetra,
# harrison, trapezoidal, lune, mecca, homing, sp\_mercator, sp\_albers
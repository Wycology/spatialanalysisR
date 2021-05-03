# Working with shapefile in R as df
library(rgdal) # library for working with shapefile in R

mySHP <- file.choose() # Choosing the shapefile from directory    

myFILE <- readOGR(mySHP)

mode(myFILE)
length(myFILE)

library(ggplot2)

myDF <- fortify(myFILE)

class(myDF)

head(myDF)
nrow(myDF)

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
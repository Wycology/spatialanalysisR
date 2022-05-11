# Working with shapefile in R as data frame
# Created on 30th August 2021
# Last edited on 11th May 2022 by the author

library(rgdal)        # Package version 1.5.31
library(tidyverse)    # Package version 1.3.1
library(broom)        # Package version 0.8.0

myFILE <- readOGR('shp/dummy.shp')

mode(myFILE) 
length(myFILE) 

myDF <- broom::tidy(myFILE)

class(myDF) 

head(myDF) 

nrow(myDF)

summary(myDF)

glimpse(myDF)

# Making a simple map

myDF |> ggplot(aes(long, lat, group = group, col = group)) +
  coord_map('mercator') +
  geom_path()

# All done
# sinusoidal, cylequalarea, cylindrical, rectangular, gall, mollweide, gilbert,
# azequidistant, gnomonic, perspective, orthographic, stereographic, laue, fisheye,
# newyorker, conic, simpleconic, lambert, albers, bonne, polyconic, aitoff, lagrange,
# bicentric, elliptic, globular, vandergrinten, eisenlohr, guyou, square, tetra,
# harrison, trapezoidal, lune, mecca, homing, sp\_mercator, sp\_albers
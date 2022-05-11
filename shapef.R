# Working with shapefile in R as df
# Created on 30th August 2021
# Last edited on 30th August 2021 by the author

base::library(rgdal) 

mySHP <- base::file.choose()    

myFILE <- rgdal::readOGR(mySHP) 

base::mode(myFILE) 
base::length(myFILE) 

library(ggplot2) 

myDF <- ggplot2::fortify(myFILE) 

base::class(myDF) 

utils::head(myDF) 

base::nrow(myDF)

base::summary(myDF)

base::library(tidyverse)
tibble::glimpse(myDF)

ggplot2::ggplot(myDF, ggplot2::aes(long, lat, group = group)) +
  ggplot2::coord_map('mercator') +
  ggplot2::geom_path()

# All done
# sinusoidal, cylequalarea, cylindrical, rectangular, gall, mollweide, gilbert,
# azequidistant, gnomonic, perspective, orthographic, stereographic, laue, fisheye,
# newyorker, conic, simpleconic, lambert, albers, bonne, polyconic, aitoff, lagrange,
# bicentric, elliptic, globular, vandergrinten, eisenlohr, guyou, square, tetra,
# harrison, trapezoidal, lune, mecca, homing, sp\_mercator, sp\_albers
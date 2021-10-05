# Working with shapefile in R as df
# Created on
# Last edited on 30th August 2021 by the author

base::library(rgdal) 

mySHP <- base::file.choose()    

myFILE <- rgdal::readOGR(mySHP) 

base::mode(myFILE) 
base::length(myFILE) 

library(ggplot2) 

myDF <- ggplot2::fortify(myFILE) # Running the fortify function on the object to create dataframe

base::class(myDF) # Confirming the class of the dataframe

utils::head(myDF) # the first few rows of the dataframe

base::nrow(myDF) # checking the number of rows

base::summary(myDF) # Unnecessary summary of the dataframe

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
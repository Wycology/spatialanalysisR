# Working with shapefile in R as data frame
# Created on 30th August 2021
# Last edited on 11th May 2022 by the author

library(rgdal)        # Package version 1.5.31
library(tidyverse)    # Package version 1.3.1
library(broom)        # Package version 0.8.0

my_shapefile <- readOGR('shp/dummy.shp')

mode(my_shapefile)     # This is an S4 object, typical of GIS files
length(my_shapefile)   # Reveals the number of features in the layer...like two roads

my_tibble <- tidy(my_shapefile)

class(my_tibble) 

head(my_tibble) 

nrow(my_tibble)

summary(my_tibble)

glimpse(my_tibble)

# Making a simple map

my_tibble |> ggplot(aes(long, lat, group = group, col = group)) +
  coord_map('mercator') +
  geom_path()

# All done
# sinusoidal, cylequalarea, cylindrical, rectangular, gall, mollweide, gilbert,
# azequidistant, gnomonic, perspective, orthographic, stereographic, laue, fisheye,
# newyorker, conic, simpleconic, lambert, albers, bonne, polyconic, aitoff, lagrange,
# bicentric, elliptic, globular, vandergrinten, eisenlohr, guyou, square, tetra,
# harrison, trapezoidal, lune, mecca, homing, sp\_mercator, sp\_albers
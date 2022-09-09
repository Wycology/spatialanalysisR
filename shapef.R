# Working with shapefile in R as data a frame      #
# Created on 30th August 2021                      #
# Last edited on 8th September 2022 by the author  #

library(rgdal)        # Version 1.5.32
library(tidyverse)    # Version 1.3.2
library(broom)        # Version 1.0.1
library(patchwork)    # Version 1.1.2

my_shapefile <- readOGR('shp/dummy.shp')

mode(my_shapefile)     # This is an S4 object.
length(my_shapefile)   # Number of features in the layer...like two roads

my_tibble <- tidy(my_shapefile) # Tidy the shapefile

class(my_tibble) 

head(my_tibble) 

nrow(my_tibble)

summary(my_t

glimpse(my_tibble)

# Making map

my_tibble |> ggplot(aes(long, lat, group = group, col = group)) +
  geom_path()
my_tibble |> ggplot(aes(long, lat, group = group, col = group)) +
  coord_map('mercator') +
  geom_path()

# Plotting regardless of coordinates
a <- my_tibble |> ggplot(aes(long, lat, group = group, col = group)) +
  geom_path()

# Consideration the spatial nature of the data
b <- my_tibble |> ggplot(aes(long, lat, group = group, col = group)) +
  coord_map('mercator') +
  geom_path()

b / a

# All done
# sinusoidal, cylequalarea, cylindrical, rectangular, gall, mollweide, gilbert,
# azequidistant, gnomonic, perspective, orthographic, stereographic, laue, fisheye,
# newyorker, conic, simpleconic, lambert, albers, bonne, polyconic, aitoff, lagrange,
# bicentric, elliptic, globular, vandergrinten, eisenlohr, guyou, square, tetra,
# harrison, trapezoidal, lune, mecca, homing, sp\_mercator, sp\_albers
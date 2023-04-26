####################################################
# Working with shapefile in R as data a frame      #
# Created on 30th August April 2021                #
# Last edited on 8th April 2023 by the author      #
####################################################

library(sf)           # Version 1.0.12
library(ggplot2)      # Version 3.4.2

my_shapefile <- st_read('shp/dummy.shp')

mode(my_shapefile)     # This i.
length(my_shapefile)   # Number of features in the layer...like two roads

head(my_shapefile) 

nrow(my_shapefile)


# Making map --------------------------------------------------------------

my_shapefile |> ggplot() +
  geom_sf()
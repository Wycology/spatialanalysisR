####################################################
# Working with shapefile in R as data a frame      #
# Created on 30th August April 2021                #
# Last edited on 8th April 2023 by the author      #
###########

library(sf)           # Version 1.5.32
library(tidyverse)    # Version 1.3.2
library(broom)        # Version 1.0.1
library(patchwork)    # Version 1.1.2

my_shapefile <- st_read('shp/dummy.shp')

mode(my_shapefile)     # This is S4 object.
length(my_shapefile)   # Number of features in the layer...like two roads

head(my_shapefile) 

nrow(my_shapefile)

# Making map -------------------------------------------------------------------------------------

my_shapefile |> ggplot() +
  geom_sf()


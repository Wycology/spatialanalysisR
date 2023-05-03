####################################################
# Working with shapefile in R as data a frame      #
# Created on 30th August April 2021                #
# Last edited on 8th April 2023 by the author      #
####################################################

library(sf)           # Version 1.0.12
library(ggplot2)      # Version 3.4.2

my_shapefile <- st_read('shp/dummy.shp')

mode(my_shapefile)     #.
length(my_shapefile)   # Number of features in the layer...like two roads

head(my_shapefile) 

nrow(my_shapefile)


# Making map --------------------------------------------------------------

my_shapefile |> ggplot() +
  geom_sf()

# Just learnt of creating gridded hexagon using sf

# Let start by creating a polygon 

poly <- st_sfc(st_polygon(list(rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 0)))))

hex_grid <- st_make_grid(poly, cellsize = 0.05, square = FALSE) 

plot(hex_grid, add = TRUE)

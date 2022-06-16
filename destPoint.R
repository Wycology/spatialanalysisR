#####################################
# Created by Wyclife Agumba Oluoch  #  
# Created on: 19th Apr 2021         #
# Last edited on 28th M 2022      #
#####################################

# Task: Calculate the coordinates of destination given the starting 
# point coordinates in degrees, bearing in degrees, while distance in meters

library(geosphere) # Loading the geosphere package

start_longitude <- 35

start_latitude <- 3

bearing_to_destination <- 60

distance_to_destination <- 125

destPoint(p = c(start_longitude, start_latitude), b = bearing_to_destination, 
                     d = distance_to_destination) # Gives coordinates of a destination point.

the_mat <- data.frame(lon = c(35, 35.2, 35.6),
                      lat = c(3, 3.2, 3.6))

the_matrix <- as.matrix(the_mat)

destPoint(p = the_matrix, b = abs(rnorm(3, 100, 30)), d = rnorm(3, 150, 25))

# The end
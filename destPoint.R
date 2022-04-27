#####################################
# Created by Wyclife Agumba Oluoch  #  
# Created on: 19th Apr 2021         #
# Last edited on 15th April 2022    #
#####################################

# Task: Calculate coordinates of destination point given the starting 
# point coordina in degrees, bearing in degrees, and distance in metres

library(geosphere) # Loading the necessary library

geosphere::destPoint(p = c(35, 3), b = 60, d = 125) # Gives coords of a point.

the_mat <- base::data.frame(lon = c(35, 35.2, 35.6),
                      lat = c(3, 3.2, 3.6))
the_matrix <- base::as.matrix(the_mat)

geosphere::destPoint(p = the_matrix, b = base::abs(stats::rnorm(3, 100, 30)), d = stats::rnorm(3, 150, 25))

# The end
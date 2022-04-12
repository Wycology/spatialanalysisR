#####################################
# Created by Wyclife Agumba Oluoch  #  
# Created on: 19th Apr 2021         #
# Last edited on 12th April 2022    #
#####################################

# Task: Calculating  coordinates of the destination point given the starting 
# point in degrees, bearing in degrees, and distance in metres

library(geosphere) # Loading the library

geosphere::destPoint(p = c(35, 3), b = 60, d = 125) # Gives the coords of a point which is 
                                        # bearing and distance

the_mat <- base::data.frame(lon = c(35, 35.2, 35.6),
                      lat = c(3, 3.2, 3.6))
the_matrix <- base::as.matrix(the_mat)

geosphere::destPoint(p = the_matrix, b = base::abs(stats::rnorm(3, 100, 30)), d = stats::rnorm(3, 150, 25))

# The end
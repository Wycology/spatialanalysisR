####################################
# Created by Wyclife Agumba Oluoch    
# Created on: 19th Apr 2021
# Last edited on 14th July 2021
####################################

# Task: Calculating the coordinates of the destination point given the starting 
# point in degrees, bearing in degrees, and distance in metres

library(geosphere) # The packagefor the analyses

destPoint(p = c(35, 3), b = 60, d = 125) # Gives the coords of a point which is 
                                        # bearing and distance

the_mat <- data.frame(lon = c(35, 35.2, 35.6),
                      lat = c(3, 3.2, 3.6))
the_matrix <- as.matrix(the_mat)

destPoint(p = the_matrix, b = abs(rnorm(3, 100, 30)), d = rnorm(3, 150, 25))
# The end

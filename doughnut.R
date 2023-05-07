# Extracting pseudo-absence from doughnut buffer around presence points

library(sf) # For wrangling simple features o
library(tidyverse) # Manipulating data frames or tibbles
library(charlatan) # For generating random points

presence <- data.frame(lon = ch_lon(n = 20), lat = ch_lat(n = 20)) |> 
  st_as_sf(coords = c("lon", "lat"))

plot(presence)

large_buffer <- st_union(st_buffer(presence, 20))
plot(large_buffer, add = T)

small_buffer <- st_union(st_buffer(presence, 9))
plot(small_buffer, add = T)

doughnut <- st_difference(large_buffer, small_buffer)

plot(doughnut)

pseudo_absence <- st_sample(doughnut, size = 500)

plot(pseudo_absence, add = T)

pseudo_absence |> 
  st_coordinates() |> 
  as.data.frame() |> 
  dplyr::select(lon = X, lat = Y) |> 
  nrow()

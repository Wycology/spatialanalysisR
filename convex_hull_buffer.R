# Library ----------------------------------------------------
pacman::p_load(sf) # For spatial analysis, Version 1.0.14

# Occurrence records -----------------------------------------
species <-
  st_read(system.file("external/species.shp", 
                      package = "sdm")) 

# Subset for presences only ----------------------------------
species <- species[species$Occurrence == 1,] 

# Building convex hull around the points ---------------------
convex_hull <- st_convex_hull(st_union(species)) 

# Getting the coordinates of the hull vertices ----------------
convex_verts <- st_as_sf(
  as.data.frame(st_coordinates(convex_hull)), 
  coords = c("X", "Y"),
  crs = st_crs(convex_hull)
)

# Centroid coordinates of the convex hull----------------------
centroid <- st_centroid(convex_hull) 

# Distances from the centroid to each of the vertices -------------
distances <- st_distance(centroid, convex_verts) 

# Beyond the maximum distance by 10% ----------------------------------
hull_ext <- st_buffer(x = convex_hull, 
                      dist = max(distances) * 0.1) 
# Background points--------------------------------------------
set.seed(248) # The number can be between -2147483648 & 2147483647
bg_points <- st_sample(x = hull_ext, size = 100)

# Visualizing ----------------------------------------------------
plot(st_geometry(hull_ext), axes = TRUE)
plot(convex_hull, add = TRUE, col = "#33a02c")
plot(st_geometry(species), col = "blue", pch = 19, add = TRUE)
plot(st_geometry(bg_points), col = "#e31a1c", pch = 16, add = TRUE)
plot(st_geometry(centroid), col = "yellow", pch = 16, cex = 4, add = TRUE)



# Write a function --------------------------------------------------------

for (i in 1:20) {

conv_calc <- function(species) {
  convex_hull <- st_convex_hull(st_union(species))
  convex_verts <- st_as_sf(
    as.data.frame(st_coordinates(convex_hull)),
    coords = c("X", "Y"),
    crs = st_crs(convex_hull)
  )
  centroid <- st_centroid(convex_hull)
  distances <- st_distance(centroid, convex_verts)
  hull_ext <- st_buffer(x = convex_hull,
                        dist = max(distances) * 0.1)
  bg_points <- st_sample(x = hull_ext, size = 1000)
  outs <- list(
               hull_ex = hull_ext,
               hull = convex_hull,
               cent = centroid, 
               spp = species, 
               bgs = bg_points)
  return(outs)
}

outed <- conv_calc(species)

plot(st_geometry(outed$hull_ex), axes = TRUE)
plot(st_geometry(outed$hull), add = T, col = "#33a02c")
plot(st_geometry(outed$cent), pch = 16, col = "red", cex = 3, add = T)
plot(st_geometry(outed$bgs), pch = 16, col = "blue", add = T)
plot(st_geometry(outed$spp), pch = 16, add = T)
}


new_sf <- st_as_sf(data.frame(lon = rnorm(n = 300, mean = 33, sd = 5.263),
                              lat = rnorm(n = 300, mean = 0, sd = 5.263)),
                   coords = c("lon", "lat"),
                   crs = 4326)


new_outs <- conv_calc(species = new_sf)


pilot <- function(outs) {
  plot(st_geometry(new_outs$hull_ex), axes = TRUE, lwd = 6,)
  plot(st_geometry(new_outs$hull), col = "#33a02c", lwd = 3, add = T)
  plot(st_geometry(new_outs$bgs), pch = 16, col = "blue", add = T)
  plot(st_geometry(new_outs$spp), pch = 16, col = "red", add = T)
  plot(st_geometry(new_outs$cent), pch = 16, col = "yellow", cex = 3, add = T)
  
}

pilot(outs = new_outs)























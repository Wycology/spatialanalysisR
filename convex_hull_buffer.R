# Library ----------------------------------------------------
pacman::p_load(sf) # For spatial, version 1.0.14

# Occurrence records -----------------------------------------
species <-
  st_read(system.file("external/species.shp", 
                      package = "sdm")) 

# Subsetting for presences only -------------------------------
species <- species[species$Occurrence == 1,] 

# Building convext hull around the points ---------------------
convex_hull <- st_convex_hull(st_union(species)) 

# Getting the coordinates of the hull vertices ----------------
convex_verts <- st_as_sf(
  as.data.frame(st_coordinates(convex_hull)), 
  coords = c("X", "Y"),
  crs = st_crs(convex_hull)
)

# Centroid coordinates of the convex hull----------------------
centroid <- st_centroid(convex_hull) 

# Distances from centroid to each of the vertices -------------
distances <- st_distance(centroid, convex_verts) 

# Beyond max distance by 10% ----------------------------------
hull_ext <- st_buffer(x = convex_hull, 
                      dist = max(distances) * 0.1) 
# Background points--------------------------------------------
bg_points <- st_sample(x = hull_ext, size = 100)

# Plotting ----------------------------------------------------
plot(st_geometry(hull_ext), axes = TRUE)
plot(convex_hull, add = TRUE, col = "#33a02c")
plot(st_geometry(species), col = "blue", pch = 19, add = T)
plot(st_geometry(bg_points), col = "#e31a1c", pch = 16, add = T)
plot(st_geometry(centroid), col = "yellow", pch = 16, cex = 4, add = T)


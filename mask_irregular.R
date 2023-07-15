############################################################
# Created by Wyclife Agumba Oluoch                         #
# As a response to an email on 8th March 2021              #
# Last edited on 6th November 2021 while at Ashanzi Lodwar #
############################################################
# Clipping occurrence records data and predictor variables with our study area
# boundary:

# Libraries
# R         version 4.1.2 a.k.a "Bird Hippie"
# RStudio   version RStudio 2021.09.1+372 "Ghost Orchid"

library(dismo)   # version 1.3.5
library(raster)  # version 3.5.2
library(sf)      # version 1.0.4

# Downloading env. data

worldclim <- getData('worldclim', var = 'tmin', res = 10)

# Study area boundary  from the GADM database

kenya <- getData('GADM', country = 'KEN', level = 0)

crs(worldclim) # CRS arguments: +proj=longlat +datum=WGS84 +no_defs 
crs(kenya)     # CRS arguments: +proj=longlat +datum=WGS84 +no_defs 

# Cropping and masking worldclim using the polygon boundary 

clim_cropped <- crop(worldclim, kenya)

# Plot the cropped raster, which is one of the stacked layers

plot(clim_cropped[[1]]) # This is using the rectangular extent of the polygon 

# Now we need to mask the cropped layer using the study area boundary

clim_masked <- mask(clim_cropped, kenya)

# Plot the cropped and now masked raster

plot(clim_masked[[1]]) # This now takes the irregular shape of your polygon and 
# covers only the cropped region

# Occurrence records---- 
# I will download it from gbif. Replace with your species of interest.
# You can also add other parameters within gbif function like ext.

occurrence <- gbif('Salvadora', 'persica', download = T, ntries = 5, sp = T)

# Let's see the distribution of occurrence globally
plot(worldclim[[1]])
plot(occurrence, add = T) # A lot of points largely in Africa but also Arabia

# Cropping occurrence records

occurrence_cropped <- crop(occurrence, kenya) # Only occurrence within kenya extent

# Plot the cropped occurrence on the Kenyan masked worldclim

plot(clim_masked[[1]])
plot(occurrence_cropped, add = T) # Yes cropped, but with kenyan extent.


# Start by converting occurrence_cropped to sf object, since my Kenya is an sf
# object

occurrence_cropped_sf <- st_as_sf(occurrence_cropped)
kenya_sf <- st_as_sf(kenya)
crs(occurrence_cropped_sf) # This confirms it has NO crs yet, so we assign it crs 

occurrence_cropped_sf_geo <- st_set_crs(occurrence_cropped_sf, 4326)

crs(occurrence_cropped_sf_geo) # WGS84
crs(clim_masked)               # WGS84
crs(kenya_sf)                  # WGS84

# Transformed to projected CRS 32636

occurrence_projected <- st_transform(occurrence_cropped_sf_geo, 32636)
kenya_sf_projected <- st_transform(kenya_sf, 32636)
clim_masked_projected <- projectRaster(from = clim_masked, crs = 32636)

crs(occurrence_projected)
crs(kenya_sf_projected)
crs(clim_masked_projected)

# Now, to clip points just within our study area use st_intersection()

occurrence_kenya = st_intersection(occurrence_projected, kenya_sf_projected)

plot(clim_masked_projected[[1]])  
plot(kenya_sf_projected, add = T) # Kenya shapefile, study area
plot(occurrence_kenya,add = T)    # Only points within the study area

# Now, our two data-sets which are within our study area are:

# 1. clim_masked_projected .... for the predictors
# 2. occurrence_kenya ... for the occurrence records.

# All the best!

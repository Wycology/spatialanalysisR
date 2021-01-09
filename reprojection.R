#########################################
# Created by Wyclife Agumba Oluoch      #
# Contacts: wyclifeoluoch@gmail.com     #
# +254729371248                         #
# TASK: Chapter 6 Geocomputation with R #
# Created on 9th January 2021           #
# Last modified on 9th January 2021     #
#########################################

# Loading the necessary libraries 

library(sf); library(raster); library(dplyr); library(spData); library(spDataLarge)

# I need to install the last two packages in the list

# install.packages("spData") # Done seamlessly
# remotes::install_github("Nowosad/spDataLarge") # Took some time, then done

# st_is_longlat is one great function for checking if geographic CRS is in use

# Here we start by a dummy dataframe and check for crs using the function

london <- data.frame(lon = -0.1,
                     lat = 51.5) %>% # Creating a dataframe 
  st_as_sf(coords = c('lon', 'lat')) # Setting the variables lon and lat as coords

st_is_longlat(london) # Returns NA

# Knowledge gained from the above failure is that one cannot just pick a some
# variables and set them as coordinates and that is enough. Setting coordinates
# in sf package needs the use of st_set_crs whereby epsg code is also given. That
# is when the full metadata will be availed to fully define the crs.
# So we can overcome the above problem by:

london_geo <- london %>% # Feeding in the above information
  st_set_crs(4326) # Setting the crs

# To show the whole process in one group of code:

london_geo <- data.frame(lon = -0.1,
                         lat = 51.5) %>%
  st_as_sf(coords = c('lon', 'lat')) %>% 
  st_set_crs(4326)

st_is_longlat(london_geo) # This is now returning TRUE because we have set crs

# If a dataset has no crs, a number of issues can arise. Let us see this when 
# trying to create buffer around the point we created for londin and london_geo

london_buff_nocrs <- st_buffer(london, dist = 1)
london_buff <- st_buffer(london_geo, dist = 1) # Returns a warning, the warning 
# implies that we should reproject the point to a projected coordinate system not 
# the longlat geographic datum.

geosphere::distGeo(c(0,0), c(1,0)) # Returns distance between two longitudes

geosphere::distGeo(c(-0.1,51.5), c(0.9,51.5))/1000 # Cool, distance between two 
# meridians in London is slightly below 70kms (69.43998).

# So here we do the reprojection of our london_geo to a projected CRS 

london_proj <- data.frame(x = 530000,
                          y = 180000) %>% 
  st_as_sf(coords = 1:2, crs = 27700) # Setting projection to OSGB 1936 British
                                      # National Grid

st_is_longlat(london_proj) # Returning FALSE because this is not longlat but a
# projected CRS. 

# Being a projected CRS, the units are now in m other than degrees

# We can check the CRS of our datasets using st_crs as follows:

st_crs(london) # Returns NA
st_crs(london_geo) # Returns EPSG 4326 which is GEOGCRS WGS 84
st_crs(london_proj)# Returns EPSG 4326 which is PROJCRS OSGB 1936 

# I cannot overemphasize the difference between geographic and projected CRSs
# Master the concept please and life will be smoother henceforth.

# Now we can do the buffer easily, using the distance between two meridians at 
# the equator which we got earlier as 111319.5m

london_proj_buff <- st_buffer(london_proj, dist = 111319.5)
str(london_proj_buff)

plot(london_proj_buff, col = 'purple')
plot(london_proj, col = 'yellow', add = TRUE)

# I know the distance between Kisumu and Kericho on google map to be 64.89 kms

kisumu <- data.frame(lon = 34.767439,
                     lat =  -0.090972) %>% 
  st_as_sf(coords = c('lon', 'lat')) %>% 
  st_set_crs(4326)

kericho <- data.frame(lon = 35.286047,
                      lat = -0.370453) %>% 
  st_as_sf(coords = c('lon', 'lat')) %>% 
  st_set_crs(4326)

st_distance(kisumu, kericho)/1000 # This is returning 65.48168 kms. Cool



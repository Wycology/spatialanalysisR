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

st_distance(kisumu, kericho)/1000 # This is returning 65.48168 kms. Cool. Even 
# though it is saying Units: [m], I know very well that it is kms due to division
# by 1000
st_distance(london_geo,london_proj) # Very informative error, both objects have 
# different crs.

# Let me transform the one in GEOCRS to PROJCRS

london2 <- st_transform(london_geo, 27700)
# Now this is cool, let me try calculating the distance again

st_distance(london2, london_proj) # Good, the two are about kms apart.

# This is an interesting one, getting length of met station at TRI. Manually in 
# Google Earth, I got it as 30.86 m.

north <- data.frame(lon = 35.349008,     # Upper corner near sunshine recorder
                    lat =  -0.371926) %>% 
  st_as_sf(coords = c('lon', 'lat')) %>% 
  st_set_crs(4326) # The GEOCRS since I have longitude and latitude

south <- data.frame(lon = 761462.23,     # Upper southern corner. Coords from 
                    lat = 9958825.90) %>% # Google Earth after changing display
  st_as_sf(coords = c('lon', 'lat')) %>%  # to UTM under Tools ==> Options
  st_set_crs(32736) # This is the EPSG code for UTM Zone 36 S. Google it.

st_distance(north, south) # Error, the two have different crs, EPSG 4326 and 
                          # EPSG 32736

north2 <- st_transform(north, 32736) # setting the GEOCRS to PROJCRS
st_distance(north2, south) # Super cool 31.37634 m, so close to what I got (30.86 m)
# I actually missed by half a metre.

31.37634 - 30.86

# Worried about EPSG code associated with any point on the planet, here is a tailor
# made function:

lonlat2UTM <- function(lonlat){
  utm <- (floor((lonlat[1] + 180)/6) %% 60) + 1
  if(lonlat[2] > 0){
    utm + 32600
  } else{
    utm + 32700
  }
}

# We can now get the utm zone and associated EPSG code for Auckland, London,
# Nairobi, TRI, and Lodwar

epsg_utm_aukland <- lonlat2UTM(c(174.7, -36.9)) # Returns 32760
epsg_utm_london <- lonlat2UTM(c(st_coordinates(london))) # Returns 32630
epsg_utm_Nairobi <- lonlat2UTM(c(36.822004,  -1.293444)) # Returns 32737
epsg_utm_tri <- lonlat2UTM(c(35.348981, -0.372091)) # Returns 32736
epsg_utm_lodwar <- lonlat2UTM(c(35.597457, 3.118125)) # Returns 32636
epsg_utm_tannenbusch <- lonlat2UTM(c(7.043103, 50.749732)) # Returns 32632
epsg_utm_honduras <- lonlat2UTM(c(-86.241367, 15.200515)) # Returns 32616

# Now, to know the zones associated with the above, it is simple just extract it
# with $ accessor symbol

st_crs(epsg_utm_lodwar)$proj4string # Zone 36 and the units are in m
st_crs(epsg_utm_honduras)$proj4string # Zone 16 and the units are in m

# Cool, can go on and on and on and on and on and on...







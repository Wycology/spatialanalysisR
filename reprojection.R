#########################################
# Created by Wyclife Agumba Oluoch      #
# Contacts: wyclifeoluoch@gmail.com     #
# +254729371248                         #
# TASK: Chapter 6 Geocomputation with R #
# Created on 9th January 2021           #
# Last modified on 9th January 2021     #
#########################################
# Chapter Six Reprojecting Geographic Data----
# Loading the necessary libraries 

library(sf)
library(raster)
library(dplyr)
library(spData)
library(spDataLarge)

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
london_buff <-
  st_buffer(london_geo, dist = 1) # Returns a warning, the warning
# implies that we should reproject the point to a projected coordinate system not 
# the longlat geographic datum.

geosphere::distGeo(c(0, 0), c(1, 0)) # Returns distance between two longitudes

geosphere::distGeo(c(-0.1, 51.5), c(0.9, 51.5)) / 1000 # Cool, distance between two
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
# Reprojecting vectors

crs_lnd <- st_crs(cycle_hire_osm)
class(crs_lnd)

crs_lnd$epsg # To get the epsg code for the object

crs_lnd$proj4string # This is giving the lengthy format of the code

# Now we know that the crs of our dataset is in GEOCRS, we can change this to
# PROJCRS as follows

cycle_hire_osm_projected <- st_transform(cycle_hire_osm, 27700)

# We can check this and confirm the change
st_crs(cycle_hire_osm_projected)$proj4string # Good. To know the name of the crs
# Just google it. Something like EPSG 27700 or seek from https://epsg.io/27700
# The website is cool even for https://epsg.io/32736 which was for tri
# There is also http://spatialreference.org to check
# which is the same as:
st_crs(cycle_hire_osm_projected)$epsg

# Cool, even just printing the spatial object in the console alone is enough to
# reveal the crs.

print(cycle_hire_osm_projected) # Projected CRS: OSGB 1936 / Britich National Grid

world_mollweide <- st_transform(world, crs = '+proj=moll') # Leaving space on
# either sides of equal sign will throw an error.

plot(world_mollweide) # Plotting all the variables in the dataset

plot(world_mollweide['type'], main = 'Mollweide Projection') # Plot only the 
# type column in the dataset

# There is a projection called Winkel triple projection which can be very useful 
# when trying to get the best in plotting the whole world to preserve area, direction
# and distance

world_wintri <-
  lwgeom::st_transform_proj(world, crs = '+proj=wintri')
plot(world_wintri['type'], main = 'Wintri Projection')

# Another interesting one is to transform the a projection to meet one's need
world_laea1 <-
  st_transform(world, crs = '+proj=laea +x_0=0 +y_0=0 +lat_0=0')
plot(world_laea1['type'], main = 'Lambert Azimuthal Equal Area')

# I can enjoy centering the map on New York City using the +lon_0 and +lat_0 parameters

world_laea2 <-
  st_transform(world, crs = '+proj=laea +x_0=0 +y_0=0 +lon_0=-74 +lat_0=40')
plot(world_laea2['type'], main = 'Lambert Azimuthal Equal Area, New York City')

# Why not center it on Lodwar, cool and interesting, at least to me.
world_laea3 <- st_transform(
  world, crs = '+proj=laea +x_0=0 +y_0=0 +lon_0=35.597456 +lat_0=3.118122')
plot(
  world_laea3['type'],
  graticule = TRUE,
  key.pos = NULL,
  axes = TRUE,
  main = 'Lambert Azimuthal Equal Area, Lodwar'
)

# Reprojecting raster files is a bit tricky. Many things can change including 
# number of rows and columns and cell values. Therefore, care must be taken

# Lets load some raster data which contains land use classes

cat_raster <-
  raster(system.file('raster/nlcd2011.tif', package = 'spDataLarge'))
crs(cat_raster)     # Cool
st_crs(cat_raster)  # Cooler
cat_raster          # Coolest

# To know the number of land cover classes in the dataset, we can have a look at:
unique(cat_raster) # Good
length(unique(cat_raster)) # Ooh, they are 14 classes

# Now we can reproject this raster of land cover classes to WGS84 using ngb method.
# ngb is nearest neighbor so that the program will look for neighboring cells to 
# assign cell values in the new file. It is good for catgorical raster like this\
# case. Otherwise for numeric we would use billinear which is averaging values of
# four neighboring cells. By the way, reprojecting a raster leads
# to creation of a new raster which can have different extents, number of cells and 
# resolution as the original one.

# It normally requires that we give the full text name of the crs not the epsg code
# However, epsg code can be used by tweaking as "+init=epsg:MY_NUMBER", let me try both

wgs84 <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'

cat_raster_wgs84 <-
  projectRaster(cat_raster, crs = wgs84, method = 'ngb')
plot(cat_raster)       # Clear
plot(cat_raster_wgs84) # Difference noted

cat_raster             # Checking the original raster attributes
cat_raster_wgs84       # Checking the new raster attributes

# There are notable differences in number of rows and columns.

# Let me try reprojecting using the tweaking +init blah blah method

cat_raster_wgs84_init <-
  projectRaster(cat_raster, crs = '+init=epsg:4326',
                method = 'ngb')
cat_raster_wgs84_init # Running this raster returns attributes of the raster quite
# similar to those of cat_raster_wgs84

# To compare the two rasters we can use compareRaster function as follows

compareRaster(cat_raster_wgs84_init, cat_raster_wgs84) # This returns TRUE, so the
# tweaking works fine.

# Now let us finish the job by reprojecting numeric raster. That is raster with
# double numbers (e.g. 12.3526) as opposed to categorical one which only had the
# categories as integers in the cell values. This will be an SRTM raster

con_raster <-
  raster(system.file('raster/srtm.tif', package = 'spDataLarge'))
# con here means continuous and the cat in the former meant ctegorical

crs(con_raster)    # Again Cool
st_crs(con_raster) # Cooler
con_raster         # Coolest

# Here we will reproject the GEOCRS to PROJCRS called equalarea

equalarea <- '+proj=laea +lat_0=37.32 +lon_0=-113.04'
con_raster_ea <-
  projectRaster(con_raster, crs = equalarea, method = 'bilinear')

crs(con_raster_ea)
st_crs(con_raster_ea)
con_raster_ea

# Again I can clearly see that the rasters have different number of rows and columns
# Even the minimum and maximum cell values changed.

# Okay, that is enough for the day.


# Chapter Seven Geographic data I/O----
 
# Here, we will use same libraries loaded at the beginning of Chapter 6 above

# Learning to import and export data into R. It is good to know both so that 
# you can load data from others and share out data which others can also utilize

# There are so many geographic data formats with pros and cons

# Retrieving open data from a number of geoportals such as:
# GEOS Portal http://www.geoportal.org/
# Copernicus Open Access Hub: https://scihub.copernicus.eu/
# NASA SEDAC Portal: http://sedac.ciesin.columbia.edu/
# INSPIRE: http://inspire-geoportal.ec.europa.eu/
# USGS EarthExplorer: https://earthexplorer.usgs.gov/

# It is rewarding to explore these data on the web to know what is available. 
# However, downloading data should be done by code to enable reproducibility.

# Data hosted online can be downloaded via url (download.file()) or specific 
# API such as Sentinel API: 
# https://scihub.copernicus.eu/twiki/do/view/SciHubWebPortal/APIHubDescription

# A simple code to access US National Parks from https://catalog.data.gov/dataset/national-parks
# is

download.file(url = "http://nrdata.nps.gov/programs/lands/nps_boundary.zip",
              destfile = "nps_boundary.zip") # Unfortunately couldn't open.
unzip(zipfile = "nps_boundary.zip")
usa_parks = st_read(dsn = "nps_boundary.shp")

# Some (few) of the packages which can help with accessing geodata include:

# getlandsat ==> Accessing Landsat 8 data
# osmdata ==> Download and import of OpenStreetMap data
# raster ==> getData() imports administrative, elevation, and WorldClim data
# rnaturalearth ==> Access to Natural Earth vector and raster data
# rnoaa ==> imports National Oceanographic and Atmospheric Administration climate data
# rWBclimate ==> Access World Bank climate data
# GSODR ==> Accessing climate data for given stations

# I will use the package GSODR to play around with the data. Here is how to access
# data for the region around Tea Research Institute, Kericho Kenya

tri_station <- nearest_stations(LAT = -0.36667, LON = 35.35, distance = 1)
# This can help to return a station code which is applicable in accessing data
# Returned a code of 637100-99999 hence used as follows 

tri_climate <- get_GSOD(years = 2010, station = "637100-99999")# Takes too long
# I escaped this to stop running it before completion.

# I will use three of the packages listed above:

# 1. rnaturalearth

library(rnaturalearth)
usa <- ne_countries(country = 'United States of America')
class(usa)

# By default, this is a Spatial feature. To convert to sf is simple.

usa_sf <- st_as_sf(usa)
class(usa_sf)
head(usa_sf) # This is now better.

# Why not run one for my beloved country
kenya <- ne_countries(country = 'Kenya')
class(kenya)

kenya_sf <- st_as_sf(kenya)
class(kenya_sf)
head(kenya)
plot(kenya_sf) # Of course this maps every attribute in Kenyan dataset
plot(kenya_sf['type']) # Picks only one of the attributes
plot(kenya_sf[c('type', 'adm0_dif')]) # Picks more than one attribute to map

#2. raster

library(raster)
worldclim_prec <- getData(name = 'worldclim', var = 'prec', res = 10)
class(worldclim_prec) # This is a RasterStack class object.

plot(worldclim_prec) # Plots the 12 months data
plot(worldclim_prec[[1]]) # Plots only the first layer.
plot(worldclim_prec[[1:2]]) # Plots two of the layers
plot(worldclim_prec[[c(1,4,6)]]) # Plots selected layers in the stack


#3. osmdata

library(osmdata) # Loading the library

parks <- opq(bbox = 'leeds uk') %>% 
  add_osm_feature(key = 'leisure', value = 'park') %>% 
  osmdata_sf()

plot(parks$osm_polygons) # There are 57 attributes, only first nine are plotted
plot(parks$osm_polygons['osm_id']) # This is the same as indexing by number
plot(parks$osm_polygons[1]) # Exactly

#plotting more

plot(parks$osm_polygons[1:4]) # Only four
plot(parks$osm_polygons[3:6]) # Picking third to sixth
plot(parks$osm_polygons[c(8,3,40,57)]) # Picking random across the attributes

dim(parks$osm_polygons) # So here we have 745 rows and 58 columns

# Data which come with packages can be accessed by attaching the dataset, by 
# running data(dataset), by calling the data from package name pkg::dataset, or 
# with system.file(). We can ccess data which come with spData by:

world2 <- spData::world
# or
world3 <- st_read(system.file('shapes/world.gpkg', package = 'spData'))

# File formats ----

# So many formats exist, am falling in love with geopackage, and GeoTiff .tif/.tiff

# I did not know limitations of shapefile before: has multifiles, at least 3; 
# no more than 255 columns; no more than 10 characters on column names; 
# maximum file size of 2GB and does not distinguish between polygons and 
# multipolygons (it therefore does not support all possible geometries).

# Nothing feels better than geopackage... such a cool format burrying shapefile

# Data input ----

# Reading vector ==> sf::st_read()
# Reading raster ==> raster::raster()
# Reading csv ==> dplyr::read_csv()

# These help to assign objects to your workspace, stored in RAM accessible from 
# the GlobalEnvironment of the R session

# Vector data
# sf supports so many vector file formats. These can be assigned by st_read() 
# To know the whole list:
st_drivers() # Quite a long list. Almost every format you can think of

# In reading files using st_read(), there is no need to provide driver name since
# the function can tell the driver to use based on the extension on the file name
 
# For example:

vector_file_path <- system.file('shapes/world.gpkg', package = 'spData')
world <- st_read(vector_file_path) # Simple way to specify path and read vector data

# Adding layer argument can be good if there are many layers in the destination 

# Reading in csv data with x and y coordinates columns

cycle_hire.txt <- system.file('misc/cycle_hire_xy.csv', package = 'spData')
cycle_hire.xy <- st_read(cycle_hire.txt, options = c('X_POSSIBLE_NAMES=X',
                                                     'Y_POSSIBLE_NAMES=Y'))
cycle_hire.xy # Will imported but no crs yet.

testingcrs <- st_set_crs(cycle_hire.xy, 4326) # Joking with assignment here
crs(testingcrs) # Confirming the crs set
plot(testingcrs) # Plotting every attribute of the dataset
plot(testingcrs['nbikes']) # Picking what to plot

# In some instances, for example with polygons, several coordinates can be stored
# in a single cell just like geometry column of sf package. These are well known
# text or binaries. These can be read as follows:

world_txt <- system.file('misc/world_wkt.csv', package = 'spData')
world_wkt <- read_sf(world_txt, options = 'GEOM_POSSIBLE_NAMES=WKT')

# This will do the same thing as the one which follows:

world_wkt2 <- st_read(world_txt, options = 'GEOM_POSSIBLE_NAMES=WKT',
                     quiet = TRUE, stringsAsFactors = FALSE, as_tibble = TRUE)
world_wkt == world_wkt2 # Everything TRUE...

plot(world_wkt2) # Now 10 are plotted
plot(world_wkt2[5]) # The same applies as has been with the other vectors with 
# very many attributes

# The next is to get KML file and read it

u <- 'https://developers.google.com/kml/documentation/KML_Samples.kml'
download.file(u, 'KML_Samples.kml')
st_layers('KML_Samples.kml') # Checking the number of layers in the file downloaded

kml <- read_sf('KML_Samples.kml', layer = 'Placemarks')
plot(kml$geometry) #Plots the three points

# For reading geojson files, one better goes for the package geojsonsf

# Raster Data ----

raster_filepath = system.file("raster/srtm.tif", package = "spDataLarge")
single_layer = raster(raster_filepath)

# In case the raster is multiband and only specific band is need, use band

multilayer_filepath = system.file("raster/landsat.tif", package = "spDataLarge")
band3 = raster(multilayer_filepath, band = 3)

# However, in case of reading all the available bands one can use stack or brick

multilayer_brick = brick(multilayer_filepath)
multilayer_stack = stack(multilayer_filepath)

# Working almost in similar manner.

# Data Output ----
# Vector data ==> To write vector, use st_write(object, dsn)

st_write(obj = world, dsn = 'world.gpkg') # As simple as such

# Writing to the same file is not possible, throws an error

st_write(obj = world, dsn = 'world.gpkg') # Error

# Since we are writing gpkg, it can contain many layers. To achieve this, we set
# append = TRUE
st_write(obj = world, dsn = 'world.gpkg', append = TRUE) # This works fine.

# To overwrite the former layer, just set append to FALSE, simple.

# Instead of st_write, there is write_sf which has append set to FALSE and quiet 
# set to TRUE by default. This is therefore overwriting the older layers whenever
# additional layers of same name are added.

write_sf(obj = world, dsn = 'world.gpkg') # OOooh, this takes longer time.

# In order to write the vector data out to text file such as csv or txt then 
# layer_options argument can be used as follows. This drops the geometry column

st_write(cycle_hire.xy, 'cycle_hire_xy.csv', layer_options = 'GEOMETRY=AS_XY')
st_write(world_wkt, 'world.wkt.csv', layer_options = 'GEOMETRY=AS_WKT')

# Raster data
# Adding desired extension is important. By default, writeRaster() will store .grd
writeRaster(single_layer, filename = "my_raster.tif", datatype = "INT2U")

# To know the supported formats just run:
writeFormats() # Quite a long list

# Other parameters can be included accordingly
writeRaster(x = single_layer,
            filename = "my_raster.tif",
            datatype = "INT2U",
            options = c("COMPRESS=DEFLATE"),
            overwrite = TRUE) # overwrite here works quite similar to append = TRUE

# Visual outputs ----
# Static graphics
# This should start by opening graphic device, then create plot, lastly close it.

png(filename = 'lifeExp.png', width = 500, height = 350) # Open graphic device
plot(world['lifeExp'])                                   # The plot
dev.off()                                                # Closes

# Apart from png, one can use pdf(), bmp(), jpeg(), and tiff()
# width, height, and resolution can be specified

# In case of tmap package, tmap_save() can be used.

library(tmap) # Just to run an example
tmap_obj <- tm_shape(world) + tm_polygons(col = 'lifeExp')
tmap_save(tm = tmap_obj, filename = 'lifeExp_tmap.png')

# Interactive maps using mapview 
library(mapview)
mapview_obj <- mapview(world, zcol = 'lifeExp', legend = TRUE)
mapshot(mapview_obj, file = 'my_interactive_map.html') # This threw an eror

mapview(world, zcol = 'lifeExp', legend = TRUE) # This worked pretty fine

# Great enough for today, tomorrow it will be Chapter 8.

# Chapter 8: Making maps with R ----
# This is one of the most loved chapters of the technology
# It will also rely on the packages loaded at Chapter Six as well as:
library(tmap)
library(leaflet)
library(ggplot2)
library(shiny)

# Map making - the art of cartography - is an ancient skill that involves 
# communication, intuition, and an element of creativity.
# Static maps
# By large, tmap package will be used

# Add fill layer to nz shape
tm_shape(nz) + 
  tm_fill()   # Gray plot of the map with bounding box

# Add border layer to nz map

tm_shape(nz) +
  tm_borders() # Draws the borders of the polygon

# Combine the above two

tm_shape(nz) +
  tm_fill() +  # Fills the shape
  tm_borders() # Adds border lines to the filled plot

# Apart from tm_fill and tm_borders, there are also others like
help('tmap-element') # tm_polygons, tm_symbols, tm_raster etc, for example:

tm_shape(nz) +
  tm_polygons() # This achieves the fill and borders in one

# Rapid one would be to use qtm (quick thematic maps)

qtm(nz) # Also achieves the same result as above, can also be added to 

qtm(nz) +
  qtm(nz_height) # Adding this layer on the earlier

# Storing map in a tmap object 
map_nz <- tm_shape(nz) +
  tm_polygons()

class(map_nz) # Confirms the class of the object created.

print(map_nz) # Plots the map
map_nz # Also plots the map    

# Additional map layers can be added on this map object in a simple manner

map_nz1 <- map_nz +    # The original map object 
  tm_shape(nz_elev) +  # Adding the raster elevation object on top
  tm_raster()          # For plotting the raster

map_nz1 # Returns the new object on the plot















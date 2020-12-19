# Working with raster and vector data in R----
# Calculating vegetation canopy heights across NEON sites in Harvard
# The following libraries need to be loaded:

library(raster) # Reading and working with raster files (e.g. dsm and dtm)
library(rgdal)  # Reading and manipulating vector data (points)
library(tidyverse) # Wrangling and plotting (visualizing) the data
library(maps) # Loading us map and having it as dataframe
library(spocc) # for gathering species records from gbif (https://www.gbif.org/)

# Loading the dsm (digital surface model) data - heights of top physical points
# Data were captured by lidar flyover in Harvard
# Typically trees but could be anything else projecting above the surface
# When I checked the study area using Google Earth (42.536910°, -72.17265°) 
# (https://www.neonscience.org/field-sites/harv), it is a forested site 
# and the sensor is clearly visible!!!

dsm_harvard <- raster("NEON-airborne/HARV_dsmCrop.tif") # Reads the raster data

# Have a quick visual appeal of the data (raster).

raster::plot(dsm_harvard) # Specifying raster:: helps to pick the right plot.

# Check the metadata of the dsm. Done by running the raster object itself

dsm_harvard # Looking closely at the source indicates it is a .tif file.
# The resolution is 1 by 1 meaning 1m by 1m very high spatial resolution!!
# Am sure the unit is m because it is indicated under the crs:....+units=m
# Convert the dsm raster (GeoTiff) to dataframe to plot using ggplot2 package
# This can be a very big file depending on size of the raster.

dsm_harvard_df <- as.data.frame(dsm_harvard, 
                                xy = TRUE) # xy true ensures that coordinates
                                           # are also stored in the created 
                                           # dataframe

# Check the head of the dataframe created to confirm the data is as expected

head(dsm_harvard_df) # In deed the coordinates are also stored in the
# newly created dataframe.

# We can create another one with xy set to FALSE.

no_xy_df <- as.data.frame(dsm_harvard, xy = FALSE)

head(no_xy_df) # This is only giving back the dsm values without xy coordinates.

# This is a very lengthy dataframe, check the length of one column

length(dsm_harvard_df$x) # This (2319799) is the same value as ncell shown by:

dsm_harvard # Weare now quite sure about how long the dataset is.

# Create the ggplot "HARV_dsmCrop" is the name of the column with cell values

ggplot(data = dsm_harvard_df, 
       aes(x = x, y = y,
           fill = HARV_dsmCrop))+
  geom_raster() + # This use of geom_raster() is important to note
  labs(x = "Longitude (m)", 
       y = "Latitude (m)", # Checking the axes, the values are projected utm meters
       title = "NEON-dsm of Harvard") 

# Importing dtm data, note that this is d t m, not d s m.

dtm_harvard <- raster("NEON-airborne/HARV_dtmCrop.tif") # Note the difference in 

# the name of the file as dtm not dsm, I repeat.

# Subtracting dtm from dsm to get the true canopy height. Difference in the 
# heights of the two layers.

canopy_height_harvard <- dsm_harvard - dtm_harvard # Simple raster subtraction

# Convert canopy height raster to dataframe for plotting

canopy_height_harvard_df <- as.data.frame(canopy_height_harvard, 
                                          xy = TRUE)
# Checking the head of the canopy dataframe.

head(canopy_height_harvard_df) # x, y, and layer well displayed. Layer here 
# indicates the values for the heights (m) of the trees within Harvard.

# Generating the plot of the canopy data

ggplot(data = canopy_height_harvard_df,
              aes(x = x, y = y, 
                  fill = layer )) +
  geom_raster() +
  labs(x = "Longitude (m)", 
       y = "Latitude (m)",
       title = "NEON-canopy height of Harvard")

# It is interesting that some of the trees grow above 30 m, quite tall!!

# Checking the height of the tallest tree (38.16998 m). Wow!

# Some few summaries of the data. 
max(canopy_height_harvard_df$layer, na.rm = TRUE)
min(canopy_height_harvard_df$layer, na.rm = TRUE)
mean(canopy_height_harvard_df$layer, na.rm = TRUE)
median(canopy_height_harvard_df$layer, na.rm = TRUE)
quantile(canopy_height_harvard_df$layer, 0.25, na.rm = TRUE) # Checking the 
# heights of trees at the 25th percentile.

summary(canopy_height_harvard_df$layer, na.rm = TRUE) # Simple and clear.

ggplot(data = canopy_height_harvard_df,
       aes(layer)) +
  geom_histogram(col = "blue", fill = "purple") +
  labs(title = "Distribution of tree heights in Harvard",
       x = "Tree heights (m)",
       y = "Frequency") +
  theme_classic() # Generating histogram to show the distribution of tree heights

hist(canopy_height_harvard_df$layer, col = "purple") # Assigning my favorite color

# Working with Vector data----
# Loading the points data of some plots in Harvard
# You do not need to specify the file extension .shp in the readOGR function
plots_harvard <- readOGR(dsn = "NEON-airborne/plot_locations", 
                         layer = "HARV_plots")
class(plots_harvard) # This is a SpatialPointsDataFrame class object

plot(plots_harvard) # Plots the five points 

head(plots_harvard) # The coordinates are hidden but plot ids displayed.

# It is time to plot the points on top of the canopy height layer.

# Both layers must be in the same coordinate reference system

crs(canopy_height_harvard) # Checking the crs of the canopy height raster
crs(plots_harvard) # Checking the crs of the points layer

# I have realized that they are different and I want the crs of canopy height
# to be the one I use because it is a projected utm and has units in m.
# That is better for small localized mapping as is in this case.
# I use spTransform function to accomplish the transformation

plots_harvard_utm <- spTransform(x = plots_harvard, 
                                 CRSobj = crs(canopy_height_harvard))

# The spTransform() is just too powerful. Can't imagine doing that in QGIS

crs(canopy_height_harvard)
crs(plots_harvard_utm)

# Now they have the same crs, this means they can now be overlaid precisely

# The next step is to convert the plots_harvard_utm to dataframe

plots_harvard_utm_df <- as.data.frame(plots_harvard_utm) # Here we don't add xy = TRUE
# because the data is point shapefile and coordinates are included by default 
# (they are what is being converted to dataframe).

head(plots_harvard_utm_df) # Checking the first 6 rows of the dataframe

# Finally drawing the plots of both raster and points
# Using geom_raster and geom_point

ggplot() +
  geom_raster(data = canopy_height_harvard_df,
              aes(x = x, y = y, fill = layer)) + # For the raster layer
  geom_point(data = plots_harvard_utm_df,
              aes(x = coords.x1, y = coords.x2), # For the points layer
              col = "yellow", cex = 4)
  
# Extracting raster values at point locations ----

plots_canopy_height <- raster::extract(x = canopy_height_harvard, 
                               plots_harvard_utm) # could also provide buffer and 
# fun as mean to get values at some pixels around the points. But for now the 
# extraction is done at the exact points.

# Extracted points are in the same order as the id of the plots

head(plots_canopy_height) # Provides a list of extracted values for the five 
                          # points

plots_harvard_utm$plot_id # Checking the presence of plot_id

# To build a dataframe for the extracted heights.

plots_canopy_height_df <- data.frame(plot_num = plots_harvard_utm$plot_id,
                                     plot_value = plots_canopy_height)

head(plots_canopy_height_df)

# In order to get average of cells around a point as the extracted value we use
# buffer argument

plots_canopy_height <- raster::extract(canopy_height_harvard, plots_harvard_utm, 
                                       buffer = 5, 
                                       fun = mean)

head(plots_canopy_height) # The second point showed a big change 
# from 15.94 to 9.93 when buffer is used. Probably it was a lonely tree 
# in a rather bare background.

# Mapping points data, kind of species occurrence ----

us_map <- ggplot2::map_data(map = "usa") # map_data is a function within ggplot2.

# The name "usa" is provided by the maps package.

ggplot() +
  geom_polygon(data = us_map,
               aes(x = long, y = lat,
                   group = group),
               fill = 'grey') # This is plotting a crazy map

ggplot() +
  geom_polygon(data = us_map,
               aes(x = long, y = lat,
                   group = group),
               fill = 'grey') +
  coord_quickmap() # Quickmap is helping to better project the map

# Adding species occurrence data on the map.

species_gbif <- spocc::occ(query = "Dipodomys ordii",
                    from = "gbif",
                    limit = 1000,
                    has_coords = TRUE) # Returns a max of 1k records with coords.

species_gbif_df <- data.frame(species_gbif$gbif$data)
head(species_gbif_df$Dipodomys_ordii.longitude)


ggplot() +
  geom_polygon(data = us_map,
               aes(x = long, y = lat,
                   group = group),
               fill = 'orange') +
  geom_point(data = species_gbif_df,
             aes(x = Dipodomys_ordii.longitude,
                 y = Dipodomys_ordii.latitude),
             col = "purple") +
  coord_quickmap()



# Working with raster data in R----

########################################
# Created by Wyclife Agumba Oluoch     #
# Contacts: wyclifeoluoch@gmail.com    #
# +254729371248                        #
# TASK: Working with spatial data in R #
# Created on 28th Dec 2020 4.24 pm     #
# Last modified on 15th Oct 2024       #
########################################

# Calculating plant canopy heights across NEON sites in Harvard using dsm & dtm.

# Here are a bunch of necessary libraries.

library(raster) # Working with raster data files (e.g. dsm and dtm)
library(rgdal)  # Reading and manipulating vector data (points, lines, polygons)
library(maps) # Loading required maps as SpatialDataframe
library(spocc) # Gathering species records from gbif (https://www.gbif.org/)
library(tidyverse) # Wrangling and visualizing especially tabular data
library(mapview) # Viewing interactive map with a range of background layers
library(leaflet) # Almost similar functionalities as mapview
library(leaflet.extras2) # For leaflet functionalities in R
# Versions of packages, Rstudio and R used in this project:

# R               = 4.0.5 "Shake and Throw", run the word version 
# RStudio         = 1.4.1103 (Go to Help ==> About RStudio)
# raster          = 3.4.5  Run packageVersion("raster")
# rgdal           = 1.5.23 Run packageVersion("rgdal")
# tidyverse       = 1.3.0  Run packageVersion("tidyverse")
# maps            = 3.3.0  Run packageVersion("maps")
# spocc           = 1.2.0  Run packageVersion("spocc")
# leaflet.extras2 = 1.1.0  Run packageVersion('leaflet.extras2')

# Loading the dsm (digital surface model) data - heights of top physical points.
# The data were captured by lidar flyover in Harvard NEON site.
# Typically trees but could be anything else projecting above earth's surface.
# When I checked the study area using Google Earth (42.536910°, -72.17265°) 
# (https://www.neonscience.org/field-sites/harv), it is a forested site 
# and even the NEON sensor pylon is clearly visible at the site!!!

dsm_harvard <- raster("NEON-airborne/HARV_dsmCrop.tif") # Reads the dsm data.

# Have a quick visual appeal of the loaded raster data by plotting it.

raster::plot(dsm_harvard) # Specifying raster:: helps to pick the right plot().

# Metadata of the dsm. Accessed by running the raster object itself.

dsm_harvard # Displays the attributes or metadata of the raster dsm
ncol(dsm_harvard) # Checking the number of columns in the data.frame.

# Looking at the 'source' shows it is a .tif data file.
# The resolution is 1 by 1 meaning 1 m by 1 m, very high spatial resolution!!
# Am sure the unit is m because it is indicated under the crs:....+units=m.
# Convert the dsm raster (GeoTiff) to data.frame to plot using ggplot2 package
# This can be a very big file depending on size of the raster. Only do it if
# the raster is small enough to allow and/or the processing capacity of your
# machine.

dsm_harvard_df <- as.data.frame(dsm_harvard, 
                                xy = TRUE) # xy = TRUE ensures that coordinates
                                           # are also stored in the created 
                                           # data.frame.

# Check the head of the data.frame created to confirm the data is well created.

head(dsm_harvard_df) # Indeed the coordinates are also stored in the newly 
# created data.frame.

# We can create another one with xy set to FALSE, for fun.

no_xy_df <- as.data.frame(dsm_harvard, xy = FALSE)

head(no_xy_df) # This is only giving back the dsm values without xy coordinates.

dsm_harvard_df %>% # Pick the data-set, and then.
  dplyr::select(y) %>% # pick the column called y, and then
  nrow() # Give the number of rows here. This can also be x or HARV_dsmCrop if
# they are accordingly included in the select() too.

View(dsm_harvard_df) # To have a look at the whole data.frame.

# This (2319799) is the same value as ncell shown by:

dsm_harvard # We are now sure about how long (number of rows) the dataset is.
# @ symbol can be used to extract only specific attribute of the layer e.g
dsm_harvard@crs # coordinate reference system of the dsm raster
dsm_harvard@extent # extent of the dsm layer
dsm_harvard@ncols # number of columns in the raster
dsm_harvard@nrows # number of rows in the raster layer

# However, for the number of cells, I use:

ncell(dsm_harvard) # Returning the number of cells in the raster dataset. Same as

dsm_harvard@ncols * dsm_harvard@nrows # Multiplying columns by rows

# Create the ggplot; "HARV_dsmCrop" is the name of the column with cell values.

ggplot(dsm_harvard_df, 
       aes(x, y,
           fill = HARV_dsmCrop))+ # fill takes the extracted cell values.
  geom_raster() + # This use of geom_raster() is important to note.
  labs(x = "Longitude (m)", 
       y = "Latitude (m)", # Checking the axes, the values are projected utm meters
       title = "NEON-dsm of Harvard") 

# Importing dtm data, note that this is d t m, not d s m.

dtm_harvard <- raster("NEON-airborne/HARV_dtmCrop.tif") # Note the difference in 

# The name of the file as dtm not dsm, I repeat, unfortunately.

# Subtracting dtm from dsm to get the true canopy height. Difference in the 
# heights of the two layers. It is dtm being subtracted FROM dsm.

# However, before subtraction is done, we shall check whether the two rasters
# are of the same extent, resolution, and projection

compareRaster(dtm_harvard, dsm_harvard) # This has returned TRUE, meaning they
# are of the same extent, resolution and projection. Otherwise, it would return
# a self-explanatory error message.

canopy_height_harvard <- dsm_harvard - dtm_harvard # Simple raster subtraction.

# The subtraction is being done cell by cell hence ease of subtraction.
# Convert canopy height raster to data.frame for plotting.

canopy_height_harvard_df <- as.data.frame(canopy_height_harvard, 
                                          xy = TRUE) # Again returning xy.
# Checking the head of the canopy data.frame.

head(canopy_height_harvard_df) # x, y, and layer well displayed. Layer here 
# indicates the values of the heights (m) of the trees within Harvard.
# Here the values are saved under column called layer.

# Generating the plot of the canopy data.

summary(canopy_height_harvard_df) # Just the normal summary of a base R.

newcanopy <- canopy_height_harvard_df %>% 
  mutate(heights = case_when(layer == 0 ~ 'Bare',
                             layer < 10.64 ~ 'Grass',
                             layer < 16.65 ~ 'Shrubs',
                             layer >= 16.65 ~ 'Trees'))

ggplot(data = newcanopy,
              aes(x = x, y = y, 
                  fill = layer )) + # Here we use layer column for fill.
  geom_raster(aes(fill = heights)) +
  labs(x = "Longitude (m)", # Adding the labels on the plot for easy interpretation.
       y = "Latitude (m)",
       title = "NEON-canopy height of Harvard")

ggplot(data = newcanopy, mapping = aes(x = heights, y = layer, fill = heights)) +
  geom_boxplot()

# It is interesting that some of the trees grow above 30 m, quite tall!!
# Let me see how many cells have values not less than 30 m tall.

canopy_height_harvard_df %>% # Picking the dataframe with canopy heights, and then.
  dplyr::select(layer) %>%  # Picking the layer column, and then.
  filter(layer >= 30) %>%  # Filtering rows with layer not less than 30 m, and then,
  nrow() # Counting their row numbers, actually the number of cells.
# It is unfortunate that my select in dplyr has been masked by select in raster.
# So I have to call it all through.
# Checking the height of the tallest tree (38.16998 m). Wow!

canopy_height_harvard_df %>% # Picks the dataframe, and then.
  dplyr::select(layer) %>%  # Picks the layer column to check values from, and then.
  summarize(maximum = max(layer, na.rm = TRUE)) # Returns the max of the value.

# Checking for the shortest tree in the dataset.

canopy_height_harvard_df %>% # Picks the dataframe, and then.
  dplyr::select(layer) %>%  # Picks the layer column to check values from, and then.
  summarize(minimum = min(layer, na.rm = TRUE)) # Returns the min of the values.

# Interesting the minimum value/shortest tree is 0 m. This is where both dsm
# and dtm have same value. Being a forested area, we can say these are flat zones
# which may need planting of trees.The question is, how many cells are "bare"?

canopy_height_harvard_df %>% # Picking the data.frame with canopy heights & then.
  dplyr::select(layer) %>%  # Picking the layer column, and then.
  filter(layer == 0) %>%  # Selecting those cells which are 0 m, and then,
  nrow() # Counting their row numbers, actually the number of cells.

# Since each cell is 1m by 1m, we can say bare ground is 51,135 m^2 which is 
# equivalent to 0.51135 ha of the region being bare. Please plan for raising
# seedlings to cover the bare area :).

# On average, we can estimate the mean canopy height at the site;

canopy_height_harvard_df %>% # Picks the dataframe, and then.
  dplyr::select(layer) %>%  # Picks the layer column to check values from & then.
  summarize(average = mean(layer, na.rm = TRUE)) # Returns the mean canopy height.

# Still we can get the median height of the trees;

canopy_height_harvard_df %>% # Picks the dataframe & then.
  dplyr::select(layer) %>%  # Picks the layer column to check values from & then.
  summarize(median = median(layer, na.rm = TRUE)) # Returns the median value.

# To calculate quantiles of the dataframe: Better summary with named rows;

canopy_height_harvard_df %>% # Picks the dataframe, and then.
  dplyr::select(layer) %>%  # Picks the layer column to check values from & then.
  summarize(quantiles = quantile(layer, na.rm = TRUE)) %>% 
  `row.names<-`(c("min", "lower quartile", "median", "upper quartile", "max"))

# Trying to have a summary of the cells

canopy_height_harvard_df %>% # A cool way to get the values as opposed to earlier.
  dplyr::select(layer) %>% 
  summary()

# Let us have a cool representation of the distribution of the canopy height values.

ggplot(data = canopy_height_harvard_df,
       aes(layer)) +
  geom_histogram(col = "blue", fill = "purple") +
  labs(title = "Distribution of tree heights in Harvard", # Adds title to the plot.
       x = "Tree heights (m)",
       y = "Frequency") +
  theme_classic() # Generating histogram to show the distribution of tree heights

# Trying out mapview

mapview(dsm_harvard)
mapview(plots_harvard) # The points can be added to the mapview too.
mapview(dsm_harvard) + mapview(plots_harvard_utm)

dsm_harvard_latlong <- projectRaster(dsm_harvard, crs = 4326)

m <- mapview(dsm_harvard_latlong) 
m1 <- mapview(plots_harvard)
m | m1

# I must get this up and running before sleeping tonight.

# Working with Vector data in R----
# Loading the points data of some plots in Harvard where dtm and dsm came from.
# You do not need to specify the file extension .shp in the readOGR function

plots_harvard <- readOGR(dsn = "NEON-airborne/plot_locations", 
                         layer = "HARV_plots")

class(plots_harvard) # This is a SpatialPointsDataFrame class object.

plot(plots_harvard) # Plots the five points. 

head(plots_harvard) # The coordinates are hidden but plot ids displayed.

# It is time to plot the points on top of the canopy height layer.

# Both layers must be in the same coordinate reference system.

crs(canopy_height_harvard) # Checking the crs of the canopy height raster.
crs(plots_harvard) # Checking the crs of the points layer.

# I have realized that they are different and I want the crs of canopy height
# to be the one I use because it is a projected utm and has units in m.
# That is better for small localized mapping as is in this case.
# I use spTransform function to accomplish the transformation.

plots_harvard_utm <- spTransform(x = plots_harvard, 
                                 CRSobj = crs(canopy_height_harvard))

# The spTransform() is just too powerful. Can't imagine doing that in QGIS.

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
# buffer argument.

plots_canopy_height_buffer <- raster::extract(canopy_height_harvard, 
                                              plots_harvard_utm, 
                                              buffer = 9,
                                              fun = mean)

head(plots_canopy_height_buffer) # The second point showed a big change 
# from 15.94m to 10.22m when buffer of 9 is used. Probably it was a lonely tree 
# in a rather bare background. Buffer of 9 computes mean height of 9 pixels
# around the central point.

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
                    has_coords = TRUE) # Returns a maximum of 1k records with coords.

species_gbif_df <- data.frame(species_gbif$gbif$data)
head(species_gbif_df$Dipodomys_ordii.longitude)


ggplot() +
  geom_polygon(data = us_map,
               aes(x = long, y = lat, # Plotting the USA map
                   group = group),
               fill = 'orange') +
  geom_point(data = species_gbif_df,
             aes(x = Dipodomys_ordii.longitude, # Plotting the species points
                 y = Dipodomys_ordii.latitude),
             col = "purple") +
  coord_quickmap() # Supporting better projection

# This is the end of the cool show of playing around with spatial data in R.


# In order to have everything in html, just click Ctrl + Shift + K
# Cool!!!

# Just learned how to mosaic my Sentinel data in SAGA when QGIS failed.
# Working around the bug in QGIS when masking or clipping is to save the masked
# or clipped image to file and not temporary folder. It worked after some struggles.


leaflet() %>% addTiles() %>% setView(35.348978,-0.372071,12)

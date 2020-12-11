# Working with raster data in R----
# Calculating vegetation canopy heights across NEON sites in Harvard
# The following libraries need to be loaded:

library(raster) # Reading and working with raster files (e.g. dsm and dtm)
library(rgdal)  # Reading and manipulating vector data (points)
library(tidyverse) # Wrangling and plotting (visualizing) the data
library(maps) # Loading us map and having it as dataframe
library(spocc) # for gathering species records from gbif

# Loading the dsm (digital surface model) data - height of top physical points
# Data were captured by lidar flyover in Harvard
# Typically trees but could be anything else projecting above the surface
# When I checked the study area using Google Earth, it is a forested site

dsm_harvard <- raster("NEON-airborne/HARV_dsmCrop.tif") # Reads the raster data

# Have a quick visual appeal of the data (raster).

raster::plot(dsm_harvard)

# Check the metadata of the dsm. Done by running the raster object

dsm_harvard

# Convert the dsm raster to dataframe to plot using ggplot2 package

dsm_harvard_df <- as.data.frame(dsm_harvard, 
                                xy = TRUE) # xy true ensures that coordinates
                                          # are also stored in the created dataframe

# Check the head of the dataframe created

head(dsm_harvard_df) # In deed the coordinates are also stored.

# Let me create another one with xy set to FALSE.

no_xy_df <- as.data.frame(dsm_harvard, xy = FALSE)
head(no_xy_df) # This is only givin back the values without xy coordinates.

# This is a very lengthy dataframe, check the length of one column

length(dsm_harvard_df$x) # This is the same value as ncell shown by

dsm_harvard

# Create the ggplot "HARV_dsmCrop" isthe name of the column with cell values

ggplot(data = dsm_harvard_df, 
       aes(x = x, y = y,
           fill = HARV_dsmCrop))+
  geom_raster() +
  labs(x = "Longitude (m)", 
       y = "Latitude (m)",
       title = "NEON-dsm of Harvard") 

# Importing the dtm data

dtm_harvard <- raster("NEON-airborne/HARV_dtmCrop.tif")

# Subtracting dtm from dsm to get the true canopy height. 
# Simple raster subtraction

canopy_height_harvard <- dsm_harvard - dtm_harvard

# Convert canopy height raster to dataframe for plotting

canopy_height_harvard_df <- as.data.frame(canopy_height_harvard, 
                                          xy = TRUE)
# Checking the head of the canopy dataframe

head(canopy_height_harvard_df)

# Generating the plot of the canopy data

ggplot() +
  geom_raster(data = canopy_height_harvard_df,
              aes(x = x, y = y, fill = layer ))

# Interesting some of the trees grow above 30 m

# Checking the height of the tallest tree (38.16998 m). Wow!

max(canopy_height_harvard_df$layer, na.rm = TRUE)

quantile(canopy_height_harvard_df$layer, 0.25, na.rm = TRUE) # Checking the 
# heights of trees at the 25th percentile.

ggplot(data = canopy_height_harvard_df,
       aes(layer)) +
  geom_histogram(col = "red", fill = "cyan") +
  labs(title = "Distribution of tree heights in Harvard",
       x = "Tree heights (m)",
       y = "Frequency") +
  theme_classic() # Generating histogram to show the distribution of tree heights



hist(canopy_height_harvard_df$layer, col = "purple") # Assigning my favorite color

# Working with Vector data----
# Loading the points data of some plots in Harvard
# You do not need to specify the file extension .shp in the readOGR function
plots_harvard <- readOGR("NEON-airborne/plot_locations", 
                         "HARV_plots")
class(plots_harvard) # This is a SpatialPointsDataFrame class object

plot(plots_harvard) # Plots the five points 

# It is time to plot the points on top of the canopy height layer.

# Both layers must be in the same coordinate reference system

crs(canopy_height_harvard) # Checking the crs of the canopy height raster
crs(plots_harvard) # Checking the crs of the points layer

# I have realized that they are different and I want the crs of canopy height
# to be the one I use because it has units in m.
# I use spTransform function to accomplish the transformation

plots_harvard_utm <- spTransform(plots_harvard, 
                                 crs(canopy_height_harvard))

crs(canopy_height_harvard)
crs(plots_harvard_utm)

# Now they have the same crs, this means they can now be overlaid precisely

# The next step is to convert the plots_harvard_utm to dataframe

plots_harvard_utm_df <- as.data.frame(plots_harvard_utm)
head(plots_harvard_utm_df) # Checking the first 6 rows of the dataframe

# Finally drawing the plots of both raster and points
# Using geom_raster and geom_point

ggplot() +
  geom_raster(data = canopy_height_harvard_df,
              aes(x = x, y = y, fill = layer)) +
  geom_point(data = plots_harvard_utm_df,
              aes(x = coords.x1, y = coords.x2),
              col = "yellow", cex = 3)
  
# Extracting raster values at point locations ----

plots_canopy_height <- raster::extract(canopy_height_harvard, 
                               plots_harvard_utm)
# Extracted points are in the same order as the id of the plots

head(plots_canopy_height)
plots_harvard_utm$plot_id # Checking the presence of plot_id

plots_canopy_height_df <- data.frame(plot_num = plots_harvard_utm$plot_id,
                                     plot_value = plots_canopy_height)

head(plots_canopy_height_df)

# In order to get average of cells around a point as the extracted value we use
# buffer argument

plots_canopy_height <- raster::extract(canopy_height_harvard, 
                                       plots_harvard_utm, buffer = 5, 
                                       fun = mean)

head(plots_canopy_height)

# Mapping points data, kind of species occurrence ----

us_map <- map_data("usa")

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

# Adding species occurrence data on the map

species_gbif <- occ(query = "Dipodomys ordii",
                    from = "gbif",
                    limit = 1000,
                    has_coords = TRUE)

species_gbif_df <- data.frame(species_gbif$gbif$data)
head(species_gbif_df$Dipodomys_ordii.longitude)


ggplot() +
  geom_polygon(data = us_map,
               aes(x = long, y = lat,
                   group = group),
               fill = 'grey') +
  geom_point(data = species_gbif_df,
             aes(x = Dipodomys_ordii.longitude,
                 y = Dipodomys_ordii.latitude),
             col = "red") +
  coord_quickmap()
devtools::install_github("babaknaimi/sdm")

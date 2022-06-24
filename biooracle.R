#######################################################################
# Created by Wyclife Agumba Oluoch                                    #
# Date 22nd Apr 2021                                                  #
# Working with bathymetry data to predict marine species distribution #
# Last edited 24th June 2022                                            #
#######################################################################

library(sdmpredictors)  # version 0.2.12
library(leaflet)        # version 2.1.1
library(tidyverse)      # version 1.3.1
library(raster)         # version 3.5.15

sdmpredictors::list_datasets() 
sdmpredictors::list_layers() 

bathy <- sdmpredictors::load_layers(c("BO_bathymin", "BO_bathymean", "BO_bathymax"))

sdmpredictors::layer_stats()

sdmpredictors::layers_correlation()

temp.max.bottom <- load_layers("BO2_tempmax_bdmax")

ne.atlantic.ext <- extent(-100, 45, 30.75, 72.5)
madagascar.ext <- extent(42.5, 50.8, -26, -11)
temp.max.bottom.crop <- crop(temp.max.bottom, madagascar.ext)

my.colors = colorRampPalette(c("#5E85B8","#EDF0C0","#C13127"))
plot(temp.max.bottom.crop,col=my.colors(1000),axes=FALSE, box=FALSE)
title(cex.sub = 1.5, sub = "Maximum temperature at the sea bottom (ºC)")

layers.bio2 <- list_layers(datasets="Bio-ORACLE")
layers.bio2

environment.bottom <- load_layers( layercodes = c("BO2_tempmax_bdmean" , 
                                                  "BO2_salinitymin_bdmean", 
                                                  "BO2_nitratemin_bdmean") , 
                                   equalarea=FALSE, rasterstack=TRUE)


# Download mean bathymetry data

bathymetry <- load_layers("BO_bathymean")

# Generate data.frame with the sites of interest

my.sites <- data.frame(Name = c("Faro, Portugal, NE Atlantic" , 
                              "Maspalomas, Spain, NE Atlantic" , 
                              "Guadeloupe, France, Caribbean Sea" , 
                              "Havana, Cuba, Caribbean Sea") , 
                       Lon = c(-7.873,-15.539,-61.208,-82.537) , 
                       Lat = c(37.047, 27.794,15.957,23.040 ) )
my.sites

# Visualize on Google map

m <- leaflet()
m <- addTiles(m)
m <- addMarkers(m, lng = my.sites$Lon, lat = my.sites$Lat, popup = my.sites$Name)
m

# Extract environmental values from layers
my.sites.environment <- data.frame(Name = my.sites$Name , 
                                   depth = raster::extract(bathymetry,my.sites[,2:3]) , 
                                   raster::extract(environment.bottom,my.sites[,2:3]))
my.sites.environment

sites_tibble <- as_tibble(my.sites.environment)

sites_tibble %>% 
  ggplot(aes(BO2_tempmax_bdmean, BO2_salinitymin_bdmean)) +
  geom_point()


my.countries <- data.frame(Name = c("Kenya_Home", "Uganda_Moroto", "Egypt_Cairo" , 
                                "Tanzania_Dar_Es_Salaam", "Germany_Tannenbusch_Mitte", 
                                "Belgium_Brussels", "France_Paris", 
                                "Netherlands_Enschede"), 
                       Lon = c(35.010684, 34.656301, 31.405672, 39.208862, 7.043222, 4.341395, 2.444176, 6.885655) , 
                       Lat = c(-0.237274, 2.527216, 30.108552, -6.872027, 50.749639, 50.894881, 48.967121, 52.223587))
my.countries

# Visualize on Google map

m <- leaflet()
m <- addTiles(m)
m <- addMarkers(m, lng = my.countries$Lon, lat = my.countries$Lat, popup = my.countries$Name)
m

####################################################
# Created by Wyclife Agumba Oluoch                 #     
# Created on 22nd April 2021                       #
# Last edited on 7th February 2023                 #
# Task: Mapping Watering points in Turkana County  #
####################################################

library(leaflet) # Version 2.1.1
library(dplyr)   # Version 1.1.0

dams <- read.csv('data/dams.csv')

dams <- dams %>% 
  mutate(popup_info = paste('County: ', County, "<br/>",
                            'Number: ',listing, "<br/>",
                            'Longitude: ', longitude, "<br/>",
                            'Latitude: ', latitude))

leaflet() %>%
  addTiles() %>%
  addCircleMarkers(data = dams,
                   lat = ~latitude, 
                   lng = ~longitude,
                   radius = ~1.8, 
                   popup = ~popup_info, 
                   color = 'blue',
                   weight = 14, 
                   opacity = 0.2,
                   fillColor = 'red',
                   fill = TRUE,
                   fillOpacity = 0.2)
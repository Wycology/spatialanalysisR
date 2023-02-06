####################################################
# Created by Wyclife Agumba Oluoch                 #     
# Created on 22nd April 2021                       #
# Last edited on 5th Febr 2023                         #
# Task: Mapping Watering points in Turkana County  #
####################################################

library(leaflet)
library(dplyr)

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

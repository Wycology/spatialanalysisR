####################################################
# Created by Wyclife Agumba Oluoch                 #     
# Created on 22nd April 2021                       #
# Last edited on th June 2022                    #
# Task: Mapping open water points in Turkana County#
####################################################

library(leaflet)
library(dplyr)

dams <- read.csv('data/dams.csv')

dams <- dams %>% dplyr::mutate(popup_info = base::paste('County: ', County, "<br/>", 
                                           'Number: ',listing, "<br/>",
                                           'Longitude: ', longitude, "<br/>", 
                                           'Latitude: ', latitude))

leaflet::leaflet() %>% 
  leaflet::addTiles() %>% 
  leaflet::addCircleMarkers(data = dams,
                            lat = ~latitude, lng = ~longitude,
                            radius = ~1.8, popup = ~popup_info)

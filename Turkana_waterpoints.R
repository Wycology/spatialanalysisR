library(leaflet)
library(tidyverse)

dams <- data.frame(longitude = c(35.582391, 35.581789, 
                                 35.579651, 35.579189, 
                                 35.577966, 35.577676,
                                 35.577371, 35.355020,
                                 35.247397, 35.254445,
                                 35.425808, 35.454905,
                                 35.398920),
                   latitude = c(3.222909, 3.222657, 
                                3.224472, 3.224581,
                                3.225090, 3.225225,
                                3.225368, 3.148446,
                                3.144302, 3.144950,
                                3.110151, 4.626087,
                                4.650788),
                   listing = c(1:13),
                   Sub_County = c('Loima', 'Loima', 'Loima', 'Loima',
                                  'Loima', 'Loima', 'Loima', 'Loima',
                                  'Loima', 'Loima', 'Loima', 'NA', 'NA'))




dams <- dams %>% mutate(popup_info = paste('Sub_County: ', Sub_County, "<br/>", 
                                           'Number: ',listing, "<br/>",
                                           'Longitude: ', longitude, "<br/>", 
                                           'Latitude: ', latitude))

leaflet() %>% 
  addTiles() %>% 
  addCircleMarkers(data = dams, 
                   lat = ~latitude,
                   lng = ~longitude, 
                   radius = ~3,
                   popup = ~popup_info)


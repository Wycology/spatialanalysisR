library(leaflet)
library(tidyverse)

dams <- data.frame(longitude = c(35.582391, 35.581789, 
                                 35.579651, 35.579189, 
                                 35.577966, 35.577676,
                                 35.577371),
                   latitude = c(3.222909, 3.222657, 
                                3.224472, 3.224581,
                                3.225090, 3.225225,
                                3.225368),
                   listing = c(1:7),
                   Sub_County = c('Loima', 'Loima', 'Loima', 'Loima',
                                  'Loima', 'Loima', 'Loima'))


leaflet() %>% 
  addTiles() %>% 
  addCircleMarkers(data = dams, 
                   lat = ~latitude,
                   lng = ~longitude, 
                   radius = ~3)

dams <- dams %>% 

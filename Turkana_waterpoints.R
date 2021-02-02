library(leaflet)
library(dplyr)

dams <- data.frame(longitude = c(35.582391, 35.581789, 35.579651, 35.579189, 
                                 35.577966, 35.577676, 35.577371, 35.355020,
                                 35.247397, 35.254445, 35.425808, 35.454905,
                                 35.398920, 35.512867, 35.639222, 35.697148,
                                 35.783488, 35.767541, 35.632926, 35.793876,
                                 35.795057, 35.795472, 35.688907, 34.988862,
                                 34.337086, 34.371066, 34.379100),
                   latitude = c(3.222909, 3.222657, 3.224472, 3.224581,
                                3.225090, 3.225225, 3.225368, 3.148446,
                                3.144302, 3.144950, 3.110151, 4.626087,
                                4.650788, 4.659210, 3.192118, 4.633935,
                                4.209630, 4.268374, 4.313147, 4.115460,
                                4.114035, 4.113375, 4.079765, 4.460585,
                                4.244939, 4.187053, 4.162871),
                   listing = c(1:27),
                   Sub_County = c('Loima', 'Loima', 'Loima', 'Loima',
                                  'Loima', 'Loima', 'Loima', 'Loima',
                                  'Loima', 'Loima', 'Loima', 'NA', 'NA', 
                                  'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 
                                  'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA'))

dams <- dams %>% mutate(popup_info = paste('Sub_County: ', Sub_County, "<br/>", 
                                           'Number: ',listing, "<br/>",
                                           'Longitude: ', longitude, "<br/>", 
                                           'Latitude: ', latitude))

leaflet() %>% addTiles() %>% 
  addCircleMarkers(data = dams, lat = ~latitude, lng = ~longitude, 
                   radius = ~4, popup = ~popup_info)
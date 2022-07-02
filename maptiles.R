# remotes::install_github('riatelab/maptiles')
library(maptiles) # For loading the tiles
library(sf) 

nc_raw <- st_read(system.file("shape/nc.shp", package = 'sf'), 
                      quiet = TRUE)

# re-project to EPSG:3857

nc <- st_transform(nc_raw, 'EPSG:3857') # picking a more appropriate crs

nc_osm <- get_tiles(nc, provider = 'CartoDB.Positron', crop = TRUE)

# Other providers

# "OpenStreetMap.MapnikBW", "OpenStreetMap", "OpenStreetMap.DE", "OpenStreetMap.France", "OpenStreetMap.HOT", 
# "Stamen.Toner", "Stamen.TonerBackground", "Stamen.TonerHybrid", "Stamen.TonerLines", "Stamen.TonerLabels", "Stamen.TonerLite", "Stamen.Watercolor", "Stamen.Terrain", "Stamen.TerrainBackground", "Stamen.TerrainLabels",
# "Esri.WorldStreetMap", "Esri.DeLorme", "Esri.WorldTopoMap", "Esri.WorldImagery", "Esri.WorldTerrain", "Esri.WorldShadedRelief", "Esri.OceanBasemap", "Esri.NatGeoWorldMap", "Esri.WorldGrayCanvas",
# "CartoDB.Positron", "CartoDB.PositronNoLabels", "CartoDB.PositronOnlyLabels", "CartoDB.DarkMatter", "CartoDB.DarkMatterNoLabels", "CartoDB.DarkMatterOnlyLabels", "CartoDB.Voyager", "CartoDB.VoyagerNoLabels", "CartoDB.VoyagerOnlyLabels",
# "Thunderforest.OpenCycleMap", "Thunderforest.Transport", "Thunderforest.TransportDark", "Thunderforest.SpinalMap", "Thunderforest.Landscape", "Thunderforest.Outdoors", "Thunderforest.Pioneer", "Thunderforest.MobileAtlas", "Thunderforest.Neighbourhood",
# "OpenTopoMap",
# "HikeBike", 
# "Wikimedia",

# -8595797, 4333842

plot_tiles(nc_osm)

plot(st_geometry(nc), col = NA, add = TRUE, axes = TRUE)

mtext(text = get_credit('CartoDB.Positron'),
      side = 1, 
      line = -1,
      adj = 0.2,
      padj = -20, 
      cex = 0.9, 
      font = 3,
      at = c(-8595797, 4333842))

# define the tile server parameters

osmnolbl <- list(src = 'osmnolabel',
                  q = 'https://{s}.tiles.wmflabs.org/osm-no-labels/{z}/{x}/{y}.png',
                  sub = c('a','b', 'c'), 
                  cit = '© OpenStreetMap contributors.')

nc_osmnolbl <- get_tiles(x = nc, provider = osmnolbl, crop = TRUE, 
                      cachedir = tempdir(), verbose = TRUE)

plot_tiles(nc_osmnolbl) # plotting the tile

mtext(text = osmnolbl$cit, side = 1, line = -1, 
      adj = 1, cex = .9, font = 3)

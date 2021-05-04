# remotes::install_github('riatelab/maptiles')
library(maptiles)
library(sf)

nc_raw <- st_read(system.file("shape/nc.shp", package = 'sf'), quiet = TRUE)

# reproject to EPSG:3857

nc <- st_transform(nc_raw, 'EPSG:3857')

nc_osm <- get_tiles(nc, crop = TRUE)

# -8595797, 4333842

plot_tiles(nc_osm)

plot(st_geometry(nc), col = NA, add = TRUE, axes = TRUE)

mtext(text = get_credit('OpenStreetMap'),
      side = 1, 
      line = -1,
      adj = 0.2,
      padj = -20, 
      cex = 0.9, 
      font = 3,
      at = c(-8595797, 4333842))

# define the tile server parameters



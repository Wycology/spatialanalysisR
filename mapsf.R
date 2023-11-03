# Task: Creating maps and inset maps
# Author: Wyclife Agumba Oluoch
# Contact: https://github.com/Wycology
# Last edited 16th May 2023 
# Last edited 18th October 2023 

pacman::p_load(mapsf, mapedit, tidyverse, sf) # Load the main map creation package

mf_theme(x = 'candy', bg = 'cyan', fg = 'purple',
  pos = 'center', inner = TRUE, line = 2,
  cex = 3, font = 16, mar = c(0, 0, 0, 0),
)

mtq <- mf_get_mtq() # Get the dataset

mtq$geom 

mtq_target <- mtq[26,] 

mtq_target <- mtq |> 
  dplyr::filter(LIBGEO == 'Sainte-Anne') 

mf_choro(mtq, var = 'POP', nbreaks = 5)

mf_map(mtq)

mf_map(mtq_target, add = TRUE, col = 'maroon')

mf_inset_on(x = mtq_target, pos = 'topright', cex = 0.4)

mf_scale(size = 0.1, pos = c(x = 702229, y = 1598119),
  lwd = 2, cex = 1.2,
  col = 'red', unit = 'km'
)

mf_map(mtq_target, add = FALSE, col = 'tomato')

mf_inset_off()

# Setting map theme ------------------------------------------------------------

mf_theme('dark')

mtq_target <- mtq[26, ] 

mf_map(mtq)

mf_map(mtq_target, add = TRUE, col = 'tomato') 

mf_inset_on(x = mtq_target, pos = 'topright', cex = .4) 

mf_init(mtq_target)

mf_shadow(mtq_target, add = TRUE) 

mf_map(mtq_target, add = TRUE, col = 'tomato')

mf_title('Saint-Anne',
  pos = 'left',
  tab = TRUE,
  cex = .9,
  line = 1,
  inner = TRUE
)

mf_scale(size = 2) # Scale is changeable

box(col = 'cyan', lwd = 3, lty = 2) # Add bounding box around inset map

mf_inset_off() # Removing the inset mapping mode

mf_title('Martinique Municipalities') # Adding title to the map

mf_scale(size = 5) # Adding a map scale

mf_credits(txt = 'T. Giraud, 2021') # Accreditation to the authors

# Inset on lower left ----------------------------------------------------------

mf_theme('default', mar = c(0, 0, 0, 0))

mf_map(mtq) # Displaying the whole map

mf_inset_on(fig = c(0, 0.25, 0, 0.25))

mf_map(mtq[9,], col = 'blue')

box(col = "green", lwd = 3, lty = 3)

mf_inset_off()

# World map as inset -----------------------------------------------------------

mf_map(mtq)

mf_inset_on(x = 'worldmap')

mf_worldmap(mtq)

box(col = "gray", lwd = 3, lty = 4) # Adding bounding box around world map inset

mf_inset_off()

# Graphs as inset --------------------------------------------------------------

bks <- mf_get_breaks(x = mtq$MED,
                     nbreaks = 5,
                     breaks = 'quantile')
bks <- round(bks, digits = -2)

pal <- hcl.colors(n = 5, palette = 'Dark Mint', rev = TRUE)

mf_theme('candy')

fg <- mf_theme()$fg

mf_map(x = mtq, var = 'MED', type = 'choro',
  pal = pal, breaks = bks, leg_pos = NA
)

mf_inset_on(fig = c(0.75, .95, .84, .99))

par(mar = c(0, 0, 1.7, 0))

hist(mtq$MED, breaks = bks, col = pal, border = fg,
  axes = F, labels = '', xlab = '', ylab = '', main = ''
)

axis(
  side = 1,
  at = bks,
  las = 2,
  tick = FALSE,
  line = -.9,
  cex.axis = .7,
  col.axis = fg
)

title(
  'Median Income\nin euros',
  cex.main = .8,
  col.main = fg,
  font.main = 1,
  adj = 0
)
box()

mf_inset_off()

mf_title('Wealth in Martinique, 2015', pos = 'left')

mf_scale(5)

mf_credits(paste0(
  'Sources: Insee and IGN, 2018\n',
  'mapsf ',
  packageVersion('mapsf')
))

# Exporting map for publication ------------------------------------------------

mtq <- mf_get_mtq()

mf_export(x = mtq,
          filename = 'fixed_width.png',
          width = 500) # The code for export

mf_map(mtq, add = TRUE) # Adding the map

mf_title(txt = 'PNG export: width = 500px, height = 605px (deducted)')

dev.off()

mf_export(
  x = mtq,
  filename = 'fixed_width_theme.png',
  width = 500,
  theme = 'dark'
)
mf_map(mtq, add = TRUE)

mf_title(txt = 'PNG export: width = 500px, height = 611px (deducted)')

dev.off()
mf_export(
  x = mtq,
  filename = 'fixed_width_expand.png',
  width = 500,
  expandBB = c(0, 0, 0, .5),
  theme = 'default'
)

# Export with space around the bounds ------------------------------------------

mf_map(mtq, add = TRUE)

mf_title(txt = 'PNG export: width=500px, height = 427px (deducted)')

dev.off()
target <- mtq[5,]

# Center map at specific region ------------------------------------------------

mf_export(x = target,
          filename = 'fixed_height_centered.png',
          height = 600)

mf_init(mtq, theme = 'candy')

mf_map(mtq, add = T)

mf_shadow(target, add = TRUE)

mf_map(target, add = TRUE)

mf_title(txt = 'PNG export: height = 600px, width = 433px (deducted)')

mf_scale(1, pos = 'bottomleft')

dev.off()
mf_export(
  x = mtq,
  export = 'svg',
  filename = 'fixed_width.svg',
  width = 5,
  theme = 'nevermind',
  bg = 'black'
)

# Additional parameters --------------------------------------------------------

mf_map(mtq, add = TRUE)

mf_title(txt = "SVG export: bg = 'black'")

dev.off()

mf_theme('agolalight')

mf_export(x = mtq,
          filename = 'theme_before.png',
          width = 200)
mf_map(mtq, add = TRUE)

mf_title(txt = 'agolalight')
dev.off()

mf_export(
  x = mtq,
  filename = 'theme_within.png',
  width = 200,
  theme = 'darkula'
)

mf_init(mtq)
mf_map(mtq, add = TRUE) # Add another layer onto the map.

mf_title(txt = 'darkula') # Assigning title to the plotted map

dev.off() # Terminating the inset

mf_map(mtq) # Starting a map again.

mf_title("Still 'darkula'")

# When nothing plots then you may need to run dev.off() several times.

# Joining sf object with a dataframe --------------------------------------

ma <- mapedit::drawFeatures()

names(ma) <- c("id", "feature_type", "geometry")

plot(st_geometry(ma), 
     axes = TRUE,
     xlab = "Longitude",
     ylab = "Latitude",
     main = "A couple of houses in Awasi",
     arrow = TRUE)

df <- data.frame(id = c( 1052,1036, 1104),
                 area = c(5, 6, 18),
                 name = c("House_1", "House_2", "House_3"))

me <- left_join(ma, df, by = "id")

library(leaflet)
leaflet() |> addTiles() |> addPolygons(data = me)

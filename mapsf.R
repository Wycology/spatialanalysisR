# Task: Creating maps and inset maps
# Author: Wyclife Agumba Oluoch
# Contact: https://github.com/Wycology

library(mapsf) # Loading the mapsf library for mapping

mf_theme(x = 'candy',
         bg = 'purple',
         fg = 'purple',
         pos = 'center',
         inner = TRUE,
         line = 2,
         cex = 3,
         font = 16,
         mar = c(0, 0, 0, 0),
         ) # Creating theme for plotting the map

mtq <- mf_get_mtq() # Loading Martinique sf data to the R environment

mtq$geom # Looking at the geometry column of the data. Projected coords   

mtq_target <- mtq[c(26), ] # Sub-setting single polygon and storing to object

mf_choro(mtq, var = 'POP', nbreaks = 5) # plotting pop data

colnames(mtq) # column names of the data-set

mf_map(mtq) # Displaying the map of the mtq object

mf_map(mtq_target, add = TRUE, col = 'yellow') # Awesome tomato color

mf_inset_on(x = mtq_target, pos = 'topright', cex = .4)

mf_scale(size = 2, 
         pos = c(x = 702229, y = 1598119), 
         lwd = 2, 
         cex = 1.2, 
         col = 'red', 
         unit = 'km')

mf_map(mtq_target, add = FALSE, col = 'tomato') # mapping the inset

mf_inset_off() # stopping the inset function from running

# set theme for the project

mf_theme('dark') # picking the dark theme of mapsf to do the mapping

mtq_target <- mtq[c(26),] # picking the right map to plot out

mf_map(mtq) # mapping the sf object to the plots window

mf_map(mtq_target, add = TRUE, col = 'tomato') # tomato is the color, fun

mf_inset_on(x = mtq_target, pos = 'topright', cex = .4) # location of the inset

mf_init(mtq_target) # initializing the map building

mf_map(mtq, add = TRUE) # Adding more maps

mf_shadow(mtq_target, add = TRUE) # adding map shadow 

mf_map(mtq_target, add = TRUE, col = 'tomato')

mf_title('Saint-Anne', pos = 'left', tab = TRUE, cex = .9, line = 1, inner = TRUE)

mf_scale(size = 2) # Scale is changeable

box()

mf_inset_off() # removing the inset mapping mode

mf_title('Martinique Municipalities') # Adding title to the map

mf_scale(size = 5) # Adding map scale

mf_credits(txt = 'T. Giraud, 2021') # Accreditation to the authors

mf_theme('default', mar = c(0, 0, 0, 0))

mf_map(mtq) # Mapping the whole map

mf_inset_on(fig = c(0, 0.25, 0, 0.25))

mf_map(mtq[9, ])

box()

mf_inset_off()

# world map inset

mf_map(mtq)

mf_inset_on(x = 'worldmap')

mf_worldmap(mtq)

mf_inset_off()

# Non-geographic insets, such as histograms

bks <- mf_get_breaks(x = mtq$MED, nbreaks = 5, breaks = 'quantile')
bks <- round(bks, digits = -2)

pal <- hcl.colors(n = 5, palette = 'Dark Mint', rev = TRUE)

mf_theme('candy')

fg <- mf_theme()$fg

mf_map(x = mtq, var = 'MED', type = 'choro', 
       pal = pal, breaks = bks, leg_pos = NA)

mf_inset_on(fig = c(0.75, .95, .84, .99))

par(mar = c(0, 0, 1.7, 0))

hist(mtq$MED, breaks = bks, col = pal, border = fg, axes = F, labels = '',
     xlab = '', ylab = '', main = '')

axis(side = 1, at = bks, las = 2, tick = FALSE, line = -.9,
     cex.axis = .7, col.axis = fg)

title('Median Income\nin euros', cex.main = .8, col.main = fg,
      font.main = 1, adj = 0)

mf_inset_off()

mf_title('Wealth in Martinique, 2015', pos = 'left')
mf_scale(5)
mf_credits(paste0('Sources: Insee and IGN, 2018\n', 'mapsf ', packageVersion('mapsf')))

# Exporting map

mtq <- mf_get_mtq()

mf_export(x = mtq, filename = 'fixed_width.png', width = 500)

mf_map(mtq, add = TRUE)

mf_title(txt = 'PNG export: width = 500px, height = 605px (deducted)')
dev.off()

mf_export(x = mtq, filename = 'fixed_width_theme.png',
          width = 500, theme = 'dark')
mf_map(mtq, add = TRUE)

mf_title(txt = 'PNG export: width = 500px, height = 611px (deducted)')
dev.off()

# Export with extra space

mf_export(x = mtq, filename = 'fixed_width_expand.png',
          width = 500, expandBB = c(0, 0, 0, .5), theme = 'default')

mf_map(mtq, add = TRUE)

mf_title(txt = 'PNG export: width=500px, height = 427px (deducted)')
dev.off()

# Center the map on a specific area

target <- mtq[5, ]
mf_export(x = target, filename = 'fixed_height_centered.png', height = 600)

mf_map(mtq, add = TRUE)

mf_shadow(target, add = TRUE)

mf_map(target, add = TRUE)

mf_title(txt = 'PNG export: height = 600px, width = 433px (deducted)')

mf_scale(1, pos = 'bottomleft')
dev.off()

# Other parameters

mf_export(x = mtq, export = 'svg', filename = 'fixed_width.svg',
          width = 5, theme = 'nevermind', bg = 'black')

mf_map(mtq, add = TRUE)

mf_title(txt = "SVG export: bg = 'black'")
dev.off()

mf_theme('agolalight')

mf_export(x = mtq, filename = 'theme_before.png', width = 200)
mf_map(mtq, add = TRUE)

mf_title(txt = 'agolalight')
dev.off()

mf_export(x = mtq, filename = 'theme_within.png', width = 200,
          theme = 'darkula')

mf_map(mtq, add = TRUE) # adding another layer to the map

mf_title(txt = 'darkula') # assigning title to the map

dev.off() # stopping the inset

mf_map(mtq)

mf_title("Still 'darkula'")

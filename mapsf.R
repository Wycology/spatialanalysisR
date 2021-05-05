# Task: Creating inset maps

library(mapsf)

mf_theme('default', mar = c(0, 0, 0, 0))
mtq <- mf_get_mtq()
mtq
mtq_target <- mtq[c(26), ]
mf_map(mtq)

mf_map(mtq_target, add = TRUE, col = 'tomato')

mf_inset_on(x = mtq_target, pos = 'topright', cex = .4)

mf_map(mtq_target, add = FALSE, col = 'tomato')

mf_inset_off()

# set theme

mf_theme('dark')

mtq_target <- mtq[c(26),]

mf_map(mtq)
mf_map(mtq_target, add = TRUE, col = 'tomato')

mf_inset_on(x = mtq_target, pos = 'topright', cex = .4)

mf_init(mtq_target)

mf_map(mtq, add = TRUE)

mf_shadow(mtq_target, add = TRUE)

mf_map(mtq_target, add = TRUE, col = 'tomato')

mf_title('Saint-Anne', pos = 'left', tab = TRUE, cex = .9, line = 1, inner = TRUE)

mf_scale(size = 2)

box()

mf_inset_off()

mf_title('Martinique Municipalities')

mf_scale(size = 5)

mf_credits(txt = 'T. Giraud, 2021')

mf_theme('default', mar = c(0, 0, 0, 0))

mf_map(mtq)

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





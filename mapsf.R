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









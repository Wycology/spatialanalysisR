# Task: Creating inset maps

library(mapsf)

mf_theme('default', mar = c(0, 0, 0, 0))
mtq <- mf_get_mtq()
mtq

mf_map(mtq)

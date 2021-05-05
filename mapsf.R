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

# 










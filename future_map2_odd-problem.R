#This is addressed in issue #30 found here:
#https://github.com/DavisVaughan/furrr/issues/30
# Will be fixed in v 0.2
DATA = data.frame(SITE = rep(c("a","b","c","d","e"), 2), DATA = rnorm(10))

site_list = as.list(unique(DATA$SITE))

map2(list(DATA), site_list, function(x,y) x %>% subset(SITE == as.character(y)))

future_map2(list(DATA), site_list, function(x,y) x %>% subset(SITE == as.character(y)))

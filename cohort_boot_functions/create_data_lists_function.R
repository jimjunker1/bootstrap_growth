create_data_lists = function(DATA,...) {
  #create a list of sites in data frame
  sites_list <<- as.list(unique(DATA$SITE))
  #create list of site data subsets
  sites_data_list <<- map2(list(DATA), sites_list, sites_subset)
  #create lists of taxa in each site
  taxa_lists <<- map(sites_data_list, function(x) as.list(unique(x$TAXON)))
  #LIST = list(sites_list, sites_data_list, taxa_lists)
  #NAMES = c("sites_list","sites_data_list", "taxa_lists")
  #list2env(NAMES, envir = parent.frame())
  #lapply(seq_along(LIST), function(x) {
  #  assign(NAMES[x], LIST[[x]], envir = parent.frame())
  #  })
  #pmap(list(LIST, list(NAMES)), function(x,y) assign(y, x, envir = parent.frame()))
  #list(sites_list = sites_list, sites_data_list = sites_data_list, taxa_lists = taxa_lists)
}
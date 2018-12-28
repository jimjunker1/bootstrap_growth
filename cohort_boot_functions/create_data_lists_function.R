create_data_lists = function(DATA,...) {
  #create a list of sites in data frame
  sites_list = as.list(unique(DATA$SITE))
  #create list of site data subsets
  sites_data_list = map2(list(DATA), sites_list, sites_subset)
  #create lists of taxa in each site
  taxa = map(sites_data_list, function(x) as.list(unique(x$TAXON)))
  return(c(sites_data_list = sites_data_list, sites_data_list = sites_data_list, taxa = taxa))
}
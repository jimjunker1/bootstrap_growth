date_reorder = function(site_taxa_data_list, cohort_date_list,...){
  map2(site_taxa_data_list, cohort_date_list, function(x,y) {
    x %>% mutate(DATE = factor(x$DATE, levels = levels(y$DATE))) %>%
      left_join(y %>% select(c(DATE, day)))})
}
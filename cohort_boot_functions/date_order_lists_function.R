date_order_lists = function(site_taxa_data_list,...){
  date_order = future_map(site_taxa_data_list, function(x) { x %>%
      group_by(DATE) %>%
      summarise(mean_mass = median(MASS, na.rm=T), julian = unique(julian, na.rm = T)) %>%
      arrange(mean_mass) %>%
      mutate(DATE = reorder(DATE, mean_mass)) %>%
      mutate(day = c(diff(julian),NA), id = 1:n()) %>%
      mutate(day = ifelse(day < 0, 365-abs(day),day))})
}
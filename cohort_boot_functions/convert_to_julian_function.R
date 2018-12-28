convert_to_julian = function(data_list,...) {
  posix_vec = as.POSIXct(data_list$DATE, format = "%d-%m-%y")
  data_list %>% mutate(DATE = as.factor(posix_vec)) %>%
    mutate(julian = as.numeric(format(posix_vec, "%j")))
}
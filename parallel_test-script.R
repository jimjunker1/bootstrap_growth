#parallel test script 
library(tidyverse)
library(furrr)
source("./cohort_boot_functions/create_data_lists_function.R")
source("./cohort_boot_functions/sites_subset_function.R")
source("./cohort_boot_functions/convert_to_julian_function.R")

df = read.csv(file = "./infrequensMeasurements.csv",T)
df <- df %>% dplyr::rename(TAXON = species, SITE = Site, DATE = date, MASS = mass)
df2 = df %>% mutate(TAXON = "infrequens2")
DATA = rbind(df,df2)

#### create_data_lists(DATA) ####

#create a list of sites in data frame
sites_list = as.list(unique(DATA$SITE))
#create list of site data subsets
sites_data_list = map2(list(DATA), sites_list, sites_subset)
#create lists of taxa in each site
taxa_lists = map(sites_data_list, function(x) as.list(as.character(unique(x$TAXON))))#list of lists
#####
#### convert_to_julian(data_list) ####
posix_vec = as.POSIXct(sites_data_list[[1]]$DATE, format = "%d-%m-%y")
DATE = as.factor(posix_vec)
julian = as.numeric(format(posix_vec, "%j"))
####
#convert all dates to julian
sites_data_list = map(sites_data_list, convert_to_julian)

##### date_reorder(site_data_list, taxa_list) ####
site_list = c(rep(unique(sites_data_list[[1]]$SITE),length(taxa_lists[[1]])))
dada = map2(list(site_list),taxa_lists[[1]], function(x,y) {
  print(x);print(y)
  sites_data_list[[1]] %>%
    filter(SITE == x & TAXON == y)})
#subset taxa/cohorts and reorder dates based start of cohort
date_reorder = function(site_data_list, taxa_list,...){
#cohort_site_date_df = 
  site_list = c(rep(unique(sites_data_list[[1]]$SITE),length(taxa_lists[[1]])))
  map2(list(site_list),taxa_list, function(x,y){
    date_order = site_data_list %>% filter(SITE == x,TAXON == y) %>%
    group_by(DATE) %>%
    summarise(mean_mass = mean(MASS, na.rm=T), julian = unique(julian, na.rm = T)) %>%
    arrange(mean_mass) %>%
    mutate(DATE = reorder(DATE, mean_mass)) %>%
    mutate(day = c(diff(julian),NA), id = 1:n()) %>%
    mutate(day = ifelse(day < 0, 365-abs(day),day))
    site_data_list$DATE = factor(site_data_list$DATE, levels = levels(date_order$DATE))
  })
  #wrap_loc = which(cohort_site_date_df$day < 0);wrap_fix = 365-abs(cohort_site_date_df[which(cohort_site_date_df$day < 0),'day'])
  ##date_df[wrap_loc, 'day'] = wrap_fix
  #site_data$date = factor(site_data$date, levels = levels(date_df$date))
}
#debugonce(date_reorder)
cohort_site_data = pmap(list(sites_data_list, taxa_lists), date_reorder)
cohort_site_data[[1]]
cohort_sub
#create a list of data subset for every taxa
tax_data_lists= map2(list(df), taxa, taxa_subset)
#create lists of 
sites = map(tax_data_lists, function(x) as.list(unique(x$Site)))



lapply(taxa, function(x) unique(df[x,"SITE"]))#list(unique(x[,"SITE"])))

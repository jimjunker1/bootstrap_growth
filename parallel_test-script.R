#parallel test script 
library(tidyverse)
library(furrr)
source("./cohort_boot_functions/create_data_lists_function.R")
source("./cohort_boot_functions/sites_subset_function.R")
source("./cohort_boot_functions/convert_to_julian_function.R")
source("./cohort_boot_functions/taxa_list_split_function.R")
source("./cohort_boot_functions/date_order_lists_function.R")
source("./cohort_boot_functions/date_reorder_function.R")

df = read.csv(file = "./infrequensMeasurements.csv",T)
df <- df %>% dplyr::rename(TAXON = species, SITE = Site, DATE = date, MASS = mass)
df2 = df %>% mutate(TAXON = "infrequens2")
DATA = rbind(df,df2)

source("./parallel_boot_function.R")
parallel_cohort_boot(DATA, parallel = TRUE)
tic();create_data_lists(DATA);toc()

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
#### taxa_list_split(site_data_list, taxa_list) ####
taxa_list_split = function(site_data_list, taxa_list,...){
  site_list = c(rep(unique(site_data_list$SITE),length(taxa_list)))
  map2(list(site_list),taxa_list, function(x,y){
    site_data_list %>% filter(SITE == x, TAXON == y)
  })
}
#subset taxa/cohorts and create list of lists of data
site_taxa_data_lists = pmap(list(sites_data_list,taxa_lists), taxa_list_split)

##### date_reorder(site_taxa_data_lists) #####
date_order_lists = function(site_taxa_data_list,...){
  date_order = map(site_taxa_data_list, function(x) { x %>%
    group_by(DATE) %>%
    summarise(mean_mass = median(MASS, na.rm=T), julian = unique(julian, na.rm = T)) %>%
    arrange(mean_mass) %>%
    mutate(DATE = reorder(DATE, mean_mass)) %>%
    mutate(day = c(diff(julian),NA)) %>%#, id = 1:n()) %>%
    mutate(day = ifelse(day < 0, 365-abs(day),day))})
}
#create list of date order based start of cohort
cohort_date_lists = map(site_taxa_data_lists, date_order_lists)

##### date_reorder
date_reorder = function(site_taxa_data_list, cohort_date_list,...){
  map2(site_taxa_data_list, cohort_date_list, function(x,y) {
  x %>% mutate(DATE = factor(x$DATE, levels = levels(y$DATE))) %>%
  left_join(y %>% select(c(DATE, day)))})
}
#reorder dates based start of cohort
site_taxa_data_lists = pmap(list(site_taxa_data_lists,cohort_date_lists), date_reorder)

###### start here to integrate parallel bootstrap function after build it ######
source("./bootstrap_function.R")#function selects data 
bootsdata = boots(site_data, nboot = nboot)



ggplot(site_taxa_data_lists[[1]][[1]], aes(x = DATE, y = MASS)) + geom_point(size = 2, position = 'jitter')
make_plot = function(site_taxa_list,...){
  ggplot(site_taxa_list, aes(x = DATE, y = MASS)) + geom_point(size = 2, position = 'jitter') +
    scale_y_continuous(limits = c(0,5))
}
map(site_taxa_data_lists, function(x) map(x,make_plot))




ggplot(site5, aes(x = Pd, y = mass, colour = as.factor(SITE))) + 
  geom_point(size =2, position = "jitter")


map()


cohort_sub
#create a list of data subset for every taxa
tax_data_lists= map2(list(df), taxa, taxa_subset)
#create lists of 
sites = map(tax_data_lists, function(x) as.list(unique(x$Site)))



lapply(taxa, function(x) unique(df[x,"SITE"]))#list(unique(x[,"SITE"])))

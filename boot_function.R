#boot strap function

cohort_boot = function(DATA,nboot = NULL,...){
  #sets the number of individuals you want to sample
  if(is.null(nboot)){
    nboot = 500
  } else{ nboot = nboot }
  
  taxa = unique(levels(DATA$species))#set taxa levels
for(i in taxa){
  tax_data = DATA[which(DATA$species == as.character(i)),]#grab single taxa
  sites = unique(levels(as.factor(tax_data$Site)))#set site levels for a single taxa
for(j in sites){
  site_data = tax_data[which(as.factor(tax_data$Site) == j),]#grab all data for a single site

#convert to dates and calculate julian day
site_data$Pd = as.POSIXct(site_data$date, format = "%d-%m-%y")
site_data$date = as.factor(site_data$Pd)
site_data$julian = as.numeric(format(site_data$Pd, "%j"));site_data$Pd = NULL

#reorder dates based start of cohort
site_data %>%
  group_by(date) %>%
  summarize(mean_mass = mean(mass, na.rm=T), julian = unique(julian, na.rm = T)) %>%
  arrange(mean_mass) %>%
  mutate(date = reorder(date, mean_mass)) %>%
  mutate(day = c(diff(julian),NA), id = 1:n()) -> date_df

#fix the weird wrap around dates
wrap_loc = which(date_df$day < 0);wrap_fix = 365-abs(date_df[which(date_df$day < 0),'day'])
date_df[wrap_loc, 'day'] = wrap_fix

#reorder the site data based on 'cohort date'
site_data$date = factor(site_data$date, levels = levels(date_df$date))

source("./bootstrap_function.R")#function selects data 
bootsdata = boots(site_data, nboot = nboot)
bootsdata = bootsdata[-1,]

bootsdata %>%
  select(species, Site, date, mass) %>%
  left_join(date_df[,c('date','id')]) -> bootsdata_id

bootsdata %>%
  select(species, Site, date, mass) %>%
  mutate(id = rep(1:nboot, nrow(date_df))) %>%
  spread(date,mass) %>%
  select(-id)-> bootsdata_wide

#now work across two columns at a time to estimate growth
#1:names(bootsdata[,-1:2]))
growth_rate <- function(x,y,z) {(log(x)/log(y))/z}

igr_df = c()
for(w in 1:(max(bootsdata_id$id)-1)){
  igr_vec = mapply(growth_rate, bootsdata_id[which(bootsdata_id$id == w),'mass'], bootsdata_id[which(bootsdata_id$id == w+1),'mass'],
                   date_df[w,'day'])
  igr_df = cbind(igr_df,igr_vec)
}

colnames(igr_df) = unique(levels(date_df$date))[1:(nrow(date_df)-1)]
igr_df = data.frame(species = as.character(i), site = j, igr_df, check.names = F)
igr_df %>%
  gather(start_date, IGR, 3:(dim(igr_df)[2]), factor_key = T) -> igr_long

write.csv(igr_long, file = paste(i,"_site-",j,"_IGR.csv",sep = ""), row.names = F) 
  }  
 }  
}

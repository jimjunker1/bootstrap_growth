#boot strap function

cohort_boot = function(DATA,nboot = NULL,...){
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
  summarize(mean_mass = mean(mass, na.rm=T), julian = unique(julian, na.rm = T)) -> date_df
date_mass = within(date_df, date <- reorder(date, mean_mass))
site_data$date = factor(site_data$date, levels = levels(date_df$date))
source("./bootstrap_function.R")
bootsdata = boots(site_data, nboot = nboot)
bootsdata = bootsdata[-1,]
browser()

bootsdata %>%
  select(species, date, mass) %>%
  mutate(id = rep(1:nboot, nrow(date_df))) %>%
  spread(date,mass) %>%
  select(-id)-> bootsdata.wide


igr = c()
start_date = c()
for(i in 1:unique(levels(bootsdata$date))){
 paste(i+1)
 igr_vec = mapply(function(w,x,y,z) (log(w)/log(x))/(z-y), bootsdata[which(bootsdata$date == i),'mass'], 
                                                                         bootsdata[which(bootsdata$date == as.character(i+1)),'mass'],
                    bootsdata[which(bootsdata$date == i),'julian'],bootsdata[which(bootsdata$date == as.character(i+1)),'julian'])
 i_date = rep(as.chracter(i), nboot)
 
 igr = rbind(igr, igr_vec)
 start_date = rbind(start_date, i_date)
} 
  

mapply(function(x,y) (log(mass[i])/log(mass[i+1]))/(julian[i+1]-julian[i]))
  


  }  
 }  
}

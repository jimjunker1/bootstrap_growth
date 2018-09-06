# script to boot strap growth rates from size-frequency histograms
library(tidyverse);theme_set(theme_bw())


##### some basic QAQC of data for my benefit #####

#read in file
df = read.csv(file = "./infrequensMeasurements.csv",T)

#check levels of date
unique(levels(df$date))
#double check species levels
unique(levels(df$species))

#need to convert dates into POSIX format
df$Pd = as.POSIXct(df$date, format = "%d-%m-%y")
#just a quick visualization
ggplot(df, aes(x = Pd, y = mass, colour = as.factor(Site))) + 
  geom_point(size =2, position = "jitter")
##### done with QAQC #####

df = read.csv(file = "./infrequensMeasurements.csv",T)
#df$Pd = as.POSIXct(df$date, format = "%d-%m-%y")
#df$date = as.factor(df$Pd)
#df$julian = as.numeric(format(df$Pd, format = "%j"));df$Pd = NULL



#load in the function
source("./boot_function.R")
debugonce(cohort_boot)
x = cohort_boot(df)


df %>%
  group_by(Pd) %>%
  summarize(mean_mass = mean(mass, na.rm = T)) -> date_mass

date_mass = within(date_mass, Pd <- reorder(Pd, mean_mass))
levels(date_mass$Pd)

df$date = factor(df$date, levels = levels(date_mass$Pd))
levels(df$date)

bootsdata = subset(df, Site == 1)

del.mass = mapply(function(x,y) (log(x)/log(y)), bootsdata[which(bootsdata$date == "17-09-17"),'mass'], 
                  bootsdata[which(bootsdata$date == "15-10-17"),'mass'])

for(i in 1:unique(levels(as.factor(df$date)))){
  
} 

unique(levels(as.factor(df$date)))

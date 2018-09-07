# script to boot strap growth rates from size-frequency histograms
library(tidyverse);theme_set(theme_bw())
library(tictoc)

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

x = subset(df, Site == "1" & date == "19-11-17")
hist(x$bodySize, breaks = 30)

#quick notes: You should take a look at some of the size frequency data (see above).
#It appears you potentially have multiple cohorts going at the same time.
#you should probably split these out at some point, otherwise this code is going 
#to do some weird stuff and you will get some wack-a-doo growth estimates
##### done with QAQC #####

df = read.csv(file = "./infrequensMeasurements.csv",T)

#load in the function
source("./boot_function.R")
#debugonce(cohort_boot)
tic()
suppressWarnings(cohort_boot(df))
toc()

#take a look at one example
x = read.csv(file = "infrequens_site-1_IGR.csv",T)

x_1 = subset(x,site == 1 & start_date == "2017-08-18")

hist(x_1$IGR)
median(x_1$IGR)
quantile(x_1$IGR, c(0.05,0.95))

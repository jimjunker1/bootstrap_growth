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
##### done with QAQC #####

#run the code to estimate growth rates
df = read.csv(file = "./infrequensMeasurements.csv",T)
#load in the function
source("./boot_function.R")
#debugonce(cohort_boot)
tic()
suppressWarnings(cohort_boot(df))
toc()

#take a look at one example
x = read.csv(file = "./output/infrequens_site-1_IGR.csv",T)

x_1 = subset(x,site == 1 & start_date == "2017-09-18")

hist(x_1$IGR)
median(x_1$IGR)
quantile(x_1$IGR, c(0.05,0.95))

#quick notes: You should take a look at some of the size frequency data (see above in QAQC).
#It appears you potentially have multiple cohorts going at the same time.
#you should probably split these out at some point, otherwise this code is going 
#to do some weird stuff and you will get some wack-a-doo growth estimates

#what i do is rename the species/taxon name with a number appended. for example,
# infrequens_1, infrequens_2, etc. for each cohort. It looks like this is not an issue
#at every site, just some. Which is interesting...

#Here is some script to run len_frequency histograms
df = read.csv(file = "./infrequensMeasurements.csv", T)
#create some variables so don't have to F with len_freq code
df$Pd = as.POSIXct(df$date, format = "%d-%m-%y")
df$DATE = format(df$Pd, "%m/%d/%Y")
df$HABITAT = "COBBLE"
df$SAMPLE = "1"
df = df[,c(1,7:9,3:4)]
#change columns to fit code
colnames(df) = c("SITE", "DATE","HABITAT", "SAMPLE","TAXON","bodysize")
#widen it out
df %>%
  group_by(SITE, DATE, HABITAT, SAMPLE, TAXON, bodysize) %>%
  tally() %>%
  spread(key = bodysize, value = n, fill = 0) -> df.wide


source("./ins_julian_function.txt")
ins_julian(df.wide, "./output/df_wide")
source("./len_freq_function.R")
df.julian = read.table(file ='./output/df_wide_julian.txt', header = T, sep = "\t", quote = "", strip.white = T,
                       check.names = F, row.names = NULL)

#this is running oddly slow on my computer if you run it on just the subset, can you let me know how 
#quickly it ran, please.
tic();len_freq(df.julian, fun = sum);toc()


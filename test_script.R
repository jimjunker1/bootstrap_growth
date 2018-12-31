# script to boot strap growth rates from size-frequency histograms
library(tidyverse);theme_mod = function() {theme_bw() %+replace% theme(panel.grid = element_blank())};theme_set(theme_mod())
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
site_x = df[which(df$Site == "1"),]
ggplot(site1, aes(x = Pd, y = mass, colour = as.factor(Site))) + 
  geom_point(size =2, position = "jitter")
levels(site_x$date)
x = subset(df, Site == "1" & date == "19-12-17")
hist(x$bodySize, breaks = 30)
##### done with QAQC #####

#run the code to estimate growth rates
df = read.csv(file = "./infrequensMeasurements.csv",T)
#load in the function
source("./boot_function.R")
tic()
suppressWarnings(cohort_boot(df, nboot = 2))
toc()
#907 secs regular
#Run parallel
#debugonce(cohort_boot)


library(parallel)
num_cores <- detectCores()
parallel_cluster<-makeCluster(num_cores-1) # without the fork
tic();parLapply(parallel_cluster, df, cohort_boot);toc()
stopCluster(parallel_cluster)
tic()
parallel_cluster<-makeCluster(num_cores-1) # without the fork
df = read.csv(file = "./infrequensMeasurements.csv",T)
#load in the function
source("./boot_function.R")
clusterExport(parallel_cluster, varlist = c("cohort_boot","df"))
suppressWarnings(cohort_boot(df, nboot = 1000))
stopCluster(parallel_cluster)
toc()
#take a look at one example
x = read.csv(file = "./output/infrequens_site-1_IGR.csv",T)
unique(levels(x$start_date))
x_1 = subset(x, site == 1 & start_date == "2017-10-15")

hist(x_1$IGR)
median(x_1$IGR)
quantile(x_1$IGR, c(0.0275,0.975))

site_date = paste(x_1$site,"_",x_1$start_date, sep = "")
df_int = data.frame(site_date = site_date, IGR = x_1$IGR)

f <-  function(u){
  which.min(abs(as.numeric(u$vals) - 0.5))
}

ecdf.plot <- function(DATA)
{
  
  x = sort(unique(DATA[,2]))
  vals = cumsum(tabulate(match(DATA[,2] , unique(DATA[,2]))))/length(DATA[,2])
  df = data.frame(x, vals)
  ggplot(df, aes(x = vals, y = x)) + geom_point(size = 3, shape = 19, colour = "#999999") + geom_line(lwd = 1.2) +
    labs(x = "Cumulative Frequency", y = "IGR (d-1)") + scale_y_continuous() +
    scale_x_continuous(limits = c(0,1)) +
    geom_hline(yintercept = as.numeric(df$x[f(df)])) +
    annotate("text", x = 0.25, y = 0.2, label = as.character(levels(droplevels(DATA$site_date)))) +
    annotate("text", x = 0.25, y = 0.18, label = as.character(paste("igr50 = ",round(df$x[f(df)],4), sep = "")))
}

ecdf.plot(df_int)

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


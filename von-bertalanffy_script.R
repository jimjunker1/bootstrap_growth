# modeling growth of bugs to get get some estimates of 

df = read.csv(file = "./infrequensMeasurements.csv",T)
df$Pd = as.POSIXct(df$date, format = "%d-%m-%y")
df$julian = as.numeric(format(df$Pd, "%j"))

site5 = df[which(df$Site == "5"),]
ggplot(site5, aes(x = Pd, y = log(mass), colour = as.factor(Site))) + 
  geom_point(size =2, position = "jitter")

#reorder dates based start of cohort
site5 %>%
  group_by(date) %>%
  summarize(mean_mass = mean(mass, na.rm=T), julian = unique(julian, na.rm = T)) %>%
  arrange(mean_mass) %>%
  mutate(date = reorder(date, mean_mass)) %>%
  mutate(day = c(diff(julian),NA), id = 1:n()) -> site5_df

#fix the weird wrap around dates
wrap_loc = which(site5_df$day < 0);wrap_fix = 365-abs(site5_df[which(site5_df$day < 0),'day'])
site5_df[wrap_loc, 'day'] = wrap_fix

site5_df %>%
  mutate(cohort_day = cumsum(day)) %>%
  mutate(cohort_day = c(0,cohort_day[!is.na(cohort_day)]))->site5_df

site5 %>%
  left_join(site5_df[,c('date','cohort_day')]) -> site5
ggplot(site5, aes(x = cohort_day, log(mass))) + geom_point(size = 2, position = 'jitter') +
  geom_smooth(span = 0.9)




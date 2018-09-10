# Bootstrap function
boots = function(allData, nboot,...){
  #browser()
  # for each sampling day resample data
  boot.data = c()
  for (j in unique(levels(allData$date))){ # for each date
    
    date.sub <- subset(allData, date == j)
    #resample rows #s of date.sub
    date.samp <- sample( 1:length(date.sub[,1]), size=nboot, replace=TRUE )
    #compile resampled data
    date.samp <- date.sub[date.samp,]
    boot.data = rbind(boot.data,date.samp)
  }
  
  boot.data %>%
    select(species, Site, date, mass) %>%
    mutate(id = rep(1:nboot, length(unique(levels(allData$date))))) %>%
    spread(date,mass) %>%
    select(-id)-> boot.data_wide
  # make sure mass change is positive
  for(k in 4:dim(boot.data_wide)[2]){
    for(l in 1:nrow(boot.data_wide)){
      if(boot.data_wide[l,k] < boot.data_wide[l,(k-1)]){ 
        date.sub = subset(allData, date == as.character(names(boot.data_wide)[k]))
        x=1
          repeat {
        date.samp = date.sub[sample(1:length(date.sub[,1]),1, replace = TRUE),]
        x = x+1#date.samp = date.sub[date.samp,]
        if(date.samp[,'mass'] >= boot.data_wide[l,k-1] | x == 300){
          boot.data_wide[l,k] <- date.samp
          break
          }}
      }
    }
    
  }
  #clean up and output data
  rownames(boot.data) = NULL
  return (boot.data)
}
# Bootstrap function
boots = function(allData, nboot,...){
  #browser()
  # for each sampling day resample data
  boot.data = 0.
  for (j in unique(levels(allData$date))){ # for each date
    
    date.sub <- subset(allData, date == j)
    #resample rows #s of date.sub
    date.samp <- sample( 1:length(date.sub[,1]), size=nboot, replace=TRUE )
    #compile resampled data
    date.samp <- date.sub[date.samp,]
    if (j==1) {
      boot.data = date.samp
    } else {
      boot.data = rbind(boot.data,date.samp)
    }
  }
  #clean up and output data
  rownames(boot.data) = NULL
  return (boot.data)
}
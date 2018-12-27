#parallel test script 

library(furrr)
library(purrr)

df = read.csv(file = "./infrequensMeasurements.csv",T)

taxa= list(unique(levels(df$TAXON)))
sites = list(unique(df$SITE))


lapply(taxa, function(x) unique(df[x,"SITE"]))#list(unique(x[,"SITE"])))

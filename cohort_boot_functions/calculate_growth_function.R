calculate_growth = function(bootdata, cohort_date_list){
  growth_rate <- function(x,y,z) {(log(y/x))/z}
  map2(bootdata, cohort_date_list, function(x,y){
    igr_df = c()
    for(date_id in 1:(max(x$id)-1)){
      igr_vec = mapply(growth_rate, x[which(x$id == date_id),'MASS'], x[which(x$id == date_id+1),'MASS'],
                       y[date_id,'day'])
      igr_df = cbind(igr_df,igr_vec)}
  colnames(igr_df) = unique(levels(y$DATE))[1:(nrow(y)-1)]
  igr_df = data.frame(TAXON = as.character(unique(x$TAXON)), SITE = as.character(unique(x$SITE)), igr_df, check.names = F)
  igr_df %>%
    gather(start_date, IGR, 3:(dim(igr_df)[2]), factor_key = T) -> igr_long
  igr_fix = which(igr_long$IGR <= 0)
  igr_long[igr_fix, 'IGR'] = 0.001
  TAXON = as.character(unique(droplevels(x$TAXON)))
  SITE = as.character(unique(droplevels(as.factor(x$SITE))))
  write.csv(igr_long, file = paste("./output/",TAXON,"_site-",SITE,"_IGR.csv",sep = ""), row.names = F) 
  return(igr_long)})
}
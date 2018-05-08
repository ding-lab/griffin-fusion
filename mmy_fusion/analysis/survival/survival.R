#Primary samples only time to event and censoring status

#Event free survival
return_efs_survival <- function(clinical_df, mmrf){
  survival_time <- as.numeric(subset(clinical_df, MMRF == mmrf & Clinical_category == "EFS")$Clinical_value)
  censored <- as.numeric(as.character(subset(clinical_df, MMRF==mmrf & Clinical_category == "EFS_censor")$Clinical_value))
  return(c(mmrf, survival_time, censored))
}
efs_survival_df <- as.data.frame(t(apply(samples, 1, function(x) return_efs_survival(clinical_df, x[1]))), stringsAsFactors=FALSE)
names(efs_survival_df) <- c("MMRF","EFS_survival","EFS_censored")

#Death
return_death_survival <- function(clinical_df, mmrf){
  D_PT_deathdy <- subset(clinical_df, MMRF == mmrf & Clinical_category == "D_PT_deathdy")$Clinical_value
  D_PT_lstalive <- subset(clinical_df, MMRF == mmrf & Clinical_category == "D_PT_lstalive")$Clinical_value
  
  if( is.na(D_PT_deathdy) ){
    survival_time <- D_PT_lstalive
    censored <- 1
  } else{
    survival_time <- D_PT_deathdy
    censored <- 0
  }
  
  return(c(mmrf, survival_time, censored))  
}
death_survival_df <- as.data.frame(t(apply(samples, 1, function(x) return_death_survival(clinical_df, x[1]))), stringsAsFactors=FALSE)
names(death_survival_df) <- c("MMRF","Death_survival","Death_censored")

#Recurrent (how many?) fusions, in order of recurrence
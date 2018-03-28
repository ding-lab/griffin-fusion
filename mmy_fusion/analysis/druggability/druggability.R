#Use DEPO annotation to look at druggable fusions

library(ggplot2)

#High risk only
seqfish <- read.table(paste0(data,"seqfish_df.tsv"), header=T)
clinical <- read.table(paste0(data,"clinical_df.tsv"), header=T)
samples <- read.table(paste0(data,"sample_list.807.txt"), header=F)
names(samples) <- c("mmrf","srr")

get_high_risk_samples <- function(){
  high_risk_list <- NULL
  for(mmrf in sort(unique(samples$mmrf))){
    
    if( mmrf %in% seqfish$MMRF & any(endsWith(as.character(subset(seqfish, MMRF==mmrf)$Study_Visit_ID), "BM"))){
      del17p <- as.logical(subset(seqfish, MMRF==mmrf & seqfish_category=="CN_del_17p13" & endsWith(as.character(seqfish$Study_Visit_ID), "BM") )$has_event)[1]
      t1416 <- as.logical(subset(seqfish, MMRF==mmrf & seqfish_category=="Translocation_MAF_14_16" & endsWith(as.character(seqfish$Study_Visit_ID), "BM") )$has_event)[1]
      t1420 <- as.logical(subset(seqfish, MMRF==mmrf & seqfish_category=="Translocation_MAFB_14_20" & endsWith(as.character(seqfish$Study_Visit_ID), "BM") )$has_event)[1]
    } else{
      del17p <- NA
      t1416 <- NA
      t1420 <- NA
    }
    
    ldh <- as.numeric(subset(clinical, MMRF==mmrf & Clinical_category=="LDH")$Clinical_value)
    if( is.na(ldh) ){
      ldh_high <- NA
    } else{
      ldh_high <- (ldh > 500)
    }
    
    efs_days <- as.numeric(as.character(subset(clinical, MMRF==mmrf & Clinical_category=="EFS")$Clinical_value))
    efs_censored <- as.logical(as.numeric(subset(clinical, MMRF==mmrf & Clinical_category=="EFS_censor")$Clinical_value))
    if( is.na(efs_days) | is.na(efs_censored) ){
      efs_lt1 <- NA
    } else{
      efs_lt1 <- (efs_days < 365) & !efs_censored 
    }
    
    iss_stage <- as.numeric(as.character(subset(clinical, MMRF==mmrf & Clinical_category=="ISS_Stage")$Clinical_value))
    if( is.na(iss_stage) ){
      iss_stage3 <- NA
    } else{
      iss_stage3 <- (iss_stage == 3)
    }
    
    high_risk <- c(del17p, t1416, t1420, ldh_high, efs_lt1, iss_stage3)
    high_risk <- high_risk[!is.na(high_risk)]
    
    if( length(high_risk) == 6 | (length(high_risk) > 1 & sum(high_risk) > 1) ){
      total <- total + 1
      if( sum(high_risk) > 1 ){
        count <- count + 1
        high_risk_list <- c(high_risk_list, mmrf)
      }
    }
    
  }
  return(high_risk_list)
}

high_risk_samples <- get_high_risk_samples()

druggable_high_risk_samples <- NULL
for(this_mmrf in high_risk_samples){
  this_df <- subset(primary_df, mmrf == this_mmrf)
  if(sum(this_df$drug_fusion | this_df$drug_geneA | this_df$drug_geneB) > 0){
    druggable_high_risk_samples <- c(druggable_high_risk_samples, this_mmrf)
  }
}

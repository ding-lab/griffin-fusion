library(survival)
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

#Recurrent (more than 1) fusions, in order of recurrence
#broken down by GeneA--GeneB
fusion_pairs <- names(sort(table(as.character(primary_df$fusion)), decreasing=TRUE)[sort(table(as.character(primary_df$fusion)), decreasing=TRUE) > 1])
#all genes involved in fusions GeneA or GeneB
fusion_genes <- names(sort(table(c(as.character(primary_df$geneA), as.character(primary_df$geneB))), decreasing=TRUE)[sort(table(c(as.character(primary_df$geneA), as.character(primary_df$geneB))), decreasing=TRUE) > 1])

fusion_survival <- function(primary_df, this_fusion, efs_survival_df, death_survival_df){
  if( grepl("--", this_fusion) ){
    mmrf_with_fusion <- as.character(subset(primary_df, as.character(fusion)==this_fusion)$mmrf)
  } else{
    mmrf_with_fusion <- as.character(subset(primary_df, as.character(geneA)==this_fusion | as.character(geneB)==this_fusion )$mmrf)
  }
  fusion_indicator <- as.numeric(efs_survival_df$MMRF %in% mmrf_with_fusion)
  efs_survival_object <- Surv(as.numeric(efs_survival_df$EFS_survival), as.numeric(efs_survival_df$EFS_censored))
  death_survival_object <- Surv(as.numeric(death_survival_df$Death_survival), as.numeric(death_survival_df$Death_censored))
  efs_survival_coxph <- summary(coxph(efs_survival_object ~ fusion_indicator))
  death_survival_coxph <- summary(coxph(death_survival_object ~ fusion_indicator))
  return(c(this_fusion, length(mmrf_with_fusion), as.vector(efs_survival_coxph$coefficients), as.vector(death_survival_coxph$coefficients)))
}

seqfish_survival <- function(primary_df, this_fish, efs_survival_df, death_survival_df){
  this_index <- which(names(primary_df) == paste0("seqfish_", this_fish), arr.ind=TRUE)
  mmrf_with_fish <- unique(as.character(subset(primary_df, !is.na(primary_df[,this_index]) & primary_df[,this_index]==1 )$mmrf))
  fish_indicator <- as.numeric(efs_survival_df$MMRF %in% mmrf_with_fish)
  efs_survival_object <- Surv(as.numeric(efs_survival_df$EFS_survival), as.numeric(efs_survival_df$EFS_censored))
  death_survival_object <- Surv(as.numeric(death_survival_df$Death_survival), as.numeric(death_survival_df$Death_censored))
  efs_survival_coxph <- summary(coxph(efs_survival_object ~ fish_indicator))
  death_survival_coxph <- summary(coxph(death_survival_object ~ fish_indicator))
  return(c(this_fish, length(mmrf_with_fish), as.vector(efs_survival_coxph$coefficients), as.vector(death_survival_coxph$coefficients)))
}

fusion_pair_survival_df <- as.data.frame(t(apply(as.matrix(fusion_pairs), 1, function(x) fusion_survival(primary_df, x[1], efs_survival_df, death_survival_df))))
names(fusion_pair_survival_df) <- c("fusion_pair","number_with_fusion_pair",paste0("EFS_",c("coef","exp(coef)","se(coef)","z","p-value")), paste0("Death_",c("coef","exp(coef)","se(coef)","z","p-value")))

fusion_gene_survival_df <- as.data.frame(t(apply(as.matrix(fusion_genes), 1, function(x) fusion_survival(primary_df, x[1], efs_survival_df, death_survival_df))))
names(fusion_gene_survival_df) <- c("fusion_gene","number_with_fusion_gene",paste0("EFS_",c("coef","exp(coef)","se(coef)","z","p-value")), paste0("Death_",c("coef","exp(coef)","se(coef)","z","p-value")))

#Summarize top 10 fusions
x <- subset(fusion_pair_survival_df, as.numeric(as.character(number_with_fusion_pair)) > 4)
head(x[order(x$`EFS_p-value`),], 10)
head(x[order(x$`Death_p-value`),], 10)

y <- subset(fusion_gene_survival_df, as.numeric(as.character(number_with_fusion_gene)) > 4)
head(y[order(y$`EFS_p-value`),], 10)
head(y[order(y$`Death_p-value`),], 10)

#SeqFISH
seqfish_events <- unique(as.character(seqfish_df$seqfish_category))
fish_survival_df <- as.data.frame(t(apply(as.matrix(seqfish_events), 1, function(x) seqfish_survival(primary_df, x[1], efs_survival_df, death_survival_df))))
names(fish_survival_df) <- c("fish_event","number_with_fish_event",paste0("EFS_",c("coef","exp(coef)","se(coef)","z","p-value")), paste0("Death_",c("coef","exp(coef)","se(coef)","z","p-value")))
fish_survival_df
setwd("/Users/sfoltz/Desktop/griffin-fusion/mmy_fusion/analysis")
data="/Users/sfoltz/Desktop/griffin-fusion/mmy_fusion/data/"

fusion_df <- read.table(paste0(data,"fusion_df.txt"),header=T)
all_samples <- read.table(paste0(data,"sample_list.807.txt"), header=F)
samples <- read.table(paste0(data,"sample_list.primary.txt"), header=F)
expr <- read.table(paste0(data,"mmy_gene_expr_with_fusions.tsv"),header=T)
primary_df <- subset(fusion_df, sample_number==1)
clinical_df <- read.table(paste0(data,"clinical_df.tsv"), header=T)

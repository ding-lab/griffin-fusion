library(ggplot2)
library(ggrepel)

fusion_df <- read.table("data/fusion_df.txt", header=T, stringsAsFactors=F)
all_samples <- read.table("data/sample_list.806.txt", header=F, stringsAsFactors=F)
names(all_samples) <- c("MMRF","SRR")
samples <- read.table("data/sample_list.primary.txt", header=F, stringsAsFactors=F)
names(samples) <- c("MMRF","SRR")

expr <- read.table("data/mmy_gene_expr_with_fusions.tsv", header=T, stringsAsFactors=F)
expr <- subset(expr, srr %in% all_samples$SRR)
primary_expr <- subset(expr, srr %in% samples$SRR)

clinical_df <- read.table("data/clinical_df.tsv", header=T, stringsAsFactors=F)
seqfish_df <- read.table("data/seqfish_df.tsv", header=T, stringsAsFactors=F)

ericsript_twocallers <- c("EF", "EI", "EP", "ES")
bad_gene_list <- c("FOSB")
readthrough_fusion_list <- names(sort(table(subset(fusion_df, chrA == chrB & strandA == strandB & abs(posA-posB) < 250000)$fusion), decreasing=TRUE))

fosb_df <- subset(fusion_df, sample_number==1 & CallerN >= 2 & !( Callers %in% ericsript_twocallers) & geneA=="FOSB")
readthrough_df <- subset(fusion_df, sample_number==1 & CallerN >= 2 & !( Callers %in% ericsript_twocallers) & !(geneA %in% bad_gene_list) & !(geneB %in% bad_gene_list) & (fusion %in% readthrough_fusion_list) )
primary_df <- subset(fusion_df, sample_number==1 & CallerN >= 2 & !( Callers %in% ericsript_twocallers) & !(geneA %in% bad_gene_list) & !(geneB %in% bad_gene_list) & !(fusion %in% readthrough_fusion_list) )

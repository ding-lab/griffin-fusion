#Plot fusion + expression data

args <- commandArgs(trailingOnly=TRUE)

df <- read.table(args[1],header=T)
over_outlier_level <- paste0("overexpression",args[2])
under_outlier_level <- paste0("underexpression",args[2])
plot_file <- args[3]
title <- args[4]
ylabel <- args[5]

gene_id <- df$gene[1]

levels(df$cnv) <- c("Deep Deletion","Shallow Deletion","Neutral","Low Amplification","High Amplification","No Data")
df$cnv[df$cnv=="-2"] <- "Deep Deletion"
df$cnv[df$cnv=="-1"] <- "Shallow Deletion"
df$cnv[df$cnv=="0"] <- "Neutral"
df$cnv[df$cnv=="1"] <- "Low Amplification"
df$cnv[df$cnv=="2"] <- "High Amplification"
df$cnv[is.na(df$cnv)] <- "No Data"
df$cnv <- factor(df$cnv, levels=c("Deep Deletion","Shallow Deletion","Neutral","Low Amplification","High Amplification","No Data"))

#fusion_status <- factor((df$fusion_info == "None_reported") + 2*(df$fusion_info == "Fusion_file_NA"))
#levels(fusion_status) <- c("Fusion", "None reported", "No Data")

fusion_status_num <- 1 + as.numeric(df$fusion_info == "None_reported") + 2*as.numeric(df$fusion_info == "Fusion_file_NA")
fusion_status <- factor(c("Fusion","None reported","No Data")[fusion_status_num], levels=c("Fusion","None reported","No Data"), ordered=TRUE)


gene_partners <- rep(NA,nrow(df))

for(i in 1:nrow(df)){
  if( !(as.character(df$fusion_info[i]) %in% c("None_reported","Fusion_file_NA")) ){
    genes_a <- strsplit(as.character(df$geneA[i]),';')[[1]]
    genes_b <- strsplit(as.character(df$geneB[i]),';')[[1]]
    gene_partners[i] <- paste( unique(sort(paste(genes_a,genes_b,sep="--"))), collapse='\n')
  }
}

set.seed(10)

outlier_status <- apply(df[,c(over_outlier_level, under_outlier_level)], 1, sum)
alpha_level <- outlier_status/2 + 0.5
plot_df <- data.frame(df[,c("gene","cnv","expression_level")], fusion_status=fusion_status, fusion_status_j=jitter(as.numeric(fusion_status)), alpha_level=alpha_level, gene_partners=gene_partners)[df$fusion_info != "Fusion_file_NA",]

if( any(plot_df$fusion_status == "Fusion") ){
  library(ggplot2)
  library(ggrepel)
  #fusion_median <- median(subset(plot_df, fusion_status=="Fusion")$expression_level)
  #nonfusion_median <- median(subset(plot_df, fusion_status=="None reported")$expression_level)
  #median_df <- data.frame(fusion_status=c("Fusion","None reported"), median=c(fusion_median, nonfusion_median))

  p <- ggplot(plot_df, aes(x=fusion_status, y=expression_level, color=cnv, fill=cnv))
  
  p <- p + geom_violin(aes(fill=NULL),color="black", draw_quantiles=c(0.5))
  p <- p + geom_point(aes(x=fusion_status_j), alpha=0.75) #, alpha=alpha_level))
  p  <- p + geom_label_repel(data=subset(plot_df, fusion_status=="Fusion"), aes(x=fusion_status_j, y=expression_level, label=gene_partners), label.size=NA, color=c('white','black')[as.numeric(subset(plot_df, fusion_status=="Fusion")$cnv=="No Data")+1], box.padding=0.35, point.padding=0.5, segment.color="grey50")
  p <- p + ylim(0, max(df$expression_level))
  p <- p + theme_bw(base_size=20)
  p <- p + scale_color_brewer(palette="Set1", drop=FALSE)
  p <- p + scale_fill_brewer(palette="Set1", drop=FALSE)
  p <- p + labs(x="Fusion Status", y=ylabel, title=title, color="CNV Status")
  p <- p + guides(fill=FALSE, alpha=FALSE)

  pdf(plot_file, 10, 10, useDingbats = FALSE)
  print(p)
  shh <- dev.off()
}


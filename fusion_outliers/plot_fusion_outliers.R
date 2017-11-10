#Plot fusion + expression data

args <- commandArgs(trailingOnly=TRUE)

df <- read.table(args[1],header=T)
over_outlier_level <- paste0("overexpression",args[2])
under_outlier_level <- paste0("underexpression",args[2])
plot_file <- args[3]
title <- args[4]
ylabel <- args[5]

gene_id <- df$gene[1]

df$cnv <- factor(df$cnv)

fusion_status <- factor((df$fusion_info == "None_reported") + 2*(df$fusion_info == "Fusion_file_NA"))
levels(fusion_status) <- c("Fusion", "None reported", "NA")

five_prime_three_prime <- factor(rep(NA,nrow(df)))
levels(five_prime_three_prime) <- c("5' and 3'", "5' only", "3' only", "None", "NA")
gene_partners <- rep(NA,nrow(df))

for(i in 1:nrow(df)){
  if(as.character(df$fusion_info[i]) == "None_reported"){
    five_prime_three_prime[i] <- "None"
  } else if(as.character(df$fusion_info[i]) == "Fusion_file_NA"){
    five_prime_three_prime[i] <- "NA"
  } else{
    a <- (as.character(df$gene[i]) %in% strsplit(as.character(df$geneA[i]),';')[[1]])
    b <- (as.character(df$gene[i]) %in% strsplit(as.character(df$geneB[i]),';')[[1]])
    if(a & b){
      five_prime_three_prime[i] <- "5' and 3'"
      all_partners <- unique(sort(c(unique(sort(strsplit(as.character(df$geneA[i]),';')[[1]])), unique(sort(strsplit(as.character(df$geneB[i]),';')[[1]])))))
      gene_partners[i] <- paste(all_partners[all_partners != as.character(df$gene[i])], collapse=",")
    } else if(a & !b){
      five_prime_three_prime[i] <- "5' only"
      all_partners <- unique(sort(strsplit(as.character(df$geneB[i]),';')[[1]]))
      gene_partners[i] <- paste(all_partners[all_partners != as.character(df$gene[i])], collapse=",")
    } else if(!a & b){
      five_prime_three_prime[i] <- "3' only"
      all_partners <- unique(sort(strsplit(as.character(df$geneA[i]),';')[[1]]))
      gene_partners[i] <- paste(all_partners[all_partners != as.character(df$gene[i])], collapse=",")
    } else{
      exit("Fusion but neither 5' nor 3'")
    }
  }
}

set.seed(10)

outlier_status <- apply(df[,c(over_outlier_level, under_outlier_level)], 1, sum)
alpha_level <- outlier_status/2 + 0.5
plot_df <- data.frame(df[,c("gene","cnv","expression_level")], fusion_status=fusion_status, fusion_status_j=jitter(as.numeric(fusion_status)), alpha_level=alpha_level, five_prime_three_prime, gene_partners=gene_partners)

library(ggplot2)
library(ggrepel)

#p <- ggplot(plot_df, aes(x=fusion_status, y=expression_level, color=five_prime_three_prime))
p <- ggplot(plot_df, aes(x=fusion_status, y=expression_level, color=cnv))
p <- p + geom_violin(color="black")
p <- p + geom_point(data=subset(plot_df, fusion_status != "Fusion"), aes(x=fusion_status_j, alpha=alpha_level))
p <- p + geom_text(data=subset(plot_df, fusion_status == "Fusion"), aes(x=fusion_status_j, alpha=alpha_level, label=five_prime_three_prime))
p <- p + geom_text_repel(data=subset(plot_df, fusion_status=="Fusion"), aes(x=fusion_status_j, y=expression_level, label=gene_partners, alpha=alpha_level), show.legend=FALSE)
p <- p + ylim(0, max(df$expression_level))
p <- p + theme_bw(base_size=20)
p <- p + scale_color_brewer(palette="Set1", drop=FALSE) #, breaks=c("5' and 3'", "5' only", "3' only"))
p <- p + labs(x="Fusion Status", y=ylabel, title=title, color=paste0(gene_id,"\nPosition"))
p <- p + guides(alpha=FALSE)

pdf(plot_file, 10, 10)
p
shh <- dev.off()


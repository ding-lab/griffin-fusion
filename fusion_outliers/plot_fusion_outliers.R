#Plot fusion + expression data

args <- commandArgs(trailingOnly=TRUE)

df <- read.table(args[1],header=T)
outlier_level <- paste0("overexpression",args[2])
plot_file <- args[3]
title <- args[4]
ylabel <- args[5]

fusion_status <- factor((df$fusion_info == "None_reported") + 2*(df$fusion_info == "Fusion_file_NA"))
levels(fusion_status) <- c("Fusion", "None_reported", "NA")

five_prime_three_prime <- rep(NA,nrow(df))
for(i in 1:nrow(df)){
  if(as.character(df$fusion_info[i]) == "None_reported"){
    five_prime_three_prime[i] <- "None"
  } else if(as.character(df$fusion_info[i]) == "Fusion_file_NA"){
    five_prime_three_prime[i] <- "NA"
  } else{
    a <- (as.character(df$gene[i]) %in% strsplit(as.character(df$geneA[i]),';')[[1]])
    b <- (as.character(df$gene[i]) %in% strsplit(as.character(df$geneB[i]),';')[[1]])
    if(a & b){
      five_prime_three_prime[i] <- "Both"
    } else if(a & !b){
      five_prime_three_prime[i] <- "5 prime"
    } else if(!a & b){
      five_prime_three_prime[i] <- "3 prime"
    } else{
      exit("Fusion but neither 5' nor 3'")
    }
  }
}

#pval <- t.test(df$expression_level[fusion_status=="Fusion"], df$expression_level[fusion_status=="None_reported"])$p.value
#if(pval < 0.0005){
#  small_pval <- "p < 0.0005"
#}

alpha_level <- df[,outlier_level]/2 + 0.5
plot_df <- data.frame(df[,c("gene","expression_level",outlier_level)], fusion_status=fusion_status, alpha_level=alpha_level, five_prime_three_prime)

library(ggplot2)

set.seed(10)
p <- ggplot(plot_df, aes(x=fusion_status, y=expression_level, color=five_prime_three_prime)) #color=fusion_status))
p <- p + geom_violin(color="black")
p <- p + geom_jitter(aes(alpha=alpha_level), height=0, width=0.25)
p <- p + ylim(0, max(df$expression_level))
p <- p + theme_bw(base_size=20)
p <- p + scale_color_brewer(palette="Set1")
p <- p + labs(x="Fusion Status", y=ylabel, title=title, color="Fusion")
p <- p + guides(alpha=FALSE) #, color=FALSE)
#if(pval < 0.0005){
#  p <- p + annotate("text", x=1,y=0,hjust=0.5,vjust=-0.2,label=small_pval)
#} else {
#  p <- p + annotate("text", x=1,y=0,hjust=0.5,vjust=-0.2,label=paste0("p = ", sprintf("%.3f",round(pval,3))))
#}

pdf(plot_file, 10, 10)
p
shh <- dev.off()


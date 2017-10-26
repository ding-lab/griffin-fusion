#Plot fusion + expression data

args <- commandArgs(trailingOnly=TRUE)

df <- read.table(args[1],header=T)
outlier_level <- paste0("overexpression",args[2])
plot_file <- args[3]
title <- args[4]
ylabel <- args[5]

fusion_status <- factor((df$fusion_info == "None") + 2*(df$fusion_info == "Fusion_file_NA"))
levels(fusion_status) <- c("Fusion", "None", "NA")

pval <- t.test(df$expression_level[fusion_status=="Fusion"], df$expression_level[fusion_status=="None"])$p.value

alpha_level <- df[,outlier_level]/2 + 0.5
plot_df <- data.frame(df[,c("gene","expression_level",outlier_level)], fusion_status=fusion_status, alpha_level=alpha_level)

library(ggplot2)

set.seed(10)
p <- ggplot(plot_df, aes(x=fusion_status, y=expression_level, color=fusion_status))
p <- p + geom_violin()
p <- p + geom_jitter(aes(alpha=alpha_level), height=0, width=0.25)
p <- p + ylim(0, max(df$expression_level))
p <- p + theme_bw(base_size=20)
p <- p + scale_color_brewer(palette="Set1")
p <- p + labs(x="Fusion Status", y=ylabel, title=title)
p <- p + guides(color=FALSE, alpha=FALSE)
p <- p + annotate("text", x=1,y=0,hjust=0.5,vjust=-0.2,label=paste0("p = ", sprintf("%.3f",round(pval,3))))

pdf(plot_file, 10, 10)
p
shh <- dev.off()


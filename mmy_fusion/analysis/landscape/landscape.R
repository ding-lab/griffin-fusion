library(ggplot2)

top_5_names <- names(sort(table(primary_df$fusion), decreasing=TRUE)[1:5])
top_5_count <- as.vector(sort(table(primary_df$fusion), decreasing=TRUE)[1:5])
top_5_rate <- vector()
for(name in top_5_names){
  val <- subset(primary_df, fusion==name)$n_discordant
  rate <- mean(val >= 3 , na.rm=T)
  top_5_rate <- c(top_5_rate, rate)
}

top_callers <- names(sort(table(primary_df$Callers), decreasing=TRUE))
top_callers_count <- as.vector(sort(table(primary_df$Callers), decreasing=TRUE))
top_callers_validation_rate <- vector()
for(name in top_callers){
  val <- subset(primary_df, Callers==name)$n_discordant
  rate <- mean(val >= 3 , na.rm=T)
  top_callers_validation_rate <- c(top_callers_validation_rate, rate)
}

hyperdiploid_status <- unique(primary_df[,c("mmrf", "seqfish_Hyperdiploidy")])

#plots

##### top 5 #####
#update coverge to show N of samples with WGS
plot_df <- data.frame(fusion=top_5_names, n_samples=top_5_count, validation_rate=format(top_5_rate,digits=1))
plot_df$fusion <- factor(plot_df$fusion, levels = rev(top_5_names))
ggplot(plot_df, aes(x=fusion, weight=n_samples, label=validation_rate)) + geom_bar() + coord_flip() + labs(x="Fusion name", y="Number of samples", title="Top 5 Recurrent Fusions") + theme_bw(base_size=20) + geom_text(aes(x = fusion, y = n_samples), hjust=1, color="white", size=5)

##### callers #####
#update coverge to show N of samples with WGS
plot_df <- data.frame(caller=top_callers, fusions=top_callers_count, validation_rate=format(top_callers_validation_rate,digits=1))
plot_df$caller <- factor(plot_df$caller, levels = rev(top_callers))
ggplot(plot_df, aes(x=caller, weight=fusions, label=validation_rate)) + geom_bar() + coord_flip() + labs(x="Fusion Caller Combination", y="Number of Fusions", title="Number of Fusions by Caller") + theme_bw(base_size=20) + geom_text(aes(x = caller, y = fusions), hjust=1, color="white", size=5)

##### distribution of fusions #####
n_fusions = vector()
for(sample in unique(samples[,1])){
  n_fusions <- c(n_fusions, dim(subset(primary_df, mmrf==sample))[1])
}
plot_df <- data.frame(srr = unique(samples[,1]), n_fusions = n_fusions)
med <- median(plot_df$n_fusions)
ggplot(plot_df, aes(x=n_fusions)) + geom_histogram(breaks=seq(min(plot_df[,2])-0.5,max(plot_df[,2])+0.5,1)) + labs(x="Number of Fusions", y="Number of Samples", title="Number of Fusions per Sample") + theme_bw(base_size=20) + geom_vline(xintercept=med, lty=2, color="white", size=1) + geom_text(label=paste0("Median=",med), x=med, y=50, angle=90, vjust=-1, color="white", size=5)

##### coverage by validation #####
plot_df <- subset(primary_df, !is.na(n_discordant))[,c("depthA","depthB","n_discordant")]
plot_df <- data.frame(plot_df, mean_cov=(plot_df$depthA + plot_df$depthB)/2, val=(plot_df$n_discordant>0))
plot_df$mean_cov[is.na(plot_df$mean_cov)] <- 0
ggplot(plot_df, aes(x=val, y=mean_cov)) + geom_violin(size=2) + theme_bw(base_size=20) + labs(x="Fusion Validated by WGS", y="WGS Coverage at Fusion", title="Fusion Validation and WGS Coverage") + geom_boxplot(width=.1, outlier.alpha = 0)

##### samples per fusion #####
plot_df <- data.frame(table(primary_df$fusion))
ggplot(plot_df, aes(x=Freq)) + geom_histogram(breaks=seq(min(plot_df[,2])-0.5,max(plot_df[,2])+0.5,1)) + theme_bw(base_size=20) + labs(x="Number of Samples", y="Number of Fusions", title="Number of Samples per Fusion")
ggplot(subset(plot_df, Freq>1), aes(x=Freq)) + geom_histogram(breaks=seq(min(plot_df[,2])-0.5,max(plot_df[,2])+0.5,1)) + theme_bw(base_size=20) + labs(x="Number of Samples", y="Number of Fusions", title="Number of Samples per Fusion")

##### FGFR3 WHSC1 #####
srr_with_fgfr3_fusion <- unique(subset(primary_df, (geneB %in% c("FGFR3","IGH@") & geneA %in% c("FGFR3","IGH@")))$srr) #| (geneB %in% c("FGFR3") & geneA=="IGH@"))$srr)
srr_with_whsc1_fusion <- unique(subset(primary_df, (geneB %in% c("WHSC1","IGH@") & geneA %in% c("WHSC1","IGH@")))$srr) #| (geneB %in% c("WHSC1") & geneA=="IGH@"))$srr)
srr_with_t414 <- unique(subset(primary_df, primary_df$seqfish_Translocation_WHSC1_4_14 == 1)$srr)
fgfr3_log10tpm <- subset(expr, gene %in% c("FGFR3"))[,c("srr","log10tpm")]
whsc1_log10tpm <- subset(expr, gene %in% c("WHSC1"))[,c("srr","log10tpm")]
has_fgfr3_fusion=(fgfr3_log10tpm[,1] %in% srr_with_fgfr3_fusion)
has_whsc1_fusion=(fgfr3_log10tpm[,1] %in% srr_with_whsc1_fusion)
has_fusion = factor(c("No fusion", "FGFR3 fusion", "WHSC1 fusion", "Both fusions")[1+1*has_fgfr3_fusion + 2*has_whsc1_fusion], levels=c("No fusion", "FGFR3 fusion", "WHSC1 fusion", "Both fusions"))
plot_df <- data.frame(srr=fgfr3_log10tpm[,1], FGFR3=fgfr3_log10tpm[,2], WHSC1=whsc1_log10tpm[,2], has_fusion=has_fusion, has_t414=(fgfr3_log10tpm[,1] %in% srr_with_t414))
ggplot(plot_df, aes(x=WHSC1, y=FGFR3, color=has_fusion, shape=has_t414)) + geom_point(size=3) + theme_bw(base_size=20) + labs(color="IGH fusion status", shape="Has t(4;14)", title="Gene Expression and IGH--FGFR3/WHSC1 Fusions", x="WHSC1 Gene TPM (log10)", y="FGFR3 Gene TPM (log10)")

##### FGFR3 WHSC1 2 #####
srr_with_fgfr3_fusion <- unique(subset(primary_df, (geneB %in% c("FGFR3","IGH@") & geneA %in% c("FGFR3","IGH@")))$srr) #| (geneB %in% c("FGFR3") & geneA=="IGH@"))$srr)
srr_with_whsc1_fusion <- unique(subset(primary_df, (geneB %in% c("WHSC1","IGH@") & geneA %in% c("WHSC1","IGH@")))$srr) #| (geneB %in% c("WHSC1") & geneA=="IGH@"))$srr)
srr_with_t414 <- unique(subset(primary_df, primary_df$seqfish_Translocation_WHSC1_4_14 == 1)$srr)
srr_NA_t414 <- unique(subset(primary_df, is.na(primary_df$seqfish_Translocation_WHSC1_4_14))$srr)
fgfr3_log10tpm <- subset(expr, gene %in% c("FGFR3"))[,c("srr","log10tpm")]
whsc1_log10tpm <- subset(expr, gene %in% c("WHSC1"))[,c("srr","log10tpm")]
has_fgfr3_fusion=(fgfr3_log10tpm[,1] %in% srr_with_fgfr3_fusion)
has_whsc1_fusion=(fgfr3_log10tpm[,1] %in% srr_with_whsc1_fusion)
has_t414=(fgfr3_log10tpm[,1] %in% srr_with_t414)
has_t414_NA=(fgfr3_log10tpm[,1] %in% srr_NA_t414)
has_t414[has_t414_NA] <- "Not available"
has_fusion = factor(c("No fusion", "FGFR3 fusion", "WHSC1 fusion", "Both fusions")[1+1*has_fgfr3_fusion + 2*has_whsc1_fusion], levels=c("No fusion", "FGFR3 fusion", "WHSC1 fusion", "Both fusions"))
plot_df <- data.frame(srr=fgfr3_log10tpm[,1], FGFR3=fgfr3_log10tpm[,2], WHSC1=whsc1_log10tpm[,2], has_fusion=has_fusion, has_t414=has_t414)
ggplot(plot_df, aes(x=WHSC1, y=FGFR3, color=has_fusion, shape=has_t414)) + geom_point(alpha=0.5, size=3) + theme_bw(base_size=20) + labs(color="IGH fusion status", shape="Has t(4;14)", title="Gene Expression and IGH--FGFR3/WHSC1 Fusions", x="WHSC1 Gene TPM (log10)", y="FGFR3 Gene TPM (log10)")

library(ggplot2)

##### FGFR3 WHSC1 #####
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
has_t414 <- factor(has_t414, levels=c(TRUE, FALSE, "Not available"))
has_fusion = factor(c("No fusion", "FGFR3 fusion", "WHSC1 fusion", "Both fusions")[1+1*has_fgfr3_fusion + 2*has_whsc1_fusion], levels=c("No fusion", "FGFR3 fusion", "WHSC1 fusion", "Both fusions"))
plot_df <- data.frame(srr=fgfr3_log10tpm[,1], FGFR3=fgfr3_log10tpm[,2], WHSC1=whsc1_log10tpm[,2], has_fusion=has_fusion, has_t414=has_t414)
ggplot(plot_df, aes(x=WHSC1, y=FGFR3, color=has_fusion, shape=has_t414)) + geom_point(alpha=1, size=3) + theme_bw(base_size=20) + labs(color="IGH fusion status", shape="Has t(4;14)", title="Gene Expression and IGH--FGFR3/WHSC1 Fusions", x="WHSC1 Gene TPM (log10)", y="FGFR3 Gene TPM (log10)") + scale_shape_manual(values=c(16,17,5))


##### Breakpoints of IGH@--WHSC1 fusions #####
srr_with_whsc1_fusion_plus_wgs <- subset(primary_df, fusion=="IGH@--WHSC1" & !is.na(seqfish_Translocation_WHSC1_4_14))$srr
bpA <- subset(primary_df, fusion=="IGH@--WHSC1" & !is.na(seqfish_Translocation_WHSC1_4_14))$posA
bpB <- subset(primary_df, fusion=="IGH@--WHSC1" & !is.na(seqfish_Translocation_WHSC1_4_14))$posB
seqfish_t414 <- subset(primary_df, fusion=="IGH@--WHSC1" & !is.na(seqfish_Translocation_WHSC1_4_14))$seqfish_Translocation_WHSC1_4_14
wgs_bam <- subset(primary_df, fusion=="IGH@--WHSC1" & !is.na(seqfish_Translocation_WHSC1_4_14))$wgs_bam
wgs_reads <- subset(primary_df, fusion=="IGH@--WHSC1" & !is.na(seqfish_Translocation_WHSC1_4_14))$discordant_reads
fgfr3_expr <- NULL
fgfr3_pct <- NULL
fgfr3_cnv <- NULL
fgfr3_over <- NULL
for(this_srr in srr_with_whsc1_fusion_plus_wgs){
  this_line <- subset(expr, srr==this_srr & gene=="FGFR3")
  fgfr3_expr <- c(fgfr3_expr, this_line$log10tpm)
  fgfr3_pct <- c(fgfr3_pct, this_line$pct)
  fgfr3_cnv <- c(fgfr3_cnv, this_line$gene_avg_cnv)
  fgfr3_over <- c(fgfr3_over, this_line$outlier_over_log10tpm)
}
plot_df <- data.frame(srr=srr_with_whsc1_fusion_plus_wgs, breakpointA=bpA, breakpointB=bpB, fgfr3_expr=fgfr3_expr, fgfr3_cnv=fgfr3_cnv, fgfr3_pct=fgfr3_pct, seqfish_t414=as.factor(seqfish_t414), wgs_bam=wgs_bam, wgs_reads=wgs_reads, category=as.factor(round(fgfr3_expr)), over=as.factor(fgfr3_over))
ggplot(plot_df, aes(x=as.factor(breakpointA), y=as.factor(breakpointB), color=category)) + geom_jitter(width=0.2, height=0.2) + theme_bw(base_size=20)



ggplot(plot_df, aes(y=c("chr14","chr4"))) + geom_segment(y="chr14", yend="chr4", x=plot_df$breakpointB, xend=plot_df$breakpointA)

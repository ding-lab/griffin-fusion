library(ggplot2)

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

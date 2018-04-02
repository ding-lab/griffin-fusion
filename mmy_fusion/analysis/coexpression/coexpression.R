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
delly <- read.table(paste0(data,"Filtered_Fusions_100000_delly_20180219.txt"), header=T)
process_sv_evidence <- function(sv_ev){
  sv_ev <- as.character(sv_ev)
  if(is.na(sv_ev)){
    return(NA)
  } else{
    sv_list <- strsplit(sv_ev,"|", fixed=TRUE)[[1]][!apply(as.matrix(strsplit(sv_ev,"|", fixed=TRUE)[[1]]), 1, function(x) grepl("*",x,fixed=TRUE))]
    chr <- NULL
    pos <- NULL
    for(sv in sv_list){
      chr <- c(chr, strsplit(sv,":")[[1]][1])
      pos <- c(pos, strsplit(strsplit(sv,":")[[1]][2],"-")[[1]][1])
    }
  }
  return(data.frame(chr=chr, pos=pos))
}

igh_whsc1_df <- subset(delly, FusionName=="IGH@--WHSC1")
sv_breakpoints <- apply(as.matrix(igh_whsc1_df$DELLY_SV_EVIDENCE), 1, function(x) process_sv_evidence(x))

bp_igh <- NULL
bp_whsc1 <- NULL
ex_fgfr3 <- NULL
ex_whsc1 <- NULL
for(this_srr in igh_whsc1_df$Sample){
  sv_bp <- process_sv_evidence(subset(igh_whsc1_df, Sample==this_srr)$DELLY_SV_EVIDENCE)
  if( is.na(sv_bp[1]) ){
    this_whsc1_breakpoint <- NA
    this_fgfr3_expression <- NA
    this_igh_breakpoint <- NA
    this_whsc1_expression <- NA
  } else{
    nrow_14 <- nrow(subset(sv_bp, chr==14))
    nrow_4 <- nrow(subset(sv_bp, chr==4))
    if(nrow_14==1 & nrow_4==1){
      this_whsc1_breakpoint <- as.numeric(as.character(subset(process_sv_evidence(subset(igh_whsc1_df, Sample==this_srr)$DELLY_SV_EVIDENCE), chr==4)$pos[1]))
      this_igh_breakpoint <- as.numeric(as.character(subset(process_sv_evidence(subset(igh_whsc1_df, Sample==this_srr)$DELLY_SV_EVIDENCE), chr==14)$pos[1]))
      this_fgfr3_expression <- subset(expr, srr==this_srr & gene=="FGFR3")$log10tpm
      this_whsc1_expression <- subset(expr, srr==this_srr & gene=="WHSC1")$log10tpm
    } else{
      this_whsc1_breakpoint <- NA
      this_fgfr3_expression <- NA
      this_igh_breakpoint <- NA
      this_whsc1_expression <- NA
    }
  }
  bp_igh <- c(bp_igh, this_igh_breakpoint)
  bp_whsc1 <- c(bp_whsc1, this_whsc1_breakpoint)
  ex_fgfr3 <- c(ex_fgfr3, this_fgfr3_expression)
  ex_whsc1 <- c(ex_whsc1, this_whsc1_expression)
}
plot_df <- data.frame(bp_igh=bp_igh, bp_whsc1=bp_whsc1, ex_fgfr3=ex_fgfr3, ex_whsc1=ex_whsc1)

ggplot(plot_df, aes(x=bp_igh, y=bp_whsc1, color=as.factor(round(ex_fgfr3)))) + geom_point()

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




############# Expression of oncogenes and TSGs
library(ggplot2)
library(scales)

plot_df <- NULL
#plot_df <- rbind(plot_df, data.frame(percentile=subset(primary_df, as.logical(geneA_kinase))$geneA_pct, gene_type="Kinase", pos="geneA"))
#plot_df <- rbind(plot_df, data.frame(percentile=subset(primary_df, as.logical(geneB_kinase))$geneB_pct, gene_type="Kinase", pos="geneB"))
plot_df <- rbind(plot_df, data.frame(percentile=subset(primary_df, as.logical(geneA_oncogene))$geneA_pct, gene_type="Oncogene", pos="geneA"))
plot_df <- rbind(plot_df, data.frame(percentile=subset(primary_df, as.logical(geneB_oncogene))$geneB_pct, gene_type="Oncogene", pos="geneB"))
plot_df <- rbind(plot_df, data.frame(percentile=subset(primary_df, as.logical(geneA_tsg))$geneA_pct, gene_type="Tumor Suppressor", pos="geneA"))
plot_df <- rbind(plot_df, data.frame(percentile=subset(primary_df, as.logical(geneB_tsg))$geneB_pct, gene_type="Tumor Suppressor", pos="geneB"))
#plot_df <- rbind(plot_df, data.frame(percentile=subset(primary_df, as.logical(geneA_driver))$geneA_pct,gene_type="Driver", pos="geneA"))
#plot_df <- rbind(plot_df, data.frame(percentile=subset(primary_df, as.logical(geneB_driver))$geneB_pct,gene_type="Driver", pos="geneB"))

ggplot(plot_df, aes(x=gene_type, y=percentile)) + geom_violin(size=2) + geom_jitter(data=plot_df, aes(color=pos), size=3, alpha=0.5) + theme_bw(base_size=20) + labs(x="Gene Category", y="Expression Percentile", title="Expression Percentile of Oncogenes and Tumor Suppressors", color="Gene\norientation")


plot_df <- NULL
oncogenes <- as.character(subset(primary_df, as.logical(geneA_oncogene))$geneA)
oncogenes <- c(oncogenes, as.character(subset(primary_df, as.logical(geneB_oncogene))$geneB))
oncogenes <- unique(oncogenes)
for(g in oncogenes){
  pct <- subset(primary_df, geneA==g)$geneA_pct
  pct <- unique(c(pct, subset(primary_df, geneB==g)$geneB_pct))
  plot_df <- rbind(plot_df, data.frame(gene_type="Oncogene", gene=g, n_samples=length(pct), median=median(pct)))
}

tsgs <- as.character(subset(primary_df, as.logical(geneA_tsg))$geneA)
tsgs <- c(tsgs, as.character(subset(primary_df, as.logical(geneB_tsg))$geneB))
tsgs <- unique(tsgs)
for(g in tsgs){
  pct <- subset(primary_df, geneA==g)$geneA_pct
  pct <- unique(c(pct, subset(primary_df, geneB==g)$geneB_pct))
  plot_df <- rbind(plot_df, data.frame(gene_type="Tumor Suppressor", gene=g, n_samples=length(pct), median=median(pct)))
}

gene_type_order <- c(order(subset(plot_df, gene_type=="Oncogene")$n_samples, decreasing=TRUE), order(subset(plot_df, gene_type=="Tumor Suppressor")$n_samples, decreasing=TRUE))
plot_df <- data.frame(plot_df, gene_type_order=gene_type_order)


ggplot(subset(plot_df, gene_type_order <= 10), aes(x=gene, y=gene_type, size=n_samples, color=median)) + geom_point(alpha=0.5) + scale_colour_gradient2(low=muted("blue"), high=muted("red"), midpoint = 0.5) + facet_grid( . ~ gene_type)


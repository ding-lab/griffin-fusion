#Analysis of fusions with outlier expression

get_expr_df <- function(df, this_gene){
  srr_with_fusion <- sort(unique(subset(df, geneA == this_gene | geneB == this_gene)$srr))
  this_gene_expr <- subset(expr, gene == this_gene)
  if(nrow(this_gene_expr) > 0){
    expr_with_fusion <- data.frame(subset(this_gene_expr, srr %in% srr_with_fusion), fusion_status="Fusion")
    expr_without_fusion <- data.frame(subset(this_gene_expr, !(srr %in% srr_with_fusion)), fusion_status="None reported")
    this_gene_expr <- rbind(expr_with_fusion, expr_without_fusion)
  
    gene_partners <- NULL
    cnv <- NULL
    for(this_srr in this_gene_expr$srr){
      if(this_srr %in% srr_with_fusion){
        gene_partners <- c(gene_partners, paste( unique(sort(subset(primary_df, srr == this_srr & (geneA == this_gene | geneB == this_gene) )$fusion)) , collapse='\n'))
      } else{
        gene_partners <- c(gene_partners, "")
      }
    
      tn_cnv_ratio = 2^(subset(this_gene_expr, srr == this_srr)$gene_avg_cnv) # gene_avg_cnv is average log2 CNV ratio
      if(is.na(tn_cnv_ratio)){
        this_cnv <- "No data"
      } else if(tn_cnv_ratio > 2){
        this_cnv <- "Major amplification"
      } else if(tn_cnv_ratio > 1.25){
        this_cnv <- "Minor amplification"
      } else if(tn_cnv_ratio > 0.75){
        this_cnv <- "Neutral"
      } else if(tn_cnv_ratio > 0.5){
        this_cnv <- "Minor deletion"
      } else{
        this_cnv <- "Major deletion"
      }  
      cnv <- c(cnv, this_cnv)
    }
  
    this_gene_expr <- data.frame(this_gene_expr, gene_partners = gene_partners, cnv = factor(cnv, levels=c("Major deletion","Minor deletion","Neutral","Minor amplification","Major amplification", "No data")))
  } else{
    this_gene_expr <- NA
  }
  return(this_gene_expr)
}

get_expr_df_multiple <- function(df, gene_list){
  these_genes_expr <- NULL
  
  for(this_gene in gene_list){
      
    srr_with_fusion <- sort(unique(subset(df, geneA == this_gene | geneB == this_gene)$srr))
    this_gene_expr <- subset(expr, gene == this_gene)
    if(nrow(this_gene_expr) > 0){
      expr_with_fusion <- data.frame(subset(this_gene_expr, srr %in% srr_with_fusion), fusion_status="Fusion")
      expr_without_fusion <- data.frame(subset(this_gene_expr, !(srr %in% srr_with_fusion)), fusion_status="None reported")
      this_gene_expr <- rbind(expr_with_fusion, expr_without_fusion)
      
      gene_partners <- NULL
      cnv <- NULL
      for(this_srr in this_gene_expr$srr){
        if(this_srr %in% srr_with_fusion){
          gene_partners <- c(gene_partners, paste( unique(sort(subset(primary_df, srr == this_srr & (geneA == this_gene | geneB == this_gene) )$fusion)) , collapse='\n'))
        } else{
          gene_partners <- c(gene_partners, "")
        }
        
        tn_cnv_ratio = 2^(subset(this_gene_expr, srr == this_srr)$gene_avg_cnv)
        if(is.na(tn_cnv_ratio)){
          this_cnv <- "No data"
        } else if(tn_cnv_ratio > 2){
          this_cnv <- "Major amplification"
        } else if(tn_cnv_ratio > 1.25){
          this_cnv <- "Minor amplification"
        } else if(tn_cnv_ratio > 0.75){
          this_cnv <- "Neutral"
        } else if(tn_cnv_ratio > 0.5){
          this_cnv <- "Minor deletion"
        } else{
          this_cnv <- "Major deletion"
        }  
        cnv <- c(cnv, this_cnv)
      }
      
      this_gene_expr <- data.frame(this_gene_expr, gene_partners = gene_partners, cnv = factor(cnv, levels=c("Major deletion","Minor deletion","Neutral","Minor amplification","Major amplification", "No data")))
    } else{
      this_gene_expr <- NA
    }
    
    if( !all(is.na(this_gene_expr)) ){
      these_genes_expr <- rbind(these_genes_expr, this_gene_expr)  
    }
    
  }
  
  if(is.null(this_gene_expr)){
    these_genes_expr <- NA
  }
  return(these_genes_expr)
}

significance <- function(df, this_gene){
  this_gene_expr <- get_expr_df(df, this_gene)
  if(all(is.na(this_gene_expr))){
    return(NA)
  }
  
  expr_with_fusion <- subset(this_gene_expr, fusion_status=="Fusion")
  expr_without_fusion <- subset(this_gene_expr, fusion_status!="Fusion")
  
  median_pct <- median(expr_with_fusion$pct)
  
  if( nrow(expr_with_fusion) > 1){
    t.test.pvalue.overexpression <- t.test(expr_with_fusion$log10tpm, expr_without_fusion$log10tpm, alternative="greater")$p.value #Fusion greater than No fusion
    t.test.pvalue.underexpression <- t.test(expr_with_fusion$log10tpm, expr_without_fusion$log10tpm, alternative="less")$p.value #Fusion less than No fusion
  } else{
    t.test.pvalue.overexpression <- NA
    t.test.pvalue.underexpression <- NA
  }
  
  if(nrow(table(this_gene_expr[,c("outlier_over_tpm", "fusion_status")])) > 1){
    fisher.pvalue.overexpression <- fisher.test(table(this_gene_expr[,c("outlier_over_tpm", "fusion_status")]), alternative="t")$p.value #No fusion is less than Fusion
  } else{
    fisher.pvalue.overexpression <- NA
  }
  if(nrow(table(this_gene_expr[,c("outlier_under_tpm", "fusion_status")])) > 1){
    fisher.pvalue.underexpression <- fisher.test(table(this_gene_expr[,c("outlier_under_tpm", "fusion_status")]), alternative="t")$p.value #No fusion is greater than Fusion
  } else{
    fisher.pvalue.underexpression <- NA
  }
  
  return(data.frame(gene=this_gene, n_samples_with_fusion=nrow(expr_with_fusion), n_samples=nrow(this_gene_expr), median_expr_pct=median_pct, ttest_over=t.test.pvalue.overexpression, ttest_under=t.test.pvalue.underexpression, fisher_over=fisher.pvalue.overexpression, fisher_under=fisher.pvalue.underexpression))
}

plotting <- function(df, this_gene, pdf_path, ymax_value, labels=TRUE, seed=10){
  set.seed(seed)
  
  this_gene_expr <- get_expr_df(df, this_gene)
  
  plot_df <- data.frame(this_gene_expr, fusion_status_j=jitter(as.numeric(this_gene_expr$fusion_status)))
  
  if( any(plot_df$fusion_status == "Fusion") ){
    #library(ggplot2)
    #library(ggrepel)
    
    p <- ggplot(plot_df, aes(x=fusion_status, y=log10tpm, color=cnv, fill=cnv))
    
    p <- p + geom_violin(aes(fill=NULL),color="black", draw_quantiles=c(0.5))
    p <- p + geom_point(aes(x=fusion_status_j), alpha=0.75, shape=16)
    if(labels){
      p <- p + geom_label_repel(data=subset(plot_df, fusion_status=="Fusion"), aes(x=fusion_status_j, y=log10tpm, label=gene_partners), label.size=NA, color=c('white','black')[as.numeric(subset(plot_df, fusion_status=="Fusion")$cnv=="No data")+1], box.padding=0.35, point.padding=0.5, segment.color="grey50")
    }
    p <- p + ylim(0, ymax_value)
    p <- p + theme_bw(base_size=20)
    p <- p + scale_color_brewer(palette="Set1", drop=FALSE)
    p <- p + scale_fill_brewer(palette="Set1", drop=FALSE)
    p <- p + labs(x="Fusion Status", y="Gene Expression TPM (log10)", color="CNV Status", title=paste0("Fusion Gene Expression (", this_gene, ")"))
    p <- p + guides(fill=FALSE, alpha=FALSE)
    
    if( !(dir.exists( dirname(pdf_path) )) ){
      dir.create(dirname(pdf_path), showWarnings = FALSE, recursive = TRUE)
    }
    pdf(pdf_path, 10, 10, useDingbats = FALSE)
    print(p)
    shh <- dev.off()
  }
}

plotting_multiple <- function(df, gene_list, pdf_path, ymax_value, labels=TRUE, seed=10){
  set.seed(seed)
  
  this_gene_expr <- get_expr_df_multiple(df, gene_list)
  
  plot_df <- data.frame(this_gene_expr, gene_x=droplevels(factor(this_gene_expr$gene, levels=gene_list)), fusion_status_j=jitter(as.numeric(droplevels(factor(this_gene_expr$gene, levels=gene_list)))))
  
  if( any(plot_df$fusion_status == "Fusion") ){
    #library(ggplot2)
    #library(ggrepel)
    
    p <- ggplot(plot_df, aes(x=gene_x, y=log10tpm, color=cnv, fill=cnv))
    
    p <- p + geom_violin(aes(fill=NULL),color="black", draw_quantiles=c(0.5))
    p <- p + geom_point(aes(x=fusion_status_j), alpha=0.75, shape=16)
    if(labels){
      p <- p + geom_label_repel(data=plot_df, aes(x=fusion_status_j, y=log10tpm, label=gene_partners), label.size=NA, color=c('white','black')[as.numeric(plot_df$cnv=="No data")+1], box.padding=0.35, point.padding=0.5, segment.color="grey50")
    }
    p <- p + ylim(0, ymax_value)
    p <- p + theme_bw(base_size=20)
    p <- p + scale_color_brewer(palette="Set1", drop=FALSE)
    p <- p + scale_fill_brewer(palette="Set1", drop=FALSE)
    p <- p + labs(x="Fusion Status", y="Gene Expression TPM (log10)", color="CNV Status", title="Fusion Gene Expression (Multiple Genes)")
    p <- p + guides(fill=FALSE, alpha=FALSE)
    
    if( !(dir.exists( dirname(pdf_path) )) ){
      dir.create(dirname(pdf_path), showWarnings = FALSE, recursive = TRUE)
    }
    pdf(pdf_path, 10+2*length(gene_list), 10, useDingbats = FALSE)
    print(p)
    shh <- dev.off()
  }
}

genes_with_fusions <- sort(unique(c(as.character(primary_df$geneA), as.character(primary_df$geneB))))
ymax_value <- max(primary_df$geneA_log10tpm, primary_df$geneB_log10tpm, na.rm=T)
n_genes <- length(genes_with_fusions)

significance_df <- NULL

if(FALSE){ #warning, takes several hours for MMY dataset
  count_up <- 0
  for(this_gene in genes_with_fusions){
    significant <- significance(primary_df, this_gene)
    if( this_gene %in% c("IGH@","IGK@","IGL@","IGH@pseudo","IGK@pseudo","IGL@pseudo") | all(is.na(significant))){
      next
    }
    
    count_up <- count_up + 1
    print(count_up/n_genes)
    significance_df <- rbind(significance_df, significant)
    
    #plot all genes
    plotting(primary_df, this_gene, paste0("analysis/fusion_outliers/all_plots/",this_gene,".pdf"), ymax_value)
    #plotting(primary_df, this_gene, paste0("analysis/fusion_outliers/all_plots_withoutlabel/",this_gene,".pdf"), ymax_value, label=FALSE)
    
  }
  
  write.table(significance_df[apply(significance_df, 1, function(x) !all(is.na(x))),], "analysis/fusion_outliers/significance.tsv", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
}

sig_df <- read.table("analysis/fusion_outliers/significance.tsv", header=T)
n_genes_recurrent_fusions <- nrow(subset(sig_df, n_samples_with_fusion > 1))
genes_recurrent_fusions <- subset(sig_df, n_samples_with_fusion > 1)$gene
pvalue_threshold <- 0.01

#if significant, plot in significant folder
for(this_gene in genes_recurrent_fusions){
  significant <- subset(sig_df, gene == this_gene)
  if(significant$n_samples_with_fusion > 1 & significant$median_expr_pct > 0.50 & ( (!is.na(significant$ttest_over) & significant$ttest_over < pvalue_threshold) | (!is.na(significant$fisher_over) & significant$fisher_over < pvalue_threshold) )){
    print(paste0("Over: ",this_gene))
    plotting(primary_df, this_gene, paste0("analysis/fusion_outliers/significant_plots_over/",this_gene,".pdf"), ymax_value)
    plotting(primary_df, this_gene, paste0("analysis/fusion_outliers/significant_plots_over_withoutlabel/",this_gene,".pdf"), ymax_value, label=FALSE)
  } else if(significant$n_samples_with_fusion > 1 & significant$median_expr_pct < 0.50 & ( (!is.na(significant$ttest_under) & significant$ttest_under < pvalue_threshold) | (!is.na(significant$fisher_under) & significant$fisher_under < pvalue_threshold) )){
    print(paste0("Under: ",this_gene))
    plotting(primary_df, this_gene, paste0("analysis/fusion_outliers/significant_plots_under/",this_gene,".pdf"), ymax_value)
    plotting(primary_df, this_gene, paste0("analysis/fusion_outliers/significant_plots_under_withoutlabel/",this_gene,".pdf"), ymax_value, label=FALSE)
  }
}


#Analysis of fusions with outlier expression

get_expr_df <- function(df, this_gene){
  srr_with_fusion <- sort(unique(subset(df, geneA == gene | geneB == gene)$srr))
  this_gene_expr <- subset(expr, gene == this_gene)
  expr_with_fusion <- data.frame(subset(this_gene_expr, srr %in% srr_with_fusion), fusion_status="Fusion")
  expr_without_fusion <- data.frame(subset(this_gene_expr, !(srr %in% srr_with_fusion)), fusion_status="None reported")
  this_gene_expr <- rbind(expr_with_fusion, expr_without_fusion)
  
  gene_partners <- NULL
  cnv <- NULL
  for(this_srr in this_gene_expr$srr){
    if(this_srr %in% srr_with_fusion){
      gene_partners <- c(gene_partners, paste( unique(sort(subset(primary_df, srr == this_srr & (geneA == this_gene | geneB == this_gene) )$fusion)) , collapse='\n'))
    } else{
      gene_partners <- c(gene_partners, NA)
    }
    
    tn_cnv_ratio = exp(subset(this_gene_expr, srr == this_srr)$gene_avg_cnv)
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
  
  return(this_gene_expr)
}

significance <- function(df, this_gene){
  this_gene_expr <- get_expr_df(df, this_gene)
  
  median_pct <- median(expr_with_fusion$pct)
  
  t.test.pvalue.overexpression <- t.test(expr_with_fusion$log10tpm, expr_without_fusion$log10tpm, alternative="greater")$p.value #Fusion greater than No fusion
  t.test.pvalue.underexpression <- t.test(expr_with_fusion$log10tpm, expr_without_fusion$log10tpm, alternative="less")$p.value #Fusion less than No fusion
  
  fisher.pvalue.overexpression <- fisher.test(table(this_gene_expr[,c("outlier_over_tpm", "fusion_status")]), alternative="less")$p.value #No fusion is less than Fusion
  fisher.pvalue.underexpression <- fisher.test(table(this_gene_expr[,c("outlier_over_tpm", "fusion_status")]), alternative="greater")$p.value #No fusion is greater than Fusion
  
  return(c(this_gene, length(srr_with_fusion), median_pct, t.test.pvalue.overexpression, t.test.pvalue.underexpression, fisher.pvalue.overexpression, fisher.pvalue.underexpression))
}

plotting <- function(df, this_gene, pdf_path){
  this_gene_expr <- get_expr_df(df, this_gene)
  
  
}








genes_with_fusions <- sort(unique(c(as.character(primary_df$geneA), as.character(primary_df$geneB))))

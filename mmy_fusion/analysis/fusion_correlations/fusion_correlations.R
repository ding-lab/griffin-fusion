# ==============================================================================
# Measure associations between fusion status and expression and clinical info
# Steven Foltz (smfoltz@wustl.edu), August 2018
# ==============================================================================

# ==============================================================================
# Functions
# ==============================================================================

# Return tbl of samples with a particular fusion, either SRRs or MMRFs
get_srrs_with_gene <- function(this_gene, fusions_tbl){
  fusions_tbl %>% filter(geneA == this_gene | geneB == this_gene) %>% 
    select(srr) %>% unique()
}
get_mmrfs_with_fusion <- function(this_fusion, fusions_tbl){
  fusions_tbl %>% filter(fusion == this_fusion ) %>% select(mmrf) %>% unique()
}

# Test clinical correlations of gene involved in a gene fusion
# Clinical testing includes seq-FISH plus other clinical measures
test_fusion_seqfish <- function(this_fusion, fusions_tbl, clinical_tbl, seqfish){
  samples_with_fusion <- get_mmrfs_with_fusion(this_fusion, fusions_tbl)
  fisher_table_counts <- clinical_tbl %>%
    mutate(has_fusion = mmrf %in% samples_with_fusion$mmrf) %>%
    mutate_at(seqfish, factor) %>%
    mutate_at("has_fusion", factor) %>%
    filter(!is.na(seqfish_Study_Visit_ID)) %>%
    group_by(eval(parse(text = seqfish)), has_fusion) %>% 
    summarize(count = n()) %>%
    complete(has_fusion, fill = list(count = 0)) %>% pull(count)
  fisher_test_table <- fisher_table_counts %>% matrix(nrow = 2)
  fisher_test_result <- fisher.test(fisher_test_table)
  return(fisher_test_result)
}

# Test gene expression for both genes in a fusion pair
test_fusion_expression <- function(this_fusion, fusions_tbl, expression_tbl){
  geneA = str_split(this_fusion, pattern = "--", simplify = TRUE)[1]
  geneB = str_split(this_fusion, pattern = "--", simplify = TRUE)[2]
  return(c(
    test_gene_expression(geneA, fusions_tbl, expression_tbl),
    test_gene_expression(geneA, fusions_tbl, expression_tbl)
  ))
}

# Returns median expression percentile, t-test, and Fisher's Exact Test p-values
test_gene_expression <- function(this_gene, fusions_tbl, expression_tbl){
  # Given a gene, find out if expression differs between fusion / not fusion
  samples_with_fusion <- get_srrs_with_gene(this_gene, fusions_tbl)
  expression_with_fusion <- expression_tbl %>% semi_join(samples_with_fusion, by = "srr") %>%
    filter(gene == this_gene)
  expression_without_fusion <- expression_tbl %>% anti_join(samples_with_fusion, by = "srr") %>%
    filter(gene == this_gene)
  
  # median expression percentile of fusion samples
  median_pct <- as.double(expression_with_fusion %>% summarize(median(pct)))
  
  # t-test
  if ( median_pct >= 0.5 ) {
    t.test.pvalue <- t.test(expression_with_fusion$log10tpm, expression_without_fusion$log10tpm, alternative = "greater")$p.value
    over_table <- expression_primary %>% 
      group_by(has_fusion = srr %in% samples_with_fusion$srr) %>%
      filter(gene == this_gene) %>% select(has_fusion, outlier_over_tpm) %>%
      summarize(outlier = sum(outlier_over_tpm == 1), not_outlier = sum(outlier_over_tpm == 0)) %>%
      select(outlier, not_outlier)
    fisher.pvalue <- fisher.test(over_table, alternative = "t")$p.value
  } else{
    t.test.pvalue <- t.test(expression_with_fusion$log10tpm, expression_without_fusion$log10tpm, alternative = "less")$p.value
    under_table <- expression_primary %>% 
      group_by(has_fusion = srr %in% samples_with_fusion$srr) %>%
      filter(gene == this_gene) %>% select(has_fusion, outlier_under_tpm) %>%
      summarize(outlier = sum(outlier_under_tpm == 1), not_outlier = sum(outlier_under_tpm == 0)) %>%
      select(outlier, not_outlier)
    fisher.pvalue <- fisher.test(under_table, alternative = "t")$p.value
    
  }
  return( c(median_pct, t.test.pvalue, fisher.pvalue))
}

# ==============================================================================
# Prepare lists of fusions and genes for testing (avoid testing rare events)
# ==============================================================================

# Fusions pairs seen in at least 3 samples
fusion_pairs_gt2 <- fusions_primary %>% group_by(fusion) %>% 
  summarize(count = n()) %>% filter(count >= 3) %>% arrange(desc(count))

# Fusion genes seen in at least 3 samples
fusion_genes_gt2 <- fusions_primary %>% gather(geneA, geneB, key = "geneAB", value = "fusion_gene") %>%
  group_by(fusion_gene) %>% filter( !(fusion_gene %in% c("IGH@", "IGK@", "IGL@")) ) %>%
  summarize(count = n()) %>% filter(count >= 3) %>% arrange(desc(count))

# ==============================================================================
# ==============================================================================

for (this_gene in fusion_genes_gt2$fusion_gene) {
  print(this_gene)
  print(c(this_gene, test_gene_expression(this_gene, fusions_primary, expression_primary)))
}

for (this_fusion in fusion_pairs_gt2$fusion) {
  fusions_primary %>% filter( fusion == this_fusion)
}

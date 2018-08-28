# ==============================================================================
# Measure associations between fusion status and expression and clinical info
# Steven Foltz (smfoltz@wustl.edu), August 2018
# ==============================================================================

# ==============================================================================
# General Functions
# ==============================================================================

# Extract sample IDs with a certain fusion
get_ids_with_fusion <- function(this_fusion, fusions_tbl){
  fusions_tbl %>% filter(fusion == this_fusion ) %>% select(mmrf, srr) %>% unique()
}
get_ids_without_fusion <- function(this_fusion, fusions_tbl){
  ids_with_fusion <- get_ids_with_fusion(this_fusion, fusions_tbl)
  fusions_tbl %>% anti_join(ids_with_fusion, by = "srr") %>%
    select(mmrf, srr) %>% unique()
}

# Extract sample IDs with a certain gene involved in a fusion
get_ids_with_gene <- function(this_gene, fusions_tbl){
  fusions_tbl %>% filter(geneA == this_gene | geneB == this_gene) %>% 
    select(mmrf, srr) %>% unique()
}
get_ids_without_gene <- function(this_gene, fusions_tbl){
  ids_with_gene <- get_ids_with_gene(this_gene, fusions_tbl)
  fusions_tbl %>% anti_join(ids_with_gene, by = "srr") %>%
    select(mmrf, srr) %>% unique()
}

# Extract sample IDs with a certain seqFISH feature
get_ids_with_seqfish <- function(this_seqfish, seqfish_tbl, samples_tbl){
  seqfish_clinical_info %>% filter( !is.na(seqfish_Study_Visit_ID) ) %>% 
    filter( eval(parse(text = this_seqfish)) == 1 ) %>%
    left_join(samples_tbl, by = "mmrf") %>% select(mmrf, srr) %>% unique()
}
get_ids_without_seqfish <- function(this_seqfish, seqfish_tbl, samples_tbl){
  seqfish_clinical_info %>% filter( !is.na(seqfish_Study_Visit_ID) ) %>% 
    filter( eval(parse(text = this_seqfish)) != 1 ) %>%
    left_join(samples_tbl, by = "mmrf") %>% select(mmrf, srr) %>% unique()
}

# ==============================================================================
# Clinical Functions
# ==============================================================================

# Test event status against binary clinical variable
test_fusion_clinical_binary <- function(samples_with, samples_without,
                                        clinical_tbl, clinical_feature){
  # Fisher Exact test
  
}
  
# Test event status against categorical clinical variable
test_fusion_clinical_categorical <- function(samples_with, samples_without,
                                             clinical_tbl, clinical_feature){
  # Chi-square test
}
  
# Test event status against continuous clinical variable
test_fusion_clinical_continuous <- function(samples_with, samples_without,
                                            clinical_tbl, clinical_feature){
  # T test or Mann-Whitney U test
}

# Fisher exact test of seqFISH translocation and fusion status
# TODO need to re-write in context of samples_with, samples_without, etc.
# Make fisher table creation explicit and not rely on assumed ordering of squares
# Don't delete until other testing functions are written!
test_fusion_seqfish <- function(this_fusion, fusions_tbl, clinical_tbl, seqfish){
  # This test is problematic because you get a significant result when
  # you test translocations and fusions that are unrelated (anticorrelated).
  # Use only if there is a direct relationship between translocation and fusion
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

test_seqfish_correlations <- function(seqfish_1, seqfish_2, clinical_tbl)
test_fusion_correlations <- function(fusion_1, fusion_2, fusions_tbl)

# ==============================================================================
# Expression Functions
# ==============================================================================

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

test_seqfish_expression <- function(this_seqfish, clinical_tbl, expression_tbl){}

# ==============================================================================
# Prepare lists of fusions and genes for testing (avoid testing rare events)
# ==============================================================================

# Fusions pairs seen in at least 3 samples
fusion_pairs_gt2 <- fusions_primary %>% group_by(fusion) %>% 
  summarize(count = n()) %>% filter(count >= 3) %>% arrange(desc(count))

# Fusion genes seen in at least 3 samples
fusion_genes_gt2 <- fusions_primary %>% gather(geneA, geneB, key = "geneAB", value = "fusion_gene") %>%
  group_by(fusion_gene) %>% summarize(count = n()) %>% filter(count >= 3) %>% 
  arrange(desc(count))

# ==============================================================================
# ==============================================================================

for (this_gene in fusion_genes_gt2$fusion_gene) {
  print(this_gene)
  print(c(this_gene, test_gene_expression(this_gene, fusions_primary, expression_primary)))
}

for (this_fusion in fusion_pairs_gt2$fusion) {
  fusions_primary %>% filter( fusion == this_fusion)
}

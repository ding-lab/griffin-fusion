# ==============================================================================
# Measure associations between fusion status and expression and clinical info
# Steven Foltz (smfoltz@wustl.edu), August 2018
# ==============================================================================

# ==============================================================================
# General Functions
# ==============================================================================

# Extract sample IDs with a certain fusion
get_ids_with_fusion <- function(this_fusion, fusions_tbl){
  fusions_tbl %>% filter(fusion == this_fusion) %>% 
    select(mmrf, srr) %>% unique()
}
get_ids_without_fusion <- function(this_fusion, fusions_tbl, samples_tbl){
  ids_with_fusion <- get_ids_with_fusion(this_fusion, fusions_tbl)
  samples_tbl %>% anti_join(ids_with_fusion, by = "srr") %>% unique()
}

# Extract sample IDs with a certain gene involved in a fusion
get_ids_with_gene <- function(this_gene, fusions_tbl){
  fusions_tbl %>% filter(geneA == this_gene | geneB == this_gene) %>% 
    select(mmrf, srr) %>% unique()
}
get_ids_without_gene <- function(this_gene, fusions_tbl, samples_tbl){
  ids_with_gene <- get_ids_with_gene(this_gene, fusions_tbl)
  samples_tbl %>% anti_join(ids_with_gene, by = "srr") %>% unique()
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

# Test event status against categorical clinical variables, including
# seqFISH results, binary, and categorical variables with Fisher Exact test
test_event_clinical_discrete <- function(samples_with, samples_without,
                                         clinical_tbl, clinical_feature,
                                         fisher_test = FALSE,
                                         return_tibble = FALSE,
                                         return_table = FALSE){
  # Fisher Exact test
  if (fisher_test) {
    fisher_test_tibble <- clinical_tbl %>%
    mutate(has_event = mmrf %in% samples_with$mmrf) %>%
    mutate_at(clinical_feature, factor) %>%
    mutate_at("has_event", factor) %>%
    filter( !is.na(eval(parse(text = clinical_feature))) ) %>%
    group_by( clinical_value = eval(parse(text = clinical_feature)), has_event) %>% 
    summarize(count = n()) %>%
    complete(has_event, fill = list(count = 0))
  fisher_test_table <- fisher_test_tibble %>% pull(count) %>%
    matrix(ncol = 2, byrow = TRUE)
  fisher_test_result <- fisher.test(fisher_test_table)
  n_na <- clinical_tbl %>% 
    filter( is.na(eval(parse(text = clinical_feature))) ) %>% nrow()
  } else {
    stop("Value of fisher_test is not TRUE. Although Fisher's Exact Test is
         the only option, fisher_test must be made explicitly TRUE.")
  }
  # Return values
  if (return_tibble & return_table) {
    stop("Both return_tibble and return_table are TRUE. Only one may be TRUE.")
  } else if (return_tibble) {
    return(fisher_test_tibble)
  } else if (return_table) {
    return(fisher_test_table)
  } else {
    return(
      data.frame( p.value = fisher_test_result$p.value, 
                  n_missing = n_na)
      )
  }
  
}
  
# Test event status against continuous clinical variable
test_event_clinical_continuous <- function(samples_with, samples_without,
                                           clinical_tbl, clinical_feature,
                                           t_test = FALSE,
                                           mwu_test = FALSE,
                                           return_tibble = FALSE){
  # T test or Mann-Whitney U test
  with_tbl <- clinical_tbl %>% 
    semi_join(samples_with, by = "mmrf") %>% 
    filter( !is.na(eval(parse(text = clinical_feature))) ) %>%
    select(mmrf, clinical_feature) %>%
    mutate(has_event = 1)
  with_values <- with_tbl %>% pull(clinical_feature)
  
  without_tbl <- clinical_tbl %>% 
    anti_join(samples_with, by = "mmrf") %>% 
    filter( !is.na(eval(parse(text = clinical_feature))) ) %>%
    select(mmrf, clinical_feature) %>%
    mutate(has_event = 0)
  without_values <- without_tbl %>% pull(clinical_feature)
  
  combined_tbl <- bind_rows(with_tbl, without_tbl)
  
  n_na <- clinical_tbl %>% 
    filter( is.na(eval(parse(text = clinical_feature))) ) %>% nrow()
  
  if (t_test & mwu_test) {
    stop("Both t_test and mwu_test are TRUE. Only one may be TRUE.")
  } else if (t_test) {
    test_result <- t.test(with_values, without_values)
  } else if (mwu_test) {
    test_result <- wilcox.test(with_values, without_values)
  } else {
    stop("Neither t_test not mwu_test is TRUE. Either t_test xor mwu_test must be TRUE.")
  }
  
  if (return_tibble) {
    return(combined_tbl)
  } else {
    return(
      data.frame( p.value = test_result$p.value, 
                  n_missing = n_na)
    )
  }
}

# ==============================================================================
# Event correlations
# ==============================================================================
test_gene_correlations <- function(gene_1, gene_2, fusions_tbl, samples_tbl, 
                                   return_tibble = FALSE){
  samples_with_gene_1 <- fusions_tbl %>%
    filter(geneA == gene_1 | geneB == gene_1) %>% select(mmrf, srr) %>%
    unique()
  samples_with_gene_2 <- fusions_tbl %>%
    filter(geneA == gene_2 | geneB == gene_2) %>% select(mmrf, srr) %>%
    unique()
  
  has_gene_1 <- as.numeric(samples_tbl$srr %in% samples_with_gene_1$srr)
  has_gene_2 <- as.numeric(samples_tbl$srr %in% samples_with_gene_2$srr)
  
  if ( any(has_gene_1 + has_gene_2  == 2) ) {
    gene_correlation <- cor.test(has_gene_1, has_gene_2)
  } else{
    gene_correlation <- NA
  }
  
  gene_tbl <- samples_tbl %>% 
    mutate( gene_1 = has_gene_1, gene_2 = has_gene_2)
  
  if (return_tibble) {
    return(gene_tbl)
  } else {
    return(gene_correlation)
  }
  
}

test_fusion_correlations <- function(fusion_1, fusion_2, fusions_tbl, 
                                     samples_tbl, return_tibble = FALSE){
  samples_with_fusion_1 <- fusions_tbl %>%
    filter(fusion == fusion_1) %>% select(mmrf, srr) %>% unique()
  samples_with_fusion_2 <- fusions_tbl %>%
    filter(fusion == fusion_2) %>% select(mmrf, srr) %>% unique()
  
  has_fusion_1 <- as.numeric(samples_tbl$srr %in% samples_with_fusion_1$srr)
  has_fusion_2 <- as.numeric(samples_tbl$srr %in% samples_with_fusion_2$srr)
  
  if ( any(has_fusion_1 + has_fusion_2 == 2) ) {
    fusion_correlation <- cor.test(has_fusion_1, has_fusion_2)
  } else {
    fusion_correlation <- NA
  }
  
  fusion_tbl <- samples_tbl %>% 
    mutate( fusion_1 = has_fusion_1, fusion_2 = has_fusion_2)
  
  if (return_tibble) {
    return(fusion_tbl)
  } else {
    return(fusion_correlation)
  }
  
}

test_seqfish_correlations <- function(seqfish_1, seqfish_2, clinical_tbl,
                                      return_tibble = FALSE) {
  seqfish_tbl <- clinical_tbl %>% filter( !is.na(seqfish_Study_Visit_ID) ) %>%
    select(mmrf, seqfish_1, seqfish_2)
  seqfish_1_vector <- seqfish_tbl %>% pull(seqfish_1)
  seqfish_2_vector <- seqfish_tbl %>% pull(seqfish_2)
  seqfish_correlation <- cor.test(seqfish_1_vector, seqfish_2_vector)

  if (return_tibble) {
    return(seqfish_tbl)
  } else {
    return(seqfish_correlation)
  }
}

# ==============================================================================
# Expression Functions
# ==============================================================================

# Test gene expression for both genes in a fusion pair
test_fusion_expression <- function(samples_with, samples_without, 
                                   fusions_tbl, expression_tbl){
  geneA = str_split(this_fusion, pattern = "--", simplify = TRUE)[1]
  geneB = str_split(this_fusion, pattern = "--", simplify = TRUE)[2]
  return(c(
    test_event_expression(samples_with, samples_without, geneA, 
                          fusions_tbl, expression_tbl),
    test_event_expression(samples_with, samples_without, geneB, 
                          fusions_tbl, expression_tbl)
  ))
}

# Returns median expression percentile, t-test, and Fisher's Exact Test p-values
test_event_expression <- function(samples_with, samples_without, this_gene,
                                  fusions_tbl, expression_tbl){
  # Given a gene, find out if expression differs between fusion / not fusion
  expression_with_fusion <- expression_tbl  %>%
    filter(gene == this_gene) %>% semi_join(samples_with, by = "srr")
  expression_without_fusion <- expression_tbl %>% 
    filter(gene == this_gene) %>% anti_join(samples_with, by = "srr")
  
  # median expression percentile of fusion samples
  median_pct <- as.double(expression_with_fusion %>% summarize(median(pct)))
  
  # t-test
  if ( median_pct >= 0.5 ) {
    t.test.pvalue <- t.test(expression_with_fusion$log10tpm, 
                            expression_without_fusion$log10tpm, 
                            alternative = "greater")$p.value
    over_table <- expression_primary %>% 
      group_by(has_fusion = srr %in% samples_with$srr) %>%
      filter(gene == this_gene) %>% select(has_fusion, outlier_over_tpm) %>%
      summarize(outlier = sum(outlier_over_tpm == 1), 
                not_outlier = sum(outlier_over_tpm == 0)) %>%
      select(outlier, not_outlier)
    fisher.pvalue <- fisher.test(over_table, alternative = "t")$p.value
  } else{
    t.test.pvalue <- t.test(expression_with_fusion$log10tpm, 
                            expression_without_fusion$log10tpm, 
                            alternative = "less")$p.value
    under_table <- expression_primary %>% 
      group_by(has_fusion = srr %in% samples_with$srr) %>%
      filter(gene == this_gene) %>% select(has_fusion, outlier_under_tpm) %>%
      summarize(outlier = sum(outlier_under_tpm == 1), 
                not_outlier = sum(outlier_under_tpm == 0)) %>%
      select(outlier, not_outlier)
    fisher.pvalue <- fisher.test(under_table, alternative = "t")$p.value
    
  }
  return( c(median_pct, t.test.pvalue, fisher.pvalue))
}

# ==============================================================================
# Assign seqFISH and clinical variables to approapiate lists
# ==============================================================================

seqfish_variable_names <- c("seqfish_CN_del_13q14", 
                            "seqfish_CN_del_13q34", 
                            "seqfish_CN_del_17p13", 
                            "seqfish_CN_gain_1q21", 
                            "seqfish_Hyperdiploidy", 
                            "seqfish_Translocation_WHSC1_4_14", 
                            "seqfish_Translocation_CCND3_6_14", 
                            "seqfish_Translocation_MYC_8_14", 
                            "seqfish_Translocation_MAFA_8_14", 
                            "seqfish_Translocation_CCND1_11_14", 
                            "seqfish_Translocation_CCND2_12_14", 
                            "seqfish_Translocation_MAF_14_16", 
                            "seqfish_Translocation_MAFB_14_20")

discrete_clinical_variable_names <- c("age_ge_66", 
                                      "Female", 
                                      "Race_White", 
                                      "Race_Black", 
                                      "Race_Other", 
                                      "race", 
                                      "ECOG", 
                                      "ISS_Stage",
                                      "Bone_lesions", 
                                      "Plamacytoma")

continuous_clinical_variable_names <- c("Age", 
                                        "BM_Plasma_Cell_Percent", 
                                        "LDH")

# ==============================================================================
# Prepare lists of fusions and genes for testing (avoid testing rare events)
# ==============================================================================

# Fusions pairs seen in at least 3 samples
fusion_pairs_gt2 <- fusions_primary %>% group_by(fusion) %>% 
  summarize(count = n()) %>% filter(count >= 3) %>% arrange(desc(count))

# Fusion genes seen in at least 3 samples
fusion_genes_gt2 <- fusions_primary %>% 
  gather(geneA, geneB, key = "geneAB", value = "fusion_gene") %>%
  select(mmrf, srr, fusion_gene) %>% distinct() %>% group_by(fusion_gene) %>% 
  summarize(count = n()) %>% filter(count >= 3) %>% arrange(desc(count))

# ==============================================================================
# Business
# ==============================================================================

for (this_gene in fusion_genes_gt2$fusion_gene) {
  print(this_gene)
  if ( str_detect(this_gene, "@") ) {
    next
  } else {
    print(c(this_gene, test_event_expression(get_ids_with_gene(this_gene, fusions_primary), get_ids_without_gene(this_gene, fusions_primary, samples_primary), this_gene, fusions_primary, expression_primary)))
  }
}

for (this_fusion in fusion_pairs_gt2$fusion) {
  fusions_primary %>% filter( fusion == this_fusion)
}

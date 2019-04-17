# ==============================================================================
# Expression (MMRF Fusions)
# Steven Foltz (smfoltz@wustl.edu), April 2019
# ==============================================================================

paper_main = "paper/main/02_expression/"
paper_supp = "paper/supplemental/02_expression/"

# Create directories 
dir.create(paper_main, recursive = TRUE, showWarnings = FALSE)
dir.create(paper_supp, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# Measure associations between fusion status and expression and clinical info
# Steven Foltz (smfoltz@wustl.edu), August 2018
# ==============================================================================

if (FALSE) {
  
  # CHANGE TO TRUE ONLY IF YOU NEED TO CREATE TESTING FILE FOR THE FIRST TIME
  recreate_testing_tbl <- FALSE
  
  # ============================================================================
  # General Functions
  # ============================================================================

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
  
  # ============================================================================
  # Clinical Functions
  # ============================================================================
  
  # Test event status against categorical clinical variables, including
  # seqFISH results, binary, and categorical variables with Fisher Exact test
  test_event_clinical_discrete <- function(samples_with, samples_without,
                                           event, event_type,
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
        group_by( clinical_value = eval(parse(text = clinical_feature)), 
                  has_event) %>%
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
    event1 <- event
    event2 <- clinical_feature
    event_type <- event_type
    test_performed <- "Fisher's Exact Test"
    n_samples_with <- samples_with %>% nrow()
    n_samples_without <- samples_without %>% nrow()
    n_samples_with_tested <- sum(fisher_test_table[,2])
    n_samples_without_tested <- sum(fisher_test_table[,1])
    n_samples_na_tested <- n_na
    if ( dim(fisher_test_table)[1] == 2) {
      test_statistic <- fisher_test_result$estimate
    } else{
      test_statistic <- NA
    }
    p.value <- fisher_test_result$p.value
    median_value <- NA
    return_test_stats <- tribble(~event1, ~event2, ~event_type, ~test_performed,
                                 ~n_samples_with, ~n_samples_without,
                                 ~n_samples_with_tested, 
                                 ~n_samples_without_tested,
                                 ~n_samples_na_tested, ~median_value, 
                                 ~test_statistic, ~p.value,
                                 event1,  event2,  event_type,  test_performed,
                                 n_samples_with,  
                                 n_samples_without,
                                 n_samples_with_tested,  n_samples_without_tested,
                                 n_samples_na_tested, median_value, 
                                 test_statistic,  p.value)
    
    if (return_tibble & return_table) {
      stop("Both return_tibble and return_table are TRUE. Only one may be TRUE.")
    } else if (return_tibble) {
      return(fisher_test_tibble)
    } else if (return_table) {
      return(fisher_test_table)
    } else {
      return(return_test_stats)
    }
    }
  
  # Test event status against continuous clinical variable
  test_event_clinical_continuous <- function(samples_with, samples_without,
                                             event, event_type,
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
    
    n_na <- clinical_tbl %>% 
      filter( is.na(eval(parse(text = clinical_feature))) ) %>% nrow()
    
    combined_tbl <- bind_rows(with_tbl, without_tbl)
    
    if (return_tibble) {
      return(combined_tbl)
    }
    
    if (t_test & mwu_test) {
      stop("Both t_test and mwu_test are TRUE. Only one may be TRUE.")
    }
    
    if ( length(with_values) < 3 | length(without_values) < 3) {
      test_result <- NA
      test_performed <- NA
      test_statistic <- NA
      p.value <- NA
    } else if (t_test) {
      test_result <- t.test(with_values, without_values)
      test_performed <- "Student's t-Test"
      test_statistic <- test_result$statistic
      p.value <- test_result$p.value
    } else if (mwu_test) {
      test_result <- wilcox.test(with_values, without_values)
      test_performed <- "Mann-Whitney U Test"
      test_statistic <- test_result$statistic
      p.value <- test_result$p.value
    } else {
      stop("Neither t_test not mwu_test is TRUE. 
           Either t_test xor mwu_test must be TRUE.")
    }    
    
    # Return values
    event1 <- event
    event2 <- clinical_feature
    event_type <- event_type
    n_samples_with <- samples_with %>% nrow()
    n_samples_without <- samples_without %>% nrow()
    n_samples_with_tested <- with_values %>% length()
    n_samples_without_tested <- without_values %>% length()
    n_samples_na_tested <- n_na
    median_value <- NA
    
    return_test_stats <- tribble(~event1, ~event2, ~event_type, ~test_performed,
                                 ~n_samples_with, ~n_samples_without,
                                 ~n_samples_with_tested, 
                                 ~n_samples_without_tested,
                                 ~n_samples_na_tested, ~median_value, 
                                 ~test_statistic, ~p.value,
                                 event1,  event2,  event_type,  test_performed,
                                 n_samples_with,  n_samples_without,
                                 n_samples_with_tested,
                                 n_samples_without_tested,
                                 n_samples_na_tested, median_value, 
                                 test_statistic,  p.value)
    
    return(return_test_stats)
    }
  
  # ============================================================================
  # Event correlations
  # ============================================================================
  test_gene_correlations <- function(gene_1, gene_2, event_type,
                                     fusions_tbl, samples_tbl, 
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
      test_statistic <- gene_correlation$statistic
      p.value <- gene_correlation$p.value
    } else{
      gene_correlation <- NA
      test_statistic <- NA
      p.value <- NA
    }
    
    gene_tbl <- samples_tbl %>% 
      mutate( gene_1 = has_gene_1, gene_2 = has_gene_2)
    
    # Return values
    event1 <- gene_1
    event2 <- gene_2
    event_type <- event_type
    test_performed <- "Paired Sample Correlation"
    n_samples_with <- NA
    n_samples_without <- NA
    n_samples_with_tested <- NA
    n_samples_without_tested <- NA
    n_samples_na_tested <- NA
    median_value <- NA
    return_test_stats <- tribble(~event1, ~event2, ~event_type, ~test_performed,
                                 ~n_samples_with, ~n_samples_without,
                                 ~n_samples_with_tested, 
                                 ~n_samples_without_tested,
                                 ~n_samples_na_tested, ~median_value, 
                                 ~test_statistic, ~p.value,
                                 event1,  event2,  event_type,  test_performed,
                                 n_samples_with,  
                                 n_samples_without,
                                 n_samples_with_tested,  n_samples_without_tested,
                                 n_samples_na_tested, median_value, 
                                 test_statistic,  p.value)
    
    if (return_tibble) {
      return(gene_tbl)
    } else {
      return(return_test_stats)
    }
    
  }
  
  test_fusion_correlations <- function(fusion_1, fusion_2, event_type,
                                       fusions_tbl, samples_tbl, 
                                       return_tibble = FALSE){
    samples_with_fusion_1 <- fusions_tbl %>%
      filter(fusion == fusion_1) %>% select(mmrf, srr) %>% unique()
    samples_with_fusion_2 <- fusions_tbl %>%
      filter(fusion == fusion_2) %>% select(mmrf, srr) %>% unique()
    
    has_fusion_1 <- as.numeric(samples_tbl$srr %in% samples_with_fusion_1$srr)
    has_fusion_2 <- as.numeric(samples_tbl$srr %in% samples_with_fusion_2$srr)
    
    if ( any(has_fusion_1 + has_fusion_2 == 2) ) {
      fusion_correlation <- cor.test(has_fusion_1, has_fusion_2)
      test_statistic <- fusion_correlation$statistic
      p.value <- fusion_correlation$p.value
    } else {
      fusion_correlation <- NA
      test_statistic <- NA
      p.value <- NA
    }
    
    fusion_tbl <- samples_tbl %>% 
      mutate( fusion_1 = has_fusion_1, fusion_2 = has_fusion_2)
    
    # Return values
    event1 <- fusion_1
    event2 <- fusion_2
    event_type <- event_type
    test_performed <- "Paired Sample Correlation"
    n_samples_with <- NA
    n_samples_without <- NA
    n_samples_with_tested <- NA
    n_samples_without_tested <- NA
    n_samples_na_tested <- NA
    median_value <- NA
    return_test_stats <- tribble(~event1, ~event2, ~event_type, ~test_performed,
                                 ~n_samples_with, ~n_samples_without,
                                 ~n_samples_with_tested, 
                                 ~n_samples_without_tested,
                                 ~n_samples_na_tested, ~median_value, 
                                 ~test_statistic, ~p.value,
                                 event1,  event2,  event_type,  test_performed,
                                 n_samples_with,  
                                 n_samples_without,
                                 n_samples_with_tested,  n_samples_without_tested,
                                 n_samples_na_tested, median_value, 
                                 test_statistic,  p.value)
    
    if (return_tibble) {
      return(fusion_tbl)
    } else {
      return(return_test_stats)
    }
    
  }
  
  test_seqfish_correlations <- function(seqfish_1, seqfish_2, event_type,
                                        clinical_tbl, return_tibble = FALSE){
    seqfish_tbl <- clinical_tbl %>% filter( !is.na(seqfish_Study_Visit_ID) ) %>%
      select(mmrf, seqfish_1, seqfish_2)
    seqfish_1_vector <- seqfish_tbl %>% pull(seqfish_1)
    seqfish_2_vector <- seqfish_tbl %>% pull(seqfish_2)
    
    if ( any(seqfish_1_vector + seqfish_2_vector == 2) ) {
      seqfish_correlation <- cor.test(seqfish_1_vector, seqfish_2_vector)
      test_statistic <- seqfish_correlation$statistic
      p.value <- seqfish_correlation$p.value
    } else {
      seqfish_correlation <- NA
      test_statistic <- NA
      p.value <- NA
    }
    
    # Return values
    event1 <- seqfish_1
    event2 <- seqfish_2
    event_type <- event_type
    test_performed <- "Paired Sample Correlation"
    n_samples_with <- NA
    n_samples_without <- NA
    n_samples_with_tested <- NA
    n_samples_without_tested <- NA
    n_samples_na_tested <- NA
    median_value <- NA
    return_test_stats <- tribble(~event1, ~event2, ~event_type, ~test_performed,
                                 ~n_samples_with, ~n_samples_without,
                                 ~n_samples_with_tested, 
                                 ~n_samples_without_tested,
                                 ~n_samples_na_tested, ~median_value, 
                                 ~test_statistic, ~p.value,
                                 event1,  event2,  event_type,  test_performed,
                                 n_samples_with,  
                                 n_samples_without,
                                 n_samples_with_tested,  n_samples_without_tested,
                                 n_samples_na_tested, median_value, 
                                 test_statistic,  p.value)
    
    if (return_tibble) {
      return(seqfish_tbl)
    } else {
      return(return_test_stats)
    }
  }
  
  # ============================================================================
  # Expression Functions
  # ============================================================================
  
  # Returns median expression percentile, t-test, and Fisher's Exact Test p-values
  test_event_expression <- function(samples_with, samples_without, 
                                    this_gene, event_type, expression_tbl,
                                    t_test = FALSE, outlier = FALSE){
    # Given a gene, find out if expression differs between event / not event
    expression_with_event <- expression_tbl  %>%
      filter(gene == this_gene) %>% semi_join(samples_with, by = "srr")
    expression_without_event <- expression_tbl %>% 
      filter(gene == this_gene) %>% anti_join(samples_with, by = "srr")
    
    # Median expression percentile of event samples
    median_pct <- as.double(expression_with_event %>% summarize(median(pct)))
    
    # Testing
    if ( t_test & outlier ) {
      stop("Variables t_test and outlier are both TRUE. 
           Either t_test xor outlier must be TRUE.")
    } else if (t_test) {
      test_performed <- "Student's t-Test"
      if ( median_pct >= 0.5 ) {
        test_result <- t.test(expression_with_event$log10tpm, 
                              expression_without_event$log10tpm, 
                              alternative = "greater")
        
      } else{
        test_result <- t.test(expression_with_event$log10tpm, 
                              expression_without_event$log10tpm, 
                              alternative = "less")
      }
      test_statistic <- test_result$statistic
      p.value <- test_result$p.value
      
    } else if (outlier) {
      test_performed <- "Fisher's Exact Test"
      if ( median_pct >= 0.5 ) {
        over_table <- expression_primary %>% 
          group_by(has_event = srr %in% samples_with$srr) %>%
          filter(gene == this_gene) %>% select(has_event, outlier_over_tpm) %>%
          summarize(outlier = sum(outlier_over_tpm == 1),
                    not_outlier = sum(outlier_over_tpm == 0)) %>%
          select(outlier, not_outlier)
        test_result <- fisher.test(over_table, alternative = "t")
      } else{
        under_table <- expression_primary %>%
          group_by(has_event = srr %in% samples_with$srr) %>%
          filter(gene == this_gene) %>% select(has_event, outlier_under_tpm) %>%
          summarize(outlier = sum(outlier_under_tpm == 1), 
                    not_outlier = sum(outlier_under_tpm == 0)) %>%
          select(outlier, not_outlier)
        test_result <- fisher.test(under_table, alternative = "t")
      }
      test_statistic <- test_result$estimate
      p.value <- test_result$p.value
      
    } else{
      stop("Variables t_test and outlier are both FALSE. 
           Either t_test xor outlier must be TRUE.")
    }
    
    # Return value
    event1 <- this_gene
    event2 <- NA
    event_type <- event_type
    n_samples_with <- samples_with %>% nrow()
    n_samples_without <- samples_without %>% nrow()
    n_samples_with_tested <- NA
    n_samples_without_tested <- NA
    n_samples_na_tested <- NA
    median_value <- median_pct
    return_test_stats <- tribble(~event1, ~event2, ~event_type, ~test_performed,
                                 ~n_samples_with, ~n_samples_without,
                                 ~n_samples_with_tested, 
                                 ~n_samples_without_tested,
                                 ~n_samples_na_tested, ~median_value, 
                                 ~test_statistic, ~p.value,
                                 event1,  event2,  event_type,  test_performed,
                                 n_samples_with,  
                                 n_samples_without,
                                 n_samples_with_tested,  n_samples_without_tested,
                                 n_samples_na_tested, median_value, 
                                 test_statistic,  p.value)
    
    }
  
  # ============================================================================
  # Assign seqFISH and clinical variables to approapiate lists
  # ============================================================================
  
  seqfish_gene_names <- c("seqfish_Translocation_WHSC1_4_14", 
                          "seqfish_Translocation_CCND3_6_14", 
                          "seqfish_Translocation_MYC_8_14", 
                          "seqfish_Translocation_MAFA_8_14", 
                          "seqfish_Translocation_CCND1_11_14", 
                          "seqfish_Translocation_CCND2_12_14", 
                          "seqfish_Translocation_MAF_14_16", 
                          "seqfish_Translocation_MAFB_14_20")
  
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
                                        "Plasmacytoma")
  
  continuous_clinical_variable_names <- c("Age", 
                                          "BM_Plasma_Cell_Percent", 
                                          "LDH")
  
  # ============================================================================
  # Prepare lists of fusions and genes for testing (avoid testing rare events)
  # ============================================================================
  
  # Fusions pairs seen in at least 3 samples
  fusion_pairs_gt2 <- fusions_primary %>% group_by(fusion) %>% 
    summarize(count = n()) %>% filter(count >= 3) %>% arrange(desc(count))
  
  # Fusion genes seen in at least 3 samples
  fusion_genes_gt2 <- fusions_primary %>% 
    gather(geneA, geneB, key = "geneAB", value = "fusion_gene") %>%
    select(mmrf, srr, fusion_gene) %>% 
    distinct() %>% group_by(fusion_gene) %>% 
    summarize(count = n()) %>% filter(count >= 3) %>% arrange(desc(count))
  
  # ============================================================================
  # Business
  # ============================================================================
  
  if (recreate_testing_tbl) {
    testing_tbl <- tribble(~event1, ~event2, ~event_type, ~test_performed,
                           ~n_samples_with, ~n_samples_without,
                           ~n_samples_with_tested, ~n_samples_without_tested,
                           ~n_samples_na_tested, ~median_value, 
                           ~test_statistic, ~p.value)
    
    # Test expression of genes involved in seqFISH translocations
    for (full_seqfish_name in seqfish_gene_names) {
      print(full_seqfish_name)
      gene <- str_split(full_seqfish_name, "_", simplify = TRUE)[3]
      print(gene)
      samples_with <- get_ids_with_seqfish(full_seqfish_name, 
                                           seqfish_tbl = seqfish_clinical_info, 
                                           samples_tbl = samples_primary)
      samples_without <- get_ids_without_seqfish(full_seqfish_name, 
                                                 seqfish_tbl = seqfish_clinical_info, 
                                                 samples_tbl = samples_primary)
      
      new_ttest_row <- test_event_expression(
        samples_with, samples_without, this_gene = gene, 
        event_type = "seqFISH Expression", 
        expression_tbl = expression_primary,
        t_test = TRUE, outlier = FALSE)
      testing_tbl <- bind_rows(testing_tbl, new_ttest_row)
      
      new_outlier_row <- test_event_expression(
        samples_with, samples_without, this_gene = gene, 
        event_type = "seqFISH Expression Outlier", 
        expression_tbl = expression_primary,
        t_test = FALSE, outlier = TRUE)
      testing_tbl <- bind_rows(testing_tbl, new_outlier_row)
    }
    
    # Test genes recurrently involved in fusions
    for (gene in fusion_genes_gt2$fusion_gene) {
      if ( gene %in% c("IGH", "IGK", "IGL", 
                       "IGHpseudo", "IGKpseudo", "IGLpseudo")) {
        next
      }
      print(gene)
      samples_with <- get_ids_with_gene(gene, fusions_primary)
      samples_without <- get_ids_without_gene(gene, fusions_primary, 
                                              samples_primary)
      
      new_ttest_row <- test_event_expression(
        samples_with, samples_without, this_gene = gene, 
        event_type = "Fusion Expression", 
        expression_tbl = expression_primary,
        t_test = TRUE, outlier = FALSE)
      testing_tbl <- bind_rows(testing_tbl, new_ttest_row)
      
      new_outlier_row <- test_event_expression(
        samples_with, samples_without, this_gene = gene, 
        event_type = "Fusion Expression Outlier", 
        expression_tbl = expression_primary,
        t_test = FALSE, outlier = TRUE)
      testing_tbl <- bind_rows(testing_tbl, new_outlier_row)
      
      for (this_feature in seqfish_variable_names) {
        new_seqfish_row <- test_event_clinical_discrete(
          samples_with, samples_without, event = gene, 
          event_type = "Fusion seqFISH", 
          clinical_tbl = seqfish_clinical_info, 
          clinical_feature = this_feature, fisher_test = TRUE, 
          return_tibble = FALSE, return_table = FALSE)
        testing_tbl <- bind_rows(testing_tbl, new_seqfish_row)
      }
      
      for (this_feature in discrete_clinical_variable_names) {
        new_discrete_row <- test_event_clinical_discrete(
          samples_with, samples_without, event = gene, 
          event_type = "Fusion Clinical", 
          clinical_tbl = seqfish_clinical_info, 
          clinical_feature = this_feature, fisher_test = TRUE, 
          return_tibble = FALSE, return_table = FALSE)
        testing_tbl <- bind_rows(testing_tbl, new_discrete_row)
      }
      
      for (this_feature in continuous_clinical_variable_names) {
        new_continuous_row_ttest <- test_event_clinical_continuous(
          samples_with, samples_without, event = gene, 
          event_type = "Fusion Clinical", 
          clinical_tbl = seqfish_clinical_info, 
          clinical_feature = this_feature, t_test = TRUE, mwu_test = FALSE, 
          return_tibble = FALSE)
        testing_tbl <- bind_rows(testing_tbl, new_continuous_row_ttest)
        
        new_continuous_row_mwu <- test_event_clinical_continuous(
          samples_with, samples_without, event = gene, 
          event_type = "Fusion Clinical", 
          clinical_tbl = seqfish_clinical_info, 
          clinical_feature = this_feature, t_test = FALSE, mwu_test = TRUE, 
          return_tibble = FALSE)
        testing_tbl <- bind_rows(testing_tbl, new_continuous_row_mwu)
      }
    }
    
    # Test seqfish events
    for (seqfish in seqfish_variable_names) {
      print(seqfish)
      samples_with <- get_ids_with_seqfish(seqfish, 
                                           seqfish_tbl = seqfish_clinical_info, 
                                           samples_tbl = samples_primary)
      samples_without <- get_ids_without_seqfish(seqfish, 
                                                 seqfish_tbl = seqfish_clinical_info, 
                                                 samples_tbl = samples_primary)
      
      for (this_feature in discrete_clinical_variable_names) {
        new_discrete_row <- test_event_clinical_discrete(
          samples_with, samples_without, event = seqfish, 
          event_type = "seqFISH Clinical", 
          clinical_tbl = seqfish_clinical_info, 
          clinical_feature = this_feature, fisher_test = TRUE, 
          return_tibble = FALSE, return_table = FALSE)
        testing_tbl <- bind_rows(testing_tbl, new_discrete_row)
      }
      
      for (this_feature in continuous_clinical_variable_names) {
        new_continuous_row_ttest <- test_event_clinical_continuous(
          samples_with, samples_without, event = seqfish, 
          event_type = "seqFISH Clinical", 
          clinical_tbl = seqfish_clinical_info, 
          clinical_feature = this_feature, t_test = TRUE, mwu_test = FALSE, 
          return_tibble = FALSE)
        testing_tbl <- bind_rows(testing_tbl, new_continuous_row_ttest)
        
        new_continuous_row_mwu <- test_event_clinical_continuous(
          samples_with, samples_without, event = gene, 
          event_type = "Fusion Clinical", 
          clinical_tbl = seqfish_clinical_info, 
          clinical_feature = this_feature, t_test = FALSE, mwu_test = TRUE, 
          return_tibble = FALSE)
        testing_tbl <- bind_rows(testing_tbl, new_continuous_row_mwu)
      }
      
    }
    
    testing_tbl_pvalue_adjusted <- testing_tbl %>% 
      filter( !is.na(test_performed) ) %>%
      mutate(bonferroni = p.adjust(p.value, method = "bonferroni")) %>%
      mutate(BH = p.adjust(p.value, method = "BH")) %>%
      mutate(BY = p.adjust(p.value, method = "BY")) %>%
      mutate(fdr = p.adjust(p.value, method = "fdr"))
    
    write_tsv(testing_tbl, str_c(paper_supp, "testing_tbl.tsv"))
    write_tsv(testing_tbl_pvalue_adjusted, 
              str_c(paper_supp, "testing_tbl_pvalue_adjusted.tsv"))
    
  } else {
    testing_tbl <- read_tsv(str_c(paper_supp, "testing_tbl.tsv"))
    testing_tbl_pvalue_adjusted <- read_tsv(
      str_c(paper_supp, "testing_tbl_pvalue_adjusted.tsv"))
  }
}

# ==============================================================================
# Plot expression of samples with and without fusion or seqFISH events
# Originally written October 2018, Updated April 2019
# ==============================================================================

if (TRUE) {

  # Set seed for reproducibility
  set.seed(10)
  
  # Change to TRUE only if you want to recreate all plots from the start
  recreate_plot_df <- FALSE
  recreate_all_plots <- FALSE
  
  dir.create(paper_supp, showWarnings = FALSE)
  dir.create(str_c(paper_supp, "all_plots"), showWarnings = FALSE)
  dir.create(str_c(paper_supp, "significant_plots"), showWarnings = FALSE)
  dir.create(str_c(paper_supp, "seqFISH"), showWarnings = FALSE)
  dir.create(str_c(paper_supp, "seqFISH/fusions"), showWarnings = FALSE)
  dir.create(str_c(paper_supp, "seqFISH/translocations"), showWarnings = FALSE)
  
  testing_tbl <- read_tsv(str_c(paper_supp, "testing_tbl.tsv"))
  testing_tbl_pvalue_adjusted <- read_tsv(
    str_c(paper_supp, "testing_tbl_pvalue_adjusted.tsv"))
  
  ymax_expression_value <- plyr::round_any(
    max(expression_primary$log10tpm), 
    accuracy = .1, f = ceiling)
  
  # ============================================================================
  # General Functions
  # ============================================================================
  
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
  
  # ============================================================================
  # Important fusion genes to keep around even if not reported in > 2 samples
  # ============================================================================
  keep_genes <- unique(c(
    fusions_primary %>% filter(geneA_driver == 1) %>% pull(geneA),
    fusions_primary %>% filter(geneB_driver == 1) %>% pull(geneB),
    fusions_primary %>% filter(geneA_oncogene == 1) %>% pull(geneA),
    fusions_primary %>% filter(geneB_oncogene == 1) %>% pull(geneB),
    fusions_primary %>% filter(geneA_kinase == 1) %>% pull(geneA),
    fusions_primary %>% filter(geneB_kinase == 1) %>% pull(geneB),
    fusions_primary %>% filter(geneA_tsg == 1) %>% pull(geneA),
    fusions_primary %>% filter(geneB_tsg == 1) %>% pull(geneB),
    fusions_primary %>% filter(geneA_mmy_known == 1) %>% pull(geneA),
    fusions_primary %>% filter(geneB_mmy_known == 1) %>% pull(geneB),
    fusions_primary %>% filter(drug_geneA == 1) %>% pull(geneA),
    fusions_primary %>% filter(drug_geneB == 1) %>% pull(geneB),
    fusions_primary %>% filter(drug_fusion == 1) %>% pull(geneA),
    fusions_primary %>% filter(drug_fusion == 1) %>% pull(geneB)
  ))
  
  # ============================================================================
  # Assign seqFISH and clinical variables to approapiate lists
  # ============================================================================
  
  seqfish_gene_names <- c("seqfish_Translocation_WHSC1_4_14", 
                          "seqfish_Translocation_CCND3_6_14", 
                          "seqfish_Translocation_MYC_8_14", 
                          "seqfish_Translocation_MAFA_8_14", 
                          "seqfish_Translocation_CCND1_11_14", 
                          "seqfish_Translocation_CCND2_12_14", 
                          "seqfish_Translocation_MAF_14_16", 
                          "seqfish_Translocation_MAFB_14_20")
  
  seqfish_gene_names_formatted <- c("Translocation t(4;14) WHSC1", 
                                    "Translocation t(6;14) CCND3", 
                                    "Translocation t(8;14) MYC", 
                                    "Translocation t(8;14) MAFA", 
                                    "Translocation t(11;14) CCND1", 
                                    "Translocation t(12;14) CCND2", 
                                    "Translocation t(14;16) MAF", 
                                    "Translocation t(14;20) MAFB")
  
  seqfish_gene_names_formatted_short <- c("t(4;14) WHSC1", 
                                          "t(6;14) CCND3", 
                                          "t(8;14) MYC", 
                                          "t(8;14) MAFA", 
                                          "t(11;14) CCND1", 
                                          "t(12;14) CCND2", 
                                          "t(14;16) MAF", 
                                          "t(14;20) MAFB")
  
  seqfish_genes <- c("WHSC1", "CCND3", "MYC", "MAFA", 
                     "CCND1", "CCND2", "MAF", "MAFB")
  
  
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
                                        "Plasmacytoma")
  
  continuous_clinical_variable_names <- c("Age", 
                                          "BM_Plasma_Cell_Percent", 
                                          "LDH")
  
  # ============================================================================
  # Prepare lists of fusions and genes for testing (avoid testing rare events)
  # ============================================================================
  
  # Fusions pairs seen in at least 3 samples
  fusion_pairs_gt2 <- fusions_primary %>% group_by(fusion) %>% 
    summarize(count = n()) %>% filter(count >= 3) %>% arrange(desc(count))
  
  # Fusion genes seen in at least 3 samples
  fusion_genes_gt2 <- fusions_primary %>% 
    gather(geneA, geneB, key = "geneAB", value = "fusion_gene") %>%
    select(mmrf, srr, fusion_gene) %>% 
    distinct() %>% group_by(fusion_gene) %>% 
    summarize(count = n()) %>% filter(count >= 3) %>% arrange(desc(count))
  
  # ============================================================================
  # Plotting functions
  # ============================================================================
  
  plot_fusion_expression_1d <- function(plot_df,
                                        gene_list,
                                        labels = FALSE,
                                        translocation = NULL,
                                        translocation_formatted = NULL,
                                        ymax_value,
                                        pdf_path,
                                        pdf_width = 10,
                                        pdf_height = 10,
                                        seed = 10,
                                        nrows = 1,
                                        pretty = FALSE){
    
    set.seed(seed)
    
    plot_df <- plot_df %>% filter(gene %in% gene_list)
    if (is.null(translocation)) {
      shape_vector <- "No shape vector given"
    } else {
      shape_vector <- plot_df %>% pull(translocation) %>%
        factor(labels = c("Not detected", "Detected", "Missing"), exclude = NULL)  
    }
    plot_df <- plot_df %>% mutate(shape_factor = shape_vector)
    
    p <- ggplot(plot_df)
    
    p <- p + facet_wrap(~ gene, nrow = nrows)
    
    p <- p + geom_violin(aes(x = fusion_indicator,
                             y = log10tpm,
                             fill = fusion_status),
                         color = "black",
                         draw_quantiles = 0.5)
    
    p <- p + geom_point(aes(x = fusion_jitter,
                            y = log10tpm,
                            color = cnv_factor,
                            shape = shape_factor))
    
    if (labels & length(gene_list) == 1) {
      p <- p + geom_label_repel(data = plot_df %>% filter(!is.na(fusion_label)),
                                aes(x = fusion_jitter,
                                    y = log10tpm,
                                    label = fusion_label,
                                    color = cnv_factor),
                                point.padding = 0.5,
                                show.legend = FALSE)
    }
    
    p <- p + scale_x_continuous(breaks = c(0, 1), 
                                labels = c("No fusion", "Fusion"))
    
    p <- p + ylim(0, ymax_value)
    
    color_scale <- c("#1f78b4", "#a6cee3", # deletions
                     "#b2df8a", # neutral 
                     "#cab2d6", #missing
                     "#fb9a99", "#e31a1c") # amplifications

    p <- p + scale_color_manual(values = color_scale, drop = FALSE)
    
    p <- p + scale_fill_manual(values = c("#ffffff", "#ffffff")) # both white
    
    if (is.null(translocation)) {
      p <- p + scale_shape_manual(values = 16)
      p <- p + guides(alpha = FALSE, fill = FALSE, shape = FALSE)
      p <- p + labs(x = NULL,
                    y = "Gene Expression TPM (log10)", 
                    color = "Copy Number")
    } else {
      p <- p + scale_shape_manual(values = c(16, 17, 4))
      p <- p + guides(alpha = FALSE, fill = FALSE)
      p <- p + labs(x = NULL,
                    y = "Gene Expression TPM (log10)", 
                    color = "Copy Number",
                    shape = translocation_formatted)
    }
    
    p <- p + ggplot2_standard_additions()
    
    if (pretty) { 
      p <- p + theme(legend.position = "bottom",
                     legend.text = element_text(size = 8),
                     legend.title = element_text(size = 10),
                     axis.ticks.y = element_blank(),
                     axis.title = element_text(size = 12),
                     axis.ticks.x = element_blank(),
                     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10, color = "grey50", face = "italic"),
                     strip.background = element_blank(),
                     strip.text = element_text(size = 10,
                                               face = "italic"),
                     panel.background = element_blank(),
                     panel.border = element_blank(),
                     panel.spacing.x = unit(0.025, units = "inches"),
                     panel.grid.major = element_line(size = 0.1),
                     panel.grid.major.x = element_blank(),
                     panel.grid.minor.x = element_blank(),
                     panel.grid.minor.y = element_blank())
    }
    
    pdf(pdf_path, width = pdf_width, height = pdf_height, useDingbats = FALSE)
    print(p)
    shh <- dev.off()
    
  }
  
  plot_translocation_expression_1d <- function(plot_df,
                                               gene_list,
                                               labels = FALSE,
                                               ymax_value,
                                               pdf_path,
                                               pdf_width = 10, 
                                               pdf_height = 10,
                                               seed = 10,
                                               pretty = FALSE){
    
    set.seed(seed)
    
    plot_df <- plot_df %>% filter(gene %in% gene_list, 
                                  !is.na(translocation_indicator))
    
    p <- ggplot(plot_df)
    
    p <- p + facet_wrap(~ seqfish_gene_names_formatted, nrow = 1)
    
    p <- p + geom_violin(aes(x = translocation_indicator,
                             y = log10tpm,
                             fill = factor(translocation_indicator)),
                         color = "black",
                         draw_quantiles = 0.5)
    
    p <- p + geom_point(aes(x = translocation_jitter,
                            y = log10tpm,
                            color = cnv_factor),
                        shape = 16)
    
    if (labels & length(gene_list) == 1) {
      p <- p + geom_label_repel(data = plot_df %>% filter(!is.na(fusion_label)),
                                aes(x = translocation_jitter,
                                    y = log10tpm,
                                    label = fusion_label,
                                    color = cnv_factor),
                                point.padding = 0.5,
                                show.legend = FALSE)
    }
    
    p <- p + scale_x_continuous(breaks = c(0, 1), 
                                labels = c("No translocation", "Translocation"))
    
    p <- p + ylim(0, ymax_value)
    
    color_scale <- c("#1f78b4", "#a6cee3", # deletions
                     "#b2df8a", # neutral 
                     "#cab2d6", # missing
                     "#fb9a99", "#e31a1c") # amplifications

    p <- p + scale_color_manual(values = color_scale, drop = FALSE)
    
    p <- p + scale_fill_manual(values = c("#ffffff", "#ffffff")) # both white
    
    p <- p + guides(alpha = FALSE, fill = FALSE)
    
    p <- p + labs(x = NULL,
                  y = "Gene Expression TPM (log10)", 
                  color = "Copy Number")
    
    p <- p + ggplot2_standard_additions()
    
    if (pretty) { 
      p <- p + theme(legend.position = "bottom",
                     legend.text = element_text(size = 8),
                     legend.title = element_text(size = 10),
                     axis.ticks.y = element_blank(),
                     axis.title = element_text(size = 12),
                     axis.ticks.x = element_blank(),
                     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10, color = "grey50", face = "italic"),
                     strip.background = element_blank(),
                     strip.text = element_text(size = 10,
                                               face = "italic"),
                     panel.background = element_blank(),
                     panel.border = element_blank(),
                     panel.spacing.x = unit(0.025, units = "inches"),
                     panel.grid.major = element_line(size = 0.1),
                     panel.grid.major.x = element_blank(),
                     panel.grid.minor.x = element_blank(),
                     panel.grid.minor.y = element_blank())
    }
    
    pdf(pdf_path, width = pdf_width, height = pdf_height, useDingbats = FALSE)
    print(p)
    shh <- dev.off()
    
  }
  
  plot_expression_2d <- function(plot_df,
                                 gene_list,
                                 fusion1 = NULL,
                                 fusion2 = NULL,
                                 translocation = NULL,
                                 translocation_formatted = NULL,
                                 ymax_value,
                                 pdf_path,
                                 pdf_width = 10,
                                 pdf_height = 10,
                                 seed = 10,
                                 pretty = FALSE){
    
    set.seed(seed)
    
    if (length(gene_list) != 2) {
      stop("Gene list can only have two genes for 2D plot.")
    }
    
    plot_df <- plot_df %>% filter(gene %in% gene_list)
    shape_vector <- plot_df %>% pull(translocation) %>% 
      factor(labels = c("Not detected", "Detected", "Missing"), exclude = NULL)
    plot_df <- plot_df %>% mutate(shape_factor = shape_vector)
    
    fusion1_reverse <- str_c(rev(
      str_split(fusion1, pattern = "--", simplify = TRUE)), collapse = "--")
    gene1 <- plot_df %>% filter(gene == gene_list[1]) %>%
      select(mmrf, srr, gene, log10tpm, fusion_label, shape_factor) %>%
      mutate(has_fusion1 = str_detect(fusion_label, pattern = fusion1) | 
               str_detect(fusion_label, pattern = fusion1_reverse))
    
    fusion2_reverse <- str_c(rev(
      str_split(fusion2, pattern = "--", simplify = TRUE)), collapse = "--")
    gene2 <- plot_df %>% filter(gene == gene_list[2]) %>%
      select(mmrf, srr, gene, log10tpm, fusion_label, shape_factor) %>%
      mutate(has_fusion2 = str_detect(fusion_label, pattern = fusion2) | 
               str_detect(fusion_label, pattern = fusion2_reverse))
    
    plot_df <- left_join(gene1, gene2, by = c("mmrf", "srr", "shape_factor")) %>%
      replace_na(list(has_fusion1 = FALSE, has_fusion2 = FALSE)) %>%
      mutate(fusion_category = has_fusion1 + 2*has_fusion2) %>%
      mutate(fusion_category = factor(fusion_category,
                                      levels = c(0,1,2,3),
                                      labels = c("None reported",
                                                 fusion1,
                                                 fusion2,
                                                 "Both reported")
                                      
      )
      )
    
    if (pretty) {
      plot_df <- plot_df %>% mutate(fusion_category = case_when( fusion_category == "None reported" ~ "None Reported",
                                                                 TRUE ~ str_c(fusion1, " or\n", fusion2)))
    }
    
    p <- ggplot(plot_df)
    
    p <- p + geom_point(aes(x = log10tpm.x,
                            y = log10tpm.y,
                            color = fusion_category,
                            shape = shape_factor))
    
    p <- p + xlim(0, ymax_value)
    
    p <- p + ylim(0, ymax_value)
    
    p <- p + scale_shape_manual(values = c(16, 17, 4))
    
    p <- p + guides(alpha = FALSE, fill = FALSE)
    
    p <- p + labs(x = str_c(gene_list[1], " Expression TPM (log10)"),
                  y = str_c(gene_list[2], " Expression TPM (log10)"), 
                  color = "Fusion status",
                  shape = translocation_formatted)
    
    p <- p + ggplot2_standard_additions()
    
    p <- p + coord_fixed(ratio = 1)
    
    if (pretty) { 
      p <- p + theme(panel.grid.major = element_line(size = 0.1),
                     panel.background = element_blank(),
                     panel.border = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.ticks.y = element_blank(),
                     legend.background = element_blank(),
                     legend.box = "vertical",
                     legend.position = "bottom",
                     legend.text = element_text(size = 8),
                     legend.title = element_text(size = 10),
                     axis.title = element_text(size = 12)) +
        #guides(shape = FALSE, color = FALSE) +
        scale_color_brewer(palette = "Oranges", direction = -1)
        
    }
    
    pdf(pdf_path, width = pdf_width*2, height = pdf_height*2, useDingbats = FALSE)
    print(p)
    shh <- dev.off()
    
    pdf(str_c(pdf_path, ".no_legend.pdf"), width = pdf_width, height = pdf_height, useDingbats = FALSE)
    p <- p + guides(color = FALSE, shape = FALSE)
    print(p)
    shh <- dev.off()
    
  }
  
  # ============================================================================
  # Create plot data frame
  # ============================================================================
  
  if (recreate_plot_df) {
    return_fusions <- function(fusions_df, this_srr, gene){
      return_value <- fusions_df %>% filter(srr == this_srr,
                                            geneA == gene | geneB == gene) %>%
        pull(fusion) %>% str_c(collapse = "\n")
      if (identical(return_value, character(0))) {
        return_value <- NA
      }
      return(return_value)
    }
    
    return_translocations <- function(seqfish_df, this_mmrf, gene, 
                                      seqfish_gene_names, seqfish_genes){
      if (gene %in% seqfish_genes) {
        return_value <- seqfish_df %>% filter(mmrf == this_mmrf) %>%
          pull( seqfish_gene_names[which(seqfish_genes == gene)] )
      } else {
        return_value <- NA
      }
      return(return_value)
    }
    
    return_translocations_formatted <- function(gene, seqfish_genes, 
                                                seqfish_gene_names_formatted){
      if (gene %in% seqfish_genes) {
        return_value <- seqfish_gene_names_formatted[which(seqfish_genes == gene)]
      } else {
        return(NA)
      }
      return(return_value)
    }
    
    categorical_cnv <- function(vec) {
      sd_factor <- sd(vec, na.rm = TRUE)
      b = c(2 + sd_factor*c(-Inf, -3, -1, 1, 3, Inf))
      cnv_categories <- vec %>% cut(breaks = b)
      return(cnv_categories)
    }
    
    plot_df <- expression_primary %>% 
      filter(gene %in% fusion_genes_gt2$fusion_gene | 
               gene %in% seqfish_genes |
               gene %in% keep_genes) %>% 
      mutate(categorical_cnv = categorical_cnv(2*2^gene_avg_cnv)) %>% 
      mutate(cnv_factor = factor(categorical_cnv, 
                                 labels = c("DELETION",
                                            "Deletion",
                                            "Neutral",
                                            "Missing",
                                            "Amplification",
                                            "AMPLIFICATION"), 
                                 exclude = NULL)) %>%
      rowwise() %>% 
      mutate(translocation_indicator = return_translocations(
        seqfish_clinical_info, mmrf, gene, seqfish_gene_names, seqfish_genes)) %>%
      mutate(seqfish_gene_names_formatted = return_translocations_formatted(
        gene, seqfish_genes, seqfish_gene_names_formatted)) %>%
      mutate(fusion_label = return_fusions(fusions_primary, srr, gene)) %>%
      ungroup() %>%
      mutate(fusion_status = !is.na(fusion_label)) %>%
      mutate(fusion_indicator = as.numeric(fusion_status)) %>%
      mutate(fusion_jitter = jitter(fusion_indicator)) %>%
      mutate(translocation_jitter = jitter(translocation_indicator)) %>%
      left_join(seqfish_clinical_info, by = "mmrf")
    
    write_tsv(plot_df, str_c(paper_supp, "expression_plot_tibble.tsv"))
  } else {
    plot_df <- read_tsv(str_c(paper_supp, "expression_plot_tibble.tsv"))
    plot_df <- plot_df %>% mutate(cnv_factor = factor(categorical_cnv,
                                                      labels = c("DELETION",
                                                                 "Deletion",
                                                                 "Neutral",
                                                                 "Missing",
                                                                 "Amplification",
                                                                 "AMPLIFICATION"), 
                                                      exclude = NULL))
  }
  
  # ============================================================================
  # Business
  # ============================================================================
  
  if (recreate_all_plots) {
    
    # Plot expression of genes involved in seqFISH translocations
    for (this_gene in seqfish_genes) {
      print(this_gene)
      n_samples_with_fusion <- plot_df %>% 
        filter(gene == this_gene, !is.na(fusion_label)) %>% nrow()
      if (n_samples_with_fusion > 2) {
        t_index = which(seqfish_genes == this_gene)
        t_name = seqfish_gene_names[t_index]
        t_format = seqfish_gene_names_formatted_short[t_index]
        plot_fusion_expression_1d(plot_df,
                                  this_gene,
                                  labels = FALSE,
                                  translocation = t_name,
                                  translocation_formatted = t_format,
                                  ymax_value = ymax_expression_value,
                                  pdf_path = str_c(paper_supp, 
                                                   "seqFISH/fusions/", 
                                                   this_gene, ".pdf"),
                                  pdf_width = 10,
                                  pdf_height = 10,
                                  seed = 10)
        plot_translocation_expression_1d(plot_df,
                                         this_gene,
                                         labels = FALSE,
                                         ymax_value = ymax_expression_value,
                                         pdf_path = str_c(paper_supp, "seqFISH/",
                                                          "translocations/",
                                                          this_gene, ".pdf"),
                                         pdf_width = 10, 
                                         pdf_height = 10,
                                         seed = 10)
      }
    }
    
    plot_fusion_expression_1d(plot_df,
                              seqfish_genes,
                              labels = FALSE,
                              ymax_value = ymax_expression_value,
                              pdf_path = str_c(paper_supp, "seqFISH/fusions/",
                                               "all.pdf"),
                              pdf_width = 4*length(seqfish_gene_names),
                              pdf_height = 10,
                              seed = 10)  
    
    plot_translocation_expression_1d(plot_df,
                                     seqfish_genes,
                                     labels = FALSE,
                                     ymax_value = ymax_expression_value,
                                     pdf_path = str_c(paper_supp, 
                                                      "seqFISH/translocations/",
                                                      "all.pdf"),
                                     pdf_width = 4*length(seqfish_gene_names),
                                     pdf_height = 10,
                                     seed = 10)  
    
    # Plot expression of genes recurrently involved in fusions
    for (this_gene in fusion_genes_gt2$fusion_gene) {
      if ( this_gene %in% c("IGH", "IGK", "IGL",
                            "IGHpseudo", "IGKpseudo", "IGLpseudo")) {
        next
      }
      print(this_gene)
      plot_fusion_expression_1d(plot_df,
                                this_gene,
                                labels = FALSE,
                                ymax_value = ymax_expression_value,
                                pdf_path = str_c(paper_supp, "all_plots/", 
                                                 this_gene, ".pdf"),
                                pdf_width = 10,
                                pdf_height = 10,
                                seed = 10)  
      plot_fusion_expression_1d(plot_df,
                                this_gene,
                                labels = TRUE,
                                ymax_value = ymax_expression_value,
                                pdf_path = str_c(paper_supp, "all_plots/", 
                                                 this_gene, ".labeled.pdf"),
                                pdf_width = 10,
                                pdf_height = 10,
                                seed = 10)
    }
    
    # Plot expression of genes from driver, oncogene, tsg, and druggable lists
    for (this_gene in keep_genes) {
      if ( this_gene %in% c("IGH", "IGK", "IGL",
                            "IGHpseudo", "IGKpseudo", "IGLpseudo")) {
        next
      }
      if ( this_gene %in% fusion_genes_gt2$fusion_gene) {
        next
      }
      print(this_gene)
      plot_fusion_expression_1d(plot_df,
                                this_gene,
                                labels = FALSE,
                                ymax_value = ymax_expression_value,
                                pdf_path = str_c(paper_supp, "all_plots/", 
                                                 this_gene, ".pdf"),
                                pdf_width = 10,
                                pdf_height = 10,
                                seed = 10)  
      plot_fusion_expression_1d(plot_df,
                                this_gene,
                                labels = TRUE,
                                ymax_value = ymax_expression_value,
                                pdf_path = str_c(paper_supp, "all_plots/", 
                                                 this_gene, ".labeled.pdf"),
                                pdf_width = 10,
                                pdf_height = 10,
                                seed = 10)
    }
    
    significant_fusion_expresion_genes <- testing_tbl_pvalue_adjusted %>% 
      filter(fdr < 0.05 | median_value > 0.9, 
             event_type %in% c("Fusion Expression", 
                               "Fusion Expression Outlier")) %>% 
      pull(event1) %>% unique()
    
    for (this_gene in significant_fusion_expresion_genes) {
      print(this_gene)
      n_samples_with_fusion <- plot_df %>% 
        filter(gene == this_gene, !is.na(fusion_label)) %>% nrow()
      if (n_samples_with_fusion > 2) {
        plot_fusion_expression_1d(plot_df,
                                  this_gene, 
                                  labels = FALSE,
                                  ymax_value = ymax_expression_value,
                                  pdf_path = str_c(paper_supp, 
                                                   "significant_plots/",
                                                   this_gene, ".pdf"),
                                  pdf_width = 10,
                                  pdf_height = 10,
                                  seed = 10)  
        plot_fusion_expression_1d(plot_df,
                                  this_gene, 
                                  labels = TRUE,
                                  ymax_value = ymax_expression_value,
                                  pdf_path = str_c(paper_supp, 
                                                   "significant_plots/",
                                                   this_gene, ".labeled.pdf"),
                                  pdf_width = 10,
                                  pdf_height = 10,
                                  seed = 10)  
      }
    }
    
    plot_fusion_expression_1d(plot_df,
                              significant_fusion_expresion_genes,
                              labels = FALSE,
                              ymax_value = ymax_expression_value,
                              pdf_path = str_c(paper_supp, "significant_plots/", 
                                               "all.pdf"),
                              pdf_width = 4*length(significant_fusion_expresion_genes),
                              pdf_height = 10,
                              seed = 10)
  }
  
  hi <- testing_tbl_pvalue_adjusted %>%
    filter(fdr < 0.05 | median_value > 0.9,
           event_type %in% c("Fusion Expression",
                             "Fusion Expression Outlier")) %>%
    pull(event1) %>% unique() %>% sort()
  
  overexpressed_interesting_genes <- sort(unique(c(
    fusions_primary %>% filter(geneA %in% hi, geneA_driver == 1) %>% pull(geneA),
    fusions_primary %>% filter(geneB %in% hi, geneB_driver == 1) %>% pull(geneB),
    fusions_primary %>% filter(geneA %in% hi, geneA_oncogene == 1) %>% pull(geneA),
    fusions_primary %>% filter(geneB %in% hi, geneB_oncogene == 1) %>% pull(geneB),
    fusions_primary %>% filter(geneA %in% hi, geneA_kinase == 1) %>% pull(geneA),
    fusions_primary %>% filter(geneB %in% hi, geneB_kinase == 1) %>% pull(geneB),
    fusions_primary %>% filter(geneA %in% hi, geneA_mmy_known == 1) %>% pull(geneA),
    fusions_primary %>% filter(geneB %in% hi, geneB_mmy_known == 1) %>% pull(geneB),
    fusions_primary %>% filter(geneA %in% hi, geneA_tsg == 1) %>% pull(geneA),
    fusions_primary %>% filter(geneB %in% hi, geneB_tsg == 1) %>% pull(geneB),
    fusions_primary %>% filter(geneA %in% hi, drug_fusion == 1) %>% pull(geneA),
    fusions_primary %>% filter(geneB %in% hi, drug_fusion == 1) %>% pull(geneB),
    fusions_primary %>% filter(geneA %in% hi, drug_geneA == 1) %>% pull(geneA),
    fusions_primary %>% filter(geneB %in% hi, drug_geneB == 1) %>% pull(geneB))))
  
  gene_attributes_A <- fusions_primary %>% 
    filter(geneA %in% overexpressed_interesting_genes) %>%
    mutate(gene = geneA) %>%
    select(gene, geneA_driver, geneA_oncogene, geneA_kinase, geneA_mmy_known, 
           geneA_tsg, drug_fusion, drug_geneA) %>% unique() %>% 
    gather("geneA_driver", "geneA_oncogene", "geneA_kinase", "geneA_mmy_known", 
           "geneA_tsg", "drug_fusion", "drug_geneA", key = "category", value = "status") %>% 
    mutate(category = str_remove_all(category, pattern = "geneA_")) %>% 
    mutate(category = str_remove_all(category, pattern = "_geneA")) %>% 
    mutate(category = case_when(category == "drug_fusion" ~ "Drug",
                                category == "drug" ~ "Drug",
                                category == "driver" ~ "Driver",
                                category == "kinase" ~ "Kinase",
                                category == "mmy_known" ~ "Multiple\nMyeloma\nRelated",
                                category == "oncogene" ~ "Oncogene",
                                category == "tsg" ~ "Tumor Supp")) %>%
    unique() %>% filter(status == 1)
  
  gene_attributes_B <- fusions_primary %>% 
    filter(geneB %in% overexpressed_interesting_genes) %>% 
    mutate(gene = geneB) %>%
    select(gene, geneB_driver, geneB_oncogene, geneB_kinase, geneB_mmy_known, 
           geneB_tsg, drug_fusion, drug_geneB) %>% unique() %>% 
    gather("geneB_driver", "geneB_oncogene", "geneB_kinase", "geneB_mmy_known", 
           "geneB_tsg", "drug_fusion", "drug_geneB", key = "category", value = "status") %>% 
    mutate(category = str_remove_all(category, pattern = "geneB_")) %>% 
    mutate(category = str_remove_all(category, pattern = "_geneB")) %>% 
    mutate(category = case_when(category == "drug_fusion" ~ "Drug",
                                category == "drug" ~ "Drug",
                                category == "driver" ~ "Driver",
                                category == "kinase" ~ "Kinase",
                                category == "mmy_known" ~ "Multiple\nMyeloma\nRelated",
                                category == "oncogene" ~ "Oncogene",
                                category == "tsg" ~ "Tumor Supp")) %>%
    unique() %>% filter(status == 1)
  
  gene_attributes <- bind_rows(gene_attributes_A, gene_attributes_B) %>% unique()
  
  interesting_ymax <- plot_df %>% filter(gene %in% overexpressed_interesting_genes) %>% 
    pull(log10tpm) %>% max() %>% ceiling()
  plot_fusion_expression_1d(plot_df, overexpressed_interesting_genes, 
                            ymax_value = interesting_ymax, 
                            pdf_path = str_c(paper_main, "overexpressed_interesting.pdf"), 
                            pdf_width = 8, 
                            pdf_height = 4,
                            nrows = 1, 
                            pretty = TRUE)
  
  ggplot(data = gene_attributes, aes(x = gene, 
                                     y = fct_reorder(category, 
                                                     desc(category)))) + 
    geom_point(shape = 18, size = 3) + 
    theme_bw() +
    labs(y = NULL, x = NULL) +
    theme(panel.background = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.border = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_text(size = 10, color = "grey50"),
          panel.grid.major = element_line(size = 0.1)) +
    ggsave(str_c(paper_main, "gene_attributes.pdf"), height = 1, width = 8)
  
  fgfr3_whsc1_ymax <- plot_df %>% filter(gene %in% c("FGFR3", "WHSC1")) %>% 
    pull(log10tpm) %>% max() %>% ceiling()  
  plot_expression_2d(plot_df,
                     c("WHSC1", "FGFR3"),
                     fusion1 = "IGH--WHSC1",
                     fusion2 = "IGH--FGFR3",
                     translocation = "seqfish_Translocation_WHSC1_4_14",
                     translocation_formatted = "t(4;14)",
                     ymax_value = fgfr3_whsc1_ymax,
                     pdf_path = str_c(paper_main, "FGFR3_WHSC1_t414.pdf"),
                     pdf_width = 4,
                     pdf_height = 4,
                     seed = 10, 
                     pretty = TRUE)
}

# ==============================================================================
# Plot expression of oncogenes, tumor suppressors, kinases
# Written April 2019
# ==============================================================================

if (TRUE) {
  plot_df <- bind_rows(fusions_primary %>% filter(geneA_oncogene == 1) %>% select(geneA, geneA_pct) %>% mutate(category = "Oncogene", geneAB = "geneA") %>% rename(gene = geneA, pct = geneA_pct),
                       fusions_primary %>% filter(geneB_oncogene == 1) %>% select(geneB, geneB_pct) %>% mutate(category = "Oncogene", geneAB = "geneB") %>% rename(gene = geneB, pct = geneB_pct),
                       fusions_primary %>% filter(geneA_tsg == 1) %>% select(geneA, geneA_pct) %>% mutate(category = "Tumor\nSuppressor", geneAB = "geneA") %>% rename(gene = geneA, pct = geneA_pct),
                       fusions_primary %>% filter(geneB_tsg == 1) %>% select(geneB, geneB_pct) %>% mutate(category = "Tumor\nSuppressor", geneAB = "geneB") %>% rename(gene = geneB, pct = geneB_pct),
                       fusions_primary %>% filter(geneA_kinase == 1) %>% select(geneA, geneA_pct) %>% mutate(category = "Kinase", geneAB = "geneA") %>% rename(gene = geneA, pct = geneA_pct),
                       fusions_primary %>% filter(geneB_kinase == 1) %>% select(geneB, geneB_pct) %>% mutate(category = "Kinase", geneAB = "geneB") %>% rename(gene = geneB, pct = geneB_pct))
  
  ggplot(data = plot_df, aes(x = geneAB, y = pct)) + 
    geom_violin(size = 1, color = "grey80", scale = "width") + 
    geom_jitter(height = 0, width = 0.2, shape = 16) + 
    facet_wrap(~ category, nrow = 1) +
    theme_bw() +
    scale_y_continuous(expand = c(0,0), limits = c(-0.05, 1.02)) +
    scale_x_discrete(expand = c(0,0)) +
    labs(x = NULL, y = "Expression Percentile") +
    annotate("text", label = "Gene A", x = "geneA", y = -0.02, size = 2, color = "grey50", vjust = 1) +
    annotate("text", label = "Gene B", x = "geneB", y = -0.02, size = 2, color = "grey50", vjust = 1) +
    theme(panel.background = element_blank(),
          panel.border = element_blank(),
          plot.background = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          strip.background = element_blank(),
          axis.ticks = element_blank(),
          panel.spacing.x = unit(0.1, units = "inches"),
          axis.title = element_text(size = 12),
          axis.text.x = element_blank(),
          strip.text = element_text(size = 10)
          ) +
    ggsave(str_c(paper_main, "kinase_oncogene_tsg.pdf"),
           width = 4, height = 4, useDingbats = FALSE)
}
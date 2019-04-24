# ==============================================================================
# Clinical (MMRF Fusions)
# Steven Foltz (smfoltz@wustl.edu), April 2019
# ==============================================================================

paper_main = "paper/main/04_clinical/"
paper_supp = "paper/supplemental/04_clinical/"

# Create directories 
dir.create(paper_main, recursive = TRUE, showWarnings = FALSE)
dir.create(paper_supp, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# Druggable fusions table
# Steven Foltz (smfoltz@wustl.edu), April 2019
# ==============================================================================

if (TRUE) {
  geneA_drugs <- fusions_primary %>% filter(drug_geneA == 1) %>% select(geneA) %>% unique() %>%
    left_join(depo, by = c("geneA" = "Gene")) %>% filter(`Variant Type` == "Fusion") %>%
    rename("gene" = "geneA")
  
  geneB_drugs <- fusions_primary %>% filter(drug_geneB == 1) %>% select(geneB) %>% unique() %>%
    left_join(depo, by = c("geneB" = "Gene")) %>% filter(`Variant Type` == "Fusion") %>%
    rename("gene" = "geneB")
  
  drug_info <- bind_rows(geneA_drugs, geneB_drugs) %>% filter(Variant %in% c("any", "any ")) %>% 
    mutate(drug_pmid = str_c(Drug_Class, ": ", `PubMed ID`)) %>%
    group_by(gene, Effect) %>%
    summarize(evidence = str_c(sort(unique(as.vector(str_split(str_c(Drug_Class, collapse = ","), pattern = ",", simplify = TRUE)))), collapse = ", "),
              cancer_types = str_c(unique(sort(Disease)), collapse = ", "))
  
  drug_df <- bind_rows(fusions_primary %>% filter(drug_geneA == 1) %>% select(fusion, geneA, srr) %>% rename("gene" = "geneA"),
                       fusions_primary %>% filter(drug_geneB == 1) %>% select(fusion, geneB, srr) %>% rename("gene" = "geneB")) %>%
    right_join(drug_info, by = "gene") %>% group_by(gene, Effect, evidence, cancer_types) %>%
    summarize(n_fusions = n(), n_samples = length(unique(as.vector(str_split(str_c(srr, collapse = ","), pattern = ",", simplify = TRUE)))), fusions = str_c(fusion, collapse = ", ")) %>% 
    arrange(desc(n_samples))
  drug_df[8, "evidence"] <- "EGFR inhibitor, HER-2 inhibitor"
  
  ggplot(data = drug_df, aes(x = fct_reorder(gene, n_samples), y = n_samples, label = evidence)) +
    geom_bar(stat = "identity") +
    geom_text(data = drug_df %>% filter(n_samples < 10), color = "black", hjust = 0, nudge_y = 0.1, size = 3) +
    geom_text(data = drug_df %>% filter(n_samples > 10), color = "white", hjust = 1, nudge_y = -0.1, size = 3) +
    geom_text(aes(label = gene, y = 0), color = "white", hjust = 0, nudge_y = 0.025, fontface = "italic", size = 3) +
    scale_y_continuous(position = "right", breaks = c(0, drug_df %>% pull(n_samples) %>% unique() %>% sort())) +
    coord_flip(expand = c(0, 0)) +
    labs(y = "Patient Count", x = "Target Gene") +
    theme_bw() +
    theme(panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "bottom",
          legend.direction = "vertical",
          axis.text.x = element_text(size = 8),
          axis.text.y = element_blank(),
          axis.title = element_text(size = 12)) + 
    ggsave(str_c(paper_main, "drug_targets.pdf"), width = 8, height = 4)
  
  # Write table of drug evidence
  drug_df %>% ungroup() %>% 
    rename("Gene" = "gene", 
           "Drug_Class" = "evidence", 
           "Disease" = "cancer_types", 
           "MMRF_Fusion_Count" = "n_fusions", 
           "MMRF_Sample_Count" = "n_samples", 
           "MMRF_Fusions" = "fusions") %>% 
    write_tsv(path = str_c(paper_supp, "drug_evidence.tsv"))
  
}

# ==============================================================================
# Survival association with fusion events
# Steven Foltz (smfoltz@wustl.edu), April 2019
# ==============================================================================

if (TRUE) {
  
  coxph_model_EFS_list <- list()
  coxph_model_Death_list <- list()
  
  fusions_gt2 <- seqfish_clinical_info %>%
    filter(!is.na(ISS_Stage), !is.na(EFS_censor), !is.na(Age)) %>%
    left_join(fusions_primary, by = "mmrf") %>% 
    filter(!is.na(fusion)) %>%
    select(fusion) %>% 
    group_by(fusion) %>%
    summarize(count = n()) %>% 
    filter(count > 10) %>% 
    pull(fusion)
  
  genes_gt2 <- bind_rows(seqfish_clinical_info %>%
                           filter(!is.na(ISS_Stage), !is.na(EFS_censor), !is.na(Age)) %>%
                           left_join(fusions_primary, by = "mmrf") %>% 
                           filter(!is.na(fusion)) %>% 
                           select(geneA) %>% 
                           rename("gene" = "geneA"), 
                         seqfish_clinical_info %>%
                           filter(!is.na(ISS_Stage), !is.na(EFS_censor), !is.na(Age)) %>%
                           left_join(fusions_primary, by = "mmrf") %>% 
                           filter(!is.na(fusion)) %>% 
                           select(geneB) %>% 
                           rename("gene" = "geneB")) %>% 
    filter(!(gene %in% c("IGH", "IGK", "IGL", "IGHpseudo"))) %>%
    group_by(gene) %>% 
    summarize(count = n()) %>% 
    filter(count > 10) %>% 
    pull(gene)
  
  n_tests_fusions <- 1
  n_tests_genes <- 1
  
  for (this_fusion in fusions_gt2) {
    print(this_fusion)
    
    EFS_tibble <- seqfish_clinical_info %>%
      filter(!is.na(ISS_Stage), !is.na(EFS_censor), !is.na(Age)) %>%
      left_join(fusions_primary %>% filter(fusion == this_fusion), by = "mmrf") %>% 
      select(mmrf, Age, fusion, ISS_Stage, EFS, EFS_censor) %>%
      replace_na(list(fusion = "_None"))
    Death_tibble <- seqfish_clinical_info %>% 
      filter(!is.na(ISS_Stage), !is.na(D_PT_lstalive), !is.na(Age)) %>%
      rowwise() %>% 
      mutate( time_on_trial = max(D_PT_deathdy, D_PT_lstalive, na.rm = TRUE), 
              death = as.numeric(!is.na(D_PT_deathdy))) %>%
      left_join(fusions_primary %>% filter(fusion == this_fusion), by = "mmrf") %>% 
      select(mmrf, Age, fusion, ISS_Stage, time_on_trial, death) %>%
      replace_na(list(fusion = "_None"))
    
    if (EFS_tibble %>% pull(fusion) %>% unique() %>% length() > 1) {
      baseline_coxph_model_EFS <- coxph(Surv(EFS, EFS_censor == 0) ~ ISS_Stage + Age, data = EFS_tibble)
      coxph_model_EFS <- coxph(Surv(EFS, EFS_censor == 0) ~ ISS_Stage + Age + fusion, data = EFS_tibble)
      p_value_EFS <- summary(coxph_model_EFS)$coefficients[str_c("fusion", this_fusion), 5]
      p_value_model_comparison <- anova(baseline_coxph_model_EFS, coxph_model_EFS)["P(>|Chi|)"][2,1]
      
      if (p_value_EFS < 0.05/n_tests_fusions & p_value_model_comparison < 0.05/n_tests_fusions) {
        coxph_model_EFS_list[[this_fusion]] <- coxph_model_EFS
        
        fit <- survfit(Surv(EFS, EFS_censor == 0) ~ fusion, data = EFS_tibble)
        pdf(str_c(paper_main, this_fusion, ".EFS.with_legend.pdf"), 
            width = 4, height = 4)
        print(ggsurvplot(fit, data = EFS_tibble, conf.int = TRUE,
                         surv.median.line = "hv", pval = TRUE,
                         legend.labs = c("No Fusion", str_c(this_fusion, " Fusion")),
                         legend = "right",
                         xlab = "Time (days)", ylab = "Event-Free Survival Probability",
                         ggtheme = theme_survminer(),
                         conf.int.alpha = 0.1))
        dev.off()
        pdf(str_c(paper_main, this_fusion, ".EFS.without_legend.pdf"),
            width = 4, height = 4)
        print(ggsurvplot(fit, data = EFS_tibble, conf.int = TRUE,
                         surv.median.line = "hv", pval = TRUE,
                         legend.labs = c("No Fusion", str_c(this_fusion, " Fusion")),
                         legend = "none",
                         xlab = "Time (days)", ylab = "Event-Free Survival Probability",
                         ggtheme = theme_survminer(),
                         conf.int.alpha = 0.1))
        dev.off()
      }
    }
    
    if (Death_tibble %>% pull(fusion) %>% unique() %>% length() > 1) {
      baseline_coxph_model_Death <- coxph(Surv(time_on_trial, death == 1) ~ ISS_Stage + Age, data = Death_tibble)
      coxph_model_Death <- coxph(Surv(time_on_trial, death == 1) ~ ISS_Stage + Age + fusion, data = Death_tibble)
      p_value_Death <- summary(coxph_model_Death)$coefficients[str_c("fusion", this_fusion), 5]
      p_value_model_comparison <- anova(baseline_coxph_model_Death, coxph_model_Death)["P(>|Chi|)"][2,1]
      
      if (p_value_Death < 0.05/n_tests_fusions & p_value_model_comparison < 0.05/n_tests_fusions) {
        coxph_model_Death_list[[this_fusion]] <- coxph_model_Death
        
        fit <- survfit(Surv(time_on_trial, death == 1) ~ fusion, data = Death_tibble)
        pdf(str_c(paper_main, this_fusion, ".Death.with_legend.pdf"), 
            width = 4, height = 4)
        print(ggsurvplot(fit, data = Death_tibble, conf.int = TRUE,
                         surv.median.line = "hv", pval = TRUE,
                         legend.labs = c("No Fusion", str_c(this_fusion, " Fusion")),
                         legend = "right",
                         xlab = "Time (days)", ylab = "Overall Survival Probability",
                         ggtheme = theme_survminer(),
                         conf.int.alpha = 0.1))
        dev.off()
        pdf(str_c(paper_main, this_fusion, ".Death.without_legend.pdf"),
            width = 4, height = 4)
        print(ggsurvplot(fit, data = Death_tibble, conf.int = TRUE,
                         surv.median.line = "hv", pval = TRUE,
                         legend.labs = c("No Fusion", str_c(this_fusion, " Fusion")),
                         legend = "none",
                         xlab = "Time (days)", ylab = "Overall Survival Probability",
                         ggtheme = theme_survminer(),
                         conf.int.alpha = 0.1))
        dev.off()
      }
      }
    }
  
  for (this_gene in genes_gt2) {
    print(this_gene)
    
    EFS_tibble <- seqfish_clinical_info %>%
      filter(!is.na(ISS_Stage), !is.na(EFS_censor), !is.na(Age)) %>%
      left_join(fusions_primary %>% filter(geneA == this_gene | geneB == this_gene), by = "mmrf") %>% 
      select(mmrf, Age, fusion, ISS_Stage, EFS, EFS_censor) %>% 
      mutate(gene = case_when(is.na(fusion) ~ "_None", 
                              TRUE ~ this_gene))
    Death_tibble <- seqfish_clinical_info %>% 
      filter(!is.na(ISS_Stage), !is.na(D_PT_lstalive), !is.na(Age)) %>% 
      rowwise() %>% 
      mutate( time_on_trial = max(D_PT_deathdy, D_PT_lstalive, na.rm = TRUE), 
              death = as.numeric(!is.na(D_PT_deathdy))) %>%
      left_join(fusions_primary %>% filter(geneA == this_gene | geneB == this_gene), by = "mmrf") %>% 
      select(mmrf, Age, fusion, ISS_Stage, time_on_trial, death) %>% 
      mutate(gene = case_when(is.na(fusion) ~ "_None", 
                              TRUE ~ this_gene))
    
    if (EFS_tibble %>% pull(gene) %>% unique() %>% length() > 1) {
      baseline_coxph_model_EFS <- coxph(Surv(EFS, EFS_censor == 0) ~ ISS_Stage + Age, data = EFS_tibble)
      coxph_model_EFS <- coxph(Surv(EFS, EFS_censor == 0) ~ ISS_Stage + Age + gene, data = EFS_tibble)
      p_value_EFS <- summary(coxph_model_EFS)$coefficients[str_c("gene", this_gene), 5]
      p_value_model_comparison <- anova(baseline_coxph_model_EFS, coxph_model_EFS)["P(>|Chi|)"][2,1]
      
      if (p_value_EFS < 0.05/n_tests_fusions & p_value_model_comparison < 0.05/n_tests_fusions) {
        coxph_model_EFS_list[[this_gene]] <- coxph_model_EFS
        
        fit <- survfit(Surv(EFS, EFS_censor == 0) ~ gene, data = EFS_tibble)
        pdf(str_c(paper_main, this_gene, ".EFS.with_legend.pdf"), 
            width = 4, height = 4)
        print(ggsurvplot(fit, data = EFS_tibble, conf.int = TRUE,
                         surv.median.line = "hv", pval = TRUE,
                         legend.labs = c("No Fusion", str_c(this_gene, " Fusion")),
                         legend = "right",
                         xlab = "Time (days)", ylab = "Event-Free Survival Probability",
                         ggtheme = theme_survminer(),
                         conf.int.alpha = 0.1))
        dev.off()
        pdf(str_c(paper_main, this_gene, ".EFS.without_legend.pdf"),
            width = 4, height = 4)
        print(ggsurvplot(fit, data = EFS_tibble, conf.int = TRUE,
                         surv.median.line = "hv", pval = TRUE,
                         legend.labs = c("No Fusion", str_c(this_gene, " Fusion")),
                         legend = "none",
                         xlab = "Time (days)", ylab = "Event-Free Survival Probability",
                         ggtheme = theme_survminer(),
                         conf.int.alpha = 0.1))
        dev.off()
      }
    }
    
    if (Death_tibble %>% pull(gene) %>% unique() %>% length() > 1) {
      baseline_coxph_model_Death <- coxph(Surv(time_on_trial, death == 1) ~ ISS_Stage + Age, data = Death_tibble)
      coxph_model_Death <- coxph(Surv(time_on_trial, death == 1) ~ ISS_Stage + Age + gene, data = Death_tibble)
      p_value_Death <- summary(coxph_model_Death)$coefficients[str_c("gene", this_gene), 5]
      p_value_model_comparison <- anova(baseline_coxph_model_Death, coxph_model_Death)["P(>|Chi|)"][2,1]
      
      if (p_value_Death < 0.05/n_tests_fusions & p_value_model_comparison < 0.05/n_tests_fusions) {
        coxph_model_Death_list[[this_gene]] <- coxph_model_Death
        
        fit <- survfit(Surv(time_on_trial, death == 1) ~ gene, data = Death_tibble)
        pdf(str_c(paper_main, this_gene, ".Death.with_legend.pdf"), 
            width = 4, height = 4)
        print(ggsurvplot(fit, data = Death_tibble, conf.int = TRUE,
                         surv.median.line = "hv", pval = TRUE,
                         legend.labs = c("No Fusion", str_c(this_gene, " Fusion")),
                         legend = "right",
                         xlab = "Time (days)", ylab = "Overall Survival Probability",
                         ggtheme = theme_survminer(),
                         conf.int.alpha = 0.1))
        dev.off()
        pdf(str_c(paper_main, this_gene, ".Death.without_legend.pdf"),
            width = 4, height = 4)
        print(ggsurvplot(fit, data = Death_tibble, conf.int = TRUE,
                         surv.median.line = "hv", pval = TRUE,
                         legend.labs = c("No Fusion", str_c(this_gene, " Fusion")),
                         legend = "none",
                         xlab = "Time (days)", ylab = "Overall Survival Probability",
                         ggtheme = theme_survminer(),
                         conf.int.alpha = 0.1))
        dev.off()
      }
    }
  }
  
  print(coxph_model_EFS_list %>% names())
  print(coxph_model_Death_list %>% names())
}

################################################################################
# Look at other clinical associations
# Written April 2019
################################################################################

if (TRUE) {
  set.seed(10)

  testing_tbl_pvalue_adjusted <- read_tsv(
    "paper/supplemental/02_expression/testing_tbl_pvalue_adjusted.tsv")
  
  clinical_features <- testing_tbl_pvalue_adjusted %>% 
    filter(event_type == "Fusion Clinical", fdr < 0.05) %>% 
    pull(event2) %>% unique()
  
  genes_list <- list()
  for (clinical_feature in clinical_features) {
    genes_list[[clinical_feature]] <- testing_tbl_pvalue_adjusted %>% 
      filter(event_type == "Fusion Clinical", 
             event2 == clinical_feature, fdr < 0.05) %>% pull(event1)
  }
  
  plot_df_together <- NULL
  for (clinical_feature in clinical_features) {
    plot_df <- seqfish_clinical_info %>% 
      select(mmrf, clinical_feature) %>%
      filter(is.na(clinical_feature) | clinical_feature != "NA") %>%
      rename("clinical_interest" = clinical_feature) %>%
      left_join(fusions_primary %>% 
                  filter(geneA %in% genes_list[[clinical_feature]] | 
                           geneB %in% genes_list[[clinical_feature]]), by = "mmrf") %>%
      select(mmrf, clinical_interest, fusion, geneA, geneB) %>%
      mutate(feature = clinical_feature,
             relevant_gene = case_when(geneA %in% genes_list[[clinical_feature]] ~ geneA,
                                       geneB %in% genes_list[[clinical_feature]] ~ geneB),
             has_fusion = case_when(is.na(relevant_gene) ~ "No Fusion",
                                    TRUE ~ "Fusion"))
    
    plot_df_together <- bind_rows(plot_df_together, plot_df)
  }
  
  add_genes <- plot_df_together %>% 
    filter(has_fusion == "Fusion") %>% 
    group_by(feature) %>% 
    summarize(max_value = max(clinical_interest, na.rm = TRUE), 
              gene_names = str_c(sort(unique(relevant_gene)), collapse = ",\n"))

  
  genes_add <- plot_df_together %>% 
    filter(has_fusion == "Fusion") %>% 
    select(feature, relevant_gene) %>%
    unique() %>% arrange(relevant_gene) %>%
    group_by(feature) %>%
    mutate(gene_number = row_number()) %>%
    mutate(separator = case_when(gene_number %% 3 == 0 ~ ",\n",
                                 TRUE ~ ", ")) %>%
    summarize(gene_list = str_sub(str_c(relevant_gene, separator, collapse = ""), 1, -3))
  
  values_add <- plot_df_together %>% 
    filter(has_fusion == "Fusion") %>% 
    group_by(feature) %>% 
    summarize(max_value = max(clinical_interest, na.rm = TRUE))
  
  text_add <- genes_add %>% left_join(values_add, by = "feature")
  
  plot_df_together <- plot_df_together %>% mutate(feature = case_when(feature == "BM_Plasma_Cell_Percent" ~ "Bone Marrow Plasma Cell (%)",
                                                  feature == "LDH" ~ "LDH (units/L)"))
  text_add <- text_add %>% mutate(feature = case_when(feature == "BM_Plasma_Cell_Percent" ~ "Bone Marrow Plasma Cell (%)",
                                                      feature == "LDH" ~ "LDH (units/L)"))
  
  p <- ggplot(plot_df_together, aes(x = has_fusion, y = clinical_interest)) +
    geom_violin(fill = "black", alpha = 0.1, scale = "width", color = NA, show.legend = FALSE) +
    geom_jitter(aes(color = relevant_gene), height = 0, width = 0.1, shape = 16) +
    scale_color_brewer(palette = "Set3", na.value = "grey50", direction = -1) +
    geom_text(data = text_add, aes(x = "Fusion",
                                   y = max_value*1.1,
                                   label = gene_list),
              hjust = 0,
              fontface = "italic",
              show.legend = FALSE) +
    facet_wrap(~ feature, nrow = 1, scales = "free_x") + 
    coord_flip() +
    labs(x = NULL, y = NULL) +
    theme_bw() +
    theme(panel.background = element_blank(),
          panel.border = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = 12),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          #panel.grid.major.y = element_blank(),
          axis.ticks = element_blank())
  
    ggsave(str_c(paper_main, "clinical_associations.pdf"), p, width = 8, height = 2, useDingbats = FALSE)
    ggsave(str_c(paper_main, "clinical_associations.no_legend.pdf"), p + guides(color = FALSE), width = 8, height = 2, useDingbats = FALSE)
}

################################################################################
# Connection to TCGA cancer fusion calls
# Written April 2019
################################################################################

if (TRUE) {
  x <- fusions_primary %>% 
    filter(drug_geneA | drug_geneB) %>% 
    group_by(fusion) %>% 
    summarize(mmrf_count = n()) %>% ungroup() %>% 
    left_join(pancan_fusions, by = c("fusion" = "Fusion")) %>% 
    filter(!is.na(Cancer)) %>% 
    select(mmrf_count, fusion, Cancer, Sample) %>% 
    group_by(Cancer, fusion) %>% 
    summarize(count = n(), 
              mmrf_fusions = str_c(unique(sort(fusion)), collapse = ", "))
  ggplot(data = x, aes(x = fct_reorder(Cancer, desc(Cancer)), y = count)) +
    geom_bar(aes(fill = mmrf_fusions), stat = "identity") +
    geom_text(data = x %>% filter(count == 1), aes(label = mmrf_fusions),
              color = "black",
              hjust = 1) +
    geom_text(data = x %>% filter(count > 1), aes(label = mmrf_fusions),
              color = "white",
              hjust = 0) +
    scale_y_reverse() + 
    scale_x_discrete(position = "top") +
    coord_flip(expand = c(0, 0)) +
    labs(x = NULL, y = NULL) +
    theme_bw() +
    theme(panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          axis.ticks = element_blank()) +
    ggsave(str_c(paper_main, "TCGA_PanCan_Fusions.pdf"),
           width = 4, height = 2)
  
  x %>% ungroup() %>%
    select(Cancer, fusion, count) %>%
    rename("Fusion" = "fusion", "TCGA_Sample_Count" = "count") %>%
    write_tsv(str_c(paper_supp, "TCGA_overlap.tsv"))
}

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
    geom_col(fill = "#6e016b") +
    geom_text(data = drug_df %>% filter(n_samples < 10), color = "#000000", hjust = 0, nudge_y = 0.1, size = 3) +
    geom_text(data = drug_df %>% filter(n_samples > 10), color = "#ffffff", hjust = 1, nudge_y = -0.1, size = 3) +
    geom_text(aes(label = gene, y = 0), color = "#ffffff", hjust = 0, nudge_y = 0.025, fontface = "italic", size = 2.5) +
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
    ggsave(str_c(paper_main, "drug_targets.pdf"), width = 7.25, height = 7.25/1.618, useDingbats = FALSE)
  
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
      }
    }
    
    if (Death_tibble %>% pull(fusion) %>% unique() %>% length() > 1) {
      baseline_coxph_model_Death <- coxph(Surv(time_on_trial, death == 1) ~ ISS_Stage + Age, data = Death_tibble)
      coxph_model_Death <- coxph(Surv(time_on_trial, death == 1) ~ ISS_Stage + Age + fusion, data = Death_tibble)
      p_value_Death <- summary(coxph_model_Death)$coefficients[str_c("fusion", this_fusion), 5]
      p_value_model_comparison <- anova(baseline_coxph_model_Death, coxph_model_Death)["P(>|Chi|)"][2,1]
      
      if (p_value_Death < 0.05/n_tests_fusions & p_value_model_comparison < 0.05/n_tests_fusions) {
        coxph_model_Death_list[[this_fusion]] <- coxph_model_Death
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
      }
    }
    
    if (Death_tibble %>% pull(gene) %>% unique() %>% length() > 1) {
      baseline_coxph_model_Death <- coxph(Surv(time_on_trial, death == 1) ~ ISS_Stage + Age, data = Death_tibble)
      coxph_model_Death <- coxph(Surv(time_on_trial, death == 1) ~ ISS_Stage + Age + gene, data = Death_tibble)
      p_value_Death <- summary(coxph_model_Death)$coefficients[str_c("gene", this_gene), 5]
      p_value_model_comparison <- anova(baseline_coxph_model_Death, coxph_model_Death)["P(>|Chi|)"][2,1]
      
      if (p_value_Death < 0.05/n_tests_fusions & p_value_model_comparison < 0.05/n_tests_fusions) {
        coxph_model_Death_list[[this_gene]] <- coxph_model_Death
      }
    }
  }
  
  print(coxph_model_EFS_list %>% names())
  print(coxph_model_Death_list %>% names())
  
  # Make EFS survival figures for IGH--WHSC1 and PVT1--IGL/MYC--IGL
  plot_survival_list <- list()
  
  #IGH--WHSC1
  EFS_tibble <- seqfish_clinical_info %>%
    filter(!is.na(ISS_Stage), !is.na(EFS_censor), !is.na(Age)) %>%
    left_join(fusions_primary %>% filter(fusion == "IGH--WHSC1"), by = "mmrf") %>% 
    select(mmrf, Age, fusion, ISS_Stage, EFS, EFS_censor) %>%
    replace_na(list(fusion = "_None"))
  
  plot_survival_list[["base_EFS"]] <- coxph(formula = Surv(EFS, EFS_censor == 0) ~ ISS_Stage + Age, data = EFS_tibble)
  plot_survival_list[["WHSC1_fusion_EFS"]] <- coxph(formula = Surv(EFS, EFS_censor == 0) ~ ISS_Stage + Age + fusion, data = EFS_tibble)
  plot_survival_list[["anova_WHSC1_fusion_EFS"]] <- anova(plot_survival_list[["base_EFS"]], plot_survival_list[["WHSC1_fusion_EFS"]]) # Significant
  
  Death_tibble <- seqfish_clinical_info %>% 
    filter(!is.na(ISS_Stage), !is.na(D_PT_lstalive), !is.na(Age)) %>%
    rowwise() %>% 
    mutate( time_on_trial = max(D_PT_deathdy, D_PT_lstalive, na.rm = TRUE), 
            death = as.numeric(!is.na(D_PT_deathdy))) %>%
    left_join(fusions_primary %>% filter(fusion == "IGH--WHSC1"), by = "mmrf") %>% 
    select(mmrf, Age, fusion, ISS_Stage, time_on_trial, death) %>%
    replace_na(list(fusion = "_None"))
  
  plot_survival_list[["base_Death"]] <- coxph(formula = Surv(time_on_trial, death == 1) ~ ISS_Stage + Age, data = Death_tibble)
  plot_survival_list[["WHSC1_fusion_Death"]] <- coxph(formula = Surv(time_on_trial, death == 1) ~ ISS_Stage + Age + fusion, data = Death_tibble)
  plot_survival_list[["anova_WHSC1_fusion_Death"]] <- anova(plot_survival_list[["base_Death"]], plot_survival_list[["WHSC1_fusion_Death"]]) # Significant
  
  fit <- survfit(Surv(EFS, EFS_censor == 0) ~ fusion, data = EFS_tibble)
  pdf(str_c(paper_supp, "WHSC1.EFS.with_legend.pdf"),
      width = 3.5, height = 3.5, useDingbats = FALSE)
  print(ggsurvplot(fit, data = EFS_tibble, conf.int = TRUE,
                   surv.median.line = "hv", pval = TRUE,
                   legend.labs = c("No Fusion", str_c("IGH--WHSC1", " Fusion")),
                   legend = "right", 
                   xlab = "Time (days)", 
                   ylab = "Event-Free Survival Probability",
                   ggtheme = theme_survminer(),
                   conf.int.alpha = 0.1))
  dev.off()
  pdf(str_c(paper_supp, "WHSC1.EFS.without_legend.pdf"),
      width = 3.5, height = 3.5, useDingbats = FALSE)
  print(ggsurvplot(fit, data = EFS_tibble, conf.int = TRUE,
                   surv.median.line = "hv", pval = TRUE,
                   legend.labs = c("No Fusion", str_c("IGH--WHSC1", " Fusion")),
                   legend = "none",
                   xlab = "Time (days)", 
                   ylab = "Event-Free Survival Probability",
                   ggtheme = theme_survminer(),
                   conf.int.alpha = 0.1))
  dev.off()
  
  fit <- survfit(Surv(time_on_trial, death == 1) ~ fusion, data = Death_tibble)
  pdf(str_c(paper_supp, "WHSC1.Death.with_legend.pdf"),
      width = 3.5, height = 3.5, useDingbats = FALSE)
  print(ggsurvplot(fit, data = Death_tibble, conf.int = TRUE,
                   surv.median.line = "hv", pval = TRUE,
                   legend.labs = c("No Fusion", str_c("IGH--WHSC1", " Fusion")),
                   legend = "right", 
                   xlab = "Time (days)", 
                   ylab = "Overall Survival Probability",
                   ggtheme = theme_survminer(),
                   conf.int.alpha = 0.1))
  dev.off()
  pdf(str_c(paper_supp, "WHSC1.Death.without_legend.pdf"),
      width = 3.5, height = 3.5, useDingbats = FALSE)
  print(ggsurvplot(fit, data = Death_tibble, conf.int = TRUE,
                   surv.median.line = "hv", pval = TRUE,
                   legend.labs = c("No Fusion", str_c("IGH--WHSC1", " Fusion")),
                   legend = "none",
                   xlab = "Time (days)", 
                   ylab = "Overall Survival Probability",
                   ggtheme = theme_survminer(),
                   conf.int.alpha = 0.1))
  dev.off()
  
  #PVT1--IGL MYC--IGL
  EFS_tibble <- seqfish_clinical_info %>%
    filter(!is.na(ISS_Stage), !is.na(EFS_censor), !is.na(Age)) %>%
    left_join(fusions_primary %>% filter(fusion %in% c("PVT1--IGL", "MYC--IGL")), by = "mmrf") %>% 
    select(mmrf, Age, fusion, ISS_Stage, EFS, EFS_censor) %>%
    replace_na(list(fusion = "_None"))
  
  plot_survival_list[["PVT1_MYC_fusion_EFS"]] <- fusion_model <- coxph(formula = Surv(EFS, EFS_censor == 0) ~ ISS_Stage + Age + fusion, data = EFS_tibble)
  plot_survival_list[["anova_PVT1_MYC_fusion_EFS"]] <- anova(plot_survival_list[["base_EFS"]], plot_survival_list[["PVT1_MYC_fusion_EFS"]]) # Significant
  
  Death_tibble <- seqfish_clinical_info %>% 
    filter(!is.na(ISS_Stage), !is.na(D_PT_lstalive), !is.na(Age)) %>%
    rowwise() %>% 
    mutate( time_on_trial = max(D_PT_deathdy, D_PT_lstalive, na.rm = TRUE), 
            death = as.numeric(!is.na(D_PT_deathdy))) %>%
    left_join(fusions_primary %>% filter(fusion %in% c("PVT1--IGL", "MYC--IGL")), by = "mmrf") %>% 
    select(mmrf, Age, fusion, ISS_Stage, time_on_trial, death) %>%
    replace_na(list(fusion = "_None"))
  
  plot_survival_list[["PVT1_MYC_fusion_Death"]] <- fusion_model <- coxph(formula = Surv(time_on_trial, death == 1) ~ ISS_Stage + Age + fusion, data = Death_tibble)
  plot_survival_list[["anova_PVT1_MYC_fusion_Death"]] <- anova(plot_survival_list[["base_Death"]], plot_survival_list[["PVT1_MYC_fusion_Death"]]) # NOT Significant
  
  fit <- survfit(Surv(EFS, EFS_censor == 0) ~ fusion, data = EFS_tibble)
  pdf(str_c(paper_supp, "PVT1_MYC.EFS.with_legend.pdf"),
      width = 3.5, height = 3.5, useDingbats = FALSE)
  print(ggsurvplot(fit, data = EFS_tibble, conf.int = TRUE,
                   surv.median.line = "hv", pval = TRUE,
                   legend.labs = c("No Fusion", str_c(c("MYC--IGL", "PVT1--IGL"), " Fusion")),
                   legend = "right", 
                   xlab = "Time (days)", 
                   ylab = "Event-Free Survival Probability",
                   ggtheme = theme_survminer(),
                   conf.int.alpha = 0.1))
  dev.off()
  pdf(str_c(paper_supp, "PVT1_MYC.EFS.without_legend.pdf"),
      width = 3.5, height = 3.5, useDingbats = FALSE)
  print(ggsurvplot(fit, data = EFS_tibble, conf.int = TRUE,
                   surv.median.line = "hv", pval = TRUE,
                   legend.labs = c("No Fusion", str_c(c("MYC--IGL", "PVT1--IGL"), " Fusion")),
                   legend = "none",
                   xlab = "Time (days)", 
                   ylab = "Event-Free Survival Probability",
                   ggtheme = theme_survminer(),
                   conf.int.alpha = 0.1))
  dev.off()
  
  fit <- survfit(Surv(time_on_trial, death == 1) ~ fusion, data = Death_tibble)
  pdf(str_c(paper_supp, "PVT1_MYC.Death.with_legend.pdf"),
      width = 3.5, height = 3.5, useDingbats = FALSE)
  print(ggsurvplot(fit, data = Death_tibble, conf.int = TRUE,
                   surv.median.line = "hv", pval = TRUE,
                   legend.labs = c("No Fusion", str_c(c("MYC--IGL", "PVT1--IGL"), " Fusion")),
                   legend = "right", 
                   xlab = "Time (days)", 
                   ylab = "Overall Survival Probability",
                   ggtheme = theme_survminer(),
                   conf.int.alpha = 0.1))
  dev.off()
  pdf(str_c(paper_supp, "PVT1_MYC.Death.without_legend.pdf"),
      width = 3.5, height = 3.5, useDingbats = FALSE)
  print(ggsurvplot(fit, data = Death_tibble, conf.int = TRUE,
                   surv.median.line = "hv", pval = TRUE,
                   legend.labs = c("No Fusion", str_c(c("MYC--IGL", "PVT1--IGL"), " Fusion")),
                   legend = "none",
                   xlab = "Time (days)", 
                   ylab = "Overall Survival Probability",
                   ggtheme = theme_survminer(),
                   conf.int.alpha = 0.1))
  dev.off()
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
  
    ggsave(str_c(paper_supp, "clinical_associations.pdf"), p, width = 7.25, height = 7.25/(2*1.618), useDingbats = FALSE)
    ggsave(str_c(paper_supp, "clinical_associations.no_legend.pdf"), p + guides(color = FALSE), width = 7.25, height = 7.25/(2*1.618), useDingbats = FALSE)
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
              mmrf_fusions = str_c(unique(sort(fusion)), collapse = ", ")) %>%
    ungroup() %>%
    select(Cancer, fusion, count) %>%
    rename("Fusion" = "fusion", "TCGA_Sample_Count" = "count") %>%
    write_tsv(str_c(paper_supp, "TCGA_overlap.tsv"))
  
  p <- ggplot(x, aes(x = factor(Fusion, levels = c("SND1--BRAF", "TPM3--NTRK1", "NOTCH2--SEC22B"), ordered = TRUE), 
                y = TCGA_Sample_Count, 
                fill = Cancer)) + 
    geom_col(alpha = 0.75) + 
    scale_fill_brewer(palette = "BuPu") +
    #scale_fill_manual(values = c("#FAD2D9","#A084BD","#7E1918","#00A99D","#F9ED32","#FBE3C7")) +
    theme_bw() +
    labs(x = NULL, y = "TCGA Sample Count") +
    scale_y_continuous(position = "right", expand = c(0.01,0.01)) +
    scale_x_discrete(position = "top", expand = c(0.01,0.01)) +
    theme(panel.background = element_blank(),
          panel.grid = element_blank(), 
          panel.border = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_text(size = 10),
          axis.text.x = element_text(size = 10),
          axis.title = element_text(size = 12),
          legend.position = "bottom",
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10))
  
  ggsave(str_c(paper_main, "TCGA_connection.pdf"), p,
               width = 10, height = 5)
  ggsave(str_c(paper_main, "TCGA_connection.no_legend.pdf"), p + guides(fill = FALSE),
               width = 3, height = 2)
}

# ==============================================================================
# Types of kinases
# Written April 2019 -- Supp
# ==============================================================================

if (TRUE) {
  p <- kinases %>%
    group_by(kinase_group_full_name) %>% 
    summarize(count = n()) %>% 
    ungroup() %>%
    ggplot(aes(x = fct_reorder(kinase_group_full_name, count), y = count)) +
    geom_col(position = "dodge") +
    coord_flip(expand = c(0,0)) +
    labs(y = "Fusion Count", x = "Kinase Group") +
    scale_y_continuous(position = "right") +
    theme_bw() +
    theme(panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "bottom",
          legend.direction = "vertical",
          axis.text.x = element_text(size = 8),
          axis.text.y = element_text(size = 10),
          axis.title = element_text(size = 12)
    )
  
  ggsave(str_c(paper_supp, "kinase_groups.with_legend.pdf"), p, 
         width = 7.25, height = 7.25/1.618, useDingbats = FALSE)
  ggsave(str_c(paper_supp, "kinase_groups.without_legend.pdf"), 
         p + guides(fill = FALSE), 
         width = 7.25, height = 7.25/1.618, useDingbats = FALSE)
  
}

# ==============================================================================
# NTRK1 fusion structures -- manually create based on agFusion output
# Written April 2019 -- Main
# ==============================================================================

if (TRUE) {
  structure_tbl <- tribble(~mmrf,          ~fusion,          ~element, ~fill_color, ~exterior_color, ~start, ~stop,   ~class,
                            1232,     "TPR--NTRK1",             "TPR",           1,               1,      0,   366,   "gene",
                            1232,     "TPR--NTRK1",           "NTRK1",           2,               1,    366,   764,   "gene",
                            1232,     "TPR--NTRK1", "Tyrosine Kinase",           4,               2,    480,   749, "domain",
                            1656,    "TPM3--NTRK1",            "TPM3",           1,               1,      0,   258,   "gene",
                            1656,    "TPM3--NTRK1",           "NTRK1",           2,               1,    258,   656,   "gene",
                            1656,    "TPM3--NTRK1",     "Tropomyosin",           3,               2,     49,   258, "domain",
                            1656,    "TPM3--NTRK1", "Tyrosine Kinase",           4,               2,    372,   641, "domain",
                            2490, "ARHGEF2--NTRK1",         "ARHGEF2",           1,               1,      0,   962,   "gene",
                            2490, "ARHGEF2--NTRK1",           "NTRK1",           2,               1,    962,  1307,   "gene",
                            2490, "ARHGEF2--NTRK1",          "RhoGEF",           3,               2,    239,	  431, "domain",
                            2490, "ARHGEF2--NTRK1",              "PH",           3,               2,    474,   570, "domain",
                            2490, "ARHGEF2--NTRK1", "Tyrosine Kinase",           4,               2,   1023,	 1292, "domain") %>%
    mutate(fusion = factor(fusion, levels = c("TPM3--NTRK1", "TPR--NTRK1", "ARHGEF2--NTRK1")),
           fill_color = factor(fill_color),
           exterior_color = factor(exterior_color),
           class = factor(class))
  
  
  ggplot(structure_tbl, aes(xmin = start, 
                            xmax = stop, 
                            ymin = as.numeric(fusion), 
                            ymax = as.numeric(fusion) + 0.25, 
                            fill = fill_color, 
                            label = element)) +
    geom_rect(show.legend = FALSE) +
    scale_fill_manual(values = c("#bdd7e7", "#fcae91", "#2171b5", "#cb181d")) +
    scale_color_manual(values = c("#bdd7e7", "#fcae91", "#2171b5", "#cb181d")) +
    geom_segment(data = structure_tbl %>% filter(element == "NTRK1"), aes(x = start, xend = start, y = as.numeric(fusion) - 0.05, yend = as.numeric(fusion) + 0.25 + 0.025)) +
    geom_text(data = structure_tbl %>% filter(class == "gene"), aes(x = start + (stop - start)/2, y = as.numeric(fusion) - 0.075), color = "#000000", size = 4, show.legend = FALSE, fontface = "italic") +
    geom_text(data = structure_tbl %>% filter(class == "domain"), aes(x = start + (stop - start)/2, y = as.numeric(fusion) + 0.125), size = 4, show.legend = FALSE, color = "white") +
    geom_text(data = structure_tbl %>% filter(class == "gene", element != "NTRK1"), aes(x = start, y = as.numeric(fusion) + 0.125), label = "5'", hjust = 1, nudge_x = -10) +
    geom_text(data = structure_tbl %>% filter(class == "gene", element == "NTRK1"), aes(x = stop, y = as.numeric(fusion) + 0.125), label = "3'", hjust = 0, nudge_x = 10) +
    facet_wrap(~ fusion, strip.position = "left", scales = "free_y", ncol = 1) +
    labs(x = "Amino Acid Position", y = NULL) +
    scale_y_continuous(expand = c(0.05, 0.05)) +
    scale_x_continuous(breaks = seq(0, 1300, 100), position = "bottom", expand = c(0.01, 0.01)) +
    theme_bw() +
    theme(panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank(),
          axis.text.x = element_text(size = 8),
          axis.title.x = element_text(size = 12),
          axis.text.y = element_blank(),
          axis.ticks = element_blank()) + 
    ggsave(str_c(paper_main, "NTRK1.structures.pdf"), 
           width = 7.25, height = 7.25/(2*1.618), useDingbats = FALSE)
}

# ==============================================================================
# 3' intact kinase TCGA overlap
# Written April 2018 Supplement
# ==============================================================================

if (TRUE) {
  plot_df <- kinases %>% 
    filter(KinasePos == "3P_KINASE", KinaseDomain == "Intact") %>% 
    group_by(geneB) %>% 
    summarize(mmrf_count = n()) %>%
    mutate(geneB_mmrf_count = str_c(geneB, " (", mmrf_count, ")")) %>%
    group_by(geneB, geneB_mmrf_count) %>%
    left_join(pancan_fusions %>% 
                separate(Fusion, into = c("geneA", "geneB"), sep = "--"), 
              by = "geneB") %>% 
    group_by(geneB_mmrf_count, Cancer) %>% 
    summarize(count2 = n()) %>%
    filter(!is.na(Cancer))
  
  ggplot(data = plot_df, aes(y = geneB_mmrf_count, x = Cancer)) + 
    geom_tile(aes(fill = factor(count2))) + 
    geom_text(data = plot_df %>% filter(count2 < 6), aes(label = count2),
              color = "#000000") +
    geom_text(data = plot_df %>% filter(count2 >= 6), aes(label = count2),
              color = "#ffffff") +
    scale_fill_brewer(palette = "Blues") +
    labs(x = "TCGA Cancer Type",
         y = "MMRF 3' Intact Kinase Fusions",
         fill = "TCGA Samples\nwith Same 3' Fusion") +
    theme_bw() +
    theme(panel.background = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "bottom",
          axis.text.y = element_text(size = 10, face = "italic"),
          axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
          axis.title = element_text(size = 12)) +
    ggsave(str_c(paper_supp, "TCGA_kinase_fusions.pdf"),
           width = 7.25, height = 7)
  
}

# ==============================================================================
# 5' Partners of NTRK1 in TCGA
# Written April 2018
# ==============================================================================

if (TRUE) {
  pancan_fusions %>% 
    filter(str_detect(Fusion, pattern = "--NTRK1")) %>%
    separate(Fusion, by = "--", into = c("geneA", "geneB")) %>% 
    group_by(Cancer, geneA) %>% 
    summarize(count = n()) %>% 
    ggplot(aes(x = Cancer, y = geneA, fill = factor(count))) + 
    geom_tile() + 
    geom_text(aes(label = count)) +
    scale_fill_brewer(palette = "Blues") +
    labs(x = "TCGA Cancer Type",
         y = "5' Partner of NTRK1",
         fill = "TCGA Sample Count") +
    theme_bw() +
    theme(panel.background = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "bottom",
          axis.text.y = element_text(size = 10, face = "italic"),
          axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
          axis.title = element_text(size = 12)) +
    ggsave(str_c(paper_supp, "TCGA_NTRK1_partners.pdf"),
           width = 7.25, height = 3)
}


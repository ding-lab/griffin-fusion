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
    left_join(depo, by = c("geneA" = "Gene")) %>% filter(Variant.Type == "Fusion") %>%
    rename("gene" = "geneA")
  
  geneB_drugs <- fusions_primary %>% filter(drug_geneB == 1) %>% select(geneB) %>% unique() %>%
    left_join(depo, by = c("geneB" = "Gene")) %>% filter(Variant.Type == "Fusion") %>%
    rename("gene" = "geneB")
  
  drug_info <- bind_rows(geneA_drugs, geneB_drugs) %>% filter(Variant %in% c("any", "any ")) %>% 
    mutate(drug_pmid = str_c(Drug_Class, ": ", PubMed.ID)) %>%
    group_by(gene, Effect) %>%
    summarize(evidence = str_c(sort(unique(as.vector(str_split(str_c(Drug_Class, collapse = ","), pattern = ",", simplify = TRUE)))), collapse = ", "))
  
  drug_df <- bind_rows(fusions_primary %>% filter(drug_geneA == 1) %>% select(fusion, geneA, srr) %>% rename("gene" = "geneA"),
                       fusions_primary %>% filter(drug_geneB == 1) %>% select(fusion, geneB, srr) %>% rename("gene" = "geneB")) %>%
    right_join(drug_info, by = "gene") %>% group_by(gene, Effect, evidence) %>%
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
}


# ==============================================================================
# Survival association with fusion events
# Steven Foltz (smfoltz@wustl.edu), April 2019
# ==============================================================================

if (TRUE) {
  
  coxph_model_list <- list()
  
  fusions_gt2 <- seqfish_clinical_info %>%
    filter(!is.na(ISS_Stage), !is.na(EFS_censor), !is.na(Age),
           !is.na(seqfish_CN_del_17p13),
           !is.na(seqfish_Translocation_WHSC1_4_14)) %>% 
    left_join(fusions_primary, by = "mmrf") %>% 
    filter(!is.na(fusion)) %>%
    select(fusion) %>% 
    group_by(fusion) %>%
    summarize(count = n()) %>% 
    filter(count > 5) %>% 
    pull(fusion)
  
  genes_gt2 <- bind_rows(seqfish_clinical_info %>%
                           filter(!is.na(ISS_Stage), !is.na(EFS_censor), !is.na(Age),
                                  !is.na(seqfish_CN_del_17p13),
                                  !is.na(seqfish_Translocation_WHSC1_4_14)) %>% 
                           left_join(fusions_primary, by = "mmrf") %>% 
                           filter(!is.na(fusion)) %>% 
                           select(geneA) %>% 
                           rename("gene" = "geneA"), 
                         seqfish_clinical_info %>%
                           filter(!is.na(ISS_Stage), !is.na(EFS_censor), !is.na(Age),
                                  !is.na(seqfish_CN_del_17p13),
                                  !is.na(seqfish_Translocation_WHSC1_4_14)) %>% 
                           left_join(fusions_primary, by = "mmrf") %>% 
                           filter(!is.na(fusion)) %>% 
                           select(geneB) %>% 
                           rename("gene" = "geneB")) %>% 
    filter(!(gene %in% c("IGH", "IGK", "IGL"))) %>%
    group_by(gene) %>% 
    summarize(count = n()) %>% 
    filter(count > 5) %>% 
    pull(gene)
  
  n_tests <- length(genes_gt2) + length(fusions_gt2)
  
  for (this_fusion in fusions_gt2) {
    print(this_fusion)
    EFS_tibble <- seqfish_clinical_info %>%
      filter(!is.na(ISS_Stage), !is.na(EFS_censor), !is.na(Age),
             !is.na(seqfish_CN_del_17p13),
             !is.na(seqfish_Translocation_WHSC1_4_14)) %>% 
      left_join(fusions_primary %>% filter(fusion == this_fusion), by = "mmrf") %>% 
      select(mmrf, Age, fusion, ISS_Stage, EFS, EFS_censor,
             seqfish_CN_del_17p13, seqfish_Translocation_WHSC1_4_14) %>% 
      replace_na(list(fusion = "_None"))
    coxph_model <- coxph(Surv(EFS, EFS_censor == 0) ~ ISS_Stage + Age + seqfish_CN_del_17p13 + seqfish_Translocation_WHSC1_4_14 + fusion, data = EFS_tibble)
    #survfit_model <- survfit(Surv(EFS, EFS_censor == 0) ~ ISS_Stage + fusion + Age, data = EFS_tibble)
    p_value <- summary(coxph_model)$coefficients[str_c("fusion", this_fusion), 5]
    #if (p_value < 0.05/n_tests) {
    if (p_value < 0.01) {
      coxph_model_list[[this_fusion]] <- coxph_model
    }
  }
  
  for (this_gene in genes_gt2) {
    print(this_gene)
    EFS_tibble <- seqfish_clinical_info %>%
      filter(!is.na(ISS_Stage), !is.na(EFS_censor), !is.na(Age),
             !is.na(seqfish_CN_del_17p13),
             !is.na(seqfish_Translocation_WHSC1_4_14)) %>% 
      left_join(fusions_primary %>% filter(geneA == this_gene | geneB == this_gene), by = "mmrf") %>% 
      select(mmrf, Age, fusion, ISS_Stage, EFS, EFS_censor,
             seqfish_CN_del_17p13, seqfish_Translocation_WHSC1_4_14) %>%
      mutate(gene = case_when(is.na(fusion) ~ "_None", 
                              TRUE ~ this_gene))
    coxph_model <- coxph(Surv(EFS, EFS_censor == 0) ~ ISS_Stage + Age + seqfish_CN_del_17p13 + seqfish_Translocation_WHSC1_4_14 + gene, data = EFS_tibble)
    p_value <- summary(coxph_model)$coefficients[str_c("gene", this_gene), 5]
    #if (p_value < 0.05/n_tests) {
    if (p_value < 0.01) {
      coxph_model_list[[this_gene]] <- coxph_model
      EFS_fit <- survfit(Surv(EFS, EFS_censor == 0) ~ gene, data = EFS_tibble)
      pdf(str_c("~/Desktop/", this_gene, ".pdf"))
      print(ggsurvplot(EFS_fit, data = EFS_tibble, conf.int = TRUE,
                       surv.median.line = "hv", pval = TRUE,
                       #legend.labs = c("ISS Stage I", "ISS Stage II", "ISS Stage III"),
                       xlab = "Time (days)",
                       ggtheme = theme_bw(base_size = 20)))
      dev.off()
    }
    
  }
  
  print(coxph_model_list)
}
  
  
  
  
  
  
  
  
  # ============================================================================
  # Create EFS survival object
  # ============================================================================
  
  # Use EFS_censor == 0 because that is TRUE for death, FALSE for censored
  # Stratify by Stage
  EFS_tibble <- seqfish_clinical_info %>% 
    filter(!is.na(ISS_Stage), !is.na(EFS_censor)) %>% left_join(fusions_primary %>% filter(fusion == "IGH--WHSC1"), by = "mmrf") %>% 
    select(mmrf, Age, fusion, ISS_Stage, EFS, EFS_censor) %>% replace_na(list(fusion = "_None"))
  EFS_fit <- survfit(Surv(EFS, EFS_censor == 0) ~ ISS_Stage + fusion, 
                     data = EFS_tibble)
  
  # Plot survival curve stratified by Stage
  pdf(str_c(paper_supp, "event_free_survival.pdf"), width = 20, height = 15,
      useDingbats = FALSE)
  print(ggsurvplot(EFS_fit, data = EFS_tibble,  conf.int = TRUE,
                   surv.median.line = "hv", pval = TRUE,
                   #legend.labs = c("ISS Stage I", "ISS Stage II", "ISS Stage III"),
                   xlab = "Time (days)",
                   ggtheme = theme_bw(base_size = 20)))
  dev.off()
  
  # ============================================================================
  # Create Death survival object
  # ============================================================================
  
  # Use EFS_censor == 0 because that is TRUE for death, FALSE for censored
  # Stratify by Stage
  death_tibble <- seqfish_clinical_info %>% 
    filter(!is.na(ISS_Stage), !is.na(D_PT_lstalive)) %>% rowwise() %>% 
    mutate( time_on_trial = max(D_PT_deathdy, D_PT_lstalive, na.rm = TRUE), 
            death = as.numeric(!is.na(D_PT_deathdy)))
  death_fit <- survfit(Surv(time_on_trial, death) ~ ISS_Stage, 
                       data = death_tibble)
  
  # Plot survival curve stratified by Stage
  pdf(str_c(paper_supp, "overall_survival.pdf"), width = 20, height = 15)
  print(ggsurvplot(death_fit, data = death_tibble,  conf.int = TRUE,
                   surv.median.line = "hv", pval = TRUE, 
                   legend.labs = c("ISS Stage I", "ISS Stage II", "ISS Stage III"),
                   xlab = "Time (days)",
                   ggtheme = theme_bw(base_size = 20)))
  dev.off()
  
}

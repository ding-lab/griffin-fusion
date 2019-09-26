# ==============================================================================
# Clinical (MMRF Fusions)
# Steven Foltz (github: envest)
# ==============================================================================

paper_main = "paper/main/07_clinical/"
paper_supp = "paper/supplementary/07_clinical/"

# Create directories 
dir.create(paper_main, recursive = TRUE, showWarnings = FALSE)
dir.create(paper_supp, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# Druggable fusions table
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
    coord_flip(expand = FALSE) +
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
# ==============================================================================

if (TRUE) {
  
  coxph_model_EFS_list <- list()
  
  fusions_ge10 <- seqfish_clinical_info %>%
    filter(!is.na(ISS_Stage), !is.na(EFS_censor), !is.na(Age)) %>%
    left_join(fusions_primary, by = "mmrf") %>% 
    filter(!is.na(fusion)) %>%
    select(fusion) %>% 
    group_by(fusion) %>%
    summarize(count = n()) %>%
    filter(count >= 10) %>% 
    pull(fusion)
  
  genes_ge10 <- bind_rows(seqfish_clinical_info %>%
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
    filter(count >= 10) %>% 
    pull(gene)
  
  for (this_fusion in fusions_ge10) {
    EFS_tibble <- seqfish_clinical_info %>%
      filter(mmrf %in% mmrf_primary_pretreatment) %>%
      filter(!is.na(ISS_Stage), !is.na(EFS_censor), !is.na(Age)) %>%
      left_join(fusions_primary %>% filter(fusion == this_fusion), by = "mmrf") %>% 
      select(mmrf, Age, fusion, ISS_Stage, EFS, EFS_censor) %>%
      replace_na(list(fusion = "_None"))
    
    if (EFS_tibble %>% pull(fusion) %>% unique() %>% length() > 1) {
      baseline_coxph_model_EFS <- coxph(Surv(EFS, EFS_censor == 0) ~ ISS_Stage + Age, data = EFS_tibble)
      coxph_model_EFS <- coxph(Surv(EFS, EFS_censor == 0) ~ ISS_Stage + Age + fusion, data = EFS_tibble)
      p_value_EFS <- summary(coxph_model_EFS)$coefficients[str_c("fusion", this_fusion), 5]
      p_value_model_comparison <- anova(baseline_coxph_model_EFS, coxph_model_EFS)["P(>|Chi|)"][2,1]
      
      if (p_value_EFS < 0.05 & p_value_model_comparison < 0.05) {
        coxph_model_EFS_list[[this_fusion]] <- coxph_model_EFS
      }
    }
  }
  
  for (this_gene in genes_ge10) {
    EFS_tibble <- seqfish_clinical_info %>%
      filter(mmrf %in% mmrf_primary_pretreatment) %>%
      filter(!is.na(ISS_Stage), !is.na(EFS_censor), !is.na(Age)) %>%
      left_join(fusions_primary %>% filter(geneA == this_gene | geneB == this_gene), by = "mmrf") %>% 
      select(mmrf, Age, fusion, ISS_Stage, EFS, EFS_censor) %>% 
      mutate(gene = case_when(is.na(fusion) ~ "_None", 
                              TRUE ~ this_gene))
    
    if (EFS_tibble %>% pull(gene) %>% unique() %>% length() > 1) {
      baseline_coxph_model_EFS <- coxph(Surv(EFS, EFS_censor == 0) ~ ISS_Stage + Age, data = EFS_tibble)
      coxph_model_EFS <- coxph(Surv(EFS, EFS_censor == 0) ~ ISS_Stage + Age + gene, data = EFS_tibble)
      p_value_EFS <- summary(coxph_model_EFS)$coefficients[str_c("gene", this_gene), 5]
      p_value_model_comparison <- anova(baseline_coxph_model_EFS, coxph_model_EFS)["P(>|Chi|)"][2,1]
      
      if (p_value_EFS < 0.05 & p_value_model_comparison < 0.05) {
        coxph_model_EFS_list[[this_gene]] <- coxph_model_EFS
      }
    }
  }
  
  #### early relapse
  coxph_model_early_list <- list()
  
  for (this_fusion in fusions_ge10) {
    early_tibble <- seqfish_clinical_info %>%
      filter(mmrf %in% mmrf_primary_pretreatment) %>%
      filter(!is.na(ISS_Stage), !is.na(early_relapse_censor), !is.na(Age)) %>%
      left_join(fusions_primary %>% filter(fusion == this_fusion), by = "mmrf") %>% 
      select(mmrf, Age, fusion, ISS_Stage, early_relapse_time, early_relapse_censor) %>%
      replace_na(list(fusion = "_None"))
    
    if (early_tibble %>% pull(fusion) %>% unique() %>% length() > 1) {
      baseline_coxph_model_early <- coxph(Surv(early_relapse_time, early_relapse_censor == 0) ~ ISS_Stage + Age, data = early_tibble)
      coxph_model_early <- coxph(Surv(early_relapse_time, early_relapse_censor == 0) ~ ISS_Stage + Age + fusion, data = early_tibble)
      p_value_early <- summary(coxph_model_early)$coefficients[str_c("fusion", this_fusion), 5]
      p_value_model_comparison <- anova(baseline_coxph_model_early, coxph_model_early)["P(>|Chi|)"][2,1]
      
      if (p_value_early < 0.05 & p_value_model_comparison < 0.05) {
        coxph_model_early_list[[this_fusion]] <- coxph_model_early
      }
    }
  }
  
  for (this_gene in genes_ge10) {
    early_tibble <- seqfish_clinical_info %>%
      filter(mmrf %in% mmrf_primary_pretreatment) %>%
      filter(!is.na(ISS_Stage), !is.na(early_relapse_censor), !is.na(Age)) %>%
      left_join(fusions_primary %>% filter(geneA == this_gene | geneB == this_gene), by = "mmrf") %>% 
      select(mmrf, Age, fusion, ISS_Stage, early_relapse_time, early_relapse_censor) %>% 
      mutate(gene = case_when(is.na(fusion) ~ "_None", 
                              TRUE ~ this_gene))
    
    if (early_tibble %>% pull(gene) %>% unique() %>% length() > 1) {
      baseline_coxph_model_early <- coxph(Surv(early_relapse_time, early_relapse_censor == 0) ~ ISS_Stage + Age, data = early_tibble)
      coxph_model_early <- coxph(Surv(early_relapse_time, early_relapse_censor == 0) ~ ISS_Stage + Age + gene, data = early_tibble)
      p_value_early <- summary(coxph_model_early)$coefficients[str_c("gene", this_gene), 5]
      p_value_model_comparison <- anova(baseline_coxph_model_early, coxph_model_early)["P(>|Chi|)"][2,1]
      
      if (p_value_early < 0.05 & p_value_model_comparison < 0.05) {
        coxph_model_early_list[[this_gene]] <- coxph_model_early
      }
    }
  }
  
  # Total fusion burden and survival
  EFS_tibble <- seqfish_clinical_info %>% 
    filter(mmrf %in% mmrf_primary_pretreatment) %>%
    filter(!is.na(ISS_Stage), !is.na(EFS), !is.na(Age), !is.na(total_fusions)) %>% 
    select(mmrf, Age, total_fusions, ISS_Stage, EFS, EFS_censor)
  total_fusions_coxph_model <- coxph(Surv(EFS, EFS_censor == 0) ~ ISS_Stage + Age + total_fusions, data = EFS_tibble)
  
  early_tibble <- seqfish_clinical_info %>% 
    filter(mmrf %in% mmrf_primary_pretreatment) %>%
    filter(!is.na(ISS_Stage), !is.na(early_relapse_time), !is.na(Age), !is.na(total_fusions)) %>% 
    select(mmrf, Age, total_fusions, ISS_Stage, early_relapse_time, early_relapse_censor)
  total_fusions_coxph_model_early <- coxph(Surv(early_relapse_time, early_relapse_censor == 0) ~ ISS_Stage + Age + total_fusions, data = early_tibble)

  # ============================================================================
  # Double hit vs. Triple hit
  # ============================================================================
  
  triple_hit <- seqfish_clinical_info %>%
    filter(mmrf %in% mmrf_primary_pretreatment) %>%
    filter(!is.na(EFS), !is.na(EFS_censor), 
           !is.na(amp1q), !is.na(del17p), !is.na(updated_seqfish_t_IGH_WHSC1), 
           !is.na(Age), !is.na(ISS_Stage)) %>%
    select(EFS, EFS_censor, 
           amp1q, del17p, 
           updated_seqfish_t_IGH_WHSC1, 
           Age, ISS_Stage)
  
  fit_double <- survfit(Surv(EFS, EFS_censor == 0) ~ amp1q + del17p, data = triple_hit)
  fit_triple <- survfit(Surv(EFS, EFS_censor == 0) ~ amp1q + del17p + updated_seqfish_t_IGH_WHSC1, data = triple_hit)
  
  pdf(str_c(paper_supp, "triple_hit.EFS.with_legend.pdf"),
      width = 4.5, height = 4.5, useDingbats = FALSE)
  print(ggsurvplot(fit_triple, data = triple_hit, conf.int = TRUE,
                   surv.median.line = "hv", pval = TRUE,
                   legend.labs = c("-/-/-", "-/-/+", "-/+/-", "-/+/+", "+/-/-", "+/-/+", "+/+/-", "+/+/+"),
                   legend = "bottom",
                   legend.title = "amp(1q)/del(17p)/t(4;14)",
                   xlab = "Time (days)", 
                   ylab = "Progression-Free Survival Probability",
                   palette = "RdBu",
                   ggtheme = theme_survminer(base_size = 12,
                                             base_family = "",
                                             font.main = c(12, "plain", "black"),
                                             font.submain = c(12, "plain", "black"),
                                             font.x = c(12, "plain", "black"),
                                             font.y = c(12, "plain", "black"),
                                             font.caption = c(12, "plain", "black"),
                                             font.tickslab = c(8, "plain", "black"),
                                             legend = c("top", "bottom", "left", "right", "none"),
                                             font.legend = c(8, "plain", "black")),
                   conf.int.alpha = 0.1))
  dev.off()
  
  pdf(str_c(paper_supp, "triple_hit.EFS.without_legend.pdf"),
      width = 3.5, height = 3.5, useDingbats = FALSE)
  print(ggsurvplot(fit_triple, data = triple_hit, conf.int = TRUE,
                   surv.median.line = "hv", pval = TRUE,
                   legend.labs = c("-/-/-", "-/-/+", "-/+/-", "-/+/+", "+/-/-", "+/-/+", "+/+/-", "+/+/+"),
                   legend = "none",
                   legend.title = "amp(1q)/del(17p)/t(4;14)",
                   xlab = "Time (days)", 
                   ylab = "Progression-Free Survival Probability",
                   palette = "RdBu",
                   ggtheme = theme_survminer(base_size = 12,
                                             base_family = "",
                                             font.main = c(12, "plain", "black"),
                                             font.submain = c(12, "plain", "black"),
                                             font.x = c(12, "plain", "black"),
                                             font.y = c(12, "plain", "black"),
                                             font.caption = c(12, "plain", "black"),
                                             font.tickslab = c(8, "plain", "black"),
                                             legend = c("top", "bottom", "left", "right", "none"),
                                             font.legend = c(8, "plain", "black")),
                   conf.int.alpha = 0.1))
  dev.off()
  
  pdf(str_c(paper_supp, "double_hit.EFS.with_legend.pdf"),
      width = 4.5, height = 4.5, useDingbats = FALSE)
  print(ggsurvplot(fit_double, data = triple_hit, conf.int = TRUE,
                   surv.median.line = "hv", pval = TRUE,
                   legend.labs = c("-/-", "-/+", "+/-", "+/+"),
                   legend = "bottom",
                   legend.title = "amp(1q)/del(17p)",
                   xlab = "Time (days)", 
                   ylab = "Progression-Free Survival Probability",
                   palette = "PRGn",
                   ggtheme = theme_survminer(base_size = 12,
                                             base_family = "",
                                             font.main = c(12, "plain", "black"),
                                             font.submain = c(12, "plain", "black"),
                                             font.x = c(12, "plain", "black"),
                                             font.y = c(12, "plain", "black"),
                                             font.caption = c(12, "plain", "black"),
                                             font.tickslab = c(8, "plain", "black"),
                                             legend = c("top", "bottom", "left", "right", "none"),
                                             font.legend = c(8, "plain", "black")),
                   conf.int.alpha = 0.1))
  dev.off()
  
  pdf(str_c(paper_supp, "double_hit.EFS.without_legend.pdf"),
      width = 3.5, height = 3.5, useDingbats = FALSE)
  print(ggsurvplot(fit_double, data = triple_hit, conf.int = TRUE,
                   surv.median.line = "hv", pval = TRUE,
                   legend.labs = c("-/-", "-/+", "+/-", "+/+"),
                   legend = "none",
                   legend.title = "amp(1q)/del(17p)",
                   xlab = "Time (days)", 
                   ylab = "Progression-Free Survival Probability",
                   palette = "PRGn",
                   ggtheme = theme_survminer(base_size = 12,
                                             base_family = "",
                                             font.main = c(12, "plain", "black"),
                                             font.submain = c(12, "plain", "black"),
                                             font.x = c(12, "plain", "black"),
                                             font.y = c(12, "plain", "black"),
                                             font.caption = c(12, "plain", "black"),
                                             font.tickslab = c(8, "plain", "black"),
                                             legend = c("top", "bottom", "left", "right", "none"),
                                             font.legend = c(8, "plain", "black")),
                   conf.int.alpha = 0.1))
  dev.off()
}

# ==============================================================================
# Connection to TCGA cancer fusion calls
# ==============================================================================

if (TRUE) {
  tcga_direct_overlap <- fusions_primary %>% 
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
    write_tsv(str_c(paper_main, "TCGA_overlap.tsv"))
  
  p <- ggplot(tcga_direct_overlap, 
              aes(x = factor(Fusion, levels = c("SND1--BRAF", "TPM3--NTRK1", "NOTCH2--SEC22B"), ordered = TRUE), 
                  y = TCGA_Sample_Count, 
                  fill = Cancer)) + 
    geom_col(alpha = 0.75) + 
    scale_fill_brewer(palette = "BuPu") +
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
# NTRK1 fusion structures -- manually created based on agFusion output
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
# ==============================================================================

if (TRUE) {
  tcga_kinase_overlap <- kinases %>% 
    filter(KinasePos == "3P_KINASE" &
             KinaseDomain == "Intact" |
             (geneB == "MAP3K14" & SampleID %in% map3k14_intact_manual_review)) %>% 
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
  
  n_3p_kinase_cancer_types <- tcga_kinase_overlap %>% pull(Cancer) %>% unique() %>% length()
  
  ggplot(data = tcga_kinase_overlap, aes(y = geneB_mmrf_count, x = Cancer)) + 
    geom_tile(aes(fill = factor(count2))) + 
    geom_text(data = tcga_kinase_overlap %>% filter(count2 < 6), aes(label = count2),
              color = "#000000") +
    geom_text(data = tcga_kinase_overlap %>% filter(count2 >= 6), aes(label = count2),
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

# ==============================================================================
# APOBEC signature association with fusion events
# ==============================================================================

if (TRUE) {
  
  apobec_q25 <- mutsig %>% pull(APOBEC) %>% quantile(.25)
  apobec_q75 <- mutsig %>% pull(APOBEC) %>% quantile(.75)
  apobec_iqr <- apobec_q75 - apobec_q25
  apobec_outlier_true <- as.numeric(apobec_q75 + 1.5*apobec_iqr)
  apobec_outlier_gt <- max(as.numeric(apobec_q75 + 1.5*apobec_iqr), 0.5)
  
  genes_gt2 <- rbind(fusions_primary %>% 
                       filter(mmrf %in% mutsig$mmrf) %>% 
                       select(mmrf, visit_number, geneA) %>% 
                       select(geneA) %>% 
                       rename("gene" = "geneA"), 
                     fusions_primary %>% 
                       filter(mmrf %in% mutsig$mmrf) %>% 
                       select(mmrf, visit_number, geneB) %>% 
                       select(geneB) %>% 
                       rename("gene" = "geneB")) %>% 
    group_by(gene) %>% 
    summarize(count = n()) %>% 
    filter(count > 2) %>% pull(gene)
  
  outlier_apobec <- rbind(fusions_primary %>%
                            filter(mmrf %in% mutsig$mmrf) %>% 
                            left_join(mutsig, by = "mmrf") %>% 
                            filter(geneA %in% genes_gt2) %>%
                            rename("gene" = "geneA") %>%
                            select(mmrf, visit_number, gene, APOBEC) %>%
                            unique(),
                          fusions_primary %>%
                            filter(mmrf %in% mutsig$mmrf) %>% 
                            left_join(mutsig, by = "mmrf") %>% 
                            filter(geneB %in% genes_gt2) %>%
                            rename("gene" = "geneB") %>%
                            select(mmrf, visit_number, gene, APOBEC) %>%
                            unique()) %>%
    unique() %>%
    group_by(gene) %>%
    summarize(apobec_median = median(APOBEC),
              count = n()) %>%
    filter(apobec_median > apobec_outlier_gt,
           count >= 3)
  
  apobec_background_high <- rbind(mutsig %>% 
                                    mutate(gene = "All samples") %>% 
                                    select(mmrf, srr, visit, gene, APOBEC),
                                  rbind(fusions_primary %>% 
                                          filter(mmrf %in% mutsig$mmrf) %>%
                                          filter(geneA %in% outlier_apobec$gene) %>%
                                          rename("gene" = "geneA", "visit" = "visit_number") %>% 
                                          select(mmrf, srr, visit, gene),
                                        fusions_primary %>% 
                                          filter(mmrf %in% mutsig$mmrf) %>%
                                          filter(geneB %in% outlier_apobec$gene) %>%
                                          rename("gene" = "geneB", "visit" = "visit_number") %>% 
                                          select(mmrf, srr, visit, gene)) %>% 
                                    unique() %>% 
                                    left_join(mutsig, by = c("mmrf", "srr", "visit")) %>% 
                                    select(mmrf, srr, visit, gene, APOBEC)) %>%
    mutate(my_alpha = case_when(gene == "All samples" ~ 0.75,
                                TRUE ~ 1)) %>%
    mutate(my_size = case_when(gene == "All samples" ~ 1,
                               TRUE ~ 2))
  
  n_high_apobec <- apobec_background_high %>% filter(APOBEC >= apobec_outlier_true) %>% nrow()
  n_with_signature_score <- apobec_background_high %>% nrow()
  
  ggplot(apobec_background_high, 
         aes(x = fct_reorder(gene, APOBEC), y = APOBEC)) +
    geom_violin(scale = "width",
                color = "black",
                draw_quantiles = 0.5) +
    geom_jitter(aes(size = my_size,
                    alpha = my_alpha), 
                color = "#6e016b",
                size = 2,
                shape = 16, height = 0, width = 0.1) +
    geom_hline(yintercept = apobec_outlier_true, linetype = 2) +
    annotate("text", x = "MAF", y = apobec_outlier_true, 
             label = str_c("Outliers >\n", round(apobec_outlier_true, 3)),
             vjust = 0.5) +
    guides(alpha = FALSE, color = FALSE) +
    scale_y_continuous(expand = c(0,0), limits = c(-0.05, 1)) +
    labs(x = "Gene", y = "APOBEC Signature Score") +
    theme_bw() +
    coord_flip() +
    theme(panel.background = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_text(vjust = 0.5, size = 8),
          axis.text.y = element_text(size = 10),
          legend.position = "bottom",
          legend.direction = "vertical",
          axis.title = element_text(size = 12)) +
    ggsave(str_c(paper_main, "apobec.pdf"), width = 2.5, height = 1, useDingbats = FALSE)
}

# ==============================================================================
# Fusion druggable/clinical paragraph output
# ==============================================================================
n_patients_druggable <- fusions_primary %>% 
  filter(drug_geneA == 1 | drug_geneB == 1 | drug_fusion == 1) %>% 
  pull(mmrf) %>% unique() %>% length()

multiple_kinase_overlap <- tcga_kinase_overlap %>% ungroup() %>% filter(count2 > 1) %>% 
  separate(geneB_mmrf_count, into = c("gene", "count", "end"), by = c(" ", "(")) %>% 
  filter(count > 1) %>% pull(gene) %>% unique() %>% print()

n_tcga_ntrk1 <- pancan_fusions %>% 
  filter(str_detect(Fusion, pattern = "--NTRK1")) %>%
  pull(Cancer) %>% unique() %>% length()

print(str_c("Proportion of patients with druggable fusion: ", 
            n_patients_druggable, "/", n_samples_primary, " = ", 
            round(100*n_patients_druggable/n_samples_primary, 2), "%"))
print("Number of fusion genes in DEPO. All sensitive?")
drug_df %>% pull(Effect) %>% table() %>% print()
print("Example: BRAF fusions are druggable:")
drug_df %>% filter(gene == "BRAF") %>% print()
print("Number of TCGA cancer types with direct overlap:")
tcga_direct_overlap %>% pull(Cancer) %>% unique() %>% length() %>% print()
print("Number of TCGA cancer types with 3' kinase overlap:")
print(n_3p_kinase_cancer_types)
print("Cancer types with multiple reported overlaps:")
print(multiple_kinase_overlap)
print(str_c("Number of TCGA cancer types with 3' NTRK1 fusion: ", n_tcga_ntrk1))

print("Info for NTRK1 story: ")
fusions_all %>% filter(geneB == "NTRK1") %>% 
  select(mmrf, srr, fusion, sample_number, has_secondary, visit_number, 
         LeftBreakpoint, RightBreakpoint, PROT_FUSION_TYPE, Callers, n_discordant, 
         geneA_pct, geneB_pct, geneA_log2ratio_cnv, geneB_log2ratio_cnv) %>% print()

print(str_c("Proportion of patients with outlier APOBEC score: ", n_high_apobec, "/", n_with_signature_score, " = ", 100*round(n_high_apobec/n_with_signature_score, 4), "%"))

print("")
print("# ===== CLINICAL INFO FOR OVERVIEW SECTION ===== #")
print("")
print("HR estimate for total_fusions:")
print(total_fusions_coxph_model)
print("Confidence intervals for total_fusions HR:")
print(exp(confint(total_fusions_coxph_model)))
print("PFS time for double/triple hit:")
print(fit_double)
print(fit_triple)
# ==============================================================================
# From 01_overview.R
# ============================================================================== 

if (FALSE) {
  # ============================================================================
  # Function to create heatmap of binary seqFISH data
  # ============================================================================
  
  plot_seqfish_heatmap <- function(input_tibble, output_filename){
    
    arranged_info <- input_tibble %>% 
      arrange(desc(seqfish_Hyperdiploidy), 
              desc(seqfish_CN_del_13q14),
              desc(seqfish_CN_del_13q34),
              desc(seqfish_CN_del_17p13), 
              desc(seqfish_CN_gain_1q21),
              desc(seqfish_Translocation_WHSC1_4_14),
              desc(seqfish_Translocation_CCND3_6_14),
              desc(seqfish_Translocation_MYC_8_14), 
              desc(seqfish_Translocation_MAFA_8_14),
              desc(seqfish_Translocation_CCND1_11_14), 
              desc(seqfish_Translocation_CCND2_12_14),
              desc(seqfish_Translocation_MAF_14_16), 
              desc(seqfish_Translocation_MAFB_14_20))
    
    pheatmap_df <- data.frame(
      arranged_info %>%
        select(seqfish_Hyperdiploidy, 
               seqfish_CN_del_13q14, 
               seqfish_CN_del_13q34, 
               seqfish_CN_del_17p13, 
               seqfish_CN_gain_1q21, 
               seqfish_Translocation_WHSC1_4_14,
               seqfish_Translocation_CCND3_6_14, 
               seqfish_Translocation_MYC_8_14, 
               seqfish_Translocation_MAFA_8_14, 
               seqfish_Translocation_CCND1_11_14, 
               seqfish_Translocation_CCND2_12_14, 
               seqfish_Translocation_MAF_14_16, 
               seqfish_Translocation_MAFB_14_20) %>% 
        mutate_all(list(~replace(., is.na(.), 2))) )
    
    names(pheatmap_df) <- c("Hyperdiploidy", 
                            "CNV del(13q14)", 
                            "CNV del(13q34)", 
                            "CNV_del(17p13)", 
                            "CNV gain(1q21)", 
                            "Translocation t(4;14) (WHSC1)",
                            "Translocation t(6;14) (CCND3)", 
                            "Translocation t(8;14) (MYC)", 
                            "Translocation t(8;14) (MAFA)", 
                            "Translocation t(11;14) (CCND1)", 
                            "Translocation t(12;14) (CCND2)", 
                            "Translocation t(14;16) (MAF)", 
                            "Translocation t(14;20) (MAFB)")
    
    pheatmap_mat <- t(as.matrix(pheatmap_df))
    colnames(pheatmap_mat) <- arranged_info %>% pull(mmrf)
    
    annotation_col_df <- data.frame(
      arranged_info %>%
        select(age_ge_66, Female, race, ISS_Stage) %>%
        mutate_all(factor))
    names(annotation_col_df) <- c("Age >= 66", "Sex", "Ethnicity", "ISS Stage")
    levels(annotation_col_df$Sex) <- c("Male", "Female")
    levels(annotation_col_df$Ethnicity) <- c("White", "Black", "Other")
    levels(annotation_col_df$`Age >= 66`) <- c("No", "Yes")
    rownames(annotation_col_df) <- arranged_info %>% pull(mmrf)
    
    pheatmap(pheatmap_mat, color = brewer.pal(n = 3, name = "Blues"), 
             annotation_col = annotation_col_df, gaps_row = c(5),
             cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = FALSE, 
             cellwidth = 10, cellheight = 10, height = 10, width = 55,  
             legend = FALSE, filename = output_filename)
  }
  
  # ============================================================================
  # Function to plot seqFISH as UpSetR
  # ==============================================================================
  
  plot_seqfish_upsetr <- function(plot_tibble, output_filename){
    upsetr_df <- plot_tibble %>% 
      select(seqfish_Hyperdiploidy, 
             seqfish_CN_del_13q14, 
             seqfish_CN_del_13q34, 
             seqfish_CN_del_17p13, 
             seqfish_CN_gain_1q21, 
             seqfish_Translocation_WHSC1_4_14, 
             seqfish_Translocation_CCND3_6_14, 
             seqfish_Translocation_MYC_8_14, 
             seqfish_Translocation_MAFA_8_14, 
             seqfish_Translocation_CCND1_11_14, 
             seqfish_Translocation_CCND2_12_14,
             seqfish_Translocation_MAF_14_16, 
             seqfish_Translocation_MAFB_14_20) %>% 
      rowwise() %>% mutate(none = as.numeric(
        sum(seqfish_Hyperdiploidy, 
            seqfish_CN_del_13q14, 
            seqfish_CN_del_13q34, 
            seqfish_CN_del_17p13, 
            seqfish_CN_gain_1q21, 
            seqfish_Translocation_WHSC1_4_14,
            seqfish_Translocation_CCND3_6_14, 
            seqfish_Translocation_MYC_8_14, 
            seqfish_Translocation_MAFA_8_14, 
            seqfish_Translocation_CCND1_11_14,
            seqfish_Translocation_CCND2_12_14, 
            seqfish_Translocation_MAF_14_16, 
            seqfish_Translocation_MAFB_14_20) 
        == 0)) %>% as.data.frame()
    names(upsetr_df) <- c("Hyperdiploidy", 
                          "CNV del(13q14)", 
                          "CNV del(13q34)", 
                          "CNV_del(17p13)", 
                          "CNV gain(1q21)", 
                          "Translocation t(4;14) (WHSC1)",
                          "Translocation t(6;14) (CCND3)", 
                          "Translocation t(8;14) (MYC)", 
                          "Translocation t(8;14) (MAFA)", 
                          "Translocation t(11;14) (CCND1)", 
                          "Translocation t(12;14) (CCND2)", 
                          "Translocation t(14;16) (MAF)", 
                          "Translocation t(14;20) (MAFB)", 
                          "No seqFISH events")
    
    pdf(file = output_filename, width = 40, height = 20, useDingbats = FALSE)
    upset(upsetr_df, 
          nsets = ncol(upsetr_df), 
          nintersects = NA, 
          order.by = "freq",
          #text.scale = c(intersection size title, 
          #               intersection size tick labels, 
          #               set size title, 
          #               set size tick labels, 
          #               set names, 
          #               numbers above bars),
          text.scale = c(4, 
                         3, 
                         4, 
                         3, 
                         3, 
                         3),
          point.size = 4, 
          line.size = 2)
    dev.off()
  }
  
  # ============================================================================
  # Plot two different heatmaps and upsetrs: hyperdiploid, non-hyperdiploid
  # ============================================================================
  
  #hyperdiploid_status <- c(1, 0)
  #hyperdiploid_string <- c("hyperdiploid", "non-hyperdiploid")
  #for (i in 1:2) {
  #  h_status <- hyperdiploid_status[i]
  #  h_string <- hyperdiploid_string[i]
  #  plot_tibble <- seqfish_clinical_info %>% 
  #    filter(HRD == h_status)
  #  plot_seqfish_heatmap(plot_tibble, 
  #                       str_c(paper_supp, "heatmap.", h_string, ".pdf"))
  #  plot_seqfish_upsetr(plot_tibble, 
  #                      str_c(paper_supp, "upsetr.", h_string, ".pdf")) 
  #  
  #}
  
  # ============================================================================
  # Plot one upsetr for all samples
  # ============================================================================
  
  #plot_tibble <- seqfish_clinical_info %>% filter(!is.na(seqfish_Hyperdiploidy))
  #plot_seqfish_upsetr(plot_tibble, str_c(paper_supp, "upsetr.both.pdf"))
  
  # ============================================================================
  # Plot continuous clinical variables
  # ============================================================================
  
  pdf(str_c(paper_supp, "clinical_variables.age.pdf"), height = 10, width = 15,
      useDingbats = FALSE)
  n_missing <- seqfish_clinical_info %>% filter(is.na(Age)) %>% nrow()
  p <- seqfish_clinical_info %>% 
    ggplot(aes(x = Age)) + geom_histogram() + 
    labs(x = "Age at onset (years)", y = "Number of patients", 
         caption = str_c("Number missing = ", n_missing)) +
    theme_bw(base_size = 20)
  print(p)
  dev.off()
  
  pdf(str_c(paper_supp, "clinical_variables.bm_pct.pdf"), height = 10, 
      width = 15, useDingbats = FALSE)
  n_missing <- seqfish_clinical_info %>% 
    filter(is.na(BM_Plasma_Cell_Percent)) %>% nrow()
  p <- seqfish_clinical_info %>% 
    ggplot(aes(x = BM_Plasma_Cell_Percent)) + geom_histogram() + 
    labs(x = "Bone marrow plasma cell (%)", y = "Number of patients",
         caption = str_c("Number missing = ", n_missing)) +
    theme_bw(base_size = 20)
  print(p)
  dev.off()
  
  pdf(str_c(paper_supp, "clinical_variables.ldh.pdf"), height = 10, width = 15,
      useDingbats = FALSE)
  n_missing <- seqfish_clinical_info %>% filter(is.na(LDH)) %>% nrow()
  p <- seqfish_clinical_info %>% 
    ggplot(aes(x = LDH)) + geom_histogram() +
    labs(x = "Lactate dehydrogenase (LDH) (U/L)", y = "Number of patients",
         caption = str_c("Number missing = ", n_missing)) +
    theme_bw(base_size = 20)
  print(p)
  dev.off()
  
  
}

# ==============================================================================
# From 02_expression.R
# ==============================================================================


#BiocManager::install("ConsensusClusterPlus")
library(ConsensusClusterPlus)

all_expression <- read_tsv("data/mmy_gene_tpm_table.tsv")

most_variable_long <- all_expression %>% 
  group_by(gene) %>% 
  summarize(variability = var(log10tpm)) %>% 
  filter(!str_detect(gene, "^IG")) %>% 
  arrange(desc(variability)) %>% 
  head(5000) %>% 
  select(gene) %>% 
  left_join(all_expression, by = "gene") %>% 
  select(gene, srr, log10tpm)

most_variable_wide <- most_variable_long %>% spread(key = srr, value = log10tpm)
most_variable_matrix <- as.matrix(data.frame(most_variable_wide[,-1], 
                                             row.names = most_variable_wide$gene))

results <- ConsensusClusterPlus(d = most_variable_matrix,
                                maxK = 15,
                                reps = 1000,
                                pItem = 1,
                                pFeature = 0.8,
                                title = "xxx",
                                clusterAlg = "pam",
                                #clusterAlg = "hc",
                                distance = "pearson",
                                #distance = "euclidean",
                                seed = 1,
                                plot = "pdf")

icl <- calcICL(results, title = "yyy", plot = "pdf")

clinical_with_cluster <- as_tibble(data.frame(srr = as.character(names(results[[12]][["consensusClass"]])), 
                                              cluster = as.vector(results[[12]][["consensusClass"]]))) %>%
  mutate(srr = as.character(srr)) %>% 
  left_join(samples_primary, by = "srr") %>% 
  left_join(seqfish_clinical_info, by = "mmrf") %>% 
  filter(!is.na(seqfish_Study_Visit_ID)) %>%
  select(mmrf, cluster, del17p, amp1q, HRD, ISS_Stage, Age, 
         EFS, EFS_censor, early_relapse_time, early_relapse_censor, OS, OS_censor,
         starts_with("updated"))


coxph(Surv(early_relapse_time, early_relapse_censor == 0) ~ as.factor(cluster == 11) + as.factor(updated_seqfish_t_IGH_WHSC1), data = clinical_with_cluster)
coxph(Surv(EFS, EFS_censor == 0) ~ as.factor(cluster == 11) + as.factor(updated_seqfish_t_IGH_WHSC1), data = clinical_with_cluster)

#fit <- survfit(Surv(early_relapse_time, early_relapse_censor == 0) ~ as.factor(cluster == 11) + as.factor(updated_seqfish_t_IGH_WHSC1), data = clinical_with_cluster)
#fit <- survfit(Surv(EFS, EFS_censor == 0) ~ as.factor(cluster == 11) + as.factor(updated_seqfish_t_IGH_WHSC1), data = clinical_with_cluster)



################################################################################
# Look at other clinical associations
# Written April 2019
################################################################################

if (FALSE) { # No longer relevant after removing t-test due to bad assumptions
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

# ==============================================================================
# Types of kinases
# ==============================================================================

if (TRUE) {
  p <- kinases %>%
    group_by(kinase_group_full_name) %>% 
    summarize(count = n()) %>% 
    ungroup() %>%
    ggplot(aes(x = fct_reorder(kinase_group_full_name, count), y = count)) +
    geom_col(position = "dodge") +
    coord_flip(expand = FALSE) +
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
  
  ggsave(str_c(paper_supp, "kinase_groups.without_legend.pdf"), 
         p + guides(fill = FALSE), 
         width = 7.25, height = 3.5, useDingbats = FALSE)
  
}

# ==============================================================================
# Expression correlation of CBX7--CSNK1E the only recurrent fusion with
# 3' intact kinase
# ==============================================================================

if (FALSE) {
  geneA <- "CBX7"
  geneB <- "CSNK1E"
  this_fusion = str_c(geneA, "--", geneB)
  
  fusion_samples <- kinases %>% 
    filter(KinasePos == "3P_KINASE", KinaseDomain == "Intact") %>% 
    filter(Fusion == this_fusion) %>% 
    pull(SampleID)
  
  geneA_expr <- expression_primary %>% filter(gene == geneA) %>%
    select(srr, gene, log10tpm, pct) %>%
    rename(geneA = gene, geneA_log10tpm = log10tpm, geneA_pct = pct)
  
  geneB_expr <- expression_primary %>% filter(gene == geneB) %>%
    select(srr, gene, log10tpm, pct) %>%
    rename(geneB = gene, geneB_log10tpm = log10tpm, geneB_pct = pct)
  
  geneA_geneB_expr <- geneA_expr %>% left_join(geneB_expr, by = "srr") %>%
    mutate(has_fusion = srr %in% fusion_samples) %>%
    arrange(has_fusion)
  
  ggplot(geneA_geneB_expr, 
         aes(x = geneA_pct, y = geneB_pct, color = has_fusion)) +
    geom_point(shape = 16, size = 2) + 
    coord_fixed() + 
    geom_segment(x = 0, xend = 1, y = 0, yend = 1, 
                 linetype = 2, show.legend = FALSE,
                 color = "grey90") +
    scale_x_continuous(expand = c(0.01, 0.01), limits = c(0,1)) +
    scale_y_continuous(expand = c(0.01, 0.01), limits = c(0,1)) +
    scale_color_manual(values = c("grey90", "black")) +
    labs(x = str_c(geneA, " Expression Percentile"),
         y = str_c(geneB, " Expression Percentile"),
         color = str_c(this_fusion, "\nFusion Reported")) +
    theme_bw() +
    theme(panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_text(size = 8),
          axis.title = element_text(size = 10),
          legend.position = "bottom"
    ) +
    ggsave(str_c(paper_supp, "CBX7--CSNK1E.expression.pdf"),
           width = 3.5, height = 3.5, useDingbats = FALSE)  
}
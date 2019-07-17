# ==============================================================================
# Multiple time points
# Steven Foltz (smfoltz@wustl.edu), April 2019
# ==============================================================================

paper_main = "paper/main/05_multiple/"
paper_supp = "paper/supplemental/05_multiple/"

# Create directories 
dir.create(paper_main, recursive = TRUE, showWarnings = FALSE)
dir.create(paper_supp, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# Explore the landscape of fusions in samples with multiple timepoints
# Originally written September 2018, Updated April 2019
# ==============================================================================

if (TRUE) {
  
  # Samples with multiple timepoints
  
  samples_with_multiple_timepoints <- fusions_all %>% 
    filter(has_secondary == 1) %>% pull(mmrf) %>% unique()

  # samples with multiple bone marrow samples
  mtp_bm_samples <- samples_all %>% 
    filter(tissue_source == "BM") %>%
    group_by(mmrf) %>%
    summarize(n_bm_samples = n(),
              bm_visits = str_c(visits, sep = ":"),
              bm_srr = str_c(srr, sep = ":")) %>%
    filter(n_bm_samples > 1) %>%
    mutate(first_srr = strsplit(bm_srr, split = ":")[1],
           second_srr = strsplit(bm_srr, split = ":")[2]) %>%
    select(mmrf, first_srr, second_srr)
  
  # samples with same time points BMs and PBs
  stp_bmpb_samples <- samples_all %>%
    group_by(mmrf, visit) %>%
    summarize(n_samples = n(),
              tissue_sources = str_c(tissue_source, sep = ":"),
              source_srr = str_c(srr, sep = ":")) %>%
    filter(n_samples > 1,
           str_detect(tissue_sources, "BM"),
           str_detect(tissue_sources, "PB")) %>%
    mutate(first_srr = strsplit(source_srr, split = ":")[1],
           second_srr = strsplit(source_srr, split = ":")[2]) %>%
    select(mmrf, first_srr, second_srr)
  
  # method to overlap fusion calls from two SRRs
  get_fusion_overlap <- function(fusion_df, srr1, srr2){
    list_of_fusions_srr1 <- fusions_df %>% 
      filter(srr == srr1) %>% 
      pull(fusion) %>% 
      sort()
    list_of_fusions_srr2 <- fusions_df %>% 
      filter(srr == srr2) %>% 
      pull(fusion) %>% 
      sort()
    n_srr1 <- length(list_of_fusions_srr1)
    n_srr2 <- length(list_of_fusions_srr2)
    overlap_srr12 <- intersect(list_of_fusions_srr1,
                               list_of_fusions_srr2) %>% length()
    return(str_c(n_srr1, n_srr2, overlap_srr12, sep = ":"))
  }
  
  mtp_bm_samples %>% mutate(overlaps = get_fusion_overlap(fusions_all, srr1, srr2)) %>%
    separate(overlaps, into = c("n_fusions_srr1", "n_fusions_srr2", "overlap_srr12"), sep = ":") %>%
    ggplot(aes(x = n_fusions_srr1, y = n_fusions_srr2, size = overlap_srr12, color = overlap_srr12)) %>%
    geom_point(shape = 16) +
    geom_point(shape = 3, size = 2, color = "#000000") +
    geom_smooth(method = "lm") +
    scale_size_area(breaks = c(0, 5, 10), range = c(0,15)) +
    scale_color_brewer(palette = "BuGn") +
    theme_bw() +
    theme(plot.background = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_text(size = 8),
          axis.ticks = element_blank(),
          axis.title = element_text(size = 12),
          legend.background = element_blank(),
          legend.text = element_text(size = 10)) +
    ggsave(str_c(paper_main, "multiple_timepoints.bm.pdf"),
           device = "pdf", width = 3.5, height = 3.5, useDingbats = FALSE)
  
  # Fusions involving important genes
  
  fusions_with_important_genes <- fusions_all %>%
    filter( str_detect(fusion, "WHSC1") |
              str_detect(fusion, "FGFR3") |
              str_detect(fusion, "MYC") |
              str_detect(fusion, "PVT1") |
              !(geneA %in% c("IGH", "IGK", "IGL") |
                  geneB %in% c("IGH", "IGK", "IGL"))) %>%
    filter( geneA_oncogene | geneA_tsg | geneA_kinase | 
              geneA_mmy_known | geneA_driver |
              geneB_oncogene | geneB_tsg | geneB_kinase |
              geneB_mmy_known | geneB_driver |
              drug_fusion | drug_geneA | drug_geneB) %>% 
    select(mmrf, fusion) %>% unique() %>% group_by(fusion) %>% 
    summarize(n()) %>%
    pull(fusion)
  
  # plot it 
  keep_these_mmrfs <- fusions_all %>% 
    filter(fusion %in% fusions_with_important_genes, has_secondary) %>% 
    pull(mmrf) %>% unique()
  
  keep_these_srrs <- samples_all %>% 
    group_by(mmrf) %>% 
    summarize(count = n()) %>% 
    filter(count > 1) %>% 
    left_join(samples_all, by = "mmrf") %>% 
    mutate(mmrf_srr = str_c(mmrf, srr, sep = ": ")) %>% 
    group_by(mmrf) %>% 
    mutate(ticker = row_number()) %>% 
    select(mmrf, count, srr, visit, mmrf_srr, ticker) %>%
    ungroup() %>% 
    filter(mmrf %in% keep_these_mmrfs)
  
  keep_these_srrs %>% 
    left_join(fusions_all, by = c("mmrf", "srr")) %>%
    filter(fusion %in% fusions_with_important_genes) %>%
    right_join(keep_these_srrs, by = c("mmrf", "srr")) %>% 
    mutate(mmrf_number_only = str_remove_all(mmrf, "MMRF_")) %>%
    replace_na(list(fusion = "Zero detected")) %>%
    mutate(updated_srr = str_c(srr, ": ", visit.y)) %>%
    left_join(tumor_purity %>% filter(Type == "RNA-Seq"), by = c("srr" = "SRR_Tumor")) %>%
    ggplot(aes(y = factor(updated_srr), x = fusion, fill = log10(FFPM + 1))) + 
    geom_tile(color = "black") +
    geom_text(aes(label = visit.y), size = 1.75, vjust = 0.5) +
    geom_text(aes(label = TumorPurity, color = TumorPurity), x = Inf) +
    theme_bw() +
    facet_wrap(~ mmrf_number_only , ncol = 1, strip.position = "right", dir = "h", scales = "free_y") +
    theme(panel.grid.major = element_line(size = 0.1),
          panel.border = element_rect(size = 0.1),
          axis.ticks = element_line(size = 0.1),
          strip.background = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5, face = "italic"),
          axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 5),
          strip.text.y = element_text(angle = 0),
          legend.position = "bottom",
          panel.spacing.y = unit(0.025, units = "inches")) +
    scale_fill_gradient(low = "#ffedbc", high = "#ed4264",
                        limits = c(0, 1.05), 
                        breaks = seq(0, 1, 0.25),
                        labels = seq(0, 1, 0.25)) + 
    scale_x_discrete(drop = FALSE) +
    labs(y = "Sample Number", x = NULL, fill = "Scaled FFPM") +
    ggsave(str_c(paper_supp, "multiple_timepoints.ffpm.pdf"),
           device = "pdf", width = 7.5, height = 9, useDingbats = FALSE)
  
  keep_these_srrs %>% left_join(fusions_all, by = c("mmrf", "srr")) %>%
    filter(fusion %in% fusions_with_important_genes) %>%
    right_join(keep_these_srrs, by = c("mmrf", "srr")) %>% 
    mutate(mmrf_number_only = str_remove_all(mmrf, "MMRF_")) %>%
    replace_na(list(fusion = "Zero detected")) %>%
    mutate(tpm_to_use = case_when((geneA %in% c("IGH", "IGK", "IGL") | geneB_oncogene | geneB_tsg | geneB_kinase | geneB_driver | drug_geneB) ~ geneB_pct,
                                  (geneB %in% c("IGH", "IGK", "IGL") | geneA_oncogene | geneA_tsg | geneA_kinase | geneA_driver | drug_geneA) ~ geneA_pct)) %>%
    mutate(updated_fusion = case_when((geneA %in% c("IGH", "IGK", "IGL") | geneB_oncogene | geneB_tsg | geneB_kinase | geneB_driver | drug_geneB) ~ str_c(fusion, "*"),
                                      (geneB %in% c("IGH", "IGK", "IGL") | geneA_oncogene | geneA_tsg | geneA_kinase | geneA_driver | drug_geneA) ~ str_c(geneA, "*", "--", geneB),
                                      TRUE ~ fusion)) %>%
    mutate(updated_srr = str_c(srr, ": ", ticker.y)) %>%
    ggplot(aes(y = factor(updated_srr), x = updated_fusion, fill = tpm_to_use)) + 
    geom_tile(color = "black") +
    geom_text(aes(label = sample_number), size = 1.75, vjust = 0.5) +
    theme_bw() +
    facet_wrap(~ mmrf_number_only , ncol = 1, strip.position = "right", dir = "h", scales = "free_y") +
    theme(panel.grid.major = element_line(size = 0.1),
          panel.border = element_rect(size = 0.1),
          axis.ticks = element_line(size = 0.1),
          strip.background = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5, face = "italic"),
          axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 5),
          strip.text.y = element_text(angle = 0),
          legend.position = "bottom",
          panel.spacing.y = unit(0.025, units = "inches")) +
    scale_fill_gradient(low = "#ffedbc", high = "#ed4264",
                        limits = c(0, 1.05), 
                        breaks = seq(0, 1, 0.25),
                        labels = seq(0, 1, 0.25)
                        ) + 
    scale_x_discrete(drop = FALSE) +
    labs(y = "Sample Number", x = NULL, fill = "TPM (log10)") +
    ggsave(str_c(paper_supp, "multiple_timepoints.tpm.pdf"),
           device = "pdf", width = 7.5, height = 9, useDingbats = FALSE)
}

# ==============================================================================
# Multiple time point fusions along with multiple exome seq mutations data
# Updated April 2019
# ==============================================================================

if (FALSE) {
  
  fusion_samples_mtp <- samples_all %>% 
    group_by(mmrf) %>% 
    summarize(count = n()) %>% 
    filter(count > 1) %>% 
    left_join(samples_all, by = "mmrf") %>% 
    arrange(mmrf, srr) %>% 
    group_by(mmrf) %>% 
    mutate(sample_n = row_number()) %>% 
    select(mmrf, srr, sample_n) %>% 
    mutate(data_type = "RNA", variant_type = "Fusion")
  
  snpindel_samples_mtp <- mutation_calls %>% 
    arrange(Tumor_Sample_Barcode) %>% 
    select(Tumor_Sample_Barcode, Matched_Norm_Sample_Barcode) %>% 
    unique() %>%
    separate(Tumor_Sample_Barcode, into = c("mmrf1", "mmrf2", "sample_n"), by = "_") %>% 
    separate(Matched_Norm_Sample_Barcode, into = c("mmrf11", "mmrf22", "srr_normal", "srr", "N")) %>% 
    mutate(sample_n = as.numeric(sample_n)) %>%
    mutate(mmrf = str_c(mmrf1, mmrf2, sep = "_"), data_type = "WES", variant_type = "SNP/INDEL") %>% 
    select(mmrf, srr, sample_n, data_type, variant_type)
  
  overlap_once <- fusion_samples_mtp %>% 
    left_join(snpindel_samples_mtp, by = c("mmrf", "sample_n")) %>% 
    filter(!is.na(srr.y))
  overlap_mtp <- overlap_once %>% 
    group_by(mmrf) %>% 
    summarize(count = n()) %>% 
    filter(count > 1) %>% 
    left_join(overlap_once, by = "mmrf")
  
  mtp_fusions <- overlap_mtp %>% 
    select(mmrf, srr.x, sample_n, data_type.x, variant_type.x) %>% 
    rename("srr" = "srr.x", "data_type" = "data_type.x", "variant_type" = "variant_type.x")
  
  mtp_mutations <- overlap_mtp %>% 
    select(mmrf, srr.y, sample_n, data_type.y, variant_type.y) %>% 
    rename("srr" = "srr.y", "data_type" = "data_type.y", "variant_type" = "variant_type.y")
  
  
  mtp_fusions_with_fusions <- fusions_all %>%
    filter( fusion == "IGH--WHSC1" | 
              !(geneA %in% c("IGH", "IGK", "IGL") |
                  geneB %in% c("IGH", "IGK", "IGL"))) %>%
    filter( geneA_oncogene | geneA_tsg | geneA_kinase | geneA_mmy_known | geneA_driver | 
              geneB_oncogene | geneB_tsg | geneB_kinase | geneB_mmy_known | geneB_driver | 
              drug_fusion | drug_geneA | drug_geneB) %>%
    right_join(mtp_fusions, by = c("mmrf", "srr")) %>%
    select(mmrf, sample_n, srr, data_type, variant_type, fusion, FFPM) %>%
    rename("variant" = "fusion", "detection_level" = "FFPM") %>%
    replace_na(list(variant = "Zero Detected"))
  
  mtp_mutations_with_mutations <- mutation_calls %>% 
    filter(Hugo_Symbol %in% drivers_kinases_oncogenes_mmy_genes) %>% 
    filter(Hugo_Symbol != "TTN") %>% # false positive
    separate(Tumor_Sample_Barcode, into = c("mmrf1", "mmrf2", "sample_n"), by = "_") %>% 
    separate(Matched_Norm_Sample_Barcode, into = c("mmrf11", "mmrf22", "srr_normal", "srr", "N")) %>% 
    mutate(mmrf = str_c(mmrf1, mmrf2, sep = "_")) %>% 
    right_join(mtp_mutations, by = c("mmrf", "srr")) %>%
    mutate(detection_level = t_alt_count/t_depth) %>%
    select(mmrf, sample_n.y, srr, data_type, variant_type, Hugo_Symbol, detection_level) %>%
    rename("sample_n" = "sample_n.y", "variant" = "Hugo_Symbol")
  
  
  mmrf_to_plot <- c("MMRF_1433", "MMRF_1783", "MMRF_2401", "MMRF_2490")
  
  fusion_plot_df <- mtp_fusions_with_fusions %>% 
    filter(mmrf %in% mmrf_to_plot) %>% 
    group_by(mmrf, variant) %>% 
    summarize(count = n()) %>% 
    filter(count > 1) %>% 
    left_join(mtp_fusions_with_fusions, by = c("mmrf", "variant"))
  
  mutation_plot_df_tp1 <- mtp_mutations_with_mutations %>%
    filter(mmrf %in% mmrf_to_plot, sample_n == 1) %>%
    group_by(mmrf, sample_n, variant) %>% 
    summarize(detection_level = mean(detection_level)) %>%
    ungroup()
  
  mutation_plot_df_tp2 <- mtp_mutations_with_mutations %>%
    filter(mmrf %in% mmrf_to_plot, sample_n == 2) %>%
    group_by(mmrf, sample_n, variant) %>% 
    summarize(detection_level = mean(detection_level)) %>%
    ungroup()
  
  mutation_plot_df <- mutation_plot_df_tp1 %>% 
    full_join(mutation_plot_df_tp2, by = c("mmrf", "variant")) %>%
    replace_na(list(detection_level.x = 0, detection_level.y = 0)) %>%
    rowwise() %>%
    mutate(label_text = case_when(variant %in% mmy21 ~ variant))
  
  p <- ggplot(fusion_plot_df, 
              aes(x = as.numeric(factor(sample_n)), y = variant)) +
    geom_point(aes(size = log10(detection_level + 1), 
                   color = log10(detection_level + 1)), shape = 16) +
    geom_text(data = fusion_plot_df %>% select(mmrf, variant) %>% unique(),
              aes(label = variant, x = 1.5), nudge_y = 0.15,
              show.legend = FALSE, fontface = "italic") +
    facet_wrap(~mmrf, scales = "free", nrow = 1) +
    labs(y = NULL, x = NULL, color = "FFPM") +
    guides(size = FALSE) +
    scale_color_gradient(low = "#ffedbc", high = "#ed4264",
                         limits = c(0,1.05), 
                         breaks = seq(0, 1, 0.25),
                         labels = seq(0, 1, 0.25)) + 
    theme_bw() +
    theme(panel.background = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.border = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = 10))
  
  q <- ggplot(mutation_plot_df,
              aes(x = detection_level.x,
                  y = detection_level.y,
                  color = label_text,
                  label = label_text)) +
    geom_point(shape = 16) +
    geom_label_repel(fontface = "italic") +
    facet_wrap(~mmrf, nrow = 1) +
    labs(x = "Time Point 1 Variant Allele Frequency (VAF)", 
         y = "Time Point 2 VAF") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), 
                       labels = c("0", "0.25", "0.50", "0.75", "1.00")) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), 
                       labels = c("0", "0.25", "0.50", "0.75", "1.00")) +
    guides(color = FALSE) +
    theme_bw() +
    theme(panel.background = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_blank(),
          strip.text = element_blank(),
          strip.background = element_blank(),
          axis.text = element_text(size = 8),
          axis.title = element_text(size = 12))
  
  ggsave(str_c(paper_main, "fusions.pdf"), p + guides(color = FALSE),
         width = 7, height = 7.25/8, useDingbats = FALSE)
  
  ggsave(str_c(paper_main, "mutations.pdf"), q,
         width = 7.25, height = 7.25/4, useDingbats = FALSE)
  
  x <- ggplot(mtp_fusions_with_fusions %>% filter(mmrf %in% mmrf_to_plot | mmrf == "MMRF_1496"),
              aes(y = factor(srr), x = variant, fill = log10(detection_level + 1))) + 
    geom_tile(color = "black") +
    geom_text(aes(label = sample_n), size = 1.75, vjust = 0.5) +
    theme_bw() +
    facet_wrap(~ mmrf, ncol = 1, strip.position = "right", dir = "h", scales = "free_y") +
    theme(panel.grid.major = element_line(size = 0.1),
          panel.border = element_rect(size = 0.1),
          axis.ticks = element_line(size = 0.1),
          strip.background = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10, face = "italic"),
          axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 10),
          strip.text.y = element_text(angle = 0, size = 10),
          legend.text = element_text(size = 10),
          legend.position = "bottom",
          legend.direction = "vertical",
          panel.spacing.y = unit(0.025, units = "inches")) +
    scale_fill_gradient(low = "#ffedbc", high = "#ed4264",
                        limits = c(0,1.05), 
                        breaks = seq(0, 1, 0.25),
                        labels = seq(0, 1, 0.25)) + 
    scale_x_discrete(drop = FALSE) +
    labs(y = NULL, x = NULL, fill = "Scaled FFPM")
  
  ggsave(str_c(paper_main, "multiple_timepoints.pdf"), x,
         device = "pdf", width = 7.25, height = 3.9, useDingbats = FALSE)
  ggsave(str_c(paper_main, "multiple_timepoints.no_legend.pdf"), x + guides(fill = FALSE),
         device = "pdf", width = 7.25, height = 3.9, useDingbats = FALSE)
}
# ==============================================================================
# MYC PVT1 story
# Steven Foltz (github: envest)
# ==============================================================================

paper_main = "paper/main/06_myc_pvt1/"
paper_supp = "paper/supplementary/06_myc_pvt1/"

# Create directories 
dir.create(paper_main, recursive = TRUE, showWarnings = FALSE)
dir.create(paper_supp, recursive = TRUE, showWarnings = FALSE)

# Expression plot data

expr_plot_df <- read_tsv("paper/supplementary/02_expression/expression_plot_tibble.tsv")
expr_plot_df <- expr_plot_df %>% mutate(cnv_factor = factor(categorical_cnv,
                                                            labels = c("DELETION",
                                                                       "Deletion",
                                                                       "Neutral",
                                                                       "Missing",
                                                                       "Amplification",
                                                                       "AMPLIFICATION"), 
                                                            exclude = NULL))

# ==============================================================================
# Plot MYC expression and fusion breakpoints
# ==============================================================================
if (TRUE) {
  
  get_t822 <- function(delly){
    chr_pos_list <- NULL
    vector_of_events <- str_split(delly, pattern = "\\*")[[1]]
    for (event in vector_of_events) {
      vector_of_sides = str_split(event, pattern = "\\|")[[1]]
      chr_pos_list <- list()
      for (side in vector_of_sides) {
        chr = str_split(side, pattern = "\\:")[[1]][1]
        pos = str_split(str_split(side, pattern = "\\:")[[1]][2], pattern = "-")[[1]][1]
        chr_pos_list[[chr]] <- as.numeric(pos)
      }
      if (all(sort(names(chr_pos_list)) == c("22", "8")) ) {
        return(chr_pos_list)
        break()
      }
    }
    return(list("8" = NA, "22" = NA))
  }
  
  mmrf_with_primary_mutation_calls <- mutation_calls %>% 
    separate(Tumor_Sample_Barcode, into = c("MMRF", "NUM", "VISIT"), by = "_") %>% 
    mutate(mmrf = str_c(MMRF, NUM, sep = "_")) %>% 
    select(mmrf, VISIT) %>%
    unique() %>% mutate(VISIT = as.numeric(VISIT)) %>%
    left_join(samples_all, by = c("mmrf" = "mmrf", "VISIT" = "visit")) %>% 
    filter(!is.na(tissue_source)) %>% 
    left_join(samples_primary, by = "srr") %>%
    filter(!is.na(mmrf.y)) %>%
    pull(mmrf.x)
  
  # MYC and PVT expression
  mmrf_with_myc_mutation <- mutation_calls %>% # There are none with PVT1
    filter(Hugo_Symbol %in% c("MYC")) %>%
    separate(Tumor_Sample_Barcode, into = c("MMRF", "NUM", "VISIT"), by = "_") %>% 
    mutate(mmrf = str_c(MMRF, NUM, sep = "_")) %>% 
    select(mmrf, VISIT) %>%
    unique() %>% mutate(VISIT = as.numeric(VISIT)) %>%
    left_join(samples_all, by = c("mmrf" = "mmrf", "VISIT" = "visit")) %>% 
    filter(!is.na(tissue_source)) %>% 
    left_join(samples_primary, by = "srr") %>%
    filter(!is.na(mmrf.y)) %>%
    select(mmrf.x) %>%
    rename("mmrf" = "mmrf.x") %>%
    mutate(fusion_labels = "MYC mutation")
  
  n_with_myc_mutation <- mmrf_with_myc_mutation %>% nrow()
  
  myc_pvt1 <- fusions_primary %>% 
    filter((geneA %in% c("PVT1", "MYC") & geneB %in% c("IGH", "IGK", "IGL")) |
             (geneB %in% c("PVT1", "MYC") & geneA %in% c("IGH", "IGK", "IGL"))) %>%
    mutate(fusion_labels = fusion) %>%
    select(mmrf, fusion_labels) %>%
    rbind(mmrf_with_myc_mutation) %>%
    right_join(expr_plot_df %>% filter(gene == "MYC"), by = "mmrf") %>% 
    select(fusion_labels, log10tpm) %>%
    mutate(ig_gene = case_when(str_detect(fusion_labels, pattern = "IGH") ~ "IGH",
                               str_detect(fusion_labels, pattern = "IGK") ~ "IGK",
                               str_detect(fusion_labels, pattern = "IGL") ~ "IGL",
                               TRUE ~ "Neither\nReported"),
           other_gene = case_when(fusion_labels == "MYC mutation" ~ "MYC\nmutation",
                                  str_detect(fusion_labels, pattern = "MYC") ~ "MYC\nfusion",
                                  str_detect(fusion_labels, pattern = "PVT1") ~ "PVT1\nfusion",
                                  TRUE ~ "No mutation\nor fusion"),
           my_alpha = case_when(fusion_labels == "MYC mutation" ~ 1,
                                str_detect(fusion_labels, pattern = "MYC") ~ 1,
                                str_detect(fusion_labels, pattern = "PVT1") ~ 1,
                                TRUE ~ 0.25)) %>%
    mutate(fusion_labels = replace_na(fusion_labels, "Neither\nReported")) %>%
    mutate(other_gene = factor(other_gene, levels = c("No mutation\nor fusion",
                                                      "MYC\nmutation",
                                                      "MYC\nfusion",
                                                      "PVT1\nfusion"),
                               ordered = TRUE),
           ig_gene = factor(ig_gene, levels = c("Neither\nReported",
                                                "IGH",
                                                "IGK",
                                                "IGL"),
                            ordered = TRUE))
  
  max_myc_expr <- ceiling(max(myc_pvt1$log10tpm))
  
  p <- ggplot(myc_pvt1) +
    geom_violin(aes(x = other_gene, y = log10tpm),
                color = "black",
                draw_quantiles = 0.5,
                scale = "width") + 
    geom_jitter(aes(x = other_gene,
                    y = log10tpm,
                    color = other_gene,
                    shape = ig_gene,
                    alpha = my_alpha),
                width = 0.1,
                height = 0) +
    scale_shape_manual(values = c(1, 17, 15, 16)) +
    scale_y_continuous(limits = c(0, max_myc_expr)) +
    scale_alpha(limits = c(0,1) ) +
    labs(x = NULL,
         y = "MYC Expression TPM (log10)",
         color = "Fusion Gene",
         shape = "IG partner gene") +
    theme_bw() +
    scale_color_manual(values = viridis(4, direction = -1)) +
    theme(panel.background = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_text(vjust = 0.5, size = 10),
          axis.text.y = element_text(size = 8),
          legend.position = "bottom",
          legend.direction = "vertical",
          axis.title = element_text(size = 12))
  
  ggsave(str_c(paper_main, "PVT1_MYC.pdf"), 
         p,
         width = 2.5, height = 4.5, useDingbats = FALSE)
  
  ggsave(str_c(paper_main, "PVT1_MYC.no_legend.pdf"),
         p + guides(shape = FALSE, color = FALSE, alpha = FALSE),
         width = 3.5, height = 3.5, useDingbats = FALSE)
  
  # Breakpoints
  plot_df <- fusions_primary %>% 
    filter((geneA %in% c("PVT1", "MYC") & geneB %in% c("IGH", "IGK", "IGL")) |
             (geneB %in% c("PVT1", "MYC") & geneA %in% c("IGH", "IGK", "IGL"))) %>%
    mutate(ig_chr = case_when(geneA %in% c("IGH", "IGK", "IGL") ~ chrA,
                              TRUE ~ chrB),
           ig_pos = case_when(geneA %in% c("IGH", "IGK", "IGL") ~ posA,
                              TRUE ~ posB),
           other_chr = case_when(geneA %in% c("MYC", "PVT1") ~ chrA,
                                 TRUE ~ chrB),
           other_pos = case_when(geneA %in% c("MYC", "PVT1") ~ posA,
                                 TRUE ~ posB)) %>%
    mutate(ig_gene = case_when(str_detect(fusion, pattern = "IGH") ~ "IGH\nchr14 Position (Mb)",
                               str_detect(fusion, pattern = "IGK") ~ "IGK\nchr2 Position (Mb)",
                               str_detect(fusion, pattern = "IGL") ~ "IGL\nchr22 Position (Mb)",
                               TRUE ~ "Neither\nReported"),
           other_gene = case_when(str_detect(fusion, pattern = "MYC") ~ "MYC",
                                  str_detect(fusion, pattern = "PVT1") ~ "PVT1",
                                  TRUE ~ "Neither\nReported")) %>%
    mutate(other_gene = factor(other_gene, levels = c("Neither\nReported",
                                                      "MYC",
                                                      "PVT1"),
                               ordered = TRUE),
           ig_gene = factor(ig_gene, levels = c("Neither\nReported",
                                                "IGH\nchr14 Position (Mb)",
                                                "IGK\nchr2 Position (Mb)",
                                                "IGL\nchr22 Position (Mb)"),
                            ordered = TRUE)) %>%
    select(mmrf, srr, fusion, ig_chr, ig_pos, other_chr, other_pos, ig_gene, other_gene)
  
  chr8_min <- plot_df %>% pull(other_pos) %>% min()/1e6
  chr8_max <- plot_df %>% pull(other_pos) %>% max()/1e6
  
  myc_pvt1_gene_bounds <- tribble(~gene, ~start, ~stop, ~ig_gene,
                                  "MYC", 128747680/1e6, 128753680/1e6, "Neither\nReported",
                                  "PVT1", 128806779/1e6, 129113499/1e6, "Neither\nReported")
  
  plot_df %>% ggplot(aes(x = other_pos/1e6, y = ig_pos/1e6)) +
    coord_cartesian(xlim = c(chr8_min, chr8_max)) +
    facet_wrap(~ ig_gene, ncol = 1, 
               scales = "free_y",
               strip.position = "left") +
    geom_point(aes(color = other_gene, shape = ig_gene), 
               alpha = 0.5,
               show.legend = FALSE) +
    geom_rect(data = myc_pvt1_gene_bounds,
              aes(xmin = start, xmax = stop, ymin = 0.5, ymax = 1),
              inherit.aes = FALSE,
              color = NA,
              fill = "#bdbdbd") +
    geom_rect(data = myc_pvt1_gene_bounds,
              aes(xmin = start, xmax = stop, ymin = -1, ymax = -0.5),
              inherit.aes = FALSE,
              color = NA,
              fill = "#bdbdbd") +
    geom_text(data = myc_pvt1_gene_bounds,
              aes(x = start,
                  y = 0.75,
                  label = gene),
              fontface = "italic",
              hjust = 0,
              nudge_x = 0.001) +
    geom_text(data = myc_pvt1_gene_bounds,
              aes(x = start,
                  y = -0.75,
                  label = gene),
              fontface = "italic",
              hjust = 0,
              nudge_x = 0.001) +
    scale_color_manual(values = viridis(4, direction = -1)[-2], drop = FALSE) +
    scale_shape_manual(values = c(1, 17, 15, 16), drop = FALSE) +
    labs(y = NULL, x = "chr8 Position (Mb)") +
    theme_bw() +
    theme(plot.background = element_blank(),
          panel.background = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.placement = "outside",
          axis.ticks = element_blank(),
          axis.text = element_text(size = 8)) +
    ggsave(str_c(paper_supp, "PVT1_MYC.breakpoints.pdf"),
           width = 7.25, height = 4, useDingbats = FALSE)
}

# ==============================================================================
# Survival MYC--IGL vs. PVT1--IGL
# ==============================================================================
if (TRUE) {
  
  EFS_tibble <- seqfish_clinical_info %>%
    filter(mmrf %in% mmrf_primary_pretreatment) %>%
    filter(!is.na(ISS_Stage), !is.na(EFS_censor), !is.na(Age)) %>%
    left_join(fusions_primary %>% filter(fusion %in% c("PVT1--IGL", "MYC--IGL")), by = "mmrf") %>% 
    select(mmrf, Age, fusion, ISS_Stage, EFS, EFS_censor) %>%
    replace_na(list(fusion = "_None")) %>% 
    mutate(myc_mutant = mmrf %in% mmrf_with_myc_mutation$mmrf,
           ISS_Stage = factor(ISS_Stage, labels = c("I", "II", "III"))) %>%
    filter(mmrf %in% mmrf_with_primary_mutation_calls)
  
  early_tibble <- seqfish_clinical_info %>%
    filter(mmrf %in% mmrf_primary_pretreatment) %>%
    filter(!is.na(ISS_Stage), !is.na(early_relapse_censor), !is.na(Age)) %>%
    left_join(fusions_primary %>% filter(fusion %in% c("PVT1--IGL", "MYC--IGL")), by = "mmrf") %>% 
    select(mmrf, Age, fusion, ISS_Stage, early_relapse_time, early_relapse_censor) %>%
    replace_na(list(fusion = "_None")) %>% 
    mutate(myc_mutant = mmrf %in% mmrf_with_myc_mutation$mmrf,
           ISS_Stage = factor(ISS_Stage, labels = c("I", "II", "III"))) %>%
    filter(mmrf %in% mmrf_with_primary_mutation_calls)
  
  # There are no overlaps in samples with all EFS parameters but could be in general
  plot_survival_list <- list()
  plot_survival_list[["base_EFS"]] <- coxph(formula = Surv(EFS, EFS_censor == 0) ~ ISS_Stage + Age, data = EFS_tibble)
  plot_survival_list[["PVT1_MYC_fusion_EFS"]] <- coxph(formula = Surv(EFS, EFS_censor == 0) ~ ISS_Stage + Age + fusion, data = EFS_tibble)
  plot_survival_list[["anova_PVT1_MYC_fusion_EFS"]] <- anova(plot_survival_list[["base_EFS"]], plot_survival_list[["PVT1_MYC_fusion_EFS"]])
  
  fit_fusion_mutant <- survfit(Surv(EFS, EFS_censor == 0) ~ fusion + myc_mutant, data = EFS_tibble)
  pdf(str_c(paper_main, "PVT1_MYC.EFS.with_legend.pdf"),
      width = 3.5, height = 3.5, useDingbats = FALSE)
  print(ggsurvplot(fit_fusion_mutant, data = EFS_tibble, conf.int = TRUE,
                   surv.median.line = "hv", pval = TRUE,
                   legend.labs = c("Neither Reported", "MYC mutation", "MYC--IGL", "PVT1--IGL"), 
                   legend = "right", 
                   xlab = "Time (days)", 
                   ylab = "Progression-Free Survival Probability",
                   palette = viridis(4, direction = -1),
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
  pdf(str_c(paper_main, "PVT1_MYC.EFS.without_legend.pdf"),
      width = 3.5, height = 3.5, useDingbats = FALSE)
  print(ggsurvplot(fit_fusion_mutant, data = EFS_tibble, conf.int = TRUE,
                   surv.median.line = "hv", pval = TRUE,
                   legend.labs = c("Neither Reported", "MYC mutation", "MYC--IGL", "PVT1--IGL"), 
                   legend = "none",
                   xlab = "Time (days)", 
                   ylab = "Progression-Free Survival Probability",
                   palette = viridis(4, direction = -1),
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
  
  pdf(str_c(paper_supp, "survival.forest.PVT1_MYC.pdf"), width = 3.5, height = 2)
  print(ggforest(plot_survival_list[["PVT1_MYC_fusion_EFS"]], data = EFS_tibble, noDigits = 3))
  dev.off()
}

# ==============================================================================
# Manuscript paragraph info
# ==============================================================================
n_myc_igl_any_stage <- fusions_primary %>% 
  filter(fusion == "MYC--IGL") %>% 
  select(mmrf, fusion) %>% 
  left_join(seqfish_clinical_info %>% 
              select(mmrf, ISS_Stage), by = "mmrf") %>% 
  filter(!is.na(ISS_Stage)) %>% 
  nrow()
n_myc_igl_stageI <- fusions_primary %>% filter(fusion == "MYC--IGL") %>% 
  select(mmrf, fusion) %>% 
  left_join(seqfish_clinical_info %>% 
              select(mmrf, ISS_Stage), by = "mmrf") %>% 
  filter(!is.na(ISS_Stage)) %>% 
  filter(ISS_Stage == "1") %>%
  nrow()

n_pvt1_igl_any_stage <- fusions_primary %>% filter(fusion == "PVT1--IGL") %>% 
  select(mmrf, fusion) %>% 
  left_join(seqfish_clinical_info %>% 
              select(mmrf, ISS_Stage), by = "mmrf") %>% 
  filter(!is.na(ISS_Stage)) %>% 
  nrow()
n_pvt1_igl_stageI <- fusions_primary %>% filter(fusion == "PVT1--IGL") %>% 
  select(mmrf, fusion) %>% 
  left_join(seqfish_clinical_info %>% 
              select(mmrf, ISS_Stage), by = "mmrf") %>% 
  filter(!is.na(ISS_Stage)) %>% 
  filter(ISS_Stage == "1") %>%
  nrow()

print(str_c("Number of samples with MYC mutation: ", mmrf_with_myc_mutation %>% nrow()))

print("Kaplan-Meier estimates for MYC/PVT1/IGL fusion:")
print(fit_fusion_mutant)
print(str_c("Proportion of PVT1--IGL Stage I: ", 
            n_pvt1_igl_stageI, "/", n_pvt1_igl_any_stage, " = ", 
            round(100*n_pvt1_igl_stageI/n_pvt1_igl_any_stage, 2), "%"))
print(str_c("Proportion of MYC--IGL Stage I: ", 
            n_myc_igl_stageI, "/", n_myc_igl_any_stage, " = ", 
            round(100*n_myc_igl_stageI/n_myc_igl_any_stage, 2), "%"))

print(plot_survival_list[["PVT1_MYC_fusion_EFS"]])
print(exp(confint(plot_survival_list[["PVT1_MYC_fusion_EFS"]])))

# ==============================================================================
# Additional investigation based on reviewer response
# ==============================================================================

x <- seqfish_clinical_info %>%
  filter(mmrf %in% mmrf_primary_pretreatment) %>%
  filter(!is.na(ISS_Stage), !is.na(EFS_censor), !is.na(Age)) %>%
  left_join(fusions_primary %>% filter(fusion %in% c("PVT1--IGL", "MYC--IGL")), by = "mmrf") %>% 
  select(mmrf, Age, fusion, ISS_Stage, EFS, EFS_censor) %>%
  replace_na(list(fusion = "_None")) %>% 
  mutate(myc_mutant = mmrf %in% mmrf_with_myc_mutation$mmrf,
               ISS_Stage = factor(ISS_Stage, labels = c("I", "II", "III"))) %>%
  filter(mmrf %in% mmrf_with_primary_mutation_calls) %>%
  left_join(expr_plot_df %>% filter(gene == "MYC") %>% select(mmrf, log10tpm, gene_avg_cnv, cnv_factor), by = "mmrf")
y <- x[complete.cases(x),]
summary(lm(log10tpm ~ fusion + myc_mutant + gene_avg_cnv, data = y))
coxph(formula = Surv(EFS, EFS_censor == 0) ~ ISS_Stage + Age + fusion, data = y)
coxph(formula = Surv(EFS, EFS_censor == 0) ~ ISS_Stage + Age + fusion + gene_avg_cnv, data = y)
anova(coxph(formula = Surv(EFS, EFS_censor == 0) ~ ISS_Stage + Age + fusion, data = y), coxph(formula = Surv(EFS, EFS_censor == 0) ~ ISS_Stage + Age + fusion + gene_avg_cnv, data = y))
anova(coxph(formula = Surv(EFS, EFS_censor == 0) ~ ISS_Stage + Age + fusion +  myc_mutant, data = y), coxph(formula = Surv(EFS, EFS_censor == 0) ~ ISS_Stage + Age + fusion + myc_mutant + gene_avg_cnv, data = y))
lm(log10tpm ~ fusion + myc_mutant + gene_avg_cnv, data = y)

# 2 samples with MYC/PVT1--IGL fusion and IGLL5 mutation
fusions_primary %>% filter(fusion %in% c("MYC--IGL", "PVT1--IGL")) %>% select(mmrf, fusion) %>% left_join(mutation_calls %>% filter(Hugo_Symbol == "IGLL5") %>% mutate(mmrf = str_sub(Tumor_Sample_Barcode, 1, 9)), by = "mmrf") %>% filter(!is.na(Hugo_Symbol))

# Exploration of fusion association with other translocation events
z <- seqfish_clinical_info %>%
  filter(mmrf %in% mmrf_primary_pretreatment) %>%
  filter(!is.na(ISS_Stage), !is.na(EFS_censor), !is.na(EFS), !is.na(Age)) %>%
  left_join(fusions_primary %>% filter(fusion %in% c("PVT1--IGL", "MYC--IGL")), by = "mmrf") %>% 
  select(mmrf, Age, fusion, ISS_Stage, EFS, EFS_censor, del17p, amp1q, updated_seqfish_t_IGH_WHSC1, updated_seqfish_t_IGL_MYC, updated_seqfish_t_IGH_MAF, updated_seqfish_t_IGH_MAFB, SeqWGS_Cp_Hyperdiploid_Call) %>%
  replace_na(list(fusion = "_None")) %>% 
  mutate(myc_mutant = mmrf %in% mmrf_with_myc_mutation$mmrf,
         ISS_Stage = factor(ISS_Stage, labels = c("I", "II", "III"))) %>%
  filter(mmrf %in% mmrf_with_primary_mutation_calls) %>%
  filter(!is.na(del17p) & !is.na(updated_seqfish_t_IGH_WHSC1))

z %>% group_by(fusion, myc_mutant) %>% summarize(n_fusion = n(), total_del17p = sum(del17p, na.rm = TRUE), total_amp1q = sum(amp1q, na.rm = TRUE), total_t414 = sum(updated_seqfish_t_IGH_WHSC1, na.rm = TRUE), total_t822 = sum(updated_seqfish_t_IGL_MYC, na.rm = TRUE), total_t1416 = sum(updated_seqfish_t_IGH_MAF, na.rm = TRUE), total_t1420 = sum(updated_seqfish_t_IGH_MAFB, na.rm = TRUE), total_hyperdiploid = sum(SeqWGS_Cp_Hyperdiploid_Call, na.rm = TRUE))

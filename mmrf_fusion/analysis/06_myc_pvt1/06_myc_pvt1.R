# ==============================================================================
# MYC PVT1 story
# ==============================================================================

paper_main = "paper/main/06_myc_pvt1/"
paper_supp = "paper/supplemental/06_myc_pvt1/"

# Create directories 
dir.create(paper_main, recursive = TRUE, showWarnings = FALSE)
dir.create(paper_supp, recursive = TRUE, showWarnings = FALSE)

# Expression plot data

expr_plot_df <- read_tsv("paper/supplemental/02_expression/expression_plot_tibble.tsv")
expr_plot_df <- expr_plot_df %>% mutate(cnv_factor = factor(categorical_cnv,
                                                            labels = c("DELETION",
                                                                       "Deletion",
                                                                       "Neutral",
                                                                       "Missing",
                                                                       "Amplification",
                                                                       "AMPLIFICATION"), 
                                                            exclude = NULL))

# START HERE

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
           other_gene = case_when(fusion_labels == "MYC mutation" ~ "MYC mutation",
                                  str_detect(fusion_labels, pattern = "MYC") ~ "MYC",
                                  str_detect(fusion_labels, pattern = "PVT1") ~ "PVT1",
                                  TRUE ~ "Neither\nReported"),
           my_alpha = case_when(fusion_labels == "MYC mutation" ~ 1,
                                str_detect(fusion_labels, pattern = "MYC") ~ 1,
                                str_detect(fusion_labels, pattern = "PVT1") ~ 1,
                                TRUE ~ 0.25)) %>%
    mutate(fusion_labels = replace_na(fusion_labels, "Neither\nReported")) %>%
    mutate(other_gene = factor(other_gene, levels = c("Neither\nReported",
                                                      "MYC mutation",
                                                      "MYC",
                                                      "PVT1"),
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
    ggsave(str_c(paper_main, "PVT1_MYC.breakpoints.pdf"),
           width = 7.5, height = 4, useDingbats = FALSE)
}

# SURVIVAL MYC--IGL vs. PVT1--IGL

if (TRUE) {
  
  EFS_tibble <- seqfish_clinical_info %>%
    filter(!is.na(ISS_Stage), !is.na(EFS_censor), !is.na(Age)) %>%
    left_join(fusions_primary %>% filter(fusion %in% c("PVT1--IGL", "MYC--IGL")), by = "mmrf") %>% 
    select(mmrf, Age, fusion, ISS_Stage, EFS, EFS_censor) %>%
    replace_na(list(fusion = "_None")) %>% 
    mutate(myc_mutant = mmrf %in% mmrf_with_myc_mutation$mmrf,
           ISS_Stage = factor(ISS_Stage, labels = c("I", "II", "III"))) %>%
    #mutate(myc_mutant = case_when(myc_mutant == TRUE & fusion != "_None" ~ fusion)) %>%
    filter(mmrf %in% mmrf_with_primary_mutation_calls)
  
  # There are no overlaps in samples with all EFS parameters but could be in general
  
  plot_survival_list[["PVT1_MYC_fusion_EFS"]] <- coxph(formula = Surv(EFS, EFS_censor == 0) ~ ISS_Stage + Age + fusion, data = EFS_tibble)
  plot_survival_list[["anova_PVT1_MYC_fusion_EFS"]] <- anova(plot_survival_list[["base_EFS"]], plot_survival_list[["PVT1_MYC_fusion_EFS"]]) # Significant
  
  fit <- survfit(Surv(EFS, EFS_censor == 0) ~ fusion + myc_mutant, data = EFS_tibble)
  pdf(str_c(paper_main, "PVT1_MYC.EFS.with_legend.pdf"),
      width = 3.5, height = 3.5, useDingbats = FALSE)
  print(ggsurvplot(fit, data = EFS_tibble, conf.int = TRUE,
                   surv.median.line = "hv", pval = TRUE,
                   legend.labs = c("Neither Reported", "MYC mutation", "MYC--IGL", "PVT1--IGL"), 
                   legend = "right", 
                   xlab = "Time (days)", 
                   ylab = "Event-Free Survival Probability",
                   palette = viridis(4, direction = -1),
                   ggtheme = theme_survminer(),
                   conf.int.alpha = 0.1))
  dev.off()
  pdf(str_c(paper_main, "PVT1_MYC.EFS.without_legend.pdf"),
      width = 3.5, height = 3.5, useDingbats = FALSE)
  print(ggsurvplot(fit, data = EFS_tibble, conf.int = TRUE,
                   surv.median.line = "hv", pval = TRUE,
                   legend.labs = c("Neither Reported", "MYC mutation", "MYC--IGL", "PVT1--IGL"), 
                   legend = "none",
                   xlab = "Time (days)", 
                   ylab = "Event-Free Survival Probability",
                   palette = viridis(4, direction = -1),
                   ggtheme = theme_survminer(),
                   conf.int.alpha = 0.1))
  dev.off()
}



# Paragraph info
n_myc_fusions_with_wgs <- fusions_primary %>% 
  filter((geneA %in% c("PVT1", "MYC") & geneB %in% c("IGH", "IGK", "IGL")) |
           (geneB %in% c("PVT1", "MYC") & geneA %in% c("IGH", "IGK", "IGL"))) %>%
  filter(!is.na(n_discordant)) %>% nrow()
n_myc_fusions_validated <- fusions_primary %>%
  filter((geneA %in% c("PVT1", "MYC") & geneB %in% c("IGH", "IGK", "IGL")) |
           (geneB %in% c("PVT1", "MYC") & geneA %in% c("IGH", "IGK", "IGL"))) %>%
  filter(!is.na(n_discordant), n_discordant >= 3) %>% nrow()


print(str_c("MYC/PVT1 IG fusions validated: ", n_myc_fusions_validated, "/", n_myc_fusions_with_wgs, " = ", round(100*n_myc_fusions_validated/n_myc_fusions_with_wgs, 2), "%"))

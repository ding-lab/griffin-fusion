# ==============================================================================
# WHSC1/FGFR3 (MMRF Fusions)
# Steven Foltz (github: envest)
# ==============================================================================

paper_main = "paper/main/04_whsc1_fgfr3/"
paper_supp = "paper/supplementary/04_whsc1_fgfr3/"

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

# WHSC1 transcript info
whsc1_tx <- read_tsv("data/WHSC1_transcripts_UCSC_GenomeBrowser.txt")
whsc1_tx <- whsc1_tx %>% select(start, stop) %>% unique()
min_whsc1 <- whsc1_tx %>% pull(start) %>% min()/1e6
max_whsc1 <- whsc1_tx %>% pull(stop) %>% max()/1e6

# ==============================================================================
# FGFR3 WHSC1 IGH story
# ==============================================================================

if (TRUE) {
  
  get_t414 <- function(delly){
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
      if (all(sort(names(chr_pos_list)) == c("14", "4")) ) {
        return(chr_pos_list)
        break()
      }
    }
    return(list("4" = NA, "14" = NA))
  }
  
  # Co-expression of FGFR3 and WHSC1
  fgfr3_whsc1_ymax <- expr_plot_df %>% filter(gene %in% c("FGFR3", "WHSC1")) %>% 
    pull(log10tpm) %>% max() %>% plyr::round_any(accuracy = 0.1, f = ceiling) 
  plot_expression_2d(expr_plot_df,
                     c("WHSC1", "FGFR3"),
                     fusion1 = "IGH--WHSC1",
                     fusion2 = "IGH--FGFR3",
                     translocation = "updated_seqfish_t_IGH_WHSC1",
                     translocation_formatted = "t(4;14)",
                     ymax_value = fgfr3_whsc1_ymax,
                     pdf_path = str_c(paper_main, "FGFR3_WHSC1_t414.pdf"),
                     pdf_width = 3.5,
                     pdf_height = 3.5,
                     seed = 10, 
                     pretty = TRUE)
  
  # WHSC1 FGFR3 CNV
  p <- bind_cols(expression_primary %>% filter(gene == "FGFR3") %>% 
                   select(mmrf, log10tpm, pct, gene_avg_cnv), 
                 expression_primary %>% filter(gene == "WHSC1") %>% 
                   select(mmrf, log10tpm, pct, gene_avg_cnv)) %>% 
    left_join(fusions_primary %>%
                filter(str_detect(fusion, "IGH") & 
                         (str_detect(fusion, "WHSC1") | 
                            str_detect(fusion, "FGFR3"))), 
              by = "mmrf") %>% 
    filter(!is.na(fusion)) %>%
    mutate(expression_level = case_when(log10tpm < 1 ~ "FGFR3 Low",
                                        TRUE ~ "FGFR3 High")) %>%
    arrange(expression_level) %>%
    ggplot(aes(x = gene_avg_cnv1, y = gene_avg_cnv)) + 
    geom_abline(linetype = 2, color = "grey50") + 
    geom_point(aes(color = expression_level), shape = 16) + 
    coord_equal() +
    theme_bw() +
    labs(color = "FGFR3 Expression", 
         x = "WHSC1 CNV (log2 ratio)",
         y = "FGFR3 CNV (log2 ratio)") +
    scale_x_continuous(limits = c(-1.25, 0.75)) +
    scale_y_continuous(limits = c(-1.25, 0.75)) +
    theme(panel.background = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          plot.background = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_text(size = 8),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 12),
          legend.position = "bottom")
  ggsave(str_c(paper_supp, "fgfr3_whsc1_cnv.pdf"), p,
         width = 5, height = 5, useDingbats = FALSE)
  ggsave(str_c(paper_supp, "fgfr3_whsc1_cnv.no_legend.pdf"), 
         p + guides(color = FALSE),
         width = 3.5, height = 3.5, useDingbats = FALSE)
  
  q <- bind_cols(expression_primary %>% filter(gene == "FGFR3") %>% 
                   select(mmrf, log10tpm, pct, gene_avg_cnv), 
                 expression_primary %>% filter(gene == "WHSC1") %>% 
                   select(mmrf, log10tpm, pct, gene_avg_cnv)) %>% 
    left_join(fusions_primary %>%
                filter(str_detect(fusion, "IGH") & 
                         (str_detect(fusion, "WHSC1") | 
                            str_detect(fusion, "FGFR3"))), 
              by = "mmrf") %>%
    mutate(expression_level = case_when(is.na(fusion) ~ "No fusion\ndetected\n(background)",
                                        log10tpm < 1 ~ "Fusion\nwith low\nFGFR3\nexpression",
                                        TRUE ~ "Fusion\nwith high\nFGFR3\nexpression")) %>%
    arrange(expression_level) %>%
    ggplot(aes(x = fct_rev(expression_level), y = gene_avg_cnv)) + 
    geom_violin(scale = "width",
                draw_quantiles = 0.5) + 
    geom_jitter(aes(color = expression_level), height = 0, width = 0.1, show.legend = FALSE, shape = 18, size = 2) + 
    scale_color_manual(values = c("#D9565C", "#D9565C", "#EDA9AB")) +
    coord_flip() +
    theme_bw() +
    labs(x = NULL, y = "FGFR3 CNV (log2 ratio)") +
    theme(panel.background = element_blank(),
          panel.border = element_blank(),
          plot.background = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_text(size = 8),
          axis.text.y = element_text(hjust = 0.5),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 12),
          legend.position = "bottom",
          panel.grid.minor = element_blank())
  ggsave(str_c(paper_main, "fgfr3_expression_cnv.pdf"), q,
         width = 3.5, height = 3.5, useDingbats = FALSE)
  
  # Fusion and genome breakpoints
  
  plot_df <- fusions_primary %>% 
    filter(fusion == "IGH--WHSC1", !is.na(delly_evidence)) %>% 
    select(mmrf, fusion, srr, LeftBreakpoint, RightBreakpoint, 
           chrA, posA, chrB, posB, delly_evidence) %>% 
    rowwise() %>%
    mutate(chr4_genome_position = get_t414(as.character(delly_evidence))[["4"]],
           chr14_genome_position = get_t414(as.character(delly_evidence))[["14"]]) %>%
    filter(!is.na(chr4_genome_position), !is.na(chr14_genome_position)) %>%
    left_join(expression_primary %>% 
                filter(gene == "FGFR3") %>% 
                select(mmrf, tpm, pct), 
              by = "mmrf") %>%
    mutate(high_expression = case_when(log10(tpm + 1) > 1 ~ "FGFR3 Expression High",
                                       TRUE ~ "FGFR3 Expression Low")) %>%
    ungroup()
  
  n_samples_with_IGHWHSC1_wgs_breakpoint <- plot_df %>% nrow()
  igh_breakpoint_max <- plot_df %>% pull(chr14_genome_position) %>% max()
  igh_breakpoint_min <- plot_df %>% pull(chr14_genome_position) %>% min()
  whsc1_breakpoint_max <- plot_df %>% pull(chr4_genome_position) %>% max()
  whsc1_breakpoint_min <- plot_df %>% pull(chr4_genome_position) %>% min()
  igh_breakpoint_range <- (igh_breakpoint_max - igh_breakpoint_min)/1e6
  whsc1_breakpoint_range <- (whsc1_breakpoint_max - whsc1_breakpoint_min)/1e6
  
  min_chr14 <- min(min(plot_df$posA), min(plot_df$chr14_genome_position))/1e6
  max_chr14 <- max(max(plot_df$posA), max(plot_df$chr14_genome_position))/1e6
  min_chr4 <- min(min(plot_df$posB), min(plot_df$chr4_genome_position))/1e6
  max_chr4 <- max(max(plot_df$posB), max(plot_df$chr4_genome_position))/1e6
  
  p <- ggplot(plot_df, aes(x = posA/1e6, y = posB/1e6)) + 
    geom_point(shape = 3, alpha = 1) +
    facet_wrap(~ high_expression, nrow = 1) +
    labs(x = "IGH Fusion Breakpoint (chr14 Mb)", 
         y = "WHSC1 Fusion Breakpoint (chr4 Mb)") +
    coord_cartesian(ylim = c(1.84, 1.92),
                    xlim = c(105.85, max_chr14)) +
    theme_bw() +
    theme(panel.background = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5),
          strip.background = element_blank(),
          strip.text = element_text(size = 10),
          axis.text = element_text(size = 8),
          axis.title = element_text(size = 10),
          legend.position = "bottom") +
    geom_segment(x = 105.9, xend = 105.9, y = min_whsc1, yend = Inf, color = "#9ecae1") +
    geom_rect(data = whsc1_tx,
              aes(ymin = start/1e6, 
                  ymax = stop/1e6, 
                  xmin = 105.9 - .05,
                  xmax = 105.9 + .05),
              inherit.aes = FALSE,
              color = NA, 
              fill = "#3182bd")
  
  ggsave(str_c(paper_supp, "igh_whsc1_fusion_breakpoints.pdf"), p,
         width = 7.25, height = 3, useDingbats = FALSE)
  
  
  q <- ggplot(plot_df, aes(x = chr14_genome_position/1e6, y = chr4_genome_position/1e6)) +
    geom_point(shape = 3, alpha = 1) +
    facet_wrap(~ high_expression, nrow = 1) +
    labs(x = "IGH Genome Breakpoint (chr14 Mb)", 
         y = "WHSC1 Genome Breakpoint (chr4 Mb)") +
    coord_cartesian(ylim = c(1.84, 1.92),
                    xlim = c(105.85, max_chr14)) +
    theme_bw() +
    theme(panel.background = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5),
          strip.background = element_blank(),
          strip.text = element_text(size = 10),
          axis.text = element_text(size = 8),
          axis.title = element_text(size = 10),
          legend.position = "bottom") +
    geom_segment(x = 105.9, xend = 105.9, y = min_whsc1, yend = Inf, color = "#9ecae1") +
    geom_rect(data = whsc1_tx,
              aes(ymin = start/1e6, 
                  ymax = stop/1e6, 
                  xmin = 105.9 - .05,
                  xmax = 105.9 + .05),
              inherit.aes = FALSE,
              color = NA, 
              fill = "#3182bd")

  ggsave(str_c(paper_supp, "igh_whsc1_genome_breakpoints.pdf"), q,
         width = 7.25, height = 3, useDingbats = FALSE)
  
  r <- ggplot(plot_df, aes(x = chr4_genome_position/1e6, y = posB/1e6)) +
    geom_abline(linetype = 2, color = "grey50") +
    geom_point(shape = 3, alpha = 1) +
    labs(x = "WHSC1 Genome Breakpoint (chr4 Mb)", 
         y = "WHSC1 Fusion Breakpoint (chr4 Mb)") +
    coord_cartesian(ylim = c(1.83, 1.92),
                    xlim = c(1.83, 1.92)) +
    theme_bw() +
    theme(panel.background = element_blank(),
          panel.border = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = 10),
          axis.text = element_text(size = 8),
          axis.title = element_text(size = 10),
          legend.position = "bottom") +
    geom_segment(x = 1.84, xend = 1.84, y = min_whsc1, yend = Inf, color = "#9ecae1") +
    geom_rect(data = whsc1_tx,
              aes(ymin = start/1e6, 
                  ymax = stop/1e6, 
                  xmin = 1.84 - .0025,
                  xmax = 1.84 + .0025),
              inherit.aes = FALSE,
              color = NA, 
              fill = "#3182bd") +
    geom_segment(y = 1.84, yend = 1.84, x = min_whsc1, xend = Inf, color = "#9ecae1") +
    geom_rect(data = whsc1_tx,
              aes(xmin = start/1e6, 
                  xmax = stop/1e6, 
                  ymin = 1.84 - .0025,
                  ymax = 1.84 + .0025),
              inherit.aes = FALSE,
              color = NA, 
              fill = "#3182bd")
  
  ggsave(str_c(paper_supp, "whsc1_genome_fusion_breakpoints.pdf"), r,
         width = 3.5, height = 3.5, useDingbats = FALSE)
  
}

# ==============================================================================
# Special cases of IGH--WHSC1 survival
# ==============================================================================

if (TRUE) {
  
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
  
  plot_survival_list <- list()
  
  FGFR3_pathogenic_mutation <- samples_primary %>% 
    left_join(samples_all, by = "srr") %>% 
    mutate(Tumor_Sample_Barcode = str_c(mmrf.x, "_", visit)) %>% 
    left_join(mutation_calls %>% 
                select(Tumor_Sample_Barcode, Matched_Norm_Sample_Barcode) %>% 
                unique(), 
              by = "Tumor_Sample_Barcode") %>% 
    left_join(mutation_calls %>% 
                filter(Hugo_Symbol == "FGFR3", CLIN_SIG == "pathogenic") %>%
                select(Hugo_Symbol, Tumor_Sample_Barcode) %>% 
                unique(), 
              by = "Tumor_Sample_Barcode") %>%
    filter(Hugo_Symbol == "FGFR3") %>%
    pull(mmrf.x)
  
  FGFR3_expression_high <- expression_primary %>% 
    filter(gene == "FGFR3", log10tpm > 1) %>%
    pull(mmrf)
  
  EFS_tibble <- seqfish_clinical_info %>%
    filter(mmrf %in% mmrf_primary_pretreatment) %>%
    filter(!is.na(ISS_Stage), !is.na(EFS_censor), !is.na(Age)) %>%
    left_join(fusions_primary %>% 
                filter(fusion %in% c("IGH--WHSC1", "IGH--FGFR3", 
                                     "WHSC1--IGH", "FGFR3--IGH")) %>% 
                mutate(fusion = "IGH--WHSC1") %>% 
                select(mmrf, fusion) %>% 
                unique(), 
              by = "mmrf") %>% 
    select(mmrf, Age, fusion, ISS_Stage, EFS, EFS_censor) %>%
    mutate(fgfr3_pathogenic_mutation = mmrf %in% FGFR3_pathogenic_mutation,
           fgfr3_expression_high = mmrf %in% FGFR3_expression_high,
           ISS_Stage = factor(ISS_Stage, labels = c("I", "II", "III"))) %>%
    replace_na(list(fusion = "_None")) %>%
    filter(mmrf %in% mmrf_with_primary_mutation_calls)
  
  n_samples_WHSC1_survival <- EFS_tibble %>% nrow()
  n_samples_WHSC1_survival_fusion <- EFS_tibble %>% 
    filter(fusion == "IGH--WHSC1") %>% nrow()
  n_samples_WHSC1_survival_fusion_expression <- EFS_tibble %>% 
    filter(fusion == "IGH--WHSC1", fgfr3_expression_high == TRUE) %>% nrow()
  n_samples_WHSC1_survival_fusion_expression_pathogenic <- EFS_tibble %>% 
    filter(fusion == "IGH--WHSC1", 
           fgfr3_expression_high == TRUE, 
           fgfr3_pathogenic_mutation == TRUE) %>% nrow()
  
  plot_survival_list[["base_EFS"]] <- coxph(formula = Surv(EFS, EFS_censor == 0) ~ ISS_Stage + Age, data = EFS_tibble)
  plot_survival_list[["WHSC1_fusion_EFS"]] <- coxph(formula = Surv(EFS, EFS_censor == 0) ~ ISS_Stage + Age + fusion, data = EFS_tibble)
  
  plot_survival_list[["with_fusion_base_EFS"]] <- coxph(formula = Surv(EFS, EFS_censor == 0) ~ ISS_Stage + Age, data = EFS_tibble %>% filter(fusion == "IGH--WHSC1"))
  plot_survival_list[["with_fusion_mutation_expression_EFS"]] <- coxph(formula = Surv(EFS, EFS_censor == 0) ~ ISS_Stage + Age + fgfr3_pathogenic_mutation + fgfr3_expression_high, data = EFS_tibble %>% filter(fusion == "IGH--WHSC1"))
  
  plot_survival_list[["anova_WHSC1_fusion_EFS"]] <- anova(plot_survival_list[["base_EFS"]], plot_survival_list[["WHSC1_fusion_EFS"]])
  plot_survival_list[["anova_with_fusion_mutation_expression_EFS"]] <- anova(plot_survival_list[["with_fusion_base_EFS"]], plot_survival_list[["with_fusion_mutation_expression_EFS"]])
  
  fit_km_fusion <- survfit(Surv(EFS, EFS_censor == 0) ~ fusion, data = EFS_tibble)
  
  pdf(str_c(paper_main, "WHSC1.EFS.with_legend.pdf"),
      width = 4.25, height = 3.5, useDingbats = FALSE)
  print(ggsurvplot(fit_km_fusion, data = EFS_tibble, conf.int = TRUE,
                   surv.median.line = "hv", pval = TRUE,
                   legend.labs = c("No Fusion", "IGH--WHSC1 Fusion"),
                   legend = "right", 
                   xlab = "Time (days)", 
                   ylab = "Progression-Free Survival Probability",
                   palette = c("#EDA9AB", "#D9565C"),
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
  
  pdf(str_c(paper_main, "WHSC1.EFS.without_legend.pdf"),
      width = 3.5, height = 3.5, useDingbats = FALSE)
  print(ggsurvplot(fit_km_fusion, data = EFS_tibble, conf.int = TRUE,
                   surv.median.line = "hv", pval = TRUE,
                   legend.labs = c("No Fusion", "IGH--WHSC1 Fusion"),
                   legend = "none",
                   xlab = "Time (days)", 
                   ylab = "Progression-Free Survival Probability",
                   palette = c("#EDA9AB", "#D9565C"),
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
  
  fit_fusion_expr_mut <- survfit(Surv(EFS, EFS_censor == 0) ~ fgfr3_expression_high + fgfr3_pathogenic_mutation, data = EFS_tibble %>% filter(fusion == "IGH--WHSC1"))
  pdf(str_c(paper_main, "WHSC1.FGFR3_expression_mutation.EFS.with_legend.pdf"),
      width = 4.25, height = 3.5, useDingbats = FALSE)
  print(ggsurvplot(fit_fusion_expr_mut, data = EFS_tibble %>% filter(fusion == "IGH--WHSC1"), conf.int = TRUE,
                   surv.median.line = "hv", pval = TRUE,
                   legend.labs = c("FGFR3 Expression Low", "FGFR3 Expression High\nNo Pathogenic Mutation", "FGFR3 Expression High\nPathogenic Mutation"),
                   legend = "right", 
                   xlab = "Time (days)", 
                   ylab = "Progression-Free Survival Probability",
                   palette = rev(c("#1BB6AF", "#088BBE", "#172869")),
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
  
  pdf(str_c(paper_main, "WHSC1.FGFR3_expression_mutation.EFS.without_legend.pdf"),
      width = 3.5, height = 3.5, useDingbats = FALSE)
  print(ggsurvplot(fit_fusion_expr_mut, data = EFS_tibble %>% filter(fusion == "IGH--WHSC1"), conf.int = TRUE,
                   surv.median.line = "hv", pval = TRUE,
                   legend.labs = c("FGFR3 Expression Low", "FGFR3 Expression High\nNo Pathogenic Mutation", "FGFR3 Expression High\nPathogenic Mutation"),
                   legend = "none", 
                   xlab = "Time (days)", 
                   ylab = "Progression-Free Survival Probability",
                   palette = rev(c("#1BB6AF", "#088BBE", "#172869")),
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
# Fusion WHSC1 paragraph output
# ==============================================================================
n_whsc1_fgfr3_fusion <- fusions_primary %>% 
  filter(fusion %in% c("IGH--WHSC1", "IGH--FGFR3", 
                       "WHSC1--IGH", "FGFR3--IGH")) %>% 
  select(mmrf, srr) %>% 
  unique() %>% 
  nrow()
n_whsc1_fgfr3_high <- fusions_primary %>% 
  filter(fusion %in% c("IGH--WHSC1", "IGH--FGFR3", 
                       "WHSC1--IGH", "FGFR3--IGH")) %>% 
  select(mmrf, srr) %>% 
  unique() %>% 
  left_join(expression_primary %>% filter(gene == "FGFR3"), by = "mmrf") %>% 
  filter(log10tpm > 1) %>% 
  nrow()
n_with_high_fgfr3 <- expression_primary %>% 
  filter(gene == "FGFR3", log10tpm > 1) %>% 
  nrow()

n_with_fgfr3_pathogenic_mutation <- length(FGFR3_pathogenic_mutation)
n_with_high_fgfr3_and_mutation_calls <-  sum(FGFR3_expression_high %in% mmrf_with_primary_mutation_calls)

print(str_c("WHSC1 fusions with high FGFR3 expression: ",
            n_whsc1_fgfr3_high, "/", 
            n_whsc1_fgfr3_fusion, " = ", 
            round(100*n_whsc1_fgfr3_high/n_whsc1_fgfr3_fusion, 2), "%"))
print(str_c("Proportion of samples with high FGFR3 expression and pathogenic mutation: ", 
            n_with_fgfr3_pathogenic_mutation, "/", 
            n_with_high_fgfr3_and_mutation_calls, " = ", 
            100*round(n_with_fgfr3_pathogenic_mutation/n_with_high_fgfr3_and_mutation_calls, 4), "%"))

print(str_c("Number of samples with IGH--WHSC1 fusion and WGS breakpoint: ", n_samples_with_IGHWHSC1_wgs_breakpoint))
print(str_c("IGH chr14 genomic breakpoint range: ", round(igh_breakpoint_range, 2)))
print(str_c("WHSC1 chr4 genomic breakpoint range: ", round(whsc1_breakpoint_range, 2)))

print(str_c("Number of samples with complete survival data (ISS_Stage, Age, EFS status): ", n_samples_WHSC1_survival))
print(str_c("Number of EFS samples with IGH--WHSC1 fusion: ", n_samples_WHSC1_survival_fusion))
print(str_c("Number of EFS samples with IGH--WHSC1 fusion and high FGFR3 expression: ", n_samples_WHSC1_survival_fusion_expression))
print(str_c("Number of EFS samples with IGH--WHSC1 fusion and high FGFR3 expression, and pathogenic FGFR3 mutation: ", n_samples_WHSC1_survival_fusion_expression_pathogenic))

print("Cox PH model for fusion:")
print(plot_survival_list[["WHSC1_fusion_EFS"]])
print("Cox PH model for fusion confidence intervals:")
print(exp(confint(plot_survival_list[["WHSC1_fusion_EFS"]])))

print("Kaplan-Meier estimates for WHSC1/FGFR3/IGH fusion:")
print(fit_km_fusion)
print("Kaplan-Meier estimates for low/high expression + mutation among WHSC1/FGFR3/IGH fusion:")
print(fit_fusion_expr_mut)
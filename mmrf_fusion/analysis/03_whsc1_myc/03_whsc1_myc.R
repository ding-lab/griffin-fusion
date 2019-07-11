# ==============================================================================
# WHSC1/FGFR3 and MYC/PVT1 (MMRF Fusions)
# Steven Foltz (smfoltz@wustl.edu), April 2019, copied here July 2019
# ==============================================================================

paper_main = "paper/main/03_whsc1_myc/"
paper_supp = "paper/supplemental/03_whsc1_myc/"

# Create directories 
dir.create(paper_main, recursive = TRUE, showWarnings = FALSE)
dir.create(paper_supp, recursive = TRUE, showWarnings = FALSE)

# Expression plot data

plot_df <- read_tsv("paper/supplemental/02_expression/expression_plot_tibble.tsv")
plot_df <- plot_df %>% mutate(cnv_factor = factor(categorical_cnv,
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
# Updated April 2019
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
  
  fgfr3_mutations <- mutation_calls %>% 
    filter(Hugo_Symbol == "FGFR3") %>% 
    separate(Tumor_Sample_Barcode, 
             by = "_", 
             into = c("MMRF", "mmrf", "sample_n")) %>% 
    filter(sample_n == 1) %>% 
    mutate(mmrf = str_c(MMRF, "_", mmrf)) %>% 
    mutate(vaf = t_alt_count/t_depth)
  
  # Panel A co-expression of FGFR3 and WHSC1
  fgfr3_whsc1_ymax <- plot_df %>% filter(gene %in% c("FGFR3", "WHSC1")) %>% 
    pull(log10tpm) %>% max() %>% plyr::round_any(accuracy = 0.1, f = ceiling) 
  plot_expression_2d(plot_df,
                     c("WHSC1", "FGFR3"),
                     fusion1 = "IGH--WHSC1",
                     fusion2 = "IGH--FGFR3",
                     translocation = "seqfish_Translocation_WHSC1_4_14",
                     translocation_formatted = "t(4;14)",
                     ymax_value = fgfr3_whsc1_ymax,
                     pdf_path = str_c(paper_main, "FGFR3_WHSC1_t414.pdf"),
                     pdf_width = 3.5,
                     pdf_height = 3.5,
                     seed = 10, 
                     pretty = TRUE)
  
  # Panel B look at FGFR3 CNV
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
    geom_violin(scale = "width") + 
    geom_jitter(aes(color = expression_level), height = 0, width = 0.1, show.legend = FALSE, shape = 18, size = 2) + #, alpha = 0.5) +
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
         width = 7.25, height = 3)
  
  
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
         width = 7.25, height = 3)
  
}


# ==============================================================================
# MYC PVT1 story
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
  
  # MYC and PVT expression
  myc_pvt1 <- fusions_primary %>% 
    filter(geneA %in% c("PVT1", "MYC"), geneB == "IGL") %>% 
    group_by(srr) %>% 
    summarize(fusion_labels = str_c(fusion, collapse = "\n")) %>% 
    right_join(plot_df %>% filter(gene == "MYC"), by = "srr") %>% 
    select(fusion_labels, log10tpm) %>%
    mutate(jitter_fusion_labels = 
             jitter(as.numeric(!is.na(fusion_labels)) + 1, factor = 1),
           fusion_labels = replace_na(fusion_labels, "Neither Reported")) %>%
    mutate(fusion_labels = factor(fusion_labels, levels = c("MYC--IGL", 
                                                            "PVT1--IGL", 
                                                            "Neither Reported"), 
                                  ordered = TRUE))
  
  max_myc_expr <- ceiling(max(myc_pvt1$log10tpm))
  
  p <- ggplot(myc_pvt1) +
    coord_flip(expand = c(0.01, 0.01)) +
    #geom_violin(aes(x = str_detect(fusion_labels, "IGL"), y = log10tpm),
    #            color = NA,
    #            fill = "black",
    #            alpha = 0.25) +
    geom_violin(aes(x = str_detect(fusion_labels, "IGL"), y = log10tpm),
                color = "black",
                draw_quantiles = 0.5,
                scale = "width") + 
    geom_point(aes(x = jitter_fusion_labels,
                   y = log10tpm,
                   color = fusion_labels),
               #alpha = 0.25,
               shape = 16) +
    scale_color_manual(values = c("#1BB6AF", "#088BBE", "#172869")) +
    scale_y_continuous(limits = c(0, max_myc_expr), position = "right") +
    labs(x = NULL,
         y = "MYC Expression TPM (log10)",
         color = "Fusion Gene") +
    theme_bw() +
    theme(panel.background = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_blank(),
          legend.position = "bottom",
          legend.direction = "vertical",
          axis.text.x = element_text(size = 8),
          axis.title = element_text(size = 12))
  
  ggsave(str_c(paper_main, "PVT1_MYC.pdf"), p, width = 8, height = 2, useDingbats = FALSE)
  ggsave(str_c(paper_main, "PVT1_MYC.no_legend.pdf"), 
         p + guides(shape = FALSE, color = FALSE),
         width = 3.5, height = 1.75, useDingbats = FALSE)
  
  # Breakpoints
  plot_df <- fusions_primary %>% 
    filter(geneA %in% c("MYC", "PVT1"),
           geneB %in% c("IGH", "IGL", "IGK")) %>% #fusion %in% c("MYC--IGL", "PVT1--IGL")) %>% #, !is.na(delly_evidence)) %>% 
    select(mmrf, fusion, srr, LeftBreakpoint, RightBreakpoint, 
           chrA, posA, chrB, posB) #%>% #, delly_evidence) %>% 
    #rowwise() %>%
    #mutate(chr8_genome_position = get_t822(as.character(delly_evidence))[["8"]],
    #       chr22_genome_position = get_t822(as.character(delly_evidence))[["22"]]) %>%
    #filter(!is.na(chr8_genome_position), !is.na(chr22_genome_position)) %>%
    #left_join(expression_primary %>% 
    #            filter(gene == "FGFR3") %>% 
    #            select(mmrf, tpm, pct), 
    #          by = "mmrf") %>%
    #mutate(high_expression = case_when(log10(tpm + 1) > 1 ~ "FGFR3 Expression High",
    #                                   TRUE ~ "FGFR3 Expression Low")) %>%
    #ungroup()
  
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
         width = 7.25, height = 3)
  
  
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
         width = 7.25, height = 3)
}

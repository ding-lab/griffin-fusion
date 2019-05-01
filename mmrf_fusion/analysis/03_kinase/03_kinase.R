# ==============================================================================
# Kinases (MMRF Fusions)
# Steven Foltz (smfoltz@wustl.edu), April 2019
# ==============================================================================

paper_main = "paper/main/03_kinase/"
paper_supp = "paper/supplemental/03_kinase/"

# Create directories 
dir.create(paper_main, recursive = TRUE, showWarnings = FALSE)
dir.create(paper_supp, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# Types of kinases
# Written April 2019 -- Main
# ==============================================================================

if (TRUE) {
  p <- kinases %>% group_by(kinase_group_full_name, KinasePos) %>% summarize(count = n()) %>% ungroup() %>%
    ggplot(aes(x = fct_reorder(kinase_group_full_name, count), fill = KinasePos, y = count)) +
    #geom_bar(stat = "identity", position = "dodge") +
    geom_col(position = "dodge") +
    coord_flip(expand = c(0,0)) +
    labs(y = "Fusion Count", x = "Kinase Group", fill = "Kinase Position") +
    scale_y_continuous(position = "right") +
    scale_fill_manual(values = c("#addd8e", "#31a354"),
                      breaks = c("5P_KINASE", "3P_KINASE"),
                      labels = c("5' Kinase", "3' Kinase")) +
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
  
  ggsave(str_c(paper_main, "kinase_groups.with_legend.pdf"), p, 
         width = 7.25, height = 7.25/1.618, useDingbats = FALSE)
  ggsave(str_c(paper_main, "kinase_groups.without_legend.pdf"), 
         p + guides(fill = FALSE), 
         width = 7.25, height = 7.25/1.618, useDingbats = FALSE)
  
}

# ==============================================================================
# Expression correlation of kinases by structure (3' or 5', in-frame or not)
# Written April 2019
# ==============================================================================

if (FALSE) {
  
  cor_tibble <- kinases %>% 
    mutate(kinase_expression = case_when(KinasePos == "5P_KINASE" ~ geneA_pct,
                                         KinasePos == "3P_KINASE" ~ geneB_pct)) %>%
    group_by(KinasePos, KinaseDomain) %>% 
    summarize(c = cor(geneA_pct, geneB_pct, use = "pairwise.complete.obs")) %>% 
    ungroup() %>% mutate(KinasePos = case_when(KinasePos == "5P_KINASE" ~ "5' Kinase",
                                               KinasePos == "3P_KINASE" ~ "3' Kinase")) %>%
    mutate(KinaseDomain = case_when(KinaseDomain == "Intact" ~ "Domain Intact",
                                    KinaseDomain == "Disrupted" ~ "Domain Disrupted"))
  
  kinases %>% mutate(KinasePos = case_when(KinasePos == "5P_KINASE" ~ "5' Kinase",
                                           KinasePos == "3P_KINASE" ~ "3' Kinase")) %>%
    mutate(KinaseDomain = case_when(KinaseDomain == "Intact" ~ "Domain Intact",
                                    KinaseDomain == "Disrupted" ~ "Domain Disrupted")) %>%
    ggplot(aes(x = geneA_pct, y = geneB_pct)) + 
    geom_smooth(method = "lm") + 
    geom_point(shape = 16) + 
    coord_fixed() + 
    scale_y_continuous(expand = c(0.03, 0.03)) +
    scale_x_continuous(expand = c(0.03, 0.03)) +
    geom_text(data = cor_tibble, aes(x = 0.5, y = 0, label = str_c("cor = ", round(c, 2))), vjust = 1, color = "blue", size = 2.5) +
    facet_grid(fct_reorder(KinasePos, desc(KinasePos)) ~ KinaseDomain) +
    labs(x = "5' Gene Expression Percentile", y = "3' Expression Percentile") +
    theme_bw() +
    theme(panel.background = element_blank(),
          panel.border = element_blank(),
          strip.background = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 8),
          axis.ticks = element_blank()) +
    ggsave(str_c(paper_main, "kinase_expression_correlation.pdf"),
           height = 3.25, width = 3.25, useDingbats = FALSE)
}

# ==============================================================================
# Expression correlation of kinases by structure (3' or 5', in-frame or not)
# Intact only
# Written April 2019 -- Main
# ==============================================================================

if (TRUE) {
  
  cor_tibble <- kinases %>% 
    mutate(kinase_expression = case_when(KinasePos == "5P_KINASE" ~ geneA_pct,
                                         KinasePos == "3P_KINASE" ~ geneB_pct)) %>%
    group_by(KinasePos, KinaseDomain) %>% 
    summarize(c = cor(geneA_pct, geneB_pct, use = "pairwise.complete.obs")) %>% 
    ungroup() %>% mutate(KinasePos = case_when(KinasePos == "5P_KINASE" ~ "5' Kinase",
                                               KinasePos == "3P_KINASE" ~ "3' Kinase")) %>%
    mutate(KinaseDomain = case_when(KinaseDomain == "Intact" ~ "Domain Intact",
                                    KinaseDomain == "Disrupted" ~ "Domain Disrupted")) %>%
    filter(KinaseDomain == "Domain Intact")
  
  kinases %>% mutate(KinasePos = case_when(KinasePos == "5P_KINASE" ~ "5' Kinase",
                                           KinasePos == "3P_KINASE" ~ "3' Kinase")) %>%
    mutate(KinaseDomain = case_when(KinaseDomain == "Intact" ~ "Domain Intact",
                                    KinaseDomain == "Disrupted" ~ "Domain Disrupted")) %>%
    filter(KinaseDomain == "Domain Intact") %>%
    ggplot(aes(x = geneA_pct, y = geneB_pct)) + 
    geom_smooth(method = "lm") + 
    geom_point(shape = 16) + 
    coord_fixed() + 
    scale_y_continuous(expand = c(0.03, 0.03)) +
    scale_x_continuous(expand = c(0.03, 0.03)) +
    geom_text(data = cor_tibble, aes(x = 0.5, y = 0, label = str_c("cor = ", round(c, 2))), vjust = 1, color = "blue", size = 2.5) +
    facet_grid(fct_reorder(KinasePos, desc(KinasePos)) ~ KinaseDomain) +
    labs(x = "5' Gene Expression Percentile", y = "3' Expression Percentile") +
    theme_bw() +
    theme(panel.background = element_blank(),
          panel.border = element_blank(),
          strip.background = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 8),
          axis.ticks = element_blank()) +
    ggsave(str_c(paper_main, "kinase_expression_correlation.intact.pdf"),
           height = 5, width = 3, useDingbats = FALSE)
}

# ==============================================================================
# Expression correlation of geneA and geneB for 3' intact kinase fusions
# Written April 2019 -- Supplemental
# ==============================================================================

if (TRUE) {
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
           width = 6, height = 6, useDingbats = FALSE)  
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
                           2490, "ARHGEF2--NTRK1",          "RhoGEF",           3,               2,    239,	 431, "domain",
                           2490, "ARHGEF2--NTRK1",              "PH",           3,               2,    474,   570, "domain",
                           2490, "ARHGEF2--NTRK1", "Tyrosine Kinase",           4,               2,   1023,	1292, "domain") %>%
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
    geom_text(data = structure_tbl %>% filter(class == "gene"), aes(x = start + (stop - start)/2, y = as.numeric(fusion) - 0.075, color = fill_color), size = 4, show.legend = FALSE, fontface = "italic") +
    geom_text(data = structure_tbl %>% filter(class == "domain"), aes(x = start + (stop - start)/2, y = as.numeric(fusion) + 0.125), size = 4, show.legend = FALSE, color = "white") +
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
    summarize(count = n()) %>% 
    left_join(pancan_fusions %>% 
                separate(Fusion, into = c("geneA", "geneB"), sep = "--"), 
              by = "geneB") %>% 
    group_by(geneB, Cancer) %>% 
    summarize(count2 = n()) %>%
    filter(!is.na(Cancer))
  
  ggplot(data = plot_df, aes(y = geneB, x = Cancer)) + 
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
                 width = 7.5, height = 9)
  
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

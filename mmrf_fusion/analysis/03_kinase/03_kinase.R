# ==============================================================================
# Overview (MMRF Fusions)
# Steven Foltz (smfoltz@wustl.edu), April 2019
# ==============================================================================

paper_main = "paper/main/03_kinase/"
paper_supp = "paper/supplemental/03_kinase/"

# Create directories 
dir.create(paper_main, recursive = TRUE, showWarnings = FALSE)
dir.create(paper_supp, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# Types of kinases
# Written April 2019
# ==============================================================================

if (TRUE) {
  p <- kinases %>% group_by(kinase_group_full_name, KinasePos) %>% summarize(count = n()) %>% ungroup() %>%
    ggplot(aes(x = fct_reorder(kinase_group_full_name, count), fill = KinasePos, y = count)) +
    geom_bar(stat = "identity", position = "dodge") +
    coord_flip(expand = c(0,0)) +
    labs(y = "Fusion Count", x = "Kinase Group", fill = "Kinase Position") +
    scale_y_continuous(position = "right") +
    scale_fill_brewer(palette = "Reds",
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
  
  ggsave(str_c(paper_main, "kinase_groups.with_legend.pdf"), p, width = 8, height = 4)
  ggsave(str_c(paper_main, "kinase_groups.without_legend.pdf"), p + guides(fill = FALSE), width = 8, height = 4)
  
}

# ==============================================================================
# Expression correlation of kinases by structure (3' or 5', in-frame or not)
# Written April 2019
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
           height = 4, width = 4, useDingbats = FALSE)
}

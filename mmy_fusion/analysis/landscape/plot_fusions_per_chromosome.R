# ==============================================================================
# Create a plot showing number of fusions per chromosome
# Steven Foltz (smfoltz@wustl.edu), August 2018
# ==============================================================================

# ==============================================================================
# Create plot for all chr x chr figure
# ==============================================================================
fusions_primary %>% select(mmrf, fusion, chrA, chrB) %>% 
  left_join(seqfish_clinical_info, by = "mmrf") %>% 
  ggplot(aes(x = factor(seqfish_Hyperdiploidy))) + geom_bar() + 
  facet_grid(chrA ~ chrB) + 
  labs(x = "Hyperdiploid Category", y = "Number of fusions") +
  ggplot2_standard_additions()
ggsave("analysis/fusion_summaries/fusions_per_chromosome.chr_by_chr.pdf", 
       device = "pdf", width = 20, height = 20)

# ==============================================================================
# Create plot for chrA = chrB figure
# ==============================================================================
fusions_primary %>% select(mmrf, fusion, chrA, chrB) %>% filter(chrA == chrB) %>%
  left_join(seqfish_clinical_info, by = "mmrf") %>% 
  ggplot(aes(x = factor(seqfish_Hyperdiploidy))) + geom_bar() + 
  facet_wrap(~ chrA, ncol = 6) + 
  labs(x = "Hyperdiploid Category", y = "Number of fusions") +
  ggplot2_standard_additions()
ggsave("analysis/fusion_summaries/fusions_per_chromosome.chr_eq_chr.pdf", 
       device = "pdf", width = 10, height = 10)
  
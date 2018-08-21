# ==============================================================================
# Create a plot number of fusions per sample
# Steven Foltz (smfoltz@wustl.edu), August 2018
# ==============================================================================

# ==============================================================================
# Function to obtain number of fusions per sample
# ==============================================================================
n_fusions_per_sample <- function(fusions_df){
  tibble( mmrf = names(table(fusions_df$mmrf)), n_fusions = as.integer(unname(table(fusions_df$mmrf))) )
}

# ==============================================================================
# Create data frame for plotting
# ==============================================================================
n_hdp <- seqfish_clinical_info %>% filter(seqfish_Hyperdiploidy == 1) %>% nrow()
n_not_hdp <- seqfish_clinical_info %>% filter(seqfish_Hyperdiploidy == 0) %>% nrow()
n_na_hdp <- seqfish_clinical_info %>% filter(is.na(seqfish_Hyperdiploidy)) %>% nrow()
hpd_key <- tribble(~seqfish_Hyperdiploidy, ~hyperdiploid_categories, ~count,
               0, str_c("Not Hyperdiploid (", n_not_hdp, ")"), n_not_hdp,
               1, str_c("Hyperdiploid (", n_hdp, ")"), n_hdp,
               NA, str_c("Not Available (", n_na_hdp, ")"), n_na_hdp
)

plot_df <- seqfish_clinical_info %>% select(mmrf, seqfish_Hyperdiploidy) %>% 
  left_join(n_fusions_per_sample(fusions_primary), by = "mmrf") %>%
  mutate_at(3, funs(replace(., is.na(.), 0))) %>%
  left_join(hpd_key, by = "seqfish_Hyperdiploidy")

# ==============================================================================
# Plot number of fusions per sample
# ==============================================================================
plot_df %>% ggplot(aes(x = n_fusions)) + geom_histogram(binwidth = 1, center = 0) + 
  facet_wrap(~ fct_reorder(hyperdiploid_categories, -count), ncol = 1) +
  labs(x = "Number of fusions detected", y = "Number of samples") +
  ggplot2_standard_additions()
ggsave("analysis/fusion_summaries/histogram_n_fusions_per_sample.pdf", device = "pdf", 
       width = 10, height = 10)

# ==============================================================================
# Plot frequency of number of fusions per sample
# ==============================================================================
plot_df %>% ggplot(aes(x = n_fusions, y = ..density..)) + 
  geom_freqpoly(aes(color = fct_reorder(hyperdiploid_categories, -count)), 
                binwidth = 1, center = 0) + 
  labs(x = "Number of fusions detected", y = "Proportion of samples",
       color = "Hyperdiploid Category") +
  ggplot2_standard_additions()
ggsave("analysis/fusion_summaries/freqpoly_n_fusions_per_sample.pdf", device = "pdf", 
       width = 10, height = 10)

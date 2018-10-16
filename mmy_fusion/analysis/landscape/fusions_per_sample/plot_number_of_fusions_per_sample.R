# ==============================================================================
# Create a plot of number of fusions per sample
# Steven Foltz (smfoltz@wustl.edu), August 2018
# ==============================================================================

plot_dir = "analysis/landscape/fusions_per_sample/"

# ==============================================================================
# Create data frame for plotting
# ==============================================================================
n_hdp <- seqfish_clinical_info %>% 
  filter(seqfish_Hyperdiploidy == 1) %>% nrow()
n_not_hdp <- seqfish_clinical_info %>% 
  filter(seqfish_Hyperdiploidy == 0) %>% nrow()
n_na_hdp <- seqfish_clinical_info %>% 
  filter(is.na(seqfish_Hyperdiploidy)) %>% nrow()
hpd_key <- tribble(~seqfish_Hyperdiploidy, ~hyperdiploid_categories, ~count,
               0, str_c("Not Hyperdiploid (", n_not_hdp, ")"), n_not_hdp,
               1, str_c("Hyperdiploid (", n_hdp, ")"), n_hdp,
               NA, str_c("Not Available (", n_na_hdp, ")"), n_na_hdp
)

plot_df <- seqfish_clinical_info %>% select(mmrf, seqfish_Hyperdiploidy) %>% 
  left_join(fusions_primary, by = "mmrf") %>%
  group_by(mmrf, seqfish_Hyperdiploidy, srr) %>% 
  summarize(n = n()) %>% ungroup() %>%
  mutate(n_fusions = n - is.na(srr)) %>%
  left_join(hpd_key, by = "seqfish_Hyperdiploidy")

# ==============================================================================
# Plot number of fusions per sample
# ==============================================================================
plot_df %>% 
  ggplot(aes(x = n_fusions)) + 
  geom_histogram(binwidth = 1, center = 0) + 
  facet_wrap(~ fct_reorder(hyperdiploid_categories, -count), ncol = 1) +
  labs(x = "Number of fusions detected", y = "Number of samples") +
  ggplot2_standard_additions()
ggsave(str_c(plot_dir, "histogram_n_fusions_per_sample.pdf"), device = "pdf", 
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
ggsave(str_c(plot_dir, "freqpoly_n_fusions_per_sample.pdf"), device = "pdf", 
       width = 10, height = 10)

# ==============================================================================
# Median number of fusions per category
# ==============================================================================
n_fusion_tibble <- plot_df %>% group_by(hyperdiploid_categories) %>% 
  summarize(min_n_fusions = min(n_fusions), 
            median_n_fusions = median(n_fusions), 
            mean_n_fusions = mean(n_fusions), 
            max_n_fusions = max(n_fusions))

write_tsv(n_fusion_tibble, str_c(plot_dir, "n_fusion_table.txt"), 
          na = "NA", append = FALSE, col_names = TRUE)

# ==============================================================================
# Characteristics of samples with large number of fusions
# Not really a defined pattern that separtes outliers from non-outliers
# ==============================================================================

outlier_n_fusions <- fusions_primary %>% group_by(srr) %>% 
  summarize(n = n()) %>% 
  mutate(z = (n - mean(n))/sd(n) ) %>% 
  filter(abs(z) > 2) %>% pull(n) %>% min()

mmrf_with_large_n <- plot_df %>% 
  filter(n_fusions >= outlier_n_fusions) %>% pull(mmrf)

mean_igl_nonoutliers <- mean(fusions_primary %>%
                               filter(!(mmrf %in% mmrf_with_large_n)) %>% 
                               mutate(geneB_IGL = geneB == "IGL@") %>% 
                               pull(geneB_IGL))

mean_igl_outliers <- mean(fusions_primary %>%
                            filter(mmrf %in% mmrf_with_large_n) %>%
                            mutate(geneB_IGL = geneB == "IGL@") %>% 
                            pull(geneB_IGL))

# ==============================================================================
# Characteristics of samples with zero fusions
# ==============================================================================

mmrf_zero_fusions <- samples_primary %>% 
  left_join(fusions_primary, by = "mmrf") %>%
  filter(is.na(srr.y)) %>% pull(mmrf)

# ==============================================================================
# t-test to compare number of fusions between hyperdiploidy and not
# NOT significant difference by t-test
# ==============================================================================

n_fusions_with_hyperdiploid_info <- seqfish_clinical_info %>% 
  left_join(fusions_primary, by = "mmrf") %>% 
  filter(!is.na(seqfish_Hyperdiploidy)) %>% 
  group_by(mmrf, srr, seqfish_Hyperdiploidy)  %>% 
  summarize(n = n()) %>% 
  ungroup() %>% 
  mutate(n_corrected = n - is.na(srr))

n_fusions_hyperdiploid <- n_fusions_with_hyperdiploid_info %>% 
  filter(seqfish_Hyperdiploidy == 1) %>% pull(n_corrected)
n_fusions_nonhyperdiploid <- n_fusions_with_hyperdiploid_info %>% 
  filter(seqfish_Hyperdiploidy == 0) %>% pull(n_corrected)
t.test(n_fusions_hyperdiploid, n_fusions_nonhyperdiploid)

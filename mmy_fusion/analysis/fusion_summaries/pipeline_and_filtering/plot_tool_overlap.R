# ==============================================================================
# Create a plot showing the overlap/agreement between tools
# Steven Foltz (smfoltz@wustl.edu), September 2018
# ==============================================================================

plot_dir = "analysis/fusion_summaries/pipeline_and_filtering/"

# ==============================================================================
# Functions
# ==============================================================================

tool_key = str_c("E = EricScript",
                 "F = FusionCatcher",
                 "I = INTEGRATE",
                 "P = PRADA",
                 "S = STAR-Fusion",
                 sep = "\n")

# ==============================================================================
# Plot tool overlap (bar chart)
# ==============================================================================

fusions_primary %>% 
  select(Callers) %>% group_by(Callers) %>% count() %>% 
  ggplot(aes(x = fct_reorder(Callers, n), y = n)) + 
  geom_bar(stat = "identity") +
  labs(x = "Fusion Detection Tool Combination", y = "Number of Fusions") +
  coord_flip() + 
  geom_label(aes(label = tool_key, x = 1, y = 150), hjust = 0, vjust = 0) +
  ggplot2_standard_additions()

ggsave(str_c(plot_dir, "tool_overlap.bar_chart.pdf"), device = "pdf", 
       width = 10, height = 10)

# ==============================================================================
# Plot tool overlap (upsetr)
# ==============================================================================
library(UpSetR)

upsetr_df <- data.frame(fusions_primary %>% select(starts_with("called_by")))
names(upsetr_df) <- c("EricScript", "FusionCatcher", "INTEGRATE", 
"PRADA", "STAR-Fusion")
pdf(file = str_c(plot_dir, "tool_overlap.upsetr.pdf"), 
    width = 20, height = 10)
upset(upsetr_df, nsets = ncol(upsetr_df), nintersects = NA, order.by = "freq",
      text.scale = 2, point.size = 3, line.size = 1)
dev.off()

# ==============================================================================
# Plot n_fusions before and after soft filtering
# ==============================================================================

kept_fusions_per_srr_primary <- fusions_primary %>% group_by(srr) %>% count()
plot_df <- fusions_hard_primary %>% group_by(Sample) %>% count() %>% 
  rename(srr = Sample) %>% 
  left_join(kept_fusions_per_srr_primary, by = "srr") %>% 
  replace_na(list(n.y = 0)) %>% 
  rename(n_fusions_hard = n.x, n_fusions_soft = n.y) %>%
  mutate(difference = n_fusions_hard - n_fusions_soft)
max_n_fusions <- plot_df %>% pull(n_fusions_hard) %>% max() %>% round(-1)

# Scatter plot comparing number of variants after hard and soft filtering
plot_df %>% ungroup() %>% 
  mutate(hard_round = plyr::round_any(n_fusions_hard, 5, floor), 
         soft_round = plyr::round_any(n_fusions_soft, 5, floor)) %>% 
  group_by(hard_round, soft_round) %>% count() %>% 
  ggplot(aes(x = hard_round, y = soft_round, 
             fill = factor(plyr::round_any(n,10, floor)))) + 
  geom_tile(color = "black") + 
  geom_text(aes(label = n, color = n >= 40)) +
  scale_color_manual(values = c("black", "white")) +
  scale_fill_brewer() + 
  coord_equal() +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  labs(x = "Number of Fusions After Hard Filtering",
       y = "Number of Fusions After Soft Filtering",
       title = "Number of Fusions During Filtering Process",
       fill = "Samples") +
  guides(color = FALSE) +
  ggplot2_standard_additions()
ggsave(str_c(plot_dir, "hard_soft_filtering.tile.pdf"), 
       width = 15, height = 10)

# Histogram of number of variants after hard filtering
median_value <- plot_df %>% pull(n_fusions_hard) %>% median()
plot_df %>% ggplot(aes(x = n_fusions_hard)) +
  geom_histogram(binwidth = 5, center = 2.5) + 
  xlim(0, max_n_fusions) +
  geom_vline(xintercept = median_value, color = "white", linetype = 2, size = 1) +
  annotate("text", x = median_value + 1, y = 0 + 1, angle = 90, color = "white",
           hjust = 0, vjust = 1, label = str_c("median = ", median_value), size = 5) +
  labs(x = "Number of Fusions After Hard Filtering", y = "Number of Samples") +
  ggplot2_standard_additions()
ggsave(str_c(plot_dir, "hard_filtering.bar.pdf"),
       width = 15, height = 10)

# Histogram of number of variants after soft filtering
median_value <- plot_df %>% pull(n_fusions_soft) %>% median()
plot_df %>% ggplot(aes(x = n_fusions_soft)) +
  geom_histogram(binwidth = 1, center = 0) + 
  xlim(0, max_n_fusions) +
  geom_vline(xintercept = median_value, color = "white", linetype = 2, size = 1) +
  annotate("text", x = median_value + 1, y = 0 + 1, angle = 90, color = "white",
           hjust = 0, vjust = 1, label = str_c("median = ", median_value), size = 5) +
  labs(x = "Number of Fusions After Soft Filtering", y = "Number of Samples") +
  ggplot2_standard_additions()
ggsave(str_c(plot_dir, "soft_filtering.bar.pdf"),
       width = 15, height = 10)

median_value <- plot_df %>% pull(n_fusions_soft) %>% median()
plot_df %>% ggplot(aes(x = n_fusions_soft)) +
  geom_histogram(binwidth = 1, center = 0) + 
  geom_vline(xintercept = median_value, color = "white", linetype = 2, size = 1) +
  annotate("text", x = median_value + 1, y = 0 + 1, angle = 90, color = "white",
           hjust = 0, vjust = 1, label = str_c("median = ", median_value), size = 5) +
  labs(x = "Number of Fusions After Soft Filtering", y = "Number of Samples") +
  ggplot2_standard_additions()
ggsave(str_c(plot_dir, "soft_filtering.bar2.pdf"),
       width = 15, height = 10)

# Histogram of difference in number of variants after hard and soft filtering
median_value <- plot_df %>% pull(difference) %>% median()
plot_df %>% ggplot(aes(x = difference)) +
  geom_histogram(binwidth = 5, center = 2.5) + 
  xlim(0, max_n_fusions) +
  geom_vline(xintercept = median_value, color = "white", linetype = 2, size = 1) +
  annotate("text", x = median_value + 1, y = 0 + 1, angle = 90, color = "white",
           hjust = 0, vjust = 1, label = str_c("median = ", median_value), size = 5) +
  labs(x = "Difference in Number of Fusions After Hard and Soft Filtering", 
       y = "Number of Samples") +
  ggplot2_standard_additions()
ggsave(str_c(plot_dir, "difference.bar.pdf"),
       width = 15, height = 10)

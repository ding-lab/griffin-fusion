# ==============================================================================
# Create a plot showing the overlap/agreement between tools
# Steven Foltz (smfoltz@wustl.edu), September 2018
# ==============================================================================

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

ggsave("analysis/fusion_summaries/tool_overlap.bar_chart.pdf", device = "pdf", 
       width = 10, height = 10)

# ==============================================================================
# Plot tool overlap (upsetr)
# ==============================================================================
library(UpSetR)

plot_df <- fusions_primary %>% select(starts_with("called_by")) %>% data.frame()
names(plot_df) <- c("EricScript", "FusionCatcher", "INTEGRATE", 
"PRADA", "STAR-Fusion")
pdf(file = "analysis/fusion_summaries/tool_overlap.upsetr.pdf", 
    width = 10, height = 10)
upset(plot_df, order.by = "freq")
dev.off()

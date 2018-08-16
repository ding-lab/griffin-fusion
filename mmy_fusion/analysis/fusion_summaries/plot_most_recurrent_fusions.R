# ==============================================================================
# Create a plot showing most recurrent fusions
# Steven Foltz (smfoltz@wustl.edu), August 2018
# ==============================================================================

# ==============================================================================
# Function to define IG and non-IG fusions
# ==============================================================================
return_fusion_ig <- function(fusion_vector){
  define_ig_fusion <- function(fusion_string){
    genes = str_split(fusion_string, pattern = "--", simplify = TRUE)
    if ("IGH@" %in% genes) {
      return("IGH fusion")
    } else if ("IGL@" %in% genes) {
      return("IGL fusion")
    } else if ("IGK@" %in% genes) {
      return("IGK fusion")
    } else{
      return("non-IG fusion")
    }
  }
  apply(as.matrix(fusion_vector), 1, define_ig_fusions)
}

# ==============================================================================
# Create data frame for plotting
# ==============================================================================
plot_df <- fusions_primary %>% select(mmrf:fusion, geneA, geneB) %>% 
  mutate(fusion_ig = return_fusion_ig(fusion))

# Make list of top genes for inclusion in overall plot
n_fusions <- 10 # number of fusions to keep in plot
topX_overall <- tail(names(sort(table(plot_df$fusion))), n = n_fusions)

# Make list of top genes for inclusion in each IG sub-panel
n_fusions_by_ig <- 5 # number of fusions to keep in each panel
topX_by_ig <- NULL
ig_categories <- unique(plot_df$fusion_ig)
for (ig_cat in ig_categories) {
  topX_by_ig <- c(topX_by_ig, tail(names(sort(
    table(plot_df$fusion[plot_df$fusion_ig == ig_cat]))), n = n_fusions_by_ig))
}
  
# ==============================================================================
# Plot overall most recurrent fusions
# ==============================================================================
plot_df %>% filter(fusion %in% topX_overall) %>% group_by(fusion) %>% 
  summarize(count = n()) %>% 
  ggplot(aes(x = fct_reorder(fusion, count), weight = count)) + geom_bar() + 
  coord_flip() + labs(x = "Fusion count (number of samples)", y = "Fusions") +
  ggplot2_standard_additions()
ggsave("analysis/fusion_summaries/most_recurrent_fusions.pdf", device = "pdf", 
       width = 10, height = 10)

# ==============================================================================
# Plot most recurrent using panels for each type of IG
# ==============================================================================
plot_df %>% filter(fusion %in% topX_by_ig) %>% group_by(fusion) %>% 
  summarize(count = n()) %>% mutate(fusion_ig = return_fusion_ig(fusion)) %>%
  ggplot(aes(x = fct_reorder(fusion, count), weight = count)) + geom_bar() + 
  coord_flip() + labs(y = "Fusion count (number of samples)", x = "Fusions") +
  facet_wrap( ~ fusion_ig, ncol = 1, scales = "free") +
  ggplot2_standard_additions()
ggsave("analysis/fusion_summaries/most_recurrent_fusions.by_ig.pdf", 
       device = "pdf", width = 10, height = 20)

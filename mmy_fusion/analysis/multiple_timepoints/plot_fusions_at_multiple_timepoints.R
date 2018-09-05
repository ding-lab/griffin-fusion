# ==============================================================================
# Explore the landscape of fusions in samples with multiple timepoints
# Steven Foltz (smfoltz@wustl.edu), September 2018
# ==============================================================================

# ==============================================================================
# General Functions
# ==============================================================================

# Plot labels of fusions at multiple timepoints
plot_fusions_at_multiple_timepoints <- function(fusions_tbl, this_mmrf, 
                                                output_directory,
                                                return_tbl = FALSE){
  
  # Create the output directory recursively if it does not already exist
  if (!dir.exists(output_directory)) {
    dir.create(output_directory, recursive = TRUE)
  }  

  plot_tbl <- fusions_tbl %>% filter(mmrf == this_mmrf) %>% 
    select(mmrf, srr, fusion, sample_number, FFPM) %>% 
    mutate(fusion_factor = as.numeric(factor(fusion))) %>% 
    mutate_at("sample_number", factor)
  
  if (return_tbl) {
    return(plot_tbl)
  } else {
    ggplot(plot_tbl, aes(x = sample_number, y = fusion_factor, 
               label = fusion, fill = FFPM)) + geom_label() + 
      scale_fill_continuous(low = "white", high = "red", limits = c(0,NA)) + 
      labs(x = str_c(this_mmrf, " Sample Number"), y = "Unique Fusion Number") +
      ggplot2_standard_additions()
    
    ggsave(str_c(output_directory, "/", this_mmrf, "multiple_timepoints.pdf"), 
           device = "pdf", width = 10, height = 10)
  }
}

# ==============================================================================
# Useful lists of samples
# ==============================================================================

samples_with_multiple_timepoints <- 
  fusions_all %>% filter(has_secondary == 1) %>% pull(mmrf) %>% unique()

# ==============================================================================
# Business
# ==============================================================================

for (this_mmrf in samples_with_multiple_timepoints) {
  plot_fusions_at_multiple_timepoints(fusions_all, this_mmrf, 
                  "analysis/multiple_timepoints/fusions_at_multiple_timepoints")
}
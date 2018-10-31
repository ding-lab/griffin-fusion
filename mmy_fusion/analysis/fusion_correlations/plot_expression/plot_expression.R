# ==============================================================================
# Plot expression of samples with and without fusion or seqFISH events
# Steven Foltz (smfoltz@wustl.edu), October 2018
# ==============================================================================

recreate_plot_df <- FALSE
recreate_all_plots <- FALSE

input_dir <- "analysis/fusion_correlations/event_associations/"
output_dir <- "analysis/fusion_correlations/plot_expression/"
dir.create(output_dir)
dir.create(str_c(output_dir, "all_plots"))
dir.create(str_c(output_dir, "significant_plots"))

testing_tbl <- read_tsv(str_c(input_dir, "testing_tbl.tsv"))
testing_tbl_pvalue_adjusted <- read_tsv(
  str_c(input_dir, "testing_tbl_pvalue_adjusted.tsv"))

ymax_expression_value <- plyr::round_any(
  max(expression_primary$log10tpm), 
  accuracy = .1, f = ceiling)

# ==============================================================================
# General Functions
# ==============================================================================

# Extract sample IDs with a certain fusion
get_ids_with_fusion <- function(this_fusion, fusions_tbl){
  fusions_tbl %>% filter(fusion == this_fusion) %>% 
    select(mmrf, srr) %>% unique()
}
get_ids_without_fusion <- function(this_fusion, fusions_tbl, samples_tbl){
  ids_with_fusion <- get_ids_with_fusion(this_fusion, fusions_tbl)
  samples_tbl %>% anti_join(ids_with_fusion, by = "srr") %>% unique()
}

# Extract sample IDs with a certain gene involved in a fusion
get_ids_with_gene <- function(this_gene, fusions_tbl){
  fusions_tbl %>% filter(geneA == this_gene | geneB == this_gene) %>% 
    select(mmrf, srr) %>% unique()
}
get_ids_without_gene <- function(this_gene, fusions_tbl, samples_tbl){
  ids_with_gene <- get_ids_with_gene(this_gene, fusions_tbl)
  samples_tbl %>% anti_join(ids_with_gene, by = "srr") %>% unique()
}

# Extract sample IDs with a certain seqFISH feature
get_ids_with_seqfish <- function(this_seqfish, seqfish_tbl, samples_tbl){
  seqfish_clinical_info %>% filter( !is.na(seqfish_Study_Visit_ID) ) %>% 
    filter( eval(parse(text = this_seqfish)) == 1 ) %>%
    left_join(samples_tbl, by = "mmrf") %>% select(mmrf, srr) %>% unique()
}
get_ids_without_seqfish <- function(this_seqfish, seqfish_tbl, samples_tbl){
  seqfish_clinical_info %>% filter( !is.na(seqfish_Study_Visit_ID) ) %>% 
    filter( eval(parse(text = this_seqfish)) != 1 ) %>%
    left_join(samples_tbl, by = "mmrf") %>% select(mmrf, srr) %>% unique()
}

# ==============================================================================
# Assign seqFISH and clinical variables to approapiate lists
# ==============================================================================

seqfish_gene_names <- c("seqfish_Translocation_WHSC1_4_14", 
                        "seqfish_Translocation_CCND3_6_14", 
                        "seqfish_Translocation_MYC_8_14", 
                        "seqfish_Translocation_MAFA_8_14", 
                        "seqfish_Translocation_CCND1_11_14", 
                        "seqfish_Translocation_CCND2_12_14", 
                        "seqfish_Translocation_MAF_14_16", 
                        "seqfish_Translocation_MAFB_14_20")

seqfish_variable_names <- c("seqfish_CN_del_13q14", 
                            "seqfish_CN_del_13q34", 
                            "seqfish_CN_del_17p13", 
                            "seqfish_CN_gain_1q21", 
                            "seqfish_Hyperdiploidy", 
                            "seqfish_Translocation_WHSC1_4_14", 
                            "seqfish_Translocation_CCND3_6_14", 
                            "seqfish_Translocation_MYC_8_14", 
                            "seqfish_Translocation_MAFA_8_14", 
                            "seqfish_Translocation_CCND1_11_14", 
                            "seqfish_Translocation_CCND2_12_14", 
                            "seqfish_Translocation_MAF_14_16", 
                            "seqfish_Translocation_MAFB_14_20")

discrete_clinical_variable_names <- c("age_ge_66", 
                                      "Female", 
                                      "Race_White", 
                                      "Race_Black", 
                                      "Race_Other", 
                                      "race", 
                                      "ECOG", 
                                      "ISS_Stage",
                                      "Bone_lesions", 
                                      "Plasmacytoma")

continuous_clinical_variable_names <- c("Age", 
                                        "BM_Plasma_Cell_Percent", 
                                        "LDH")

# ==============================================================================
# Prepare lists of fusions and genes for testing (avoid testing rare events)
# ==============================================================================

# Fusions pairs seen in at least 3 samples
fusion_pairs_gt2 <- fusions_primary %>% group_by(fusion) %>% 
  summarize(count = n()) %>% filter(count >= 3) %>% arrange(desc(count))

# Fusion genes seen in at least 3 samples
fusion_genes_gt2 <- fusions_primary %>% 
  gather(geneA, geneB, key = "geneAB", value = "fusion_gene") %>%
  select(mmrf, srr, fusion_gene) %>% 
  distinct() %>% group_by(fusion_gene) %>% 
  summarize(count = n()) %>% filter(count >= 3) %>% arrange(desc(count))

# ==============================================================================
# Plotting functions
# ==============================================================================

plot_expression_1d <- function(plot_df,
                               gene_list, 
                               labels = TRUE,
                               label_feature = NULL,
                               color_feature = NULL,
                               fill_feature = NULL,
                               shape_feature = NULL,
                               color_label = NULL,
                               fill_label = NULL,
                               shape_label = NULL,
                               pdf_path,
                               ymax_value, 
                               seed = 10){
  
  library(ggplot2)
  library(ggrepel)
  set.seed(seed)
  
  plot_df <- plot_df %>% filter(gene %in% gene_list) %>%
    mutate(fusion_status = !is.na(fusion_label))
  
  p <- ggplot(plot_df, aes_string(x = quote(fill_feature),
                                  y = quote(log10tpm),
                                  label = label_feature,
                                  color = color_feature,
                                  fill = fill_feature,
                                  shape = shape_feature))
  
  p <- p + facet_wrap(~ gene)
  
  p <- p + geom_violin(color = "black",
                       alpha = 0,
                       draw_quantiles = 0.5) #,
                       #position = position_dodge())
  
  # simultaneously jitter and dodge
  # https://ggplot2.tidyverse.org/reference/position_jitterdodge.html
  p <- p + geom_point(aes_string(color = fill_feature, fill = fill_feature),
                      position = position_jitterdodge(jitter.width = 0.5),
                      shape = 16,
                      alpha = 0.75)
  
  # if (labels) {
  #   p <- p + geom_label_repel(data = subset(plot_df, fusion_status == TRUE), 
  #                             aes(x = fusion_jitter,
  #                                 y = log10tpm), 
  #                             label.size = NA, 
  #                             #color = c('white','black')[as.numeric(subset(plot_df, fusion_status=="Fusion")$cnv=="No data")+1], 
  #                             box.padding = 0.35, 
  #                             point.padding = 0.5, 
  #                             segment.color = "grey50")
  # }
  
  p <- p + ylim(0, ymax_value)
  p <- p + scale_color_brewer(palette = "Set1", drop = FALSE)
  p <- p + scale_fill_brewer(palette = "Set1", drop = FALSE)
  p <- p + guides(alpha = FALSE)
  p <- p + labs(x = "Fusion Status", 
                y = "Gene Expression TPM (log10)", 
                color = color_label,
                shape = shape_label,
                fill = fill_label)
                #title = paste0("Fusion Gene Expression (", this_gene, ")"))
  p <- p + ggplot2_standard_additions()
    
  pdf(pdf_path, 10, 10, useDingbats = FALSE)
  print(p)
  shh <- dev.off()
  
}

plot_expression_2d <- function(plot_df,
                               gene_list, 
                               labels = TRUE,
                               label_feature = NULL,
                               color_feature = NULL,
                               fill_feature = NULL,
                               shape_feature = NULL,
                               color_label = NULL,
                               fill_label = NULL,
                               shape_label = NULL,
                               pdf_path,
                               ymax_value, 
                               seed = 10){
  
}

# ==============================================================================
# Create plot data frame
# ==============================================================================

if (recreate_plot_df) {
  return_fusions <- function(fusion_df, this_srr, gene){
    return_value <- fusion_df %>% filter(srr == this_srr, geneA == gene | geneB == gene) %>%
      pull(fusion) %>% str_c(collapse = "\n")
    if (identical(return_value, character(0))) {
      return_value <- NA
    }
    return(return_value)
  }
  
  categorical_cnv <- function(vec) {
    sd_factor <- sd(vec, na.rm = TRUE)
    b = c(2 + sd_factor*c(-Inf, -3, -1, 1, 3, Inf))
    cnv_categories <- vec %>% cut(breaks = b)
    return(cnv_categories)
  }
  
  plot_df <- expression_primary %>% 
    filter(gene %in% fusion_genes_gt2$fusion_gene) %>% 
    mutate(categorical_cnv = categorical_cnv(2*2^gene_avg_cnv)) %>% 
    rowwise() %>% 
    mutate(fusion_label = return_fusions(fusions_primary, srr, gene)) %>%
    left_join(seqfish_clinical_info, by = "mmrf")
  
  write_tsv(plot_df, str_c(output_dir, "expression_plot_tibble.tsv"))
} else {
  plot_df <- read_tsv(str_c(output_dir, "expression_plot_tibble.tsv"))
}

# ==============================================================================
# Business
# ==============================================================================

if (recreate_all_plots) {
  
  # Plot expression of genes involved in seqFISH translocations
  for (full_seqfish_name in seqfish_gene_names) {
    print(full_seqfish_name)
    gene <- str_split(full_seqfish_name, "_", simplify = TRUE)[3]
    print(gene)
    samples_with <- get_ids_with_seqfish(full_seqfish_name, 
                                         seqfish_tbl = seqfish_clinical_info, 
                                         samples_tbl = samples_primary)
    samples_without <- get_ids_without_seqfish(full_seqfish_name, 
                                            seqfish_tbl = seqfish_clinical_info, 
                                            samples_tbl = samples_primary)
    plot_expression_1d(gene_name = gene,
                       expr_file = expression_primary,
                       with_mmrf = samples_with,
                       without_mmrf = samples_without,
                       output_dir = str_c(output_dir, "all_plots/",
                                          gene, ".seqFISH.pdf"))
  }
  
  # Plot expression of genes recurrently involved in fusions
  for (gene in fusion_genes_gt2$fusion_gene) {
    if ( str_detect(gene, "@") ) {
      next
    }
    print(gene)
    samples_with <- get_ids_with_gene(gene, fusions_primary)
    samples_without <- get_ids_without_gene(gene, fusions_primary, 
                                            samples_primary)
    

  }
  
}
# ==============================================================================
# Plot expression of samples with and without fusion or seqFISH events
# Steven Foltz (smfoltz@wustl.edu), October 2018
# ==============================================================================

################################################################################
# Set seed for reproducibility
################################################################################
set.seed(10)

################################################################################
# Set seed for reproducibility
################################################################################

recreate_plot_df <- FALSE
recreate_all_plots <- FALSE

input_dir <- "analysis/fusion_correlations/event_associations/"
output_dir <- "analysis/fusion_correlations/plot_expression/"
dir.create(output_dir, showWarnings = FALSE)
dir.create(str_c(output_dir, "all_plots"), showWarnings = FALSE)
dir.create(str_c(output_dir, "significant_plots"), showWarnings = FALSE)
dir.create(str_c(output_dir, "seqFISH"), showWarnings = FALSE)
dir.create(str_c(output_dir, "seqFISH/fusions"), showWarnings = FALSE)
dir.create(str_c(output_dir, "seqFISH/translocations"), showWarnings = FALSE)

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

seqfish_gene_names_formatted <- c("Translocation t(4;14) WHSC1", 
                                  "Translocation t(6;14) CCND3", 
                                  "Translocation t(8;14) MYC", 
                                  "Translocation t(8;14) MAFA", 
                                  "Translocation t(11;14) CCND1", 
                                  "Translocation t(12;14) CCND2", 
                                  "Translocation t(14;16) MAF", 
                                  "Translocation t(14;20) MAFB")

seqfish_gene_names_formatted_short <- c("t(4;14) WHSC1", 
                                        "t(6;14) CCND3", 
                                        "t(8;14) MYC", 
                                        "t(8;14) MAFA", 
                                        "t(11;14) CCND1", 
                                        "t(12;14) CCND2", 
                                        "t(14;16) MAF", 
                                        "t(14;20) MAFB")

seqfish_genes <- c("WHSC1", "CCND3", "MYC", "MAFA", 
                   "CCND1", "CCND2", "MAF", "MAFB")


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

plot_fusion_expression_1d <- function(plot_df,
                                      gene_list,
                                      labels = FALSE,
                                      translocation = NULL,
                                      translocation_formatted = NULL,
                                      ymax_value,
                                      pdf_path,
                                      pdf_width = 10,
                                      pdf_height = 10,
                                      seed = 10){
  
  library(ggplot2)
  library(ggrepel)
  set.seed(seed)
  
  plot_df <- plot_df %>% filter(gene %in% gene_list)
  if (is.null(translocation)) {
    shape_vector <- "No shape vector given"
  } else {
    shape_vector <- plot_df %>% pull(translocation) %>%
      factor(labels = c("Not detected", "Detected", "Missing"), exclude = NULL)  
  }
  plot_df <- plot_df %>% mutate(shape_factor = shape_vector)
  
  p <- ggplot(plot_df)
  
  p <- p + facet_wrap(~ gene, nrow = 1)
  
  p <- p + geom_violin(aes(x = fusion_indicator,
                           y = log10tpm,
                           fill = fusion_status),
                       color = "black",
                       draw_quantiles = 0.5)
  
  p <- p + geom_point(aes(x = fusion_jitter,
                          y = log10tpm,
                          color = cnv_factor,
                          shape = shape_factor))
  
  if (labels & length(gene_list) == 1) {
    p <- p + geom_label_repel(data = plot_df %>% filter(!is.na(fusion_label)),
                              aes(x = fusion_jitter,
                                  y = log10tpm,
                                  label = fusion_label,
                                  color = cnv_factor),
                              point.padding = 0.5,
                              show.legend = FALSE)
  }
  
  p <- p + scale_x_continuous(breaks = c(0, 1), 
                              labels = c("No fusion", "Fusion"))
  
  p <- p + ylim(0, ymax_value)
  
  color_scale <- c("#1f78b4", "#a6cee3", # deletions
                   "#b2df8a", # neutral 
                   "#fb9a99", "#e31a1c", # amplifications
                   "#cab2d6") #missing
  p <- p + scale_color_manual(values = color_scale, drop = FALSE)
  
  p <- p + scale_fill_manual(values = c("#ffffff", "#ffffff")) # both white
  
  if (is.null(translocation)) {
    p <- p + scale_shape_manual(values = 16)
    p <- p + guides(alpha = FALSE, fill = FALSE, shape = FALSE)
    p <- p + labs(x = NULL,
                  y = "Gene Expression TPM (log10)", 
                  color = "Copy Number")
  } else {
    p <- p + scale_shape_manual(values = c(16, 17, 4))
    p <- p + guides(alpha = FALSE, fill = FALSE)
    p <- p + labs(x = NULL,
                  y = "Gene Expression TPM (log10)", 
                  color = "Copy Number",
                  shape = translocation_formatted)
  }
  
  p <- p + ggplot2_standard_additions()
    
  pdf(pdf_path, width = pdf_width, height = pdf_height, useDingbats = FALSE)
  print(p)
  shh <- dev.off()
  
}

plot_translocation_expression_1d <- function(plot_df,
                                             gene_list,
                                             labels = FALSE,
                                             ymax_value,
                                             pdf_path,
                                             pdf_width = 10, 
                                             pdf_height = 10,
                                             seed = 10){
  
  library(ggplot2)
  library(ggrepel)
  set.seed(seed)
  
  plot_df <- plot_df %>% filter(gene %in% gene_list, 
                                !is.na(translocation_indicator))
  
  p <- ggplot(plot_df)
  
  p <- p + facet_wrap(~ seqfish_gene_names_formatted, nrow = 1)
  
  p <- p + geom_violin(aes(x = translocation_indicator,
                           y = log10tpm,
                           fill = factor(translocation_indicator)),
                       color = "black",
                       draw_quantiles = 0.5)
  
  p <- p + geom_point(aes(x = translocation_jitter,
                          y = log10tpm,
                          color = cnv_factor),
                      shape = 16)
  
  if (labels & length(gene_list) == 1) {
    p <- p + geom_label_repel(data = plot_df %>% filter(!is.na(fusion_label)),
                              aes(x = translocation_jitter,
                                  y = log10tpm,
                                  label = fusion_label,
                                  color = cnv_factor),
                              point.padding = 0.5,
                              show.legend = FALSE)
  }
  
  p <- p + scale_x_continuous(breaks = c(0, 1), 
                              labels = c("No translocation", "Translocation"))
  
  p <- p + ylim(0, ymax_value)
  
  color_scale <- c("#1f78b4", "#a6cee3", # deletions
                   "#b2df8a", # neutral 
                   "#fb9a99", "#e31a1c", # amplifications
                   "#cab2d6") #missing
  p <- p + scale_color_manual(values = color_scale, drop = FALSE)
  
  p <- p + scale_fill_manual(values = c("#ffffff", "#ffffff")) # both white
  
  p <- p + guides(alpha = FALSE, fill = FALSE)
  
  p <- p + labs(x = NULL,
                y = "Gene Expression TPM (log10)", 
                color = "Copy Number")
  
  p <- p + ggplot2_standard_additions()
  
  pdf(pdf_path, width = pdf_width, height = pdf_height, useDingbats = FALSE)
  print(p)
  shh <- dev.off()
  
}

plot_expression_2d <- function(plot_df,
                               gene_list,
                               fusion1 = NULL,
                               fusion2 = NULL,
                               translocation = NULL,
                               translocation_formatted = NULL,
                               ymax_value,
                               pdf_path,
                               pdf_width = 10,
                               pdf_height = 10,
                               seed = 10){
  
  library(ggplot2)
  library(ggrepel)
  set.seed(seed)
  
  if (length(gene_list) != 2) {
    stop("Gene list can only have two genes for 2D plot.")
  }
  
  plot_df <- plot_df %>% filter(gene %in% gene_list)
  shape_vector <- plot_df %>% pull(translocation) %>% 
    factor(labels = c("Not detected", "Detected", "Missing"), exclude = NULL)
  plot_df <- plot_df %>% mutate(shape_factor = shape_vector)
  
  fusion1_reverse <- str_c(rev(
    str_split(fusion1, pattern = "--", simplify = TRUE)), collapse = "--")
  gene1 <- plot_df %>% filter(gene == gene_list[1]) %>%
    select(mmrf, srr, gene, log10tpm, fusion_label, shape_factor) %>%
    mutate(has_fusion1 = str_detect(fusion_label, pattern = fusion1) | 
                     str_detect(fusion_label, pattern = fusion1_reverse))
  
  fusion2_reverse <- str_c(rev(
    str_split(fusion2, pattern = "--", simplify = TRUE)), collapse = "--")
  gene2 <- plot_df %>% filter(gene == gene_list[2]) %>%
    select(mmrf, srr, gene, log10tpm, fusion_label, shape_factor) %>%
    mutate(has_fusion2 = str_detect(fusion_label, pattern = fusion2) | 
             str_detect(fusion_label, pattern = fusion2_reverse))
  
  plot_df <- left_join(gene1, gene2, by = c("mmrf", "srr", "shape_factor")) %>%
    replace_na(list(has_fusion1 = FALSE, has_fusion2 = FALSE)) %>%
    mutate(fusion_category = has_fusion1 + 2*has_fusion2) %>%
    mutate(fusion_category = factor(fusion_category,
                                    levels = c(0,1,2,3),
                                    labels = c("None reported",
                                               fusion1,
                                               fusion2,
                                               "Both reported")
                                    
                                    )
           )
  
  p <- ggplot(plot_df)
  
  p <- p + geom_point(aes(x = log10tpm.x,
                          y = log10tpm.y,
                          color = fusion_category,
                          shape = shape_factor))
  
  p <- p + xlim(0, ymax_value)
  
  p <- p + ylim(0, ymax_value)
  
  p <- p + scale_shape_manual(values = c(16, 17, 4))
  
  p <- p + guides(alpha = FALSE, fill = FALSE)
  
  p <- p + labs(x = str_c(gene_list[1], " Expression TPM (log10)"),
                y = str_c(gene_list[2], " Expression TPM (log10)"), 
                color = "Fusion status",
                shape = translocation_formatted)
  
  p <- p + ggplot2_standard_additions()
  
  p <- p + coord_fixed(ratio = 1)
  
  pdf(pdf_path, width = pdf_width, height = pdf_height, useDingbats = FALSE)
  print(p)
  shh <- dev.off()
  
}

# ==============================================================================
# Create plot data frame
# ==============================================================================

if (recreate_plot_df) {
  return_fusions <- function(fusions_df, this_srr, gene){
    return_value <- fusions_df %>% filter(srr == this_srr,
                                          geneA == gene | geneB == gene) %>%
      pull(fusion) %>% str_c(collapse = "\n")
    if (identical(return_value, character(0))) {
      return_value <- NA
    }
    return(return_value)
  }
  
  return_translocations <- function(seqfish_df, this_mmrf, gene, 
                                    seqfish_gene_names, seqfish_genes){
    if (gene %in% seqfish_genes) {
      return_value <- seqfish_df %>% filter(mmrf == this_mmrf) %>%
        pull( seqfish_gene_names[which(seqfish_genes == gene)] )
    } else {
      return_value <- NA
    }
    return(return_value)
  }
  
  return_translocations_formatted <- function(gene, seqfish_genes, 
                                              seqfish_gene_names_formatted){
    if (gene %in% seqfish_genes) {
      return_value <- seqfish_gene_names_formatted[which(seqfish_genes == gene)]
    } else {
      return(NA)
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
    filter(gene %in% fusion_genes_gt2$fusion_gene | 
             gene %in% seqfish_genes) %>% 
    mutate(categorical_cnv = categorical_cnv(2*2^gene_avg_cnv)) %>% 
    mutate(cnv_factor = factor(categorical_cnv, 
                               labels = c("DELETION",
                                          "Deletion",
                                          "Neutral",
                                          "Amplification",
                                          "AMPLIFICATION",
                                          "Missing"), 
                               exclude = NULL)) %>%
    rowwise() %>% 
    mutate(translocation_indicator = return_translocations(
      seqfish_clinical_info, mmrf, gene, seqfish_gene_names, seqfish_genes)) %>%
    mutate(seqfish_gene_names_formatted = return_translocations_formatted(
      gene, seqfish_genes, seqfish_gene_names_formatted)) %>%
    mutate(fusion_label = return_fusions(fusions_primary, srr, gene)) %>%
    ungroup() %>%
    mutate(fusion_status = !is.na(fusion_label)) %>%
    mutate(fusion_indicator = as.numeric(fusion_status)) %>%
    mutate(fusion_jitter = jitter(fusion_indicator)) %>%
    mutate(translocation_jitter = jitter(translocation_indicator)) %>%
    left_join(seqfish_clinical_info, by = "mmrf")
  
  write_tsv(plot_df, str_c(output_dir, "expression_plot_tibble.tsv"))
} else {
  plot_df <- read_tsv(str_c(output_dir, "expression_plot_tibble.tsv"))
  plot_df <- plot_df %>% mutate(cnv_factor = factor(categorical_cnv,
                                                    labels = c("DELETION",
                                                               "Deletion",
                                                               "Neutral",
                                                               "Amplification",
                                                               "AMPLIFICATION",
                                                               "Missing"), 
                                                    exclude = NULL))
}

# ==============================================================================
# Business
# ==============================================================================

if (recreate_all_plots) {
  
  # Plot expression of genes involved in seqFISH translocations
  for (this_gene in seqfish_genes) {
    print(this_gene)
    n_samples_with_fusion <- plot_df %>% 
      filter(gene == this_gene, !is.na(fusion_label)) %>% nrow()
    if (n_samples_with_fusion > 2) {
      t_index = which(seqfish_genes == this_gene)
      t_name = seqfish_gene_names[t_index]
      t_format = seqfish_gene_names_formatted_short[t_index]
      plot_fusion_expression_1d(plot_df,
                                this_gene,
                                labels = FALSE,
                                translocation = t_name,
                                translocation_formatted = t_format,
                                ymax_value = ymax_expression_value,
                                pdf_path = str_c(output_dir, 
                                                 "seqFISH/fusions/", 
                                                 this_gene, ".pdf"),
                                pdf_width = 10,
                                pdf_height = 10,
                                seed = 10)
      plot_translocation_expression_1d(plot_df,
                                       this_gene,
                                       labels = FALSE,
                                       ymax_value = ymax_expression_value,
                                       pdf_path = str_c(output_dir, "seqFISH/",
                                                        "translocations/",
                                                        this_gene, ".pdf"),
                                       pdf_width = 10, 
                                       pdf_height = 10,
                                       seed = 10)
    }
  }
  
  plot_fusion_expression_1d(plot_df,
                            seqfish_genes,
                            labels = FALSE,
                            ymax_value = ymax_expression_value,
                            pdf_path = str_c(output_dir, "seqFISH/fusions/",
                                             "all.pdf"),
                            pdf_width = 4*length(seqfish_gene_names),
                            pdf_height = 10,
                            seed = 10)  
  
  plot_translocation_expression_1d(plot_df,
                                   seqfish_genes,
                                   labels = FALSE,
                                   ymax_value = ymax_expression_value,
                                   pdf_path = str_c(output_dir, 
                                                    "seqFISH/translocations/",
                                                    "all.pdf"),
                                   pdf_width = 4*length(seqfish_gene_names),
                                   pdf_height = 10,
                                   seed = 10)  
  
  # Plot expression of genes recurrently involved in fusions
  for (this_gene in fusion_genes_gt2$fusion_gene) {
    if ( str_detect(this_gene, "@") ) {
      next
    }
    print(this_gene)
    plot_fusion_expression_1d(plot_df,
                              this_gene,
                              labels = FALSE,
                              ymax_value = ymax_expression_value,
                              pdf_path = str_c(output_dir, "all_plots/", 
                                              this_gene, ".pdf"),
                              pdf_width = 10,
                              pdf_height = 10,
                              seed = 10)  
    plot_fusion_expression_1d(plot_df,
                              this_gene,
                              labels = TRUE,
                              ymax_value = ymax_expression_value,
                              pdf_path = str_c(output_dir, "all_plots/", 
                                               this_gene, ".labeled.pdf"),
                              pdf_width = 10,
                              pdf_height = 10,
                              seed = 10)
  }
  
  significant_fusion_expresion_genes <- testing_tbl_pvalue_adjusted %>% 
    filter(fdr < 0.05 | median_value > 0.9, 
           event_type %in% c("Fusion Expression", 
                             "Fusion Expression Outlier")) %>% 
    pull(event1) %>% unique()
  
  for (this_gene in significant_fusion_expresion_genes) {
    print(this_gene)
    n_samples_with_fusion <- plot_df %>% 
      filter(gene == this_gene, !is.na(fusion_label)) %>% nrow()
    if (n_samples_with_fusion > 2) {
      plot_fusion_expression_1d(plot_df,
                                this_gene, 
                                labels = FALSE,
                                ymax_value = ymax_expression_value,
                                pdf_path = str_c(output_dir, 
                                                 "significant_plots/",
                                                this_gene, ".pdf"),
                               pdf_width = 10,
                               pdf_height = 10,
                               seed = 10)  
      plot_fusion_expression_1d(plot_df,
                                this_gene, 
                                labels = TRUE,
                                ymax_value = ymax_expression_value,
                                pdf_path = str_c(output_dir, 
                                                 "significant_plots/",
                                                 this_gene, ".labeled.pdf"),
                                pdf_width = 10,
                                pdf_height = 10,
                                seed = 10)  
    }
  }
  
  plot_fusion_expression_1d(plot_df,
                     significant_fusion_expresion_genes,
                     labels = FALSE,
                     ymax_value = ymax_expression_value,
                     pdf_path = str_c(output_dir, "significant_plots/", 
                                      "all.pdf"),
                     pdf_width = 4*length(significant_fusion_expresion_genes),
                     pdf_height = 10,
                     seed = 10)
}
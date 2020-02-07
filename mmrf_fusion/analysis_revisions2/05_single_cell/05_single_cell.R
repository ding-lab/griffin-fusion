# ==============================================================================
# scRNA analysis (in-house samples)
# Steven Foltz (github: envest)
# ==============================================================================

paper_main = "paper/main/05_single_cell/"
paper_supp = "paper/supplementary/05_single_cell/"

# Create directories
for (id in c("27522_1/", "27522_4/")) {
  dir.create(str_c(paper_main, id), recursive = TRUE, showWarnings = FALSE)
  dir.create(str_c(paper_supp, id), recursive = TRUE, showWarnings = FALSE)
}

# ==============================================================================
# Get t-SNE and UMAP mappings in same tibble
# ==============================================================================
get_tsne_umap <- function(cell_types, seurat_object){

  tsne_coords <- Embeddings(object = seurat_object,
                            reduction = "tsne") %>%
    as_tibble() %>%
    mutate(barcode = row.names(Embeddings(object = seurat_object,
                                          reduction = "tsne"))) %>%
    select(barcode, tSNE_1, tSNE_2)

  umap_coords <- Embeddings(object = seurat_object,
                            reduction = "umap") %>%
    as_tibble() %>%
    mutate(barcode = row.names(Embeddings(object = seurat_object,
                                          reduction = "umap"))) %>%
    select(barcode, UMAP_1, UMAP_2)

  cell_types %>%
    mutate(cell_type = factor(cell_type,
                              levels = c("B",
                                         "CD4+T",
                                         "CD8+T",
                                         "DC",
                                         "CD16+Mono",
                                         "CD14+Mono",
                                         "NK",
                                         "Plasma"),
                              labels = c("B Cells",
                                         "CD4+ T Cells",
                                         "CD8+ T Cells",
                                         "Dendritic Cells",
                                         "Macrophages",
                                         "Monocytes",
                                         "Natural Killer Cells",
                                         "Plasma Cells"))) %>%
    left_join(tsne_coords, by = "barcode") %>%
    left_join(umap_coords, by = "barcode") %>%
    arrange(cell_type) %>%
    return()
}

# ==============================================================================
# Get expression of one gene and match with barcode tsne and umap
# ==============================================================================
get_gene_expression <- function(seurat_object, tsne_umap, ensg){
  gene_expression <- tibble(barcode = seurat_object@assays$RNA@data@Dimnames[[2]],
                            expression = seurat_object@assays$RNA@data[ensg,] %>% as.vector())
  tsne_umap %>%
    left_join(gene_expression, by = "barcode") %>%
    return()
}

# ==============================================================================
# Get CNV of one gene and match with barcode tsne and umap
# ==============================================================================
get_gene_cnv <- function(infercnv, tsne_umap, this_gene){
  cnv_level <- tibble(barcode = infercnv %>% select(-gene) %>% colnames(),
                      cnv = infercnv %>% filter(gene == this_gene) %>%
                        select(-gene) %>% t() %>% as.vector())
  tsne_umap %>%
    left_join(cnv_level, by = "barcode") %>%
    return()
}

# ==============================================================================
# Get mean CNV of one chromosome and match with barcode tsne and umap
# ==============================================================================
get_chr_cnv <- function(chr, infercnv, tsne_umap, genes = gene_spans){

  genes %>%
    filter(chromosome == chr) %>%
    left_join(infercnv, by = c("gene_name" = "gene")) %>%
    arrange(start) %>%
    select(ends_with("-1")) %>%
    t() %>% as.data.frame() %>%
    rownames_to_column(var = "barcode") %>%
    gather(starts_with("V"), key = "gene", value = "cnv") %>%
    filter(!is.na(cnv)) %>%
    group_by(barcode) %>%
    summarize(mean_cnv = mean(cnv)) %>%
    arrange(mean_cnv) %>%
    right_join(tsne_umap, by = "barcode") %>%
    return()
  }

# ==============================================================================
# Get expression of two genes and match with barcode tsne and umap
# ==============================================================================
get_two_genes_expression <- function(seurat_object, tsne_umap, ensg1, ensg2){

  two_genes_expression <- tibble(barcode = seurat_object@assays$RNA@data@Dimnames[[2]],
                                 expression1 = seurat_object@assays$RNA@data[ensg1,] %>% as.vector(),
                                 expression2 = seurat_object@assays$RNA@data[ensg2,] %>% as.vector())
  tsne_umap %>%
    left_join(two_genes_expression, by = "barcode") %>%
    return()
}

# ==============================================================================
# Plot cell types
# ==============================================================================
plot_cell_types <- function(cell_types, seurat_object, reduction = "UMAP", id, dir = paper_main){

  if (reduction %in% c("UMAP", "t-SNE")) {
    plot_df = get_tsne_umap(cell_types, seurat_object)
    if (reduction == "UMAP") {
      p <- ggplot(data = plot_df, aes(x = UMAP_1, y = UMAP_2))
      p <- p + labs(x = "UMAP 1", y = "UMAP 2")
    } else {
      p <- ggplot(data = plot_df, aes(x = tSNE_1, y = tSNE_2))
      p <- p + labs(x = "t-SNE 1", y = "t-SNE 2")
    }
  } else {
  stop("Reduction parameter reduction must be 'UMAP' or 't-SNE'.")
  }

  cell_type_positions <- plot_df %>% group_by(cell_type) %>%
    summarize(mean_tSNE_1 = mean(tSNE_1),
              mean_tSNE_2 = mean(tSNE_2),
              mean_UMAP_1 = mean(UMAP_1),
              mean_UMAP_2 = mean(UMAP_2))

  p <- p + geom_point(aes(color = cell_type),
                      shape = 16,
                      size = 1.5,
                      show.legend = FALSE,
                      alpha = 0.25)

  if (reduction == "UMAP") {
    p <- p + geom_text(data = cell_type_positions,
                       aes(x = mean_UMAP_1, y = mean_UMAP_2,
                           label = cell_type,
                           color = cell_type),
                       size = 3.5,
                       show.legend = FALSE)
  } else {
    p <- p + geom_text(data = cell_type_positions,
                       aes(x = mean_tSNE_1, y = mean_tSNE_2,
                           label = cell_type,
                           color = cell_type),
                       size = 3.5,
                       show.legend = FALSE)
  }

  p <- p + annotate("text", label = "Cell Types", x = -Inf, y = Inf,
                vjust = 1, hjust = 0,
                size = 3.5) +
    theme_bw() +
    #coord_fixed() +
    scale_color_brewer(palette = "Set3", drop = FALSE, direction = -1) +
    scale_x_continuous(expand = c(0.01, 0.01)) +
    scale_y_continuous(expand = c(0.01, 0.01)) +
    theme(panel.background = element_blank(),
          panel.border = element_blank(),
          plot.background = element_blank(),
          panel.grid = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10))

  ggsave(str_c(dir, "/", id, ".cell_types.", reduction, ".pdf"),
         p,
         width = 3.5, height = 3.5, useDingbats = FALSE)

}

# ==============================================================================
# Plot single gene expression
# ==============================================================================
plot_gene_expression <- function(seurat_object, tsne_umap, ensg, reduction = "UMAP", id, gene, dir = paper_supp) {

  if (reduction %in% c("UMAP", "t-SNE")) {
    plot_df <- get_gene_expression(seurat_object, tsne_umap, ensg) %>% arrange(expression)
    if (reduction == "UMAP") {
      p <- ggplot(data = plot_df, aes(x = UMAP_1, y = UMAP_2))
      p <- p + labs(x = "UMAP 1", y = "UMAP 2")
    } else {
      p <- ggplot(data = plot_df, aes(x = tSNE_1, y = tSNE_2))
      p <- p + labs(x = "t-SNE 1", y = "t-SNE 2")
    }
  } else {
    stop("Reduction parameter reduction must be 'UMAP' or 't-SNE'.")
  }

  p <- p + geom_point(aes(color = expression),
                      shape = 16,
                      size = 1.5,
                      show.legend = FALSE,
                      alpha = 0.25) +
  annotate("text", label = gene,
           x = -Inf, y = Inf,
           vjust = 1, hjust = 0,
           size = 3.5,
           fontface = "italic") +
  theme_bw() +
  #coord_fixed() +
  scale_color_continuous(low = "grey", high = "red") +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        plot.background = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))

  ggsave(str_c(dir, "/", id, ".expression.", gene, ".", reduction, ".pdf"),
         p,
         width = 3.5, height = 3.5, useDingbats = FALSE)
}

# ==============================================================================
# Plot single gene expression -- violin
# ==============================================================================
plot_gene_expression_violin <- function(bulk_sc, seurat_object, tsne_umap, ensg, id, gene, dir = paper_supp, max = NULL) {
  set.seed(10)

  gene_df <- get_gene_expression(seurat_object, tsne_umap, ensg) %>% arrange(expression)

  if (identical(max, NULL)) {
    max <- gene_df %>% pull(expression) %>% max(na.rm = TRUE)
  }

  plot_df <- bulk_sc %>%
    filter(data_type == "Single Cell Chimeric Transcript") %>%
    separate(identifier,
             into = c("cell_barcode", "molecular_barcode"),
             sep = ":") %>%
    select(cell_barcode, data_type) %>%
    unique() %>%
    right_join(gene_df, by = c("cell_barcode" = "barcode")) %>%
    filter(cell_type == "Plasma Cells")

  p <- ggplot(data = plot_df, aes(x = !is.na(data_type), y = expression)) +
    geom_violin(scale = "width",
                color = "black",
                draw_quantiles = 0.5) +
    geom_jitter(height = 0,
                shape = 16,
                size = 1.5,
                show.legend = FALSE,
                alpha = 0.25) +
    annotate("text", label = gene,
             x = -Inf, y = Inf,
             vjust = 1, hjust = 0,
             size = 3.5,
             fontface = "italic") +
    theme_bw() +
    coord_equal() +
    labs(x = NULL, y = "Scaled Expression") +
    scale_color_continuous(low = "grey", high = "red") +
    scale_x_discrete(labels = c("No CT\nDetected", "CT Detected")) +
    scale_y_continuous(expand = c(0.01, 0.01), limits = c(0, max)) +
    theme(panel.background = element_blank(),
          panel.border = element_blank(),
          plot.background = element_blank(),
          panel.grid = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10))

  ggsave(str_c(dir, "/", id, ".expression_violin.", gene, ".pdf"),
         p,
         width = 3.5, height = 3.5, useDingbats = FALSE)
}

# ==============================================================================
# Plot two genes expression
# ==============================================================================
plot_two_genes_expression <- function(seurat_object, tsne_umap, ensg1, ensg2, id, gene1, gene2, dir = paper_main) {

  plot_df <- get_two_genes_expression(seurat_object, tsne_umap, ensg1, ensg2)

  max_x <- plot_df %>% pull(expression_1) %>% max()
  may_y <- plot_df %>% pull(expression_2) %>% max()
  max_xy <- max(max_x, max_y)

  p <- ggplot(data = plot_df, aes(x = expression_1, y = expression_2))

  p <- p + labs(x = str_c(gene1, " Expression"), y = str_c(gene2, " Expression"))

  p <- p + geom_point(aes(color = expression),
                      shape = 16,
                      size = 1.5,
                      show.legend = FALSE,
                      alpha = 0.25) +
    theme_bw() +
    #coord_fixed() +
    scale_x_continuous(limits = c(0, max_xy), expand = c(0.01, 0.01)) +
    scale_y_continuous(limits = c(0, max_xy), expand = c(0.01, 0.01)) +
    theme(panel.background = element_blank(),
          panel.border = element_blank(),
          plot.background = element_blank(),
          panel.grid = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10))

  ggsave(str_c(dir, "/", id, ".joint_expression.", gene1, ".", gene2, ".pdf"),
         p,
         width = 3.5, height = 3.5, useDingbats = FALSE)
}

# ==============================================================================
# Plot single gene CNV
# ==============================================================================
plot_gene_cnv <- function(infercnv, tsne_umap, reduction = "UMAP", id, gene, dir = paper_supp) {

  if (reduction %in% c("UMAP", "t-SNE")) {
    plot_df <- get_gene_cnv(infercnv, tsne_umap, this_gene = gene) %>%
      rowwise() %>%
      mutate(diff = max(1, cnv) - min(1, cnv)) %>%
      ungroup() %>%
      arrange(!is.na(cnv), diff)
    if (reduction == "UMAP") {
      p <- ggplot(data = plot_df, aes(x = UMAP_1, y = UMAP_2))
      p <- p + labs(x = "UMAP 1", y = "UMAP 2")
    } else {
      p <- ggplot(data = plot_df, aes(x = tSNE_1, y = tSNE_2))
      p <- p + labs(x = "t-SNE 1", y = "t-SNE 2")
    }
  } else {
    stop("Reduction parameter reduction must be 'UMAP' or 't-SNE'.")
  }

  p <- p + geom_point(aes(color = cnv),
                      shape = 16,
                      size = 1.5,
                      #show.legend = FALSE,
                      alpha = 0.25) +
    annotate("text", label = gene,
             x = -Inf, y = Inf,
             vjust = 1, hjust = 0,
             size = 3.5,
             fontface = "italic") +
    theme_bw() +
    #coord_fixed() +
    scale_color_gradient2(high = "#d73027", mid = "#ffffbf",
                         low = "#4575b4", na.value = "grey",
                         midpoint = 1,
                         limits = c(0,2)) +
  labs(color = "CNV Ratio") +
    scale_x_continuous(expand = c(0.01, 0.01)) +
    scale_y_continuous(expand = c(0.01, 0.01)) +
    theme(panel.background = element_blank(),
          panel.border = element_blank(),
          plot.background = element_blank(),
          panel.grid = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10))

  ggsave(str_c(dir, "/", id, ".copy_number.", gene, ".pdf"),
         p,
         width = 3.5, height = 3.5, useDingbats = FALSE)
}

# ==============================================================================
# Plot chromosome CNV
# ==============================================================================
plot_chr_cnv <- function(infercnv, tsne_umap, reduction = "UMAP", id, chr, dir = paper_supp) {

  if (reduction %in% c("UMAP", "t-SNE")) {
    plot_df <- get_chr_cnv(chr, infercnv, tsne_umap) %>%
      rowwise() %>%
      mutate(diff = max(1, mean_cnv) - min(1, mean_cnv)) %>%
      ungroup() %>%
      arrange(!is.na(mean_cnv), diff)
    if (reduction == "UMAP") {
      p <- ggplot(data = plot_df, aes(x = UMAP_1, y = UMAP_2))
      p <- p + labs(x = "UMAP 1", y = "UMAP 2")
    } else {
      p <- ggplot(data = plot_df, aes(x = tSNE_1, y = tSNE_2))
      p <- p + labs(x = "t-SNE 1", y = "t-SNE 2")
    }
  } else {
    stop("Reduction parameter reduction must be 'UMAP' or 't-SNE'.")
  }

  p <- p + geom_point(aes(color = mean_cnv),
                      shape = 16,
                      size = 1.5,
                      alpha = 0.25) +
    annotate("text", label = chr,
             x = -Inf, y = Inf,
             vjust = 1, hjust = 0,
             size = 3.5) +
    theme_bw() +
    #coord_fixed() +
    scale_color_gradient2(high = "#d73027", mid = "#ffffbf",
                          low = "#4575b4", na.value = "grey",
                          midpoint = 1,
                          limits = c(0,2)) +
    labs(color = "Mean\nCNV Ratio") +
    scale_x_continuous(expand = c(0.01, 0.01)) +
    scale_y_continuous(expand = c(0.01, 0.01)) +
    theme(panel.background = element_blank(),
          panel.border = element_blank(),
          plot.background = element_blank(),
          panel.grid = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10))

  ggsave(str_c(dir, "/", id, ".copy_number.", chr, ".pdf"),
         p,
         width = 3.5, height = 3.5, useDingbats = FALSE)
}

# ==============================================================================
# Plot read length distribution for reads supporting scRNA chimeric transcripts
# ==============================================================================
plot_sc_read_length_distribution <- function(dis_reads, cutoff = 1024, id, dir = paper_supp){

  diff <- dis_reads %>%
    mutate(diff = log2(end - start)) %>%
    pull(diff)
  max_length <- max(diff) %>% ceiling()
  line_x = log2(cutoff)
  line_y = mean(line_x > diff)


  plot_df <- tibble(log2length = seq(1, max_length)) %>%
    rowwise() %>%
    mutate(proportion = mean(log2length > diff))

  p <- ggplot(data = plot_df, aes(x = log2length, y = proportion)) +
    geom_line(size = 2, color = "grey50") +
    geom_segment(x = line_x, xend = line_x, y = 0, yend = line_y,
                 linetype = 2) +
    geom_segment(x = 0, xend = line_x, y = line_y, yend = line_y,
                 linetype = 2) +
    annotate(geom = "text", x = line_x + 0.5, y = 0.01 ,
             label = "Discard",
             vjust = 0, hjust = 0) +
    annotate(geom = "text", x = line_x - 0.5, y = 0.01 ,
             label = "Keep",
             vjust = 0, hjust = 1) +
    annotate(geom = "text", x = 1.1, y = line_y - 0.01,
             label = round(line_y, 4),
             vjust = 1, hjust = 0) +
    scale_x_continuous(breaks = seq(1, max_length),
                       labels = seq(1, max_length),
                       expand = c(0,0)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1),
                       labels = seq(0, 1, by = 0.1),
                       expand = c(0,0)) +
    labs(x = "Read length (base pairs, log2)",
         y = "Cumulative Proportion of Reads") +
    theme_bw() +
    theme(panel.background = element_blank(),
          panel.border = element_blank(),
          plot.background = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10))

  ggsave(str_c(dir, "/", id, ".read_length_distribution.pdf"),
         width = 3.5, height = 3.5, useDingbats = FALSE)
}

# ==============================================================================
# Identify cells with chimeric transcripts
# ==============================================================================
plot_cell_chimeric_transcripts <- function(bulk_sc, tsne_umap, reduction = "UMAP", id, dir = paper_main, facet = FALSE, color_yes_no = TRUE){

  if (reduction %in% c("UMAP", "t-SNE")) {

    plot_df <- bulk_sc %>%
      filter(data_type == "Single Cell Chimeric Transcript") %>%
      separate(identifier,
               into = c("cell_barcode", "molecular_barcode"),
               sep = ":") %>%
      group_by(cell_barcode) %>%
      summarize(color_column = str_c(sort(unique(color_column)), collapse = ", ")) %>%
      right_join(tsne_umap,
                 by = c("cell_barcode" = "barcode")) %>%
      replace_na(list(color_column = "No CT Detected")) %>%
      mutate(color_column_yes_no = case_when(color_column == "No CT Detected" ~ "No",
                                             TRUE ~ "Yes")) %>%
      mutate(color_column = fct_infreq(color_column)) %>%
      arrange(color_column)

    if (reduction == "UMAP") {
      p <- ggplot(data = plot_df, aes(x = UMAP_1, y = UMAP_2))
      p <- p + labs(x = "UMAP 1", y = "UMAP 2")
    } else {
      p <- ggplot(data = plot_df, aes(x = tSNE_1, y = tSNE_2))
      p <- p + labs(x = "t-SNE 1", y = "t-SNE 2")
    }
  } else {
    stop("Reduction parameter reduction must be 'UMAP' or 't-SNE'.")
  }

  n_colors <- plot_df %>% pull(color_column) %>% levels() %>% length()

  if (color_yes_no) {
    p <- p + geom_point(aes(color = color_column_yes_no),
                        shape = 16,
                        size = 1.5,
                        alpha = 1) +
      labs(color = "Chimeric\nTranscript\nDetected") +
      scale_color_manual(values = c("#d9d9d9", "#ff00ff"))
  } else {
    p <- p + geom_point(aes(color = color_column),
                        shape = 16,
                        size = 1.5,
                        alpha = 1) +
      labs(color = "Cell Contains\nChimeric Transcript") +
      scale_color_brewer(palette = "Set1")
  }

  p <- p +
    theme_bw() +
    coord_equal() +
    scale_x_continuous(expand = c(0.01, 0.01)) +
    scale_y_continuous(expand = c(0.01, 0.01)) +
    theme(panel.background = element_blank(),
          panel.border = element_blank(),
          plot.background = element_blank(),
          panel.grid = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10))

  if (facet) {
    p <- p + facet_wrap(~ color_column, nrow = 1)
    ggsave(str_c(dir, "/", id, ".cells_chimeric_transcripts.", reduction, ".pdf"),
           p,
           width = n_colors*3.5 + 3.5, height = 3.5, useDingbats = FALSE)

    ggsave(str_c(dir, "/", id, ".cells_chimeric_transcripts.no_legend.", reduction, ".pdf"),
           p + guides(color = FALSE),
           width = n_colors*3.5, height = 3.5, useDingbats = FALSE)
  } else {
    ggsave(str_c(dir, "/", id, ".cells_chimeric_transcripts.", reduction, ".pdf"),
           p,
           width = 7, height = 3.5, useDingbats = FALSE)

    ggsave(str_c(dir, "/", id, ".cells_chimeric_transcripts.no_legend.", reduction, ".pdf"),
           p + guides(color = FALSE),
           width = 3.5, height = 3.5, useDingbats = FALSE)
  }
}

# ==============================================================================
# Plot correlation of two genes expression
# ==============================================================================
plot_two_genes_correlation <- function(two_genes_expression, bulk_sc, dir = paper_main, id, gene1, gene2){

  plot_df <- bulk_sc %>%
    filter(data_type == "Single Cell Chimeric Transcript") %>%
    separate(identifier,
             into = c("cell_barcode", "molecular_barcode"),
             sep = ":") %>%
    group_by(cell_barcode, data_type, chr14_position, chr4_position) %>%
    summarize(molecular_barcodes = str_c(molecular_barcode, collapse = ", ")) %>%
    right_join(two_genes_expression, by = c("cell_barcode" = "barcode")) %>%
    filter(cell_type == "Plasma Cells") %>%
    mutate(chimeric_transcript_detected = !is.na(data_type)) %>%
    arrange(chimeric_transcript_detected) %>%
    ungroup()

  gene12_correlation <- plot_df %>%
    filter(expression1 > 0, expression2 > 0) %>%
    select(expression1, expression2) %>%
    as.matrix() %>%
    cor()

  max_expr1 <- plot_df %>% pull(expression1) %>% max()
  max_expr2 <- plot_df %>% pull(expression2) %>% max()
  overall_max <- max(max_expr1, max_expr2)

  p <- ggplot(data = plot_df, aes(x = expression1, y = expression2))
  p <- p + geom_point(aes(color = chimeric_transcript_detected),
                      shape = 16,
                      size = 1.5,
                      show.legend = FALSE,
                      alpha = 1) +
    theme_bw() +
    coord_equal() +
    scale_x_continuous(limits = c(0, overall_max), expand = c(0.01, 0.01)) +
    scale_y_continuous(limits = c(0, overall_max), expand = c(0.01, 0.01)) +
    labs(x = str_c(gene1, " Expression"), y = str_c(gene2, " Expression")) +
    theme(panel.background = element_blank(),
          plot.background = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10))

  ggsave(str_c(dir, "/", id, ".", gene1, ".", gene2, ".correlation.pdf"),
         p,
         width = 3.5, height = 3.5, useDingbats = FALSE)
}

#General analysis both samples
if (TRUE) {

  # ==============================================================================
  # Plot cell types of each sample
  # UMAP in paper_main; t-SNE in paper_supp
  # ==============================================================================
  if (TRUE) {
    plot_cell_types(cell_types_27522_1, seurat_object_27522_1, id = "27522_1", dir = str_c(paper_main, "27522_1"))
    plot_cell_types(cell_types_27522_4, seurat_object_27522_4, id = "27522_4", dir = str_c(paper_main, "27522_4"))

    plot_cell_types(cell_types_27522_1, seurat_object_27522_1, id = "27522_1", reduction = "t-SNE", dir = str_c(paper_supp, "27522_1"))
    plot_cell_types(cell_types_27522_4, seurat_object_27522_4, id = "27522_4", reduction = "t-SNE", dir = str_c(paper_supp, "27522_4"))
  }

  # ==============================================================================
  # Plot interesting gene expressions
  # 27522: t(4;14) WHSC1 (ENSG00000109685) and FGFR3 (ENSG00000068078)
  # ==============================================================================
  if (TRUE) {
    for (red in c("UMAP", "t-SNE")) {
      t414_gene_names <- c("WHSC1", "FGFR3", "SDC1", "TNFRSF17", "SLAMF7")
      t414_ensg <- c("ENSG00000109685", "ENSG00000068078", "ENSG00000115884", "ENSG00000048462", "ENSG00000026751")
      for (i in 1:length(t414_gene_names)) {
        plot_gene_expression(seurat_object_27522_1,
                             get_tsne_umap(cell_types = cell_types_27522_1,
                                           seurat_object = seurat_object_27522_1),
                             ensg = t414_ensg[i], id = "27522_1", gene = t414_gene_names[i],
                             dir = str_c(paper_supp, "27522_1"), reduction = red)
        plot_gene_expression(seurat_object_27522_4,
                             get_tsne_umap(cell_types = cell_types_27522_4,
                                           seurat_object = seurat_object_27522_4),
                             ensg = t414_ensg[i], id = "27522_4", gene = t414_gene_names[i],
                             dir = str_c(paper_supp, "27522_4"), reduction = red)
      }
    }
  }

}

# ==============================================================================
# Work with 27522
# ==============================================================================
if (TRUE) {
  get_bulk_fusion_reads_27522 <- function(bulk_reads, star_fusion = NULL, use_SF_only = FALSE){

    return_fusions <- function(this_read, star_fusion){
      if (identical(star_fusion, NULL)) {
        return("No STAR-Fusion Results Provided")
      }

      fusions <- star_fusion %>%
        filter(str_detect(JunctionReads, this_read) |
                 str_detect(SpanningFrags, this_read)) %>%
        pull(`#FusionName`) %>%
        unique() %>%
        sort() %>%
        str_c(collapse = ", ")
      if (identical(fusions, character(0))) {
        return("Not Associated With Any Fusion")
      } else {
        return(fusions)
      }
    }

    all_reads <- bulk_reads %>%
      filter(chromosome_donor %in% c("chr4", "chr14"),
             chromosome_acceptor %in% c("chr4", "chr14")) %>%
      select(read_name,
             chromosome_donor, first_base_donor,
             chromosome_acceptor, first_base_acceptor) %>%
      unique() %>%
      mutate(chr14_position = case_when(chromosome_donor == "chr14" ~ first_base_donor,
                                        TRUE ~ first_base_acceptor),
             chr4_position = case_when(chromosome_donor == "chr4" ~ first_base_donor,
                                       TRUE ~ first_base_acceptor)) %>%
      filter(chr14_position > 105500000,
             chr14_position < 107000000) %>%
      filter(chr4_position > 1800000,
             chr4_position < 2000000) %>%
      select(read_name, chr14_position, chr4_position) %>%
      rowwise() %>%
      mutate(fusion = return_fusions(read_name, star_fusion)) %>%
      ungroup()

    if (use_SF_only) {
      all_reads %>%
        filter(fusion != "Not Associated With Any Fusion") %>%
        return()
    } else {
      return(all_reads)
    }
  }
  get_sc_chimeric_transcripts_27522 <- function(dis_reads){

    sc_chr14_min_max <- dis_reads %>%
      filter(chrom == 14) %>%
      filter(end - start < 128) %>%
      group_by(cell_barcode, molecular_barcode) %>%
      summarize(min_start = min(start), max_start = max(start),
                min_end = min(end), max_end = max(end),
                n_reads = n())

    sc_chr4_min_max <- dis_reads %>%
      filter(chrom == 4) %>%
      filter(end - start < 128) %>%
      group_by(cell_barcode, molecular_barcode) %>%
      summarize(min_start = min(start), max_start = max(start),
                min_end = min(end), max_end = max(end),
                n_reads = n())

    sc_chr14_min_max %>%
      full_join(sc_chr4_min_max,
                by = c("cell_barcode", "molecular_barcode")) %>%
      filter(!any(is.na(min_start.x), is.na(min_start.y))) %>%
      return()
  }
  get_bulk_sc_plot_df_27522 <- function(bulk_fusion_reads, sc_chimeric_transcripts){

    chr4_breakpoint <- 1871964
    chr14_breakpoint <- 105858090
    between_genes <- 0.2e6

    bulk_fusion_reads <- bulk_fusion_reads %>%
      filter(chr14_position >= 105835000, chr14_position <= 105935000)

    sc_chimeric_transcripts <- sc_chimeric_transcripts %>%
      filter(min_start.x >= 105835000, max_end.x <= 105935000)

    bulk_chr14_min <- bulk_fusion_reads %>% pull(chr14_position) %>% min()
    bulk_chr14_max <- bulk_fusion_reads %>% pull(chr14_position) %>% max()
    bulk_chr4_min <- bulk_fusion_reads %>% pull(chr4_position) %>% min()
    bulk_chr4_max <- bulk_fusion_reads %>% pull(chr4_position) %>% max()

    sc_chr14_min <- sc_chimeric_transcripts %>% pull(min_start.x) %>% min()
    sc_chr14_max <- sc_chimeric_transcripts %>% pull(min_start.x) %>% max()
    sc_chr4_min <- sc_chimeric_transcripts %>% pull(min_start.y) %>% min()
    sc_chr4_max <- sc_chimeric_transcripts %>% pull(min_start.y) %>% max()

    overall_chr14_min <- min(bulk_chr14_min, sc_chr14_min)
    overall_chr14_max <- max(bulk_chr14_max, sc_chr14_max)
    overall_chr4_min <- min(bulk_chr4_min, sc_chr4_min)
    overall_chr4_max <- max(bulk_chr4_max, sc_chr4_max)

    chr4_shift = 0 - overall_chr4_min + overall_chr14_max - overall_chr14_min + between_genes
    chr14_shift = 0 - overall_chr14_min

    shifted_chr4_breakpoint = chr4_breakpoint + chr4_shift
    shifted_chr14_breakpoint = chr14_breakpoint + chr14_shift

    bulk_plot_df <- bulk_fusion_reads %>%
      mutate(category = case_when(chr4_position <= chr4_breakpoint & chr14_position <= chr14_breakpoint ~ "Supports Reciprocal FGFR3/WHSC1--IGH Fusion",
                                  chr4_position > chr4_breakpoint & chr14_position <= chr14_breakpoint ~ "Supports General IGH Fusion with WHSC1",
                                  chr4_position <= chr4_breakpoint & chr14_position > chr14_breakpoint ~ "Supports General IGH Fusion with WHSC1",
                                  chr4_position > chr4_breakpoint & chr14_position > chr14_breakpoint ~ "Supports IGH--WHSC1 Fusion")) %>%
      mutate(category = factor(category, levels = c("Supports IGH--WHSC1 Fusion", "Supports Reciprocal FGFR3/WHSC1--IGH Fusion", "Supports General IGH Fusion with WHSC1"), ordered = TRUE)) %>%
      mutate(shifted_chr4 = chr4_position + chr4_shift,
             shifted_chr14 = chr14_position + chr14_shift) %>%
      mutate(identifier = read_name,
             data_type = "Bulk Spanning/Junction Read Pair",
             y_position_offset = 0,
             curvature = -0.25) %>%
      select(identifier, data_type, chr14_position, chr4_position, shifted_chr14, shifted_chr4, category, fusion)

    sc_plot_df <- sc_chimeric_transcripts %>%
      mutate(category = case_when((min_start.x <= chr14_breakpoint & max_end.x > chr14_breakpoint) | (min_start.y < chr4_breakpoint & max_end.y > chr4_breakpoint) ~ "Supports General IGH Fusion with WHSC1",
                                  min_start.x > chr14_breakpoint & min_start.y > chr4_breakpoint ~ "Supports IGH--WHSC1 Fusion",
                                  max_end.x <= chr14_breakpoint & max_end.y <= chr4_breakpoint ~ "Supports Reciprocal FGFR3/WHSC1--IGH Fusion",
                                  max_end.x <= chr14_breakpoint & min_start.y > chr4_breakpoint ~ "Supports General IGH Fusion with WHSC1",
                                  TRUE ~ "Supports General IGH Fusion with WHSC1")) %>%
      mutate(chr4_position = case_when(category == "Supports IGH--WHSC1 Fusion" ~ min_start.y,
                                       category == "Supports Reciprocal FGFR3/WHSC1--IGH Fusion" ~ max_end.y,
                                       category == "Supports General IGH Fusion with WHSC1" ~ min_start.y),
             chr14_position = case_when(category == "Supports IGH--WHSC1 Fusion" ~ min_start.x,
                                        category == "Supports Reciprocal FGFR3/WHSC1--IGH Fusion" ~ max_end.x,
                                        category == "Supports General IGH Fusion with WHSC1" ~ min_start.x)) %>%
      mutate(category = factor(category, levels = c("Supports IGH--WHSC1 Fusion", "Supports Reciprocal FGFR3/WHSC1--IGH Fusion", "Supports General IGH Fusion with WHSC1"), ordered = TRUE)) %>%
      mutate(shifted_chr4 = chr4_position + chr4_shift,
             shifted_chr14 = chr14_position + chr14_shift) %>%
      mutate(identifier = str_c(cell_barcode, ":", molecular_barcode),
             data_type = "Single Cell Chimeric Transcript",
             y_position_offset = -0.5,
             curvature = 0.25,
             fusion = NA) %>%
      ungroup() %>%
      select(identifier, data_type, chr14_position, chr4_position, shifted_chr14, shifted_chr4, category, fusion)

    bind_rows(bulk_plot_df, sc_plot_df) %>%
      mutate(color_column = category) %>%
      filter(chr14_position >= 105835000, chr14_position <= 105935000) %>%
      return()
  }
  plot_bulk_sc_27522 <- function(bulk_sc_plot_df, genes = gene_spans, dir, id, use_fusion_names = TRUE, bulk_color = TRUE, color_indicates = "Fusion Name"){

    reported_translocation_breakpoints <- tribble(~trans, ~chrom,      ~pos,
                                                  "t(4;14)",      4,   1871964,
                                                  "t(4;14)",     14, 105858090)
    chr4_breakpoint <- 1871964
    chr14_breakpoint <- 105858090
    between_genes <- 0.2e6
    chr4_offset_y <- 0.5

    bulk_sc_plot_df <- bulk_sc_plot_df %>%
      mutate(color_column = case_when(!bulk_color ~ "black",
                                      use_fusion_names ~ str_remove_all(str_remove_all(fusion, "NSD2"), "--"),
                                      TRUE ~ as.character(str_detect(fusion, "IGH"))),
             reciprocal = str_detect(fusion, "^NSD2"))

    overall_chr14_min <- bulk_sc_plot_df %>% pull(chr14_position) %>% min()
    overall_chr14_max <- bulk_sc_plot_df %>% pull(chr14_position) %>% max()
    overall_chr4_min <- bulk_sc_plot_df %>% pull(chr4_position) %>% min()
    overall_chr4_max <- bulk_sc_plot_df %>% pull(chr4_position) %>% max()

    chr4_shift = 0 - overall_chr4_min + overall_chr14_max - overall_chr14_min + between_genes
    chr14_shift = 0 - overall_chr14_min

    shifted_chr4_breakpoint = chr4_breakpoint + chr4_shift
    shifted_chr14_breakpoint = chr14_breakpoint + chr14_shift

    chr14_gene_spans <- genes %>%
      filter(chromosome == "chr14") %>%
      filter( (start <= overall_chr14_min & end >= overall_chr14_min) |
                (start >= overall_chr14_min & end <= overall_chr14_max) |
                (start <= overall_chr14_max & end >= overall_chr14_max) ) %>%
      filter(str_detect(gene_name, "IGH")) %>%
      filter(!str_detect(gene_name, "@")) %>%
      mutate(shifted_chr14_start = start + chr14_shift,
             shifted_chr14_end = end + chr14_shift)

    chr4_gene_spans <- genes %>%
      filter(chromosome == "chr4") %>%
      filter( gene_name %in% c("FGFR3", "WHSC1") |
                (start <= overall_chr4_min & end >= overall_chr4_min) |
                (start >= overall_chr4_min & end <= overall_chr4_max) |
                (start <= overall_chr4_max & end >= overall_chr4_max) ) %>%
      filter(strand == "+", type == "protein_coding") %>%
      mutate(shifted_chr4_start = start + chr4_shift,
             shifted_chr4_end = end + chr4_shift) %>%
      mutate(gene_name = case_when(gene_name == "NSD2" ~ "WHSC1",
                                   TRUE ~ gene_name))

    shifted_gene_boundaries <- c(chr14_gene_spans %>% pull(shifted_chr14_start) %>% min(),
                                 chr14_gene_spans %>% pull(shifted_chr14_end) %>% max(),
                                 chr4_gene_spans %>% pull(shifted_chr4_start) %>% min(),
                                 chr4_gene_spans %>% pull(shifted_chr4_start) %>% max(),
                                 chr4_gene_spans %>% pull(shifted_chr4_end) %>% max())

    gene_boundaries <- c(chr14_gene_spans %>% pull(start) %>% min(),
                         chr14_gene_spans %>% pull(end) %>% max(),
                         chr4_gene_spans %>% pull(start) %>% min(),
                         chr4_gene_spans %>% pull(start) %>% max(),
                         chr4_gene_spans %>% pull(end) %>% max())

    p <- ggplot(data = bulk_sc_plot_df) +
      scale_y_continuous(limits = c(-1.5, 2)) +
      scale_x_continuous(breaks = shifted_gene_boundaries,
                         labels = gene_boundaries) +
      #chr14 genes
      geom_rect(data = chr14_gene_spans,
                aes(xmin = shifted_chr14_start,
                    xmax = shifted_chr14_end,
                    ymin = 0,
                    ymax = 1,
                    fill = "chr14"),
                #fill = "#f03b20"
                alpha = 1) +
      #chr4 genes
      geom_rect(data = chr4_gene_spans,
                aes(xmin = shifted_chr4_start,
                    xmax = shifted_chr4_end,
                    ymin = 0 - chr4_offset_y,
                    ymax = 1 - chr4_offset_y,
                    fill = "chr4"),
                #fill = "#feb24c",
                alpha = 1) +
      # draw bulk read curves
      geom_curve(data = bulk_sc_plot_df %>%
                   filter(data_type == "Bulk Spanning/Junction Read Pair") %>%
                   arrange(str_detect(fusion, "IGH")),
                 aes(x = shifted_chr14,
                     xend = shifted_chr4,
                     color = color_column,
                     linetype = reciprocal,
                     y = 1,
                     yend = 1 - chr4_offset_y),
                 curvature = -0.25,
                 ncp = 10,
                 lineend = "round",
                 alpha = 0.25) +
      # draw sc read curves
      geom_curve(data = bulk_sc_plot_df %>%
                   filter(data_type == "Single Cell Chimeric Transcript"),
                 aes(x = shifted_chr14,
                     xend = shifted_chr4,
                     y = 0,
                     yend = 0 - chr4_offset_y,
                     color = "black"),
                 curvature = 0.25,
                 ncp = 10,
                 lineend = "round",
                 show.legend = FALSE,
                 alpha = 0.25) +
      # label chr4 genes
      geom_text(data = chr4_gene_spans,
                aes(x = (shifted_chr4_end + shifted_chr4_start)/2,
                    y = 1 - chr4_offset_y - 0.5,
                    label = gene_name),
                fontface = "italic") +
      # translocation breakpoints
      geom_vline(xintercept = shifted_chr4_breakpoint,
                 linetype = 2) +
      geom_vline(xintercept = shifted_chr14_breakpoint,
                 linetype = 2) +
      annotate(geom = "text", x = shifted_chr4_breakpoint + between_genes*0.02, y = -1.5, label = chr4_breakpoint, hjust = 0, vjust = 1) +
      annotate(geom = "text", x = shifted_chr14_breakpoint - between_genes*0.02, y = -1.5, label = chr14_breakpoint, hjust = 1, vjust = 1) +
      # facets
      facet_wrap(~ category, ncol = 1) +
      # colors
      scale_color_brewer(palette = "Paired") +
      scale_fill_brewer(palette = "Accent") +
      labs(color = color_indicates,
           linetype = "Reciprocal Fusion") +
      theme_bw() +
      theme(panel.background = element_blank(),
            panel.grid.minor = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.border = element_blank(),
            plot.background = element_blank(),
            axis.title = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),
            legend.background = element_blank(),
            legend.position = "bottom")

    ggsave(str_c(dir, "/", id, ".discordant_read_positions.pdf"),
           p,
           width = 14.5, height = 12,
           useDingbats = FALSE)

    ggsave(str_c(dir, "/", id, ".discordant_read_positions.no_legend.pdf"),
           p + guides(color = FALSE, linetype = FALSE, fill = FALSE),
           width = 7.25, height = 12,
           useDingbats = FALSE)
  }
  plot_bulk_sc_27522_main_figure <- function(bulk_sc_plot_df, genes = gene_spans, dir, id, use_fusion_names = TRUE, bulk_color = TRUE, color_indicates = "Fusion Name", plot_star_fusion_reads = FALSE){

    whsc1_tx_hg38 <- read_tsv("data/WHSC1_transcripts_UCSC_GenomeBrowser.hg38.txt")
    whsc1_tx_hg38 <- whsc1_tx_hg38 %>% select(start, stop) %>% unique()
    min_whsc1_hg38 <- whsc1_tx_hg38 %>% pull(start) %>% min()
    max_whsc1_hg38 <- whsc1_tx_hg38 %>% pull(stop) %>% max()

    reported_translocation_breakpoints <- tribble(   ~trans, ~chrom,      ~pos,
                                                  "t(4;14)",      4,   1871964,
                                                  "t(4;14)",     14, 105858090)
    chr4_breakpoint <- 1871964
    chr14_breakpoint <- 105858090

    bulk_sc_plot_df <- bulk_sc_plot_df %>%
      mutate(color_column = case_when(!bulk_color ~ "black",
                                      use_fusion_names ~ str_remove_all(str_remove_all(fusion, "NSD2"), "--"),
                                      TRUE ~ as.character(str_detect(fusion, "IGH"))),
             reciprocal = str_detect(fusion, "^NSD2"))

    overall_chr14_min <- bulk_sc_plot_df %>% pull(chr14_position) %>% min()
    overall_chr14_max <- bulk_sc_plot_df %>% pull(chr14_position) %>% max()
    overall_chr4_min <- bulk_sc_plot_df %>% pull(chr4_position) %>% min()
    overall_chr4_max <- bulk_sc_plot_df %>% pull(chr4_position) %>% max()

    chr4_shift = 0 - chr4_breakpoint
    chr14_shift = 0 - chr14_breakpoint

    shifted_chr4_breakpoint = 0
    shifted_chr14_breakpoint = 0

    chr14_gene_spans <- genes %>%
      filter(chromosome == "chr14") %>%
      filter( (start <= overall_chr14_min & end >= overall_chr14_min) |
                (start >= overall_chr14_min & end <= overall_chr14_max) |
                (start <= overall_chr14_max & end >= overall_chr14_max) ) %>%
      filter(str_detect(gene_name, "IGH")) %>%
      filter(!str_detect(gene_name, "@")) %>%
      mutate(shifted_chr14_start = start + chr14_shift,
             shifted_chr14_end = end + chr14_shift)

    chr4_gene_spans <- genes %>%
      filter(chromosome == "chr4") %>%
      filter((start <= overall_chr4_min & end >= overall_chr4_min) |
               (start >= overall_chr4_min & end <= overall_chr4_max) |
               (start <= overall_chr4_max & end >= overall_chr4_max) ) %>%
      filter(strand == "+", type == "protein_coding") %>%
      mutate(shifted_chr4_start = start + chr4_shift,
             shifted_chr4_end = end + chr4_shift) %>%
      mutate(gene_name = case_when(gene_name == "NSD2" ~ "WHSC1",
                                   TRUE ~ gene_name))

    shifted_gene_boundaries <- c(chr14_gene_spans %>% pull(shifted_chr14_start) %>% min(),
                                 chr14_gene_spans %>% pull(shifted_chr14_end) %>% max(),
                                 chr4_gene_spans %>% pull(shifted_chr4_start) %>% min(),
                                 chr4_gene_spans %>% pull(shifted_chr4_start) %>% max(),
                                 chr4_gene_spans %>% pull(shifted_chr4_end) %>% max())

    gene_boundaries <- c(chr14_gene_spans %>% pull(start) %>% min(),
                         chr14_gene_spans %>% pull(end) %>% max(),
                         chr4_gene_spans %>% pull(start) %>% min(),
                         chr4_gene_spans %>% pull(start) %>% max(),
                         chr4_gene_spans %>% pull(end) %>% max())

    p <- ggplot(data = bulk_sc_plot_df %>% arrange(data_type)) +
      coord_cartesian(ylim = c(-2, 2)) +
      scale_y_continuous(limits = c(-2, 2)) +
      scale_x_continuous(breaks = shifted_gene_boundaries,
                         labels = gene_boundaries) +
      #chr14 genes
      geom_rect(data = chr14_gene_spans,
                aes(xmin = shifted_chr14_start,
                    xmax = shifted_chr14_end,
                    ymin = 0.5,
                    ymax = 1.75),
                fill = "#bdbdbd",
                alpha = 1) +
      #chr4 genes
      geom_segment(x = min_whsc1_hg38 + chr4_shift,
                   xend = max_whsc1_hg38 + chr4_shift,
                   y = -1.125, yend = -1.125,
                   size = 5,
                   color = "#d9d9d9",
                   alpha = 0.5) +
      geom_rect(data = whsc1_tx_hg38,
                aes(ymin = -1.75,
                    ymax = -0.5,
                    xmin = start + chr4_shift,
                    xmax = stop + chr4_shift),
                inherit.aes = FALSE,
                color = NA,
                fill = "#bdbdbd")

    p <- p + geom_curve(#data = bulk_sc_plot_df, # %>% arrange(data_type),
                        aes(x = chr14_position + chr14_shift,
                            xend = chr4_position + chr4_shift,
                            linetype = data_type,
                            y = 0.5,
                            yend = -0.5,
                            color = data_type),
                        curvature = 0,
                        ncp = 10,
                        #alpha = 0.5,
                        lineend = "round")

    # label chr4 genes
    p <- p + geom_text(data = chr4_gene_spans,
                       aes(x = (shifted_chr4_end + shifted_chr4_start)/2,
                           y = -1.125,
                           label = gene_name),
                       fontface = "italic") +
      # translocation breakpoints
      geom_vline(xintercept = 0, linetype = 2) +
      annotate(geom = "text", x = 0 + 1e3, y = -2, label = chr4_breakpoint, hjust = 0, vjust = 0) +
      annotate(geom = "text", x = 0 + 1e3, y = 2, label = chr14_breakpoint, hjust = 0, vjust = 1) +
      annotate(geom = "text", x = 0 - 1e3, y = -2, label = "chr4", hjust = 1, vjust = 0) +
      annotate(geom = "text", x = 0 - 1e3, y = 2, label = "chr14", hjust = 1, vjust = 1) +
      # transcription direction arrows
      geom_curve(inherit.aes = FALSE, x = 6e3, xend = 6e3, y = 0.45, yend = -0.45, curvature = 1, ncp = 10, lineend = "round", show.legend = FALSE, arrow = arrow(angle = 45, length = unit(0.05, "inches"), type = "open"), color = "#000000") +
      geom_curve(inherit.aes = FALSE, x = -6e3, xend = -6e3, y = -0.45, yend = 0.45, curvature = 1, ncp = 10, lineend = "round", show.legend = FALSE, arrow = arrow(angle = -45, length = unit(0.05, "inches"), type = "open"), color = "#000000") +
      annotate(geom = "text", x = -6.75e3, y = 0, label = "Direction of\nTranscription", size = 2, color = "#000000") +
      # colors
      scale_linetype_manual(values = c(2,1)) +
      scale_color_manual(values = c("#a6cee3", "#33a02c")) +
      scale_fill_brewer(palette = "Accent") +
      labs(y = NULL,
           x = "Genomic Coordinates (GRCh38)",
           color = "Data Source",
           linetype = "Data Source") +
      theme_bw() +
      theme(panel.background = element_blank(),
            panel.grid.minor = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.border = element_blank(),
            plot.background = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            legend.background = element_blank(),
            legend.position = "bottom")

    ggsave(str_c(dir, "/", id, ".discordant_read_positions.main.pdf"),
           p,
           width = 14.5, height = 12,
           useDingbats = FALSE)

    ggsave(str_c(dir, "/", id, ".discordant_read_positions.main.no_legend.pdf"),
           p + guides(color = FALSE, linetype = FALSE, fill = FALSE),
           width = 7.25, height = 2.75,
           useDingbats = FALSE)
  }

  # 27522_1
  if (TRUE) {

    bulk_fusion_reads_27522_1_SF_only <- get_bulk_fusion_reads_27522(bulk_reads_27522_1, star_fusion_calls_27522_1, use_SF_only = TRUE)
    bulk_fusion_reads_27522_1 <- get_bulk_fusion_reads_27522(bulk_reads_27522_1, star_fusion_calls_27522_1, use_SF_only = FALSE)
    plot_sc_read_length_distribution(dis_reads = dis_reads_27522_1_discover, cutoff = 128, id = "27522_1", dir = str_c(paper_supp, "27522_1"))
    sc_chimeric_transcripts_27522_1 <- get_sc_chimeric_transcripts_27522(dis_reads_27522_1_discover)
    bulk_sc_plot_df_27522_1_SF_only <- get_bulk_sc_plot_df_27522(bulk_fusion_reads_27522_1_SF_only, sc_chimeric_transcripts_27522_1)
    bulk_sc_plot_df_27522_1 <- get_bulk_sc_plot_df_27522(bulk_fusion_reads_27522_1, sc_chimeric_transcripts_27522_1)

    plot_bulk_sc_27522_main_figure(bulk_sc_plot_df_27522_1_SF_only, id = "27522_1", bulk_color = FALSE, dir = str_c(paper_main, "27522_1"), plot_star_fusion_reads = FALSE)
    plot_bulk_sc_27522(bulk_sc_plot_df_27522_1_SF_only, id = "27522_1", bulk_color = FALSE, dir = str_c(paper_supp, "27522_1"))

    plot_cell_chimeric_transcripts(bulk_sc = bulk_sc_plot_df_27522_1_SF_only,
                                   tsne_umap = get_tsne_umap(cell_types = cell_types_27522_1, seurat_object = seurat_object_27522_1),
                                   id = "27522_1",
                                   reduction = "UMAP",
                                   dir = str_c(paper_main, "27522_1"))
    plot_cell_chimeric_transcripts(bulk_sc = bulk_sc_plot_df_27522_1_SF_only,
                                   tsne_umap = get_tsne_umap(cell_types = cell_types_27522_1, seurat_object = seurat_object_27522_1),
                                   id = "27522_1",
                                   reduction = "t-SNE",
                                   dir = str_c(paper_supp, "27522_1"))

    plot_gene_expression_violin(bulk_sc = bulk_sc_plot_df_27522_1,
                                seurat_object = seurat_object_27522_1,
                                tsne_umap = get_tsne_umap(cell_types = cell_types_27522_1,
                                                          seurat_object = seurat_object_27522_1),
                                ensg = "ENSG00000068078",
                                id = "27522_1",
                                gene = "FGFR3",
                                dir = str_c(paper_supp, "27522_1"),
                                max = 4.4)
    plot_gene_expression_violin(bulk_sc = bulk_sc_plot_df_27522_1,
                                seurat_object = seurat_object_27522_1,
                                tsne_umap = get_tsne_umap(cell_types = cell_types_27522_1, seurat_object = seurat_object_27522_1),
                                ensg = "ENSG00000109685",
                                id = "27522_1",
                                gene = "WHSC1",
                                dir = str_c(paper_supp, "27522_1"),
                                max = 4.4)
    # https://www.nature.com/articles/nrc2780/figures/2
    plot_gene_expression_violin(bulk_sc = bulk_sc_plot_df_27522_1,
                                seurat_object = seurat_object_27522_1,
                                tsne_umap = get_tsne_umap(cell_types = cell_types_27522_1, seurat_object = seurat_object_27522_1),
                                ensg = "ENSG00000142208",
                                id = "27522_1",
                                gene = "AKT1",
                                dir = str_c(paper_supp, "27522_1"),
                                max = 4.4)
    plot_gene_expression_violin(bulk_sc = bulk_sc_plot_df_27522_1,
                                seurat_object = seurat_object_27522_1,
                                tsne_umap = get_tsne_umap(cell_types = cell_types_27522_1, seurat_object = seurat_object_27522_1),
                                ensg = "ENSG00000100030",
                                id = "27522_1",
                                gene = "MAPK1",
                                dir = str_c(paper_supp, "27522_1"),
                                max = 4.4)

    plot_two_genes_correlation(two_genes_expression = get_two_genes_expression(seurat_object = seurat_object_27522_1,
                                                                               tsne_umap = get_tsne_umap(cell_types = cell_types_27522_1, seurat_object = seurat_object_27522_1),
                                                                               ensg1 = "ENSG00000109685",
                                                                               ensg2 = "ENSG00000068078"),
                               bulk_sc = bulk_sc_plot_df_27522_1,
                               dir = str_c(paper_supp, "27522_1"),
                               id = "27522_1",
                               gene1 = "WHSC1",
                               gene2 = "FGFR3")
  }

  # 27522_4
  if (TRUE) {
    bulk_fusion_reads_27522_4 <- get_bulk_fusion_reads_27522(bulk_reads_27522_4)
    plot_sc_read_length_distribution(dis_reads = dis_reads_27522_4_discover, cutoff = 128, id = "27522_4", dir = str_c(paper_supp, "27522_4"))
    sc_chimeric_transcripts_27522_4 <- get_sc_chimeric_transcripts_27522(dis_reads_27522_4_discover)
    bulk_sc_plot_df_27522_4 <- get_bulk_sc_plot_df_27522(bulk_fusion_reads_27522_4, sc_chimeric_transcripts_27522_4)
    plot_bulk_sc_27522(bulk_sc_plot_df_27522_4, id = "27522_4", bulk_color = FALSE, dir = str_c(paper_main, "27522_4"))

    plot_cell_chimeric_transcripts(bulk_sc = bulk_sc_plot_df_27522_4,
                                   tsne_umap = get_tsne_umap(cell_types = cell_types_27522_4, seurat_object = seurat_object_27522_4),
                                   id = "27522_4",
                                   reduction = "UMAP",
                                   dir = str_c(paper_main, "27522_4"))
    plot_cell_chimeric_transcripts(bulk_sc = bulk_sc_plot_df_27522_4,
                                   tsne_umap = get_tsne_umap(cell_types = cell_types_27522_4, seurat_object = seurat_object_27522_4),
                                   id = "27522_4",
                                   reduction = "t-SNE",
                                   dir = str_c(paper_supp, "27522_4"))

    plot_gene_expression_violin(bulk_sc = bulk_sc_plot_df_27522_4,
                                seurat_object = seurat_object_27522_4,
                                tsne_umap = get_tsne_umap(cell_types = cell_types_27522_4,
                                                          seurat_object = seurat_object_27522_4),
                                ensg = "ENSG00000068078",
                                id = "27522_4",
                                gene = "FGFR3",
                                dir = str_c(paper_supp, "27522_4"),
                                max = 4.4)
    plot_gene_expression_violin(bulk_sc = bulk_sc_plot_df_27522_4,
                                seurat_object = seurat_object_27522_4,
                                tsne_umap = get_tsne_umap(cell_types = cell_types_27522_4, seurat_object = seurat_object_27522_4),
                                ensg = "ENSG00000109685",
                                id = "27522_4",
                                gene = "WHSC1",
                                dir = str_c(paper_supp, "27522_4"),
                                max = 4.4)
    # https://www.nature.com/articles/nrc2780/figures/2
    plot_gene_expression_violin(bulk_sc = bulk_sc_plot_df_27522_4,
                                seurat_object = seurat_object_27522_4,
                                tsne_umap = get_tsne_umap(cell_types = cell_types_27522_4, seurat_object = seurat_object_27522_4),
                                ensg = "ENSG00000142208",
                                id = "27522_4",
                                gene = "AKT1",
                                dir = str_c(paper_supp, "27522_4"),
                                max = 4.4)
    plot_gene_expression_violin(bulk_sc = bulk_sc_plot_df_27522_4,
                                seurat_object = seurat_object_27522_4,
                                tsne_umap = get_tsne_umap(cell_types = cell_types_27522_4, seurat_object = seurat_object_27522_4),
                                ensg = "ENSG00000100030",
                                id = "27522_4",
                                gene = "MAPK1",
                                dir = str_c(paper_supp, "27522_4"),
                                max = 4.4)

    plot_two_genes_correlation(two_genes_expression = get_two_genes_expression(seurat_object = seurat_object_27522_4,
                                                                               tsne_umap = get_tsne_umap(cell_types = cell_types_27522_4, seurat_object = seurat_object_27522_4),
                                                                               ensg1 = "ENSG00000109685",
                                                                               ensg2 = "ENSG00000068078"),
                               bulk_sc = bulk_sc_plot_df_27522_4,
                               dir = str_c(paper_supp, "27522_4"),
                               id = "27522_4",
                               gene1 = "WHSC1",
                               gene2 = "FGFR3")

    # ============================================================================
    # Plot chr3 CNV in 25722_4
    # ============================================================================
    plot_chr_cnv(infercnv_27522_4,
                 get_tsne_umap(cell_types = cell_types_27522_4,
                               seurat_object = seurat_object_27522_4),
                 id = "27522_4", chr = "chr3",
                 dir = str_c(paper_supp, "27522_4"))

    # ============================================================================
    # Plot chr13 CNV in 25722_4
    # ============================================================================
    plot_chr_cnv(infercnv_27522_4,
                 get_tsne_umap(cell_types = cell_types_27522_4,
                               seurat_object = seurat_object_27522_4),
                 id = "27522_4", chr = "chr13",
                 dir = str_c(paper_supp, "27522_4"))

    # ============================================================================
    # Plot chr16 CNV in 25722_4
    # ============================================================================
    plot_chr_cnv(infercnv_27522_4,
                 get_tsne_umap(cell_types = cell_types_27522_4,
                               seurat_object = seurat_object_27522_4),
                 id = "27522_4", chr = "chr16",
                 dir = str_c(paper_supp, "27522_4"))


  }
}

# ==============================================================================
# Explore QC of samples
# ==============================================================================
if (TRUE) {
  bad_df <- read_tsv("data/bad_qc_reads_plot_df.tsv") %>%
    filter(chromosome == 14) %>%
    group_by(cell_barcode, molecular_barcode, sample, test, chromosome) %>%
    summarize(start = min(start), end = max(end)) %>%
    ungroup() %>%
    separate(test, sep = "_", into = c("ig", "other_gene")) %>%
    mutate(ig_name = case_when(chromosome == 2 ~ "IGK",
                               chromosome == 14 ~ "IGH",
                               chromosome == 22 ~ "IGL"))

  original_df <- bind_rows(dis_reads_27522_1_discover_preQC %>% mutate(sample = "27522_1")) %>%
    filter(end - start < 128) %>%
    group_by(sample, cell_barcode, molecular_barcode, chrom) %>%
    summarize(start = min(start), end = max(end)) %>%
    mutate(other_gene = "original") %>%
    filter(chrom == 14) %>%
    mutate(ig_name = case_when(chrom == 2 ~ "IGK",
                               chrom == 14 ~ "IGH",
                               chrom == 22 ~ "IGL"))

  ggplot(original_df, aes(x = start/1e6, fill = sample)) +
    geom_histogram(binwidth = 10000/1e6) +
    facet_wrap(~ ig_name, ncol = 1, scales = "free") +
    labs(fill = "Sample ID",
         x = "Position (Mb)",
         y = "Number of 'Chimeric Transcripts'") +
    theme_bw() +
    ggsave(str_c(paper_supp, "pre_QC_discovery.pdf"), width = 10, height = 20, useDingbats = FALSE)


  ggplot(bad_df, aes(x = start/1e6, fill = other_gene)) +
    geom_histogram(binwidth = 10000/1e6) +
    facet_wrap(~ ig_name, ncol = 1, scales = "free") +
    labs(fill = "'Partner' Gene",
         x = "Position (Mb)",
         y = "Number of 'Chimeric Transcripts'") +
    theme_bw() #+

  bind_rows(bad_df %>%
              select(other_gene, start) %>%
              rename("my_color" = "other_gene") %>%
              mutate(my_facet = "IGH with Other Genes (multiple samples)"),
            original_df %>% ungroup() %>%
              select(start) %>%
              mutate(my_color = "WHSC1",
                     my_facet = "IGH with WHSC1 (27522_1 only)")) %>%
    ggplot(aes(x = start/1e6, fill = my_color)) +
    geom_vline(xintercept = 105854501/1e6, linetype = 2, color = "#bdbdbd") +
    geom_vline(xintercept = 105922264/1e6, linetype = 2, color = "#bdbdbd") +
    geom_histogram(binwidth = 10000/1e6) +
    facet_wrap( ~ my_facet, ncol = 1, scales = "free_y") +
    labs(fill = "'Partner' Gene",
         x = "IGH chr14 Genomic Coordinates (Mb) (GRCh38)",
         y = "Number of 'Chimeric Transcripts'") +
    theme_bw() +
    theme(panel.background = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.border = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title = element_text(size = 12),
          axis.text.x = element_text(size = 8),
          strip.background = element_blank(),
          strip.text = element_text(size = 12)) +
    ggsave(str_c(paper_supp, "single_cell_QC.pdf"), width = 7.25, height = 4.0, useDingbats = FALSE)

}

# ==============================================================================
# Fusion single cells paragraph output
# ==============================================================================
n_cells_27522_1 <- get_tsne_umap(cell_types = cell_types_27522_1,
                                 seurat_object = seurat_object_27522_1) %>%
  nrow()
cells_27522_1 <- get_tsne_umap(cell_types = cell_types_27522_1,
                                    seurat_object = seurat_object_27522_1) %>%
  pull(cell_type) %>% table() %>% sort()
print(str_c("Total number of cells: ", n_cells_27522_1))
print("Number of cells by cell type:")
print(cells_27522_1)
print("Percentage of cells by cell type:")
print(cells_27522_1/(n_cells_27522_1/100))

n_cells_chimeric_transcripts <- bulk_sc_plot_df_27522_1_SF_only %>%
  filter(data_type == "Single Cell Chimeric Transcript") %>%
  separate(identifier, into = c("cell", "tx"), sep = ":") %>%
  pull(cell) %>% unique() %>% length()
n_plasma_cells <- get_tsne_umap(cell_types = cell_types_27522_1,
                                seurat_object = seurat_object_27522_1) %>%
  filter(cell_type == "Plasma Cells") %>% nrow()
print(str_c("Percentage of plasma cells with chimeric transcript detected: ", n_cells_chimeric_transcripts, "/", n_plasma_cells, " = ", round(100*n_cells_chimeric_transcripts/n_plasma_cells, 2), "%"))

################################################################################
# Reviewer responses 27522 section
################################################################################
# not every read with CT is in seurat object (affects overall proportion 98/2477)
chr13_del_df <- get_chr_cnv(chr = "chr13", infercnv = infercnv_27522_1, tsne_umap = get_tsne_umap(cell_types = cell_types_27522_1, seurat_object = seurat_object_27522_1), genes = gene_spans) %>% filter(cell_type == "Plasma Cells")
chr16_del_df <- get_chr_cnv(chr = "chr16", infercnv = infercnv_27522_1, tsne_umap = get_tsne_umap(cell_types = cell_types_27522_1, seurat_object = seurat_object_27522_1), genes = gene_spans) %>% filter(cell_type == "Plasma Cells")
cells_with_chimeric_transcript <- get_sc_chimeric_transcripts_27522(dis_reads = dis_reads_27522_1_discover) %>% ungroup() %>% filter(cell_barcode %in% get_tsne_umap(cell_types = cell_types_27522_1, seurat_object = seurat_object_27522_1)$barcode) %>% pull(cell_barcode) %>% unique()
chr13_del_df %>% mutate(has_ct = barcode %in% cells_with_chimeric_transcript) %>% filter(mean_cnv <= 0.6) %>% pull(has_ct) %>% mean()


get_gene_cnv(infercnv = infercnv_27522_1, tsne_umap = get_tsne_umap(cell_types = cell_types_27522_1, seurat_object = seurat_object_27522_1), this_gene = "FGFR3")

################################################################################
# Reviewer responses do all the scRNA samples
################################################################################
#data ranges
#4:1800000-2000000
#8:127660000-128200000
#11:69640000-69660000
#14:105500000-107000000

if (TRUE) {
  data_dir <- "data/scRNA_revision/"
  output_dir <- "paper/supplementary/05_single_cell/all_samples/"

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  sample_ids <- c("27522_1",
                  "27522_4",
                  "47499_p",
                  "56203_1",
                  "56203_2",
                  "77570",
                  "81012_1",
                  "81012_2")

  translocation_events <- c("t(4;14)",
                            "t(4;14)",
                            "t(11;14)",
                            "t(8;14)",
                            "t(8;14)",
                            "t(11;14)",
                            "t(11;14)",
                            "t(11;14)")

  gene_names <- c("WHSC1",
                  "WHSC1",
                  "CCND1",
                  "MYC",
                  "MYC",
                  "CCND1",
                  "CCND1",
                  "CCND1")

  ensg_names <- c("ENSG00000109685",
                  "ENSG00000109685",
                  "ENSG00000110092",
                  "ENSG00000136997",
                  "ENSG00000136997",
                  "ENSG00000110092",
                  "ENSG00000110092",
                  "ENSG00000110092")

  chrA_vector <- c(4, 4, 11, 8, 8, 11, 11, 11)
  chrB_vector <- rep(14, length(chrA_vector))

  rds_file_paths <- str_c(data_dir,
                          "rds/",
                          c("backup_object_cell_type_in_sample_27522_1.rds",
                            "backup_object_cell_type_in_sample_27522_4.rds",
                            "backup_object_cell_type_in_sample_47499_p.rds",
                            "backup_object_cell_type_in_sample_56203_1.rds",
                            "backup_object_cell_type_in_sample_56203_2.rds",
                            "backup_object_cell_type_in_sample_77570.rds",
                            "backup_object_cell_type_in_sample_81012_1.rds",
                            "backup_object_cell_type_in_sample_81012_2.rds"))

  dis_reads_file_paths <- str_c(data_dir,
                                "dis_reads/",
                                c("27522_1/chr4chr14.discovered_discordant_reads.tsv",
                                  "27522_4/chr4chr14.discovered_discordant_reads.tsv",
                                  "47499/chr11chr14.discovered_discordant_reads.tsv",
                                  "56203_1/chr8chr14.discovered_discordant_reads.tsv",
                                  "56203_2/chr8chr14.discovered_discordant_reads.tsv",
                                  "77570/chr11chr14.discovered_discordant_reads.tsv",
                                  "81012_1/chr11chr14.discovered_discordant_reads.tsv",
                                  "81012_2/chr11chr14.discovered_discordant_reads.tsv"))


  get_sc_chimeric_transcripts <- function(dis_reads, chrA, chrB){

    sc_chrA_min_max <- dis_reads %>%
      filter(chrom == chrA) %>%
      filter(end - start < 128) %>%
      group_by(cell_barcode, molecular_barcode) %>%
      summarize(min_start = min(start), max_start = max(start),
                min_end = min(end), max_end = max(end),
                n_reads = n())

    sc_chrB_min_max <- dis_reads %>%
      filter(chrom == chrB) %>%
      filter(end - start < 128) %>%
      group_by(cell_barcode, molecular_barcode) %>%
      summarize(min_start = min(start), max_start = max(start),
                min_end = min(end), max_end = max(end),
                n_reads = n())

    sc_chrA_min_max %>%
      full_join(sc_chrB_min_max,
                by = c("cell_barcode", "molecular_barcode")) %>%
      filter(!any(is.na(min_start.x), is.na(min_start.y))) %>%
      return()
  }
  plot_cell_chimeric_transcripts <- function(dis_reads_input, chrA, chrB, tsne_umap, reduction = "UMAP", id, dir, sample_name, translocation, sens, spec){

    if (reduction %in% c("UMAP", "t-SNE")) {

      plot_df <- get_sc_chimeric_transcripts(dis_reads = dis_reads_input, chrA, chrB) %>%
        mutate(color_column = "Yes") %>%
        right_join(tsne_umap,
                   by = c("cell_barcode" = "barcode")) %>%
        replace_na(list(color_column = "No")) %>%
        mutate(color_column = fct_infreq(color_column)) %>%
        arrange(color_column)

      if (reduction == "UMAP") {
        p <- ggplot(data = plot_df, aes(x = UMAP_1, y = UMAP_2))
        p <- p + labs(x = "UMAP 1", y = "UMAP 2")
      } else {
        p <- ggplot(data = plot_df, aes(x = tSNE_1, y = tSNE_2))
        p <- p + labs(x = "t-SNE 1", y = "t-SNE 2")
      }
    } else {
      stop("Reduction parameter reduction must be 'UMAP' or 't-SNE'.")
    }

    p <- p + geom_point(aes(color = color_column),
                        shape = 16,
                        size = 1.5,
                        alpha = 0.25) +
      labs(color = "Chimeric\nTranscript\nDetected") +
      scale_color_manual(values = c("#d9d9d9", "#ff00ff"))

    p <- p + annotate("text", x = Inf, y = -Inf, label = sample_name, size = 2.5, show.legend = FALSE, hjust = 1, vjust = 0)
    p <- p + annotate("text", x = Inf, y = Inf, label = translocation, size = 2.5, show.legend = FALSE, hjust = 1, vjust = 1, color = "#ff00ff")

    p <- p +
      theme_bw() +
      scale_x_continuous(expand = c(0.01, 0.01)) +
      scale_y_continuous(expand = c(0.01, 0.01)) +
      theme(panel.background = element_blank(),
            panel.border = element_blank(),
            panel.grid = element_blank(),
            plot.background = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            axis.title = element_text(size = 8),
            legend.title = element_text(size = 8),
            legend.text = element_text(size = 6))

    ggsave(str_c(dir, "/", id, ".cells_chimeric_transcripts.", reduction, ".pdf"),
           p,
           width = 3.5, height = 1.75, useDingbats = FALSE)

    ggsave(str_c(dir, "/", id, ".cells_chimeric_transcripts.no_legend.", reduction, ".pdf"),
           p + guides(color = FALSE),
           width = 1.75, height = 1.75, useDingbats = FALSE)

  }
  plot_cell_types <- function(cell_types, seurat_object, reduction = "UMAP", id, dir){

    if (reduction %in% c("UMAP", "t-SNE")) {
      plot_df = get_tsne_umap(cell_types, seurat_object)
      if (reduction == "UMAP") {
        p <- ggplot(data = plot_df, aes(x = UMAP_1, y = UMAP_2))
        p <- p + labs(x = "UMAP 1", y = "UMAP 2")
      } else {
        p <- ggplot(data = plot_df, aes(x = tSNE_1, y = tSNE_2))
        p <- p + labs(x = "t-SNE 1", y = "t-SNE 2")
      }
    } else {
      stop("Reduction parameter reduction must be 'UMAP' or 't-SNE'.")
    }

    cell_type_positions <- plot_df %>% group_by(cell_type) %>%
      summarize(mean_tSNE_1 = mean(tSNE_1),
                mean_tSNE_2 = mean(tSNE_2),
                mean_UMAP_1 = mean(UMAP_1),
                mean_UMAP_2 = mean(UMAP_2))

    p <- p + geom_point(aes(color = cell_type),
                        shape = 16,
                        size = 1.5,
                        show.legend = FALSE,
                        alpha = 0.25)

    if (reduction == "UMAP") {
      p <- p + geom_text(data = cell_type_positions,
                         aes(x = mean_UMAP_1, y = mean_UMAP_2,
                             label = cell_type,
                             color = cell_type),
                         size = 2.5,
                         show.legend = FALSE)
    } else {
      p <- p + geom_text(data = cell_type_positions,
                         aes(x = mean_tSNE_1, y = mean_tSNE_2,
                             label = cell_type,
                             color = cell_type),
                         size = 2.5,
                         show.legend = FALSE)
    }

    p <- p +
      theme_bw() +
      scale_color_brewer(palette = "Set3", drop = FALSE, direction = -1) +
      scale_x_continuous(expand = c(0.01, 0.01)) +
      scale_y_continuous(expand = c(0.01, 0.01)) +
      theme(panel.background = element_blank(),
            panel.border = element_blank(),
            plot.background = element_blank(),
            panel.grid = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            axis.title = element_text(size = 8),
            legend.title = element_text(size = 8),
            legend.text = element_text(size = 6))

    ggsave(str_c(dir, "/", id, ".cell_types.", reduction, ".pdf"),
           p,
           width = 1.75, height = 1.75, useDingbats = FALSE)

  }

  cell_types_all <- read_tsv("data/scRNA.sample_barcode_celltype.txt")

  for (i in 1:length(sample_ids)) {

    sample <- sample_ids[i]
    print(sample)
    trans <- translocation_events[i]
    chrA_value = chrA_vector[i]
    chrB_value = chrB_vector[i]
    rds_path <- rds_file_paths[i]
    ds_path <- dis_reads_file_paths[i]

    seurat_object_sample <- RunUMAP(UpdateSeuratObject(read_rds(rds_path)), dims = 1:20)
    dis_reads_sample <- read_tsv(ds_path) %>% mutate(supports = trans)
    cell_types_sample <- cell_types_all %>% filter(sample_id == sample)

    plot_df <- get_sc_chimeric_transcripts(dis_reads = dis_reads_sample, chrA_value, chrB_value) %>%
      mutate(color_column = "Yes") %>%
      right_join(get_tsne_umap(cell_types = cell_types_sample, seurat_object = seurat_object_sample),
                 by = c("cell_barcode" = "barcode")) %>%
      replace_na(list(color_column = "No")) %>%
      mutate(color_column = fct_infreq(color_column)) %>%
      arrange(color_column) %>% ungroup()

    n_plasma_cells <- cell_types_sample %>%
      filter(cell_type == "Plasma") %>%
      nrow()
    n_plasma_cells_with_ct <- plot_df %>%
      filter(color_column == "Yes", cell_type == "Plasma Cells") %>%
      nrow()
    n_other_cells_with_ct <- plot_df %>%
      filter(color_column == "Yes", cell_type != "Plasma Cells") %>%
      nrow()

    plasma_cell_sensitivity <- str_c(round(100*n_plasma_cells_with_ct/n_plasma_cells, 2), "%")
    plasma_cell_specificity <- str_c(round(100*n_plasma_cells_with_ct/(n_plasma_cells_with_ct + n_other_cells_with_ct), 2), "%")

    print(str_c("Sensitivity: ", plasma_cell_sensitivity))
    print(str_c("Specificity: ", plasma_cell_specificity))

    plot_cell_chimeric_transcripts(dis_reads_input = dis_reads_sample,
                                   chrA = chrA_value,
                                   chrB = chrB_value,
                                   tsne_umap = get_tsne_umap(cell_types = cell_types_sample, seurat_object = seurat_object_sample),
                                   reduction = "UMAP",
                                   id = sample,
                                   dir = output_dir,
                                   sample_name = sample,
                                   translocation = trans,
                                   sens = plasma_cell_sensitivity,
                                   spec = plasma_cell_specificity)

    plot_cell_types(cell_types_sample, seurat_object_sample, id = sample, dir = output_dir)
    plot_gene_expression(seurat_object = seurat_object_sample,
                         tsne_umap = get_tsne_umap(cell_types = cell_types_sample,
                                                   seurat_object = seurat_object_sample),
                         ensg = ensg_names[i],
                         reduction = "UMAP",
                         id = sample,
                         gene = gene_names[i],
                         dir = output_dir)
  }

}

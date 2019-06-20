# ==============================================================================
# scRNA analysis
# Steven Foltz (smfoltz@wustl.edu), May 2019
# ==============================================================================

stop("This script is deprecated because it tries to run samples with bad QC results -- re-run 06_single_cell.R instead.")

paper_main = "paper/main/06_single_cell/"
paper_supp = "paper/supplemental/06_single_cell/"

# Create directories
for (id in c("27522_1", "27522_4",    #t(4;14)
             "47499",                 #t(11;14)
             "56203_1", "56203_2",    #t(8;*) 
             "77570",                 #t(11;14)
             "81012_1", "81012_2")) { #t(11;14)
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
    coord_fixed() +
    scale_color_brewer(palette = "Paired", drop = FALSE, direction = -1) +
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
  
  ggsave(str_c(dir, id, ".cell_types.", reduction, ".pdf"),
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
  coord_fixed() +
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
  
  ggsave(str_c(dir, id, ".expression.", gene, ".", reduction, ".pdf"),
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
    geom_jitter(#aes(color = expression),
                #width = 0.5,
                height = 0,
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
  
  ggsave(str_c(dir, id, ".expression_violin.", gene, ".pdf"),
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
    coord_fixed() +
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
  
  ggsave(str_c(dir, id, ".joint_expression.", gene1, ".", gene2, ".pdf"),
         p, 
         width = 3.5, height = 3.5, useDingbats = FALSE)
}

# ==============================================================================
# Plot single gene CNV
# ==============================================================================
plot_gene_cnv <- function(infercnv, tsne_umap, reduction = "UMAP", id, gene, dir = paper_supp) {
  
  if (reduction %in% c("UMAP", "t-SNE")) {
    #plot_df <- get_gene_cnv(infercnv, tsne_umap, this_gene = gene) %>% arrange(desc(cnv))
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
    coord_fixed() +
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
  
  ggsave(str_c(dir, id, ".copy_number.", gene, ".pdf"),
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
                      #show.legend = FALSE,
                      alpha = 0.25) +
    annotate("text", label = chr,
             x = -Inf, y = Inf, 
             vjust = 1, hjust = 0,
             size = 3.5) + 
    theme_bw() +
    coord_fixed() +
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
  
  ggsave(str_c(dir, id, ".copy_number.", chr, ".pdf"),
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
  
  ggsave(str_c(dir, id, ".read_length_distribution.pdf"),
         width = 3.5, height = 3.5, useDingbats = FALSE)
}

# ==============================================================================
# Identify cells with chimeric transcripts
# ==============================================================================
plot_cell_chimeric_transcripts <- function(bulk_sc, tsne_umap, reduction = "UMAP", id, dir = paper_main, facet = FALSE){
  
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
  
  p <- p + geom_point(aes(color = color_column),
                      shape = 16, 
                      size = 1.5, 
                      alpha = 1)
  
  p <- p + 
    labs(color = "Cell Contains\nChimeric Transcript") +
    theme_bw() +
    coord_equal() +
    scale_color_brewer(palette = "Set1") +
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
    ggsave(str_c(dir, id, ".cells_chimeric_transcripts.", reduction, ".pdf"),
           p,
           width = n_colors*3.5, height = 3.5, useDingbats = FALSE)
    
    ggsave(str_c(dir, id, ".cells_chimeric_transcripts.no_legend.", reduction, ".pdf"),
           p + guides(color = FALSE),
           width = n_colors*3.5, height = 3.5, useDingbats = FALSE)
  } else {
    ggsave(str_c(dir, id, ".cells_chimeric_transcripts.", reduction, ".pdf"),
           p,
           width = 3.5, height = 3.5, useDingbats = FALSE)
    
    ggsave(str_c(dir, id, ".cells_chimeric_transcripts.no_legend.", reduction, ".pdf"),
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
  
  print(str_c(id, " correlation between ", gene1, " and ", gene2, ":"))
  
  gene12_correlation <- plot_df %>%
    filter(expression1 > 0, expression2 > 0) %>%
    select(expression1, expression2) %>%
    as.matrix() %>%
    cor()
  
  print(gene12_correlation)
  
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
  
  ggsave(str_c(dir, id, ".", gene1, ".", gene2, ".correlation.pdf"),
         p,
         width = 3.5, height = 3.5, useDingbats = FALSE)
}

#General analysis (all 5 patients)
if (TRUE) {
  
  # ==============================================================================
  # ANALYSIS
  # ==============================================================================
  
  # ==============================================================================
  # Plot cell types of each sample
  # UMAP in paper_main; t-SNE in paper_supp
  # ==============================================================================
  if (TRUE) {
    plot_cell_types(cell_types_27522_1, seurat_object_27522_1, id = "27522_1", dir = str_c(paper_main, "27522_1/"))
    plot_cell_types(cell_types_27522_4, seurat_object_27522_4, id = "27522_4", dir = str_c(paper_main, "27522_4/"))
    plot_cell_types(cell_types_47499, seurat_object_47499, id = "47499", dir = str_c(paper_main, "47499/"))
    plot_cell_types(cell_types_56203_1, seurat_object_56203_1, id = "56203_1", dir = str_c(paper_main, "56203_1/"))
    plot_cell_types(cell_types_56203_2, seurat_object_56203_2, id = "56203_2", dir = str_c(paper_main, "56203_2/"))
    plot_cell_types(cell_types_77570, seurat_object_77570, id = "77570", dir = str_c(paper_main, "77570/"))
    plot_cell_types(cell_types_81012_1, seurat_object_81012_1, id = "81012_1", dir = str_c(paper_main, "81012_1/"))
    plot_cell_types(cell_types_81012_2, seurat_object_81012_2, id = "81012_2", dir = str_c(paper_main, "81012_2/"))
    
    plot_cell_types(cell_types_27522_1, seurat_object_27522_1, id = "27522_1", reduction = "t-SNE", dir = str_c(paper_supp, "27522_1/"))
    plot_cell_types(cell_types_27522_4, seurat_object_27522_4, id = "27522_4", reduction = "t-SNE", dir = str_c(paper_supp, "27522_4/"))
    plot_cell_types(cell_types_47499, seurat_object_47499, id = "47499", reduction = "t-SNE", dir = str_c(paper_supp, "47499/"))
    plot_cell_types(cell_types_56203_1, seurat_object_56203_1, id = "56203_1", reduction = "t-SNE", dir = str_c(paper_supp, "56203_1/"))
    plot_cell_types(cell_types_56203_2, seurat_object_56203_2, id = "56203_2", reduction = "t-SNE", dir = str_c(paper_supp, "56203_2/"))
    plot_cell_types(cell_types_77570, seurat_object_77570, id = "77570", reduction = "t-SNE", dir = str_c(paper_supp, "77570/"))
    plot_cell_types(cell_types_81012_1, seurat_object_81012_1, id = "81012_1", reduction = "t-SNE", dir = str_c(paper_supp, "81012_1/"))
    plot_cell_types(cell_types_81012_2, seurat_object_81012_2, id = "81012_2", reduction = "t-SNE", dir = str_c(paper_supp, "81012_2/"))
  }
  
  # ==============================================================================
  # Plot interesting gene expressions
  # 27522: t(4;14) WHSC1 (ENSG00000109685) and FGFR3 (ENSG00000068078)
  # 47499: t(11;14) CCND1 (ENSG00000110092)
  # 56203: t(8;*) MYC (ENSG00000136997) and PVT1 (ENSG00000249859)
  # 77570: t(11;14) CCND1 (ENSG00000110092)
  # 81012: t(11;14) CCND1 (ENSG00000110092)
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
                             dir = str_c(paper_supp, "27522_1/"), reduction = red)  
        plot_gene_expression(seurat_object_27522_4, 
                             get_tsne_umap(cell_types = cell_types_27522_4, 
                                           seurat_object = seurat_object_27522_4), 
                             ensg = t414_ensg[i], id = "27522_4", gene = t414_gene_names[i],
                             dir = str_c(paper_supp, "27522_4/"), reduction = red)
      }
      
      t1114_gene_names <- c("CCND1", "SDC1", "TNFRSF17", "SLAMF7")
      t1114_ensg <- c("ENSG00000110092", "ENSG00000115884", "ENSG00000048462", "ENSG00000026751")
      for (i in 1:length(t1114_gene_names)) {
        plot_gene_expression(seurat_object_47499, 
                             get_tsne_umap(cell_types = cell_types_47499, 
                                           seurat_object = seurat_object_47499), 
                             ensg = t1114_ensg[i], id = "47499", gene = t1114_gene_names[i],
                             dir = str_c(paper_supp, "47499/"), reduction = red)
        plot_gene_expression(seurat_object_77570, 
                             get_tsne_umap(cell_types = cell_types_77570, 
                                           seurat_object = seurat_object_77570), 
                             ensg = "ENSG00000110092", id = "77570", gene = "CCND1",
                             dir = str_c(paper_supp, "77570/"), reduction = red)
        plot_gene_expression(seurat_object_81012_1, 
                             get_tsne_umap(cell_types = cell_types_81012_1, 
                                           seurat_object = seurat_object_81012_1), 
                             ensg = "ENSG00000110092", id = "81012_1", gene = "CCND1",
                             dir = str_c(paper_supp, "81012_1/"), reduction = red)
        plot_gene_expression(seurat_object_81012_2, 
                             get_tsne_umap(cell_types = cell_types_81012_2, 
                                           seurat_object = seurat_object_81012_2), 
                             ensg = "ENSG00000110092", id = "81012_2", gene = "CCND1",
                             dir = str_c(paper_supp, "81012_2/"), reduction = red)
      }
      
      t814_gene_names <- c("MYC", "SDC1", "TNFRSF17", "SLAMF7")
      t814_ensg <- c("ENSG00000136997", "ENSG00000115884", "ENSG00000048462", "ENSG00000026751")
      for (i in 1:length(t814_gene_names)) {
        plot_gene_expression(seurat_object_56203_1, 
                             get_tsne_umap(cell_types = cell_types_56203_1, 
                                           seurat_object = seurat_object_56203_1), 
                             ensg = t814_ensg[i], id = "56203_1", gene = t814_gene_names[i],
                             dir = str_c(paper_supp, "56203_1/"), reduction = red)
        plot_gene_expression(seurat_object_56203_2, 
                             get_tsne_umap(cell_types = cell_types_56203_2, 
                                           seurat_object = seurat_object_56203_2), 
                             ensg = t814_ensg[i], id = "56203_2", gene = t814_gene_names[i],
                             dir = str_c(paper_supp, "56203_2/"), reduction = red)
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
      filter(end - start < 1024) %>%
      group_by(cell_barcode, molecular_barcode) %>% 
      summarize(min_start = min(start), max_start = max(start), 
                min_end = min(end), max_end = max(end), 
                n_reads = n())
    
    sc_chr4_min_max <- dis_reads %>% 
      filter(chrom == 4) %>% 
      filter(end - start < 1024) %>% 
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
    
    ggsave(str_c(dir, id, ".discordant_read_positions.pdf"),
           p,
           width = 14.5, height = 12,
           useDingbats = FALSE)
    
    ggsave(str_c(dir, id, ".discordant_read_positions.no_legend.pdf"),
           p + guides(color = FALSE, linetype = FALSE, fill = FALSE),
           width = 7.25, height = 12,
           useDingbats = FALSE)
  }
  
  # 27522_1
  if (TRUE) {
    
    bulk_fusion_reads_27522_1_SF_only <- get_bulk_fusion_reads_27522(bulk_reads_27522_1, star_fusion_calls_27522_1, use_SF_only = TRUE)
    bulk_fusion_reads_27522_1 <- get_bulk_fusion_reads_27522(bulk_reads_27522_1, star_fusion_calls_27522_1, use_SF_only = FALSE)
    plot_sc_read_length_distribution(dis_reads = dis_reads_27522_1_discover, cutoff = 1024, id = "27522_1", dir = str_c(paper_supp, "27522_1/"))
    sc_chimeric_transcripts_27522_1 <- get_sc_chimeric_transcripts_27522(dis_reads_27522_1_discover)
    bulk_sc_plot_df_27522_1_SF_only <- get_bulk_sc_plot_df_27522(bulk_fusion_reads_27522_1_SF_only, sc_chimeric_transcripts_27522_1)
    bulk_sc_plot_df_27522_1 <- get_bulk_sc_plot_df_27522(bulk_fusion_reads_27522_1, sc_chimeric_transcripts_27522_1)
    
    plot_bulk_sc_27522(bulk_sc_plot_df_27522_1_SF_only, id = "27522_1", bulk_color = FALSE, dir = str_c(paper_main, "27522_1/"))
    #plot_bulk_sc_27522(bulk_sc_plot_df_27522_1, id = "27522_1", bulk_color = FALSE, dir = str_c(paper_supp, "27522_1/"))
    
    plot_cell_chimeric_transcripts(bulk_sc = bulk_sc_plot_df_27522_1_SF_only, 
                                   tsne_umap = get_tsne_umap(cell_types = cell_types_27522_1, seurat_object = seurat_object_27522_1), 
                                   id = "27522_1", 
                                   reduction = "UMAP", 
                                   dir = str_c(paper_main, "27522_1/"))
    plot_cell_chimeric_transcripts(bulk_sc = bulk_sc_plot_df_27522_1_SF_only, 
                                   tsne_umap = get_tsne_umap(cell_types = cell_types_27522_1, seurat_object = seurat_object_27522_1), 
                                   id = "27522_1", 
                                   reduction = "t-SNE", 
                                   dir = str_c(paper_supp, "27522_1/"))
    
    plot_gene_expression_violin(bulk_sc = bulk_sc_plot_df_27522_1, 
                                seurat_object = seurat_object_27522_1, 
                                tsne_umap = get_tsne_umap(cell_types = cell_types_27522_1, 
                                                          seurat_object = seurat_object_27522_1), 
                                ensg = "ENSG00000068078", 
                                id = "27522_1", 
                                gene = "FGFR3", 
                                dir = str_c(paper_supp, "27522_1/"),
                                max = 4.4)
    plot_gene_expression_violin(bulk_sc = bulk_sc_plot_df_27522_1, 
                                seurat_object = seurat_object_27522_1, 
                                tsne_umap = get_tsne_umap(cell_types = cell_types_27522_1, seurat_object = seurat_object_27522_1), 
                                ensg = "ENSG00000109685", 
                                id = "27522_1", 
                                gene = "WHSC1", 
                                dir = str_c(paper_supp, "27522_1/"),
                                max = 4.4)
    # https://www.nature.com/articles/nrc2780/figures/2
    plot_gene_expression_violin(bulk_sc = bulk_sc_plot_df_27522_1, 
                                seurat_object = seurat_object_27522_1, 
                                tsne_umap = get_tsne_umap(cell_types = cell_types_27522_1, seurat_object = seurat_object_27522_1), 
                                ensg = "ENSG00000142208", 
                                id = "27522_1", 
                                gene = "AKT1", 
                                dir = str_c(paper_supp, "27522_1/"),
                                max = 4.4)
    plot_gene_expression_violin(bulk_sc = bulk_sc_plot_df_27522_1, 
                                seurat_object = seurat_object_27522_1, 
                                tsne_umap = get_tsne_umap(cell_types = cell_types_27522_1, seurat_object = seurat_object_27522_1), 
                                ensg = "ENSG00000100030", 
                                id = "27522_1",
                                gene = "MAPK1", 
                                dir = str_c(paper_supp, "27522_1/"),
                                max = 4.4)
    
    plot_two_genes_correlation(two_genes_expression = get_two_genes_expression(seurat_object = seurat_object_27522_1,
                                                                               tsne_umap = get_tsne_umap(cell_types = cell_types_27522_1, seurat_object = seurat_object_27522_1),
                                                                               ensg1 = "ENSG00000109685",
                                                                               ensg2 = "ENSG00000068078"),
                               bulk_sc = bulk_sc_plot_df_27522_1, 
                               dir = str_c(paper_supp, "27522_1/"),
                               id = "27522_1", 
                               gene1 = "WHSC1", 
                               gene2 = "FGFR3")
  }
  
  # 27522_4
  if (TRUE) {
    bulk_fusion_reads_27522_4 <- get_bulk_fusion_reads_27522(bulk_reads_27522_4)
    plot_sc_read_length_distribution(dis_reads = dis_reads_27522_4_discover, cutoff = 1024, id = "27522_4", dir = str_c(paper_supp, "27522_4/"))
    sc_chimeric_transcripts_27522_4 <- get_sc_chimeric_transcripts_27522(dis_reads_27522_4_discover)
    bulk_sc_plot_df_27522_4 <- get_bulk_sc_plot_df_27522(bulk_fusion_reads_27522_4, sc_chimeric_transcripts_27522_4)
    plot_bulk_sc_27522(bulk_sc_plot_df_27522_4, id = "27522_4", bulk_color = FALSE, dir = str_c(paper_main, "27522_4/"))
    
    plot_cell_chimeric_transcripts(bulk_sc = bulk_sc_plot_df_27522_4, 
                                   tsne_umap = get_tsne_umap(cell_types = cell_types_27522_4, seurat_object = seurat_object_27522_4), 
                                   id = "27522_4",
                                   reduction = "UMAP",
                                   dir = str_c(paper_main, "27522_4/"))
    plot_cell_chimeric_transcripts(bulk_sc = bulk_sc_plot_df_27522_4, 
                                   tsne_umap = get_tsne_umap(cell_types = cell_types_27522_4, seurat_object = seurat_object_27522_4), 
                                   id = "27522_4", 
                                   reduction = "t-SNE", 
                                   dir = str_c(paper_supp, "27522_4/"))
    
    plot_gene_expression_violin(bulk_sc = bulk_sc_plot_df_27522_4, 
                                seurat_object = seurat_object_27522_4, 
                                tsne_umap = get_tsne_umap(cell_types = cell_types_27522_4, 
                                                          seurat_object = seurat_object_27522_4), 
                                ensg = "ENSG00000068078", 
                                id = "27522_4", 
                                gene = "FGFR3", 
                                dir = str_c(paper_supp, "27522_4/"),
                                max = 4.4)
    plot_gene_expression_violin(bulk_sc = bulk_sc_plot_df_27522_4, 
                                seurat_object = seurat_object_27522_4, 
                                tsne_umap = get_tsne_umap(cell_types = cell_types_27522_4, seurat_object = seurat_object_27522_4), 
                                ensg = "ENSG00000109685", 
                                id = "27522_4", 
                                gene = "WHSC1", 
                                dir = str_c(paper_supp, "27522_4/"),
                                max = 4.4)
    # https://www.nature.com/articles/nrc2780/figures/2
    plot_gene_expression_violin(bulk_sc = bulk_sc_plot_df_27522_4, 
                                seurat_object = seurat_object_27522_4, 
                                tsne_umap = get_tsne_umap(cell_types = cell_types_27522_4, seurat_object = seurat_object_27522_4), 
                                ensg = "ENSG00000142208", 
                                id = "27522_4", 
                                gene = "AKT1", 
                                dir = str_c(paper_supp, "27522_4/"),
                                max = 4.4)
    plot_gene_expression_violin(bulk_sc = bulk_sc_plot_df_27522_4, 
                                seurat_object = seurat_object_27522_4, 
                                tsne_umap = get_tsne_umap(cell_types = cell_types_27522_4, seurat_object = seurat_object_27522_4), 
                                ensg = "ENSG00000100030", 
                                id = "27522_4",
                                gene = "MAPK1", 
                                dir = str_c(paper_supp, "27522_4/"),
                                max = 4.4)
    
    plot_two_genes_correlation(two_genes_expression = get_two_genes_expression(seurat_object = seurat_object_27522_4,
                                                                               tsne_umap = get_tsne_umap(cell_types = cell_types_27522_4, seurat_object = seurat_object_27522_4),
                                                                               ensg1 = "ENSG00000109685",
                                                                               ensg2 = "ENSG00000068078"),
                               bulk_sc = bulk_sc_plot_df_27522_4, 
                               dir = str_c(paper_supp, "27522_4/"), 
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
                 dir = str_c(paper_supp, "27522_4/"))
    
    # ============================================================================
    # Plot chr13 CNV in 25722_4
    # ============================================================================
    plot_chr_cnv(infercnv_27522_4, 
                 get_tsne_umap(cell_types = cell_types_27522_4, 
                               seurat_object = seurat_object_27522_4), 
                 id = "27522_4", chr = "chr13",
                 dir = str_c(paper_supp, "27522_4/"))
    
    # ============================================================================
    # Plot chr16 CNV in 25722_4
    # ============================================================================
    plot_chr_cnv(infercnv_27522_4, 
                 get_tsne_umap(cell_types = cell_types_27522_4, 
                               seurat_object = seurat_object_27522_4), 
                 id = "27522_4", chr = "chr16",
                 dir = str_c(paper_supp, "27522_4/"))
    
    
  }
}

# ==============================================================================
# Work with 47499
# ==============================================================================
if (TRUE) {
  get_bulk_fusion_reads_47499 <- function(bulk_reads, star_fusion = NULL, use_SF_only = FALSE){
    
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
      filter(chromosome_donor %in% c("chr11", "chr14"), 
             chromosome_acceptor %in% c("chr11", "chr14")) %>%
      select(read_name, 
             chromosome_donor, first_base_donor, 
             chromosome_acceptor, first_base_acceptor) %>% 
      unique() %>% 
      mutate(chr14_position = case_when(chromosome_donor == "chr14" ~ first_base_donor,
                                        TRUE ~ first_base_acceptor),
             chr11_position = case_when(chromosome_donor == "chr11" ~ first_base_donor,
                                        TRUE ~ first_base_acceptor)) %>%
      filter(chr14_position > 105500000,
             chr14_position < 107000000) %>%
      filter(chr11_position > 69640000,
             chr11_position < 69660000) %>%
      select(read_name, chr14_position, chr11_position) %>%
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
  get_sc_chimeric_transcripts_47499 <- function(dis_reads){
    
    sc_chr14_min_max <- dis_reads %>% 
      filter(chrom == 14) %>% 
      filter(end - start < 1024) %>%
      group_by(cell_barcode, molecular_barcode) %>% 
      summarize(min_start = min(start), max_start = max(start), 
                min_end = min(end), max_end = max(end), 
                n_reads = n())
    
    sc_chr11_min_max <- dis_reads %>% 
      filter(chrom == 11) %>% 
      filter(end - start < 1024) %>% 
      group_by(cell_barcode, molecular_barcode) %>% 
      summarize(min_start = min(start), max_start = max(start), 
                min_end = min(end), max_end = max(end), 
                n_reads = n())
    
    sc_chr14_min_max %>% 
      full_join(sc_chr11_min_max, 
                by = c("cell_barcode", "molecular_barcode")) %>%
      filter(!any(is.na(min_start.x), is.na(min_start.y))) %>%
      return()
  }
  get_bulk_sc_plot_df_47499 <- function(bulk_fusion_reads, sc_chimeric_transcripts){
    
    chr11_breakpoint <- 69560348
    chr14_breakpoint <- 105861673
    between_genes <- 0.2e6
    
    if (bulk_fusion_reads %>% nrow() == 0) {
      bulk_fusion_reads <- NULL
    }
    
    if (!identical(bulk_fusion_reads, NULL)) {
      bulk_chr14_min <- bulk_fusion_reads %>% pull(chr14_position) %>% min()
      bulk_chr14_max <- bulk_fusion_reads %>% pull(chr14_position) %>% max()
      bulk_chr11_min <- bulk_fusion_reads %>% pull(chr11_position) %>% min()
      bulk_chr11_max <- bulk_fusion_reads %>% pull(chr11_position) %>% max()
    } else {
      bulk_chr11_min <- NULL
      bulk_chr11_max <- NULL
      bulk_chr14_min <- NULL
      bulk_chr14_max <- NULL
    }
    
    sc_chr14_min <- sc_chimeric_transcripts %>% pull(min_start.x) %>% min()
    sc_chr14_max <- sc_chimeric_transcripts %>% pull(min_start.x) %>% max()
    sc_chr11_min <- sc_chimeric_transcripts %>% pull(min_start.y) %>% min()
    sc_chr11_max <- sc_chimeric_transcripts %>% pull(min_start.y) %>% max()
    
    overall_chr14_min <- min(bulk_chr14_min, sc_chr14_min)
    overall_chr14_max <- max(bulk_chr14_max, sc_chr14_max)
    overall_chr11_min <- min(bulk_chr11_min, sc_chr11_min)
    overall_chr11_max <- max(bulk_chr11_max, sc_chr11_max)
    
    chr11_shift = 0 - overall_chr11_min + overall_chr14_max - overall_chr14_min + between_genes
    chr14_shift = 0 - overall_chr14_min
    
    shifted_chr11_breakpoint = chr11_breakpoint + chr11_shift
    shifted_chr14_breakpoint = chr14_breakpoint + chr14_shift
    
    if (identical(bulk_fusion_reads, NULL)) {
      bulk_plot_df <- NULL
    } else {
      bulk_plot_df <- bulk_fusion_reads %>%
        mutate(category = case_when(chr11_position <= chr11_breakpoint & chr14_position <= chr14_breakpoint ~ "Supports Reciprocal CCND1--IGH Fusion",
                                    chr11_position > chr11_breakpoint & chr14_position <= chr14_breakpoint ~ "Supports General IGH Fusion with CCND1",
                                    chr11_position <= chr11_breakpoint & chr14_position > chr14_breakpoint ~ "Supports General IGH Fusion with CCND1",
                                    chr11_position > chr11_breakpoint & chr14_position > chr14_breakpoint ~ "Supports IGH--CCND1 Fusion")) %>%
        mutate(category = factor(category, levels = c("Supports IGH--CCND1 Fusion", "Supports Reciprocal CCND1--IGH Fusion", "Supports General IGH Fusion with CCND1"), ordered = TRUE)) %>%
        mutate(shifted_chr11 = chr11_position + chr11_shift,
               shifted_chr14 = chr14_position + chr14_shift) %>%
        mutate(identifier = read_name,
               data_type = "Bulk Spanning/Junction Read Pair",
               y_position_offset = 0,
               curvature = -0.25) %>%
        select(identifier, data_type, chr14_position, chr11_position, shifted_chr14, shifted_chr11, category, fusion)  
    }
    
    sc_plot_df <- sc_chimeric_transcripts %>%
      mutate(category = case_when((min_start.x <= chr14_breakpoint & max_end.x > chr14_breakpoint) | (min_start.y < chr11_breakpoint & max_end.y > chr11_breakpoint) ~ "Supports General IGH Fusion with CCND1",
                                  min_start.x > chr14_breakpoint & min_start.y > chr11_breakpoint ~ "Supports IGH--CCND1 Fusion",
                                  max_end.x <= chr14_breakpoint & max_end.y <= chr11_breakpoint ~ "Supports Reciprocal CCND1--IGH Fusion",
                                  max_end.x <= chr14_breakpoint & min_start.y > chr11_breakpoint ~ "Supports General IGH Fusion with CCND1",
                                  TRUE ~ "Supports General IGH Fusion with CCND1")) %>%
      mutate(chr11_position = case_when(category == "Supports IGH--CCND1 Fusion" ~ min_start.y,
                                        category == "Supports Reciprocal CCND1--IGH Fusion" ~ max_end.y,
                                        category == "Supports General IGH Fusion with CCND1" ~ min_start.y),
             chr14_position = case_when(category == "Supports IGH--CCND1 Fusion" ~ min_start.x,
                                        category == "Supports Reciprocal CCND1--IGH Fusion" ~ max_end.x,
                                        category == "Supports General IGH Fusion with CCND1" ~ min_start.x)) %>%
      mutate(category = factor(category, levels = c("Supports IGH--CCND1 Fusion", "Supports Reciprocal CCND1--IGH Fusion", "Supports General IGH Fusion with CCND1"), ordered = TRUE)) %>%
      mutate(shifted_chr11 = chr11_position + chr11_shift,
             shifted_chr14 = chr14_position + chr14_shift) %>%
      mutate(identifier = str_c(cell_barcode, ":", molecular_barcode),
             data_type = "Single Cell Chimeric Transcript",
             y_position_offset = -0.5,
             curvature = 0.25,
             fusion = NA) %>%
      ungroup() %>%
      select(identifier, data_type, chr14_position, chr11_position, shifted_chr14, shifted_chr11, category, fusion)
    
    bind_rows(bulk_plot_df, sc_plot_df) %>% 
      mutate(color_column = category) %>%
      return()
  }
  plot_bulk_sc_47499 <- function(bulk_sc_plot_df, genes = gene_spans, dir, id, use_fusion_names = TRUE, bulk_color = TRUE, color_indicates = "Fusion Name"){
    reported_translocation_breakpoints <- tribble(~trans, ~chrom,      ~pos,
                                                  "t(11;14)",     11,  69560348,
                                                  "t(11;14)",     14, 105861673)
    chr11_breakpoint <- 69560348
    chr14_breakpoint <- 105861673
    between_genes <- 0.2e6
    chr11_offset_y <- 0.5
    
    bulk_sc_plot_df <- bulk_sc_plot_df %>% 
      mutate(color_column = case_when(!bulk_color ~ "black",
                                      use_fusion_names ~ str_remove_all(str_remove_all(fusion, "CCND1"), "--"),
                                      TRUE ~ as.character(str_detect(fusion, "IGH"))),
             reciprocal = str_detect(fusion, "^CCND1"))
    
    overall_chr14_min <- bulk_sc_plot_df %>% pull(chr14_position) %>% min()
    overall_chr14_max <- bulk_sc_plot_df %>% pull(chr14_position) %>% max() 
    overall_chr11_min <- bulk_sc_plot_df %>% pull(chr11_position) %>% min()
    overall_chr11_max <- bulk_sc_plot_df %>% pull(chr11_position) %>% max()
    
    chr11_shift = 0 - overall_chr11_min + overall_chr14_max - overall_chr14_min + between_genes
    chr14_shift = 0 - overall_chr14_min
    
    shifted_chr11_breakpoint = chr11_breakpoint + chr11_shift
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
    
    chr11_gene_spans <- genes %>% 
      filter(chromosome == "chr11") %>% 
      filter( gene_name %in% c("CCND1") |
                (start <= overall_chr11_min & end >= overall_chr11_min) | 
                (start >= overall_chr11_min & end <= overall_chr11_max) | 
                (start <= overall_chr11_max & end >= overall_chr11_max) ) %>% 
      filter(strand == "+", type == "protein_coding") %>%
      mutate(shifted_chr11_start = start + chr11_shift,
             shifted_chr11_end = end + chr11_shift)
    
    shifted_gene_boundaries <- c(chr14_gene_spans %>% pull(shifted_chr14_start) %>% min(),
                                 chr14_gene_spans %>% pull(shifted_chr14_end) %>% max(),
                                 chr11_gene_spans %>% pull(shifted_chr11_start) %>% min(),
                                 chr11_gene_spans %>% pull(shifted_chr11_start) %>% max(),
                                 chr11_gene_spans %>% pull(shifted_chr11_end) %>% max())
    
    gene_boundaries <- c(chr14_gene_spans %>% pull(start) %>% min(),
                         chr14_gene_spans %>% pull(end) %>% max(),
                         chr11_gene_spans %>% pull(start) %>% min(),
                         chr11_gene_spans %>% pull(start) %>% max(),
                         chr11_gene_spans %>% pull(end) %>% max())
    
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
      #chr11 genes
      geom_rect(data = chr11_gene_spans,
                aes(xmin = shifted_chr11_start, 
                    xmax = shifted_chr11_end,
                    ymin = 0 - chr11_offset_y, 
                    ymax = 1 - chr11_offset_y,
                    fill = "chr11"), 
                #fill = "#feb24c",
                alpha = 1)
    
    if (bulk_sc_plot_df %>% 
        filter(data_type == "Bulk Spanning/Junction Read Pair") %>% 
        nrow() > 0) {
      # draw bulk read curves
      p <- p + geom_curve(data = bulk_sc_plot_df %>% 
                            filter(data_type == "Bulk Spanning/Junction Read Pair"),
                          aes(x = shifted_chrN, 
                              xend = shifted_chr8,
                              color = color_column,
                              linetype = reciprocal,
                              y = 1,
                              yend = 1 - chr8_offset_y),
                          curvature = -0.25,
                          ncp = 10,
                          lineend = "round",
                          alpha = 0.25)
    } 
    # draw sc read curves
    p <- p + geom_curve(data = bulk_sc_plot_df %>% 
                          filter(data_type == "Single Cell Chimeric Transcript"),
                        aes(x = shifted_chr14, 
                            xend = shifted_chr11,
                            y = 0,
                            yend = 0 - chr11_offset_y,
                            color = "black"),
                        curvature = 0.25,
                        ncp = 10,
                        lineend = "round",
                        show.legend = FALSE,
                        alpha = 0.25) +
      # label chr11 genes
      geom_text(data = chr11_gene_spans,
                aes(x = (shifted_chr11_end + shifted_chr11_start)/2,
                    y = 1 - chr11_offset_y - 0.5,
                    label = gene_name),
                fontface = "italic") +
      # translocation breakpoints 
      geom_vline(xintercept = shifted_chr11_breakpoint,
                 linetype = 2) +
      geom_vline(xintercept = shifted_chr14_breakpoint,
                 linetype = 2) +
      annotate(geom = "text", x = shifted_chr11_breakpoint + between_genes*0.02, y = -1.5, label = chr11_breakpoint, hjust = 0, vjust = 1) +
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
    
    ggsave(str_c(dir, id, ".discordant_read_positions.pdf"),
           p,
           width = 14.5, height = 12,
           useDingbats = FALSE)
    
    ggsave(str_c(dir, id, ".discordant_read_positions.no_legend.pdf"),
           p + guides(color = FALSE, linetype = FALSE, fill = FALSE),
           width = 7.25, height = 12,
           useDingbats = FALSE)
  }
  
  # 47499
  if (TRUE) {
    
    bulk_fusion_reads_47499 <- get_bulk_fusion_reads_47499(bulk_reads = bulk_reads_47499, star_fusion = NULL, use_SF_only = FALSE)
    plot_sc_read_length_distribution(dis_reads = dis_reads_47499_discover, cutoff = 1024, id = "47499", dir = str_c(paper_supp, "47499/"))
    sc_chimeric_transcripts_47499 <- get_sc_chimeric_transcripts_47499(dis_reads_47499_discover)
    bulk_sc_plot_df_47499 <- get_bulk_sc_plot_df_47499(bulk_fusion_reads_47499, sc_chimeric_transcripts_47499)
    
    plot_bulk_sc_47499(bulk_sc_plot_df_47499, id = "47499", bulk_color = FALSE, dir = str_c(paper_supp, "47499/"))
    
    plot_cell_chimeric_transcripts(bulk_sc = bulk_sc_plot_df_47499, 
                                   tsne_umap = get_tsne_umap(cell_types = cell_types_47499, seurat_object = seurat_object_47499), 
                                   id = "47499", 
                                   reduction = "UMAP", 
                                   dir = str_c(paper_main, "47499/"))
    plot_cell_chimeric_transcripts(bulk_sc = bulk_sc_plot_df_47499, 
                                   tsne_umap = get_tsne_umap(cell_types = cell_types_47499, seurat_object = seurat_object_47499), 
                                   id = "47499", 
                                   reduction = "t-SNE", 
                                   dir = str_c(paper_supp, "47499/"))
    plot_gene_expression_violin(bulk_sc = bulk_sc_plot_df_47499, 
                                seurat_object = seurat_object_47499, 
                                tsne_umap = get_tsne_umap(cell_types = cell_types_47499, 
                                                          seurat_object = seurat_object_47499), 
                                ensg = "ENSG00000110092", 
                                id = "47499", 
                                gene = "CCND1", 
                                dir = str_c(paper_supp, "47499/"),
                                max = 4.75)
  }
}
# ==============================================================================
# Work with 56203
# ==============================================================================
if (TRUE) {
  get_bulk_fusion_reads_56203 <- function(bulk_reads, star_fusion = NULL, use_SF_only = FALSE){
    
    if (identical(bulk_reads, NULL)) {
      return(NULL)
    }
    
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
      filter((chromosome_donor == "chr8" & chromosome_acceptor %in% c("chr2", "chr14", "chr22")) |
               (chromosome_acceptor == "chr8" & chromosome_donor %in% c("chr2", "chr14", "chr22"))) %>%
      select(read_name, 
             chromosome_donor, first_base_donor, 
             chromosome_acceptor, first_base_acceptor) %>% 
      unique() %>% 
      mutate(chr8_position = case_when(chromosome_donor == "chr8" ~ first_base_donor,
                                       TRUE ~ first_base_acceptor),
             chrN_position = case_when(chromosome_donor != "chr8" ~ first_base_donor,
                                       TRUE ~ first_base_acceptor),
             chrN = case_when(chromosome_donor == "chr8" ~ chromosome_acceptor,
                              TRUE ~ chromosome_donor)) %>% 
      filter(chr8_position > 127660000 & chr8_position < 128200000)
    
    if ("chr2" %in% all_reads$chrN) {
      all_reads_chr2 <- all_reads %>% filter(chrN == "chr2" & chrN_position > 88750000 & chrN_position < 90500000)
    } else {
      all_reads_chr2 <- NULL
    }
    if ("chr14" %in% all_reads$chrN) {
      all_reads_chr14 <- all_reads %>% filter(chrN == "chr14" & chrN_position > 105500000 & chrN_position < 107000000)
    } else {
      all_reads_chr14 <- NULL
    }
    if ("chr22" %in% all_reads$chrN) {
      all_reads_chr22 <- all_reads %>% filter(chrN == "chr22" & chrN_position > 22000000 & chrN_position < 24000000)
    } else {
      all_reads_chr22 <- NULL
    }
    
    all_reads <- bind_rows(all_reads_chr2, all_reads_chr14, all_reads_chr22) %>%
      select(read_name, chr8_position, chrN_position, chrN) %>%
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
  get_sc_chimeric_transcripts_56203 <- function(dis_reads){
    
    sc_chr2_min_max <- dis_reads %>% 
      filter(chrom == 2) %>% 
      filter(end - start < 1024) %>%
      group_by(cell_barcode, molecular_barcode) %>% 
      summarize(min_start = min(start), max_start = max(start), 
                min_end = min(end), max_end = max(end), 
                n_reads = n())
    
    sc_chr8_min_max <- dis_reads %>% 
      filter(chrom == 8) %>% 
      filter(end - start < 1024) %>%
      group_by(cell_barcode, molecular_barcode) %>% 
      summarize(min_start = min(start), max_start = max(start), 
                min_end = min(end), max_end = max(end), 
                n_reads = n())
    
    sc_chr14_min_max <- dis_reads %>% 
      filter(chrom == 14) %>% 
      filter(end - start < 1024) %>%
      group_by(cell_barcode, molecular_barcode) %>% 
      summarize(min_start = min(start), max_start = max(start), 
                min_end = min(end), max_end = max(end), 
                n_reads = n())
    
    sc_chr22_min_max <- dis_reads %>% 
      filter(chrom == 22) %>% 
      filter(end - start < 1024) %>%
      group_by(cell_barcode, molecular_barcode) %>% 
      summarize(min_start = min(start), max_start = max(start), 
                min_end = min(end), max_end = max(end), 
                n_reads = n())
    
    chr8chr2_combined <- sc_chr8_min_max %>%
      full_join(sc_chr2_min_max, 
                by = c("cell_barcode", "molecular_barcode")) %>%
      filter(!any(is.na(min_start.x), is.na(min_start.y))) %>%
      mutate(chr.y = "chr2")
    
    chr8chr14_combined <- sc_chr8_min_max %>% 
      full_join(sc_chr14_min_max, 
                by = c("cell_barcode", "molecular_barcode")) %>%
      filter(!any(is.na(min_start.x), is.na(min_start.y))) %>%
      mutate(chr.y = "chr14")
    
    chr8chr22_combined <- sc_chr8_min_max %>% 
      full_join(sc_chr22_min_max, 
                by = c("cell_barcode", "molecular_barcode")) %>%
      filter(!any(is.na(min_start.x), is.na(min_start.y))) %>%
      mutate(chr.y = "chr22")
    
    bind_rows(chr8chr2_combined, chr8chr14_combined, chr8chr22_combined) %>%
      ungroup() %>% 
      return()
    
  }
  get_bulk_sc_plot_df_56203 <- function(bulk_fusion_reads, sc_chimeric_transcripts){
    
    between_genes <- 0.2e6
    
    if (!identical(bulk_fusion_reads, NULL)) {
      bulk_chr8_min <- bulk_fusion_reads %>% pull(chr8_position) %>% min()
      bulk_chr8_max <- bulk_fusion_reads %>% pull(chr8_position) %>% max()
      
      if ("chr2" %in% bulk_fusion_reads$chrN) {
        bulk_chr2_min <- bulk_fusion_reads %>% filter(chrN == "chr2") %>% pull(chrN_position) %>% min()
        bulk_chr2_max <- bulk_fusion_reads %>% filter(chrN == "chr2") %>% pull(chrN_position) %>% max()
      } else {
        bulk_chr2_min <- NULL
        bulk_chr2_max <- NULL
      }
      if ("chr14" %in% bulk_fusion_reads$chrN) {
        bulk_chr14_min <- bulk_fusion_reads %>% filter(chrN == "chr14") %>% pull(chrN_position) %>% min()
        bulk_chr14_max <- bulk_fusion_reads %>% filter(chrN == "chr14") %>% pull(chrN_position) %>% max()
      } else {
        bulk_chr14_min <- NULL
        bulk_chr14_max <- NULL
      }
      if ("chr22" %in% bulk_fusion_reads$chrN) {
        bulk_chr22_min <- bulk_fusion_reads %>% filter(chrN == "chr22") %>% pull(chrN_position) %>% min()
        bulk_chr22_max <- bulk_fusion_reads %>% filter(chrN == "chr22") %>% pull(chrN_position) %>% max()  
      } else {
        bulk_chr22_min <- NULL
        bulk_chr22_max <- NULL
      }  
    } else{
      bulk_chr8_min <- NULL
      bulk_chr8_max <- NULL
      bulk_chr2_min <- NULL
      bulk_chr2_max <- NULL
      bulk_chr14_min <- NULL
      bulk_chr14_max <- NULL
      bulk_chr22_min <- NULL
      bulk_chr22_max <- NULL
    }
    
    sc_chr8_min <- sc_chimeric_transcripts %>% pull(min_start.x) %>% min()
    sc_chr8_max <- sc_chimeric_transcripts %>% pull(min_start.x) %>% max()
    
    if ("chr2" %in% sc_chimeric_transcripts$chr.y) {
      sc_chr2_min <- sc_chimeric_transcripts %>% filter(chr.y == "chr2") %>% pull(min_start.y) %>% min()
      sc_chr2_max <- sc_chimeric_transcripts %>% filter(chr.y == "chr2") %>% pull(min_start.y) %>% max()  
    } else {
      sc_chr2_min <- NULL
      sc_chr2_max <- NULL
    }
    if ("chr14" %in% sc_chimeric_transcripts$chr.y) {
      sc_chr14_min <- sc_chimeric_transcripts %>% filter(chr.y == "chr14") %>% pull(min_start.y) %>% min()
      sc_chr14_max <- sc_chimeric_transcripts %>% filter(chr.y == "chr14") %>% pull(min_start.y) %>% max()  
    } else {
      sc_chr14_min <- NULL
      sc_chr14_max <- NULL
    }
    if ("chr22" %in% sc_chimeric_transcripts$chr.y) {
      sc_chr22_min <- sc_chimeric_transcripts %>% filter(chr.y == "chr22") %>% pull(min_start.y) %>% min()
      sc_chr22_max <- sc_chimeric_transcripts %>% filter(chr.y == "chr22") %>% pull(min_start.y) %>% max()  
    } else {
      sc_chr22_min <- NULL
      sc_chr22_max <- NULL
    }
    
    if (identical(bulk_chr2_min, NULL) & identical(sc_chr2_min, NULL)) {
      overall_chr2_min <- NULL
      overall_chr2_max <- NULL
      chr2_shift <- NULL
    } else {
      overall_chr2_min <- min(bulk_chr2_min, sc_chr2_min)
      overall_chr2_max <- max(bulk_chr2_max, sc_chr2_max)
      chr2_shift <- 0 - overall_chr2_min
    }
    if (identical(bulk_chr14_min, NULL) & identical(sc_chr14_min, NULL)) {
      overall_chr14_min <- NULL
      overall_chr14_max <- NULL
      chr14_shift <- NULL
    } else {
      overall_chr14_min <- min(bulk_chr14_min, sc_chr14_min)
      overall_chr14_max <- max(bulk_chr14_max, sc_chr14_max)
      chr14_shift = 0 - overall_chr14_min
    }
    if (identical(bulk_chr22_min, NULL) & identical(sc_chr22_min, NULL)) {
      overall_chr22_min <- NULL
      overall_chr22_max <- NULL
      chr22_shift <- NULL
    } else {
      overall_chr22_min <- min(bulk_chr22_min, sc_chr22_min)
      overall_chr22_max <- max(bulk_chr22_max, sc_chr22_max)
      chr22_shift = 0 - overall_chr22_min
    }
    
    get_bulk_plot_df <- function(bulk_fusion_reads, this_chrN, chrN_shift, chr8_shift){
      
      if (identical(bulk_fusion_reads, NULL)) {
        return(NULL)
      } else {
        bulk_fusion_reads %>% 
          filter(chrN == this_chrN) %>%
          mutate(shifted_chrN = chrN_position + chrN_shift,
                 shifted_chr8 = chr8_position + chr8_shift) %>%
          mutate(category = "None") %>%
          mutate(identifier = read_name,
                 data_type = "Bulk Spanning/Junction Read Pair") %>%
          ungroup() %>%
          select(identifier, data_type, chr8_position, chrN_position, shifted_chr8, shifted_chrN, category, fusion, chrN) %>%
          return()  
      }
    } 
    
    get_sc_plot_df <- function(sc_chimeric_transcripts, this_chrN, chrN_shift, chr8_shift) {
      
      sc_chimeric_transcripts %>%
        filter(chr.y == this_chrN) %>%
        mutate(chrN_position = min_start.y,
               chr8_position = max_end.x) %>%
               #chr8_position = min_start.x) %>%
        mutate(shifted_chrN = chrN_position + chrN_shift,
               shifted_chr8 = chr8_position + chr8_shift) %>%
        mutate(category = "None") %>%
        mutate(identifier = str_c(cell_barcode, ":", molecular_barcode),
               data_type = "Single Cell Chimeric Transcript",
               fusion = NA) %>%
        ungroup() %>%
        select(identifier, data_type, chr8_position, chrN_position, shifted_chr8, shifted_chrN, category, fusion, chr.y) %>%
        rename("chrN" = "chr.y") %>%
        return()
      
    }
    
    overall_chr8_min <- min(bulk_chr8_min, sc_chr8_min)
    overall_chr8_max <- max(bulk_chr8_max, sc_chr8_max)
    
    chr2_bulk_plot_df <- get_bulk_plot_df(bulk_fusion_reads, "chr2", chr2_shift, chr8_shift = 0 - overall_chr8_min + overall_chr2_max - overall_chr2_min + between_genes)
    chr14_bulk_plot_df <- get_bulk_plot_df(bulk_fusion_reads, "chr14", chr14_shift, chr8_shift = 0 - overall_chr8_min + overall_chr14_max - overall_chr14_min + between_genes)
    chr22_bulk_plot_df <- get_bulk_plot_df(bulk_fusion_reads, "chr22", chr22_shift, chr8_shift = 0 - overall_chr8_min + overall_chr22_max - overall_chr22_min + between_genes)
    
    chr2_sc_plot_df <- get_sc_plot_df(sc_chimeric_transcripts, "chr2", chr2_shift, chr8_shift = 0 - overall_chr8_min + overall_chr2_max - overall_chr2_min + between_genes)
    chr14_sc_plot_df <- get_sc_plot_df(sc_chimeric_transcripts, "chr14", chr14_shift, chr8_shift = 0 - overall_chr8_min + overall_chr14_max - overall_chr14_min + between_genes)
    chr22_sc_plot_df <- get_sc_plot_df(sc_chimeric_transcripts, "chr22", chr22_shift, chr8_shift = 0 - overall_chr8_min + overall_chr22_max - overall_chr22_min + between_genes)
    
    chr2_bulk_sc_plot_df <- bind_rows(chr2_bulk_plot_df, chr2_sc_plot_df)
    chr14_bulk_sc_plot_df <- bind_rows(chr14_bulk_plot_df, chr14_sc_plot_df)
    chr22_bulk_sc_plot_df <- bind_rows(chr22_bulk_plot_df, chr22_sc_plot_df)
    
    bind_rows(chr2_bulk_sc_plot_df, 
              chr14_bulk_sc_plot_df, 
              chr22_bulk_sc_plot_df) %>%
      mutate(color_column = chrN) %>%
      return()
  }
  plot_bulk_sc_56203 <- function(bulk_sc_plot_df, genes = gene_spans, dir, id, use_fusion_names = TRUE, bulk_color = TRUE, color_indicates = "Fusion Name", partner_chr = NULL){
    
    if (identical(partner_chr, NULL) | !(partner_chr %in% c("chr2", "chr14", "chr22"))) {
      stop("Partner chromosome partner_chr must be one of 'chr2', 'chr14', or 'chr22'")
    }
    
    between_genes <- 0.2e6
    chr8_offset_y <- 0.5
    
    bulk_sc_plot_df <- bulk_sc_plot_df %>%
      filter(chrN == partner_chr) %>%
      mutate(color_column = "black",
             reciprocal = FALSE)
    
    overall_chr8_min <- bulk_sc_plot_df %>% pull(chr8_position) %>% min()
    overall_chr8_max <- bulk_sc_plot_df %>% pull(chr8_position) %>% max()
    overall_chrN_min <- bulk_sc_plot_df %>% pull(chrN_position) %>% min()
    overall_chrN_max <- bulk_sc_plot_df %>% pull(chrN_position) %>% max()
    
    chr8_shift <- 0 - overall_chr8_min + overall_chrN_max - overall_chrN_min + between_genes
    chrN_shift <- 0 - overall_chrN_min
    
    chr8_gene_spans <- genes %>% 
      filter(chromosome == "chr8") %>% 
      filter( gene_name %in% c("MYC", "PVT1") |
                (start <= overall_chr8_min & end >= overall_chr8_min) | 
                (start >= overall_chr8_min & end <= overall_chr8_max) | 
                (start <= overall_chr8_max & end >= overall_chr8_max) ) %>% 
      filter(strand == "+", type %in% c("protein_coding", "lincRNA")) %>%
      mutate(shifted_chr8_start = start + chr8_shift,
             shifted_chr8_end = end + chr8_shift)
    
    chrN_gene_spans <- genes %>% 
      filter(chromosome == partner_chr) %>% 
      filter( (start <= overall_chrN_min & end >= overall_chrN_min) | 
                (start >= overall_chrN_min & end <= overall_chrN_max) | 
                (start <= overall_chrN_max & end >= overall_chrN_max) ) %>% 
      filter(str_detect(gene_name, "^IG")) %>%
      filter(!str_detect(gene_name, "@")) %>%
      mutate(shifted_chrN_start = start + chrN_shift,
             shifted_chrN_end = end + chrN_shift)
    
    shifted_gene_boundaries <- c(chr8_gene_spans %>% pull(shifted_chr8_start) %>% min(),
                                 chr8_gene_spans %>% pull(shifted_chr8_start) %>% max(),
                                 chr8_gene_spans %>% pull(shifted_chr8_end) %>% max(),
                                 chrN_gene_spans %>% pull(shifted_chrN_start) %>% min(),
                                 chrN_gene_spans %>% pull(shifted_chrN_end) %>% max())
    
    gene_boundaries <- c(chr8_gene_spans %>% pull(start) %>% min(),
                         chr8_gene_spans %>% pull(start) %>% max(),
                         chr8_gene_spans %>% pull(end) %>% max(),
                         chrN_gene_spans %>% pull(start) %>% min(),
                         chrN_gene_spans %>% pull(end) %>% max())
    
    p <- ggplot(data = bulk_sc_plot_df) +
      scale_y_continuous(limits = c(-1.5, 2)) +
      scale_x_continuous(breaks = shifted_gene_boundaries,
                         labels = gene_boundaries) +
      #chr8 genes
      geom_rect(data = chr8_gene_spans,
                aes(xmin = shifted_chr8_start, 
                    xmax = shifted_chr8_end,
                    ymin = 0 - chr8_offset_y, 
                    ymax = 1 - chr8_offset_y,
                    fill = "chr8"),
                #fill = "#f03b20"
                alpha = 1) +
      #chrN genes
      geom_rect(data = chrN_gene_spans,
                aes(xmin = shifted_chrN_start, 
                    xmax = shifted_chrN_end,
                    ymin = 0, 
                    ymax = 1,
                    fill = "chrN"), 
                #fill = "#feb24c",
                alpha = 1)
    
    if (bulk_sc_plot_df %>% 
        filter(data_type == "Bulk Spanning/Junction Read Pair") %>% 
        nrow() > 0) {
      # draw bulk read curves
      p <- p + geom_curve(data = bulk_sc_plot_df %>% 
                            filter(data_type == "Bulk Spanning/Junction Read Pair"),
                          aes(x = shifted_chrN, 
                              xend = shifted_chr8,
                              color = color_column,
                              linetype = reciprocal,
                              y = 1,
                              yend = 1 - chr8_offset_y),
                          curvature = -0.25,
                          ncp = 10,
                          lineend = "round",
                          alpha = 0.25)
    }
    
    # draw sc read curves
    p <- p + geom_curve(data = bulk_sc_plot_df %>% 
                          filter(data_type == "Single Cell Chimeric Transcript"),
                        aes(x = shifted_chrN, 
                            xend = shifted_chr8,
                            y = 0,
                            yend = 0 - chr8_offset_y,
                            color = "black"),
                        curvature = 0.25,
                        ncp = 10,
                        lineend = "round",
                        show.legend = FALSE,
                        alpha = 0.25) +
      # label chr8 genes
      geom_text(data = chr8_gene_spans,
                aes(x = (shifted_chr8_end + shifted_chr8_start)/2,
                    y = 0,
                    label = gene_name),
                fontface = "italic") +
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
    
    ggsave(str_c(dir, id, ".", partner_chr, ".discordant_read_positions.pdf"),
           p,
           width = 14.5, height = 12,
           useDingbats = FALSE)
    
    ggsave(str_c(dir, id, ".", partner_chr, ".discordant_read_positions.no_legend.pdf"),
           p + guides(color = FALSE, linetype = FALSE, fill = FALSE),
           width = 7.25, height = 12,
           useDingbats = FALSE)
  }
  
  # 56203_1
  if (TRUE) {
    bulk_sc_plot_df_56203_1 <- get_bulk_sc_plot_df_56203(bulk_fusion_reads = NULL,
                                                         sc_chimeric_transcripts = get_sc_chimeric_transcripts_56203(dis_reads = dis_reads_56203_1_discover))
    
    plot_bulk_sc_56203(bulk_sc_plot_df_56203_1, dir = str_c(paper_supp, "56203_1/"), id = "56203_1", partner_chr = "chr14")
    plot_bulk_sc_56203(bulk_sc_plot_df_56203_1, dir = str_c(paper_supp, "56203_1/"), id = "56203_1", partner_chr = "chr2")
    plot_bulk_sc_56203(bulk_sc_plot_df_56203_1, dir = str_c(paper_supp, "56203_1/"), id = "56203_1", partner_chr = "chr22")
    
    plot_cell_chimeric_transcripts(bulk_sc = bulk_sc_plot_df_56203_1,
                                   tsne_umap = get_tsne_umap(cell_types = cell_types_56203_1, seurat_object = seurat_object_56203_1), 
                                   id = "56203_2", 
                                   reduction = "UMAP", 
                                   dir = str_c(paper_main, "56203_1/"), 
                                   facet = TRUE)
    plot_cell_chimeric_transcripts(bulk_sc = bulk_sc_plot_df_56203_1,
                                   tsne_umap = get_tsne_umap(cell_types = cell_types_56203_1, seurat_object = seurat_object_56203_1), 
                                   id = "56203_2", 
                                   reduction = "UMAP", 
                                   dir = str_c(paper_main, "56203_1/"), 
                                   facet = FALSE)
    
    plot_gene_expression_violin(bulk_sc = bulk_sc_plot_df_56203_1, 
                                seurat_object = seurat_object_56203_1, 
                                tsne_umap = get_tsne_umap(cell_types = cell_types_56203_1, seurat_object = seurat_object_56203_1), 
                                ensg = "ENSG00000136997", 
                                id = "56203_1", 
                                gene = "MYC", 
                                dir = str_c(paper_supp, "56203_1/"), 
                                max = 3.5)
    plot_gene_expression_violin(bulk_sc = bulk_sc_plot_df_56203_1, 
                                seurat_object = seurat_object_56203_1, 
                                tsne_umap = get_tsne_umap(cell_types = cell_types_56203_1, seurat_object = seurat_object_56203_1), 
                                ensg = "ENSG00000249859", 
                                id = "56203_1", 
                                gene = "PVT1", 
                                dir = str_c(paper_supp, "56203_1/"), 
                                max = 3.5) 
  }
  
  # 56203_2
  if (TRUE) {
    bulk_sc_plot_df_56203_2 <- get_bulk_sc_plot_df_56203(bulk_fusion_reads = get_bulk_fusion_reads_56203(bulk_reads = bulk_reads_56203_2), 
                                                         sc_chimeric_transcripts = get_sc_chimeric_transcripts_56203(dis_reads = dis_reads_56203_2_discover))
    
    plot_bulk_sc_56203(bulk_sc_plot_df_56203_2, dir = str_c(paper_supp, "56203_2/"), id = "56203_2", partner_chr = "chr14")
    plot_bulk_sc_56203(bulk_sc_plot_df_56203_2, dir = str_c(paper_supp, "56203_2/"), id = "56203_2", partner_chr = "chr2")
    plot_bulk_sc_56203(bulk_sc_plot_df_56203_2, dir = str_c(paper_supp, "56203_2/"), id = "56203_2", partner_chr = "chr22")
    
    plot_cell_chimeric_transcripts(bulk_sc = bulk_sc_plot_df_56203_2,
                                   tsne_umap = get_tsne_umap(cell_types = cell_types_56203_2, seurat_object = seurat_object_56203_2), 
                                   id = "56203_2", 
                                   reduction = "UMAP", 
                                   dir = str_c(paper_main, "56203_2/"), 
                                   facet = TRUE)
    plot_cell_chimeric_transcripts(bulk_sc = bulk_sc_plot_df_56203_2,
                                   tsne_umap = get_tsne_umap(cell_types = cell_types_56203_2, seurat_object = seurat_object_56203_2), 
                                   id = "56203_2", 
                                   reduction = "UMAP", 
                                   dir = str_c(paper_supp, "56203_2/"), 
                                   facet = FALSE)
    
    plot_gene_expression_violin(bulk_sc = bulk_sc_plot_df_56203_2, 
                                seurat_object = seurat_object_56203_2, 
                                tsne_umap = get_tsne_umap(cell_types = cell_types_56203_2, seurat_object = seurat_object_56203_2), 
                                ensg = "ENSG00000136997", 
                                id = "56203_2", 
                                gene = "MYC", 
                                dir = str_c(paper_supp, "56203_2/"), 
                                max = 3.75)
    plot_gene_expression_violin(bulk_sc = bulk_sc_plot_df_56203_2, 
                                seurat_object = seurat_object_56203_2, 
                                tsne_umap = get_tsne_umap(cell_types = cell_types_56203_2, seurat_object = seurat_object_56203_2), 
                                ensg = "ENSG00000249859", 
                                id = "56203_2", 
                                gene = "PVT1", 
                                dir = str_c(paper_supp, "56203_2/"), 
                                max = 3.75) 
  }
  
}

# ==============================================================================
# Work with 77570
# ==============================================================================
if (TRUE) {
  get_bulk_fusion_reads_77570 <- function(bulk_reads, star_fusion = NULL, use_SF_only = FALSE){
    
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
      filter(chromosome_donor %in% c("chr11", "chr14"), 
             chromosome_acceptor %in% c("chr11", "chr14")) %>%
      select(read_name, 
             chromosome_donor, first_base_donor, 
             chromosome_acceptor, first_base_acceptor) %>% 
      unique() %>% 
      mutate(chr14_position = case_when(chromosome_donor == "chr14" ~ first_base_donor,
                                        TRUE ~ first_base_acceptor),
             chr11_position = case_when(chromosome_donor == "chr11" ~ first_base_donor,
                                        TRUE ~ first_base_acceptor)) %>%
      filter(chr14_position > 105500000,
             chr14_position < 107000000) %>%
      filter(chr11_position > 69640000,
             chr11_position < 69660000) %>%
      select(read_name, chr14_position, chr11_position) %>%
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
  get_sc_chimeric_transcripts_77570 <- function(dis_reads){
    
    sc_chr14_min_max <- dis_reads %>% 
      filter(chrom == 14) %>% 
      filter(end - start < 1024) %>%
      group_by(cell_barcode, molecular_barcode) %>% 
      summarize(min_start = min(start), max_start = max(start), 
                min_end = min(end), max_end = max(end), 
                n_reads = n())
    
    sc_chr11_min_max <- dis_reads %>% 
      filter(chrom == 11) %>% 
      filter(end - start < 1024) %>% 
      group_by(cell_barcode, molecular_barcode) %>% 
      summarize(min_start = min(start), max_start = max(start), 
                min_end = min(end), max_end = max(end), 
                n_reads = n())
    
    sc_chr14_min_max %>% 
      full_join(sc_chr11_min_max, 
                by = c("cell_barcode", "molecular_barcode")) %>%
      filter(!any(is.na(min_start.x), is.na(min_start.y))) %>%
      return()
  }
  get_bulk_sc_plot_df_77570 <- function(bulk_fusion_reads, sc_chimeric_transcripts){
    
    chr11_breakpoint <- 69404964
    chr14_breakpoint <- 105741943
    between_genes <- 0.4e6
    
    if (bulk_fusion_reads %>% nrow() == 0) {
      bulk_fusion_reads <- NULL
    }
    
    if (!identical(bulk_fusion_reads, NULL)) {
      bulk_chr14_min <- bulk_fusion_reads %>% pull(chr14_position) %>% min()
      bulk_chr14_max <- bulk_fusion_reads %>% pull(chr14_position) %>% max()
      bulk_chr11_min <- bulk_fusion_reads %>% pull(chr11_position) %>% min()
      bulk_chr11_max <- bulk_fusion_reads %>% pull(chr11_position) %>% max()
    } else {
      bulk_chr11_min <- NULL
      bulk_chr11_max <- NULL
      bulk_chr14_min <- NULL
      bulk_chr14_max <- NULL
    }
    
    sc_chr14_min <- sc_chimeric_transcripts %>% pull(min_start.x) %>% min()
    sc_chr14_max <- sc_chimeric_transcripts %>% pull(min_start.x) %>% max()
    sc_chr11_min <- sc_chimeric_transcripts %>% pull(min_start.y) %>% min()
    sc_chr11_max <- sc_chimeric_transcripts %>% pull(min_start.y) %>% max()
    
    overall_chr14_min <- min(bulk_chr14_min, sc_chr14_min)
    overall_chr14_max <- max(bulk_chr14_max, sc_chr14_max)
    overall_chr11_min <- min(bulk_chr11_min, sc_chr11_min)
    overall_chr11_max <- max(bulk_chr11_max, sc_chr11_max)
    
    chr11_shift = 0 - overall_chr11_min + overall_chr14_max - overall_chr14_min + between_genes
    chr14_shift = 0 - overall_chr14_min
    
    shifted_chr11_breakpoint = chr11_breakpoint + chr11_shift
    shifted_chr14_breakpoint = chr14_breakpoint + chr14_shift
    
    if (identical(bulk_fusion_reads, NULL)) {
      bulk_plot_df <- NULL
    } else {
      bulk_plot_df <- bulk_fusion_reads %>%
        mutate(category = case_when(chr11_position <= chr11_breakpoint & chr14_position <= chr14_breakpoint ~ "Supports Reciprocal CCND1--IGH Fusion",
                                    chr11_position > chr11_breakpoint & chr14_position <= chr14_breakpoint ~ "Supports General IGH Fusion with CCND1",
                                    chr11_position <= chr11_breakpoint & chr14_position > chr14_breakpoint ~ "Supports General IGH Fusion with CCND1",
                                    chr11_position > chr11_breakpoint & chr14_position > chr14_breakpoint ~ "Supports IGH--CCND1 Fusion")) %>%
        mutate(category = factor(category, levels = c("Supports IGH--CCND1 Fusion", "Supports Reciprocal CCND1--IGH Fusion", "Supports General IGH Fusion with CCND1"), ordered = TRUE)) %>%
        mutate(shifted_chr11 = chr11_position + chr11_shift,
               shifted_chr14 = chr14_position + chr14_shift) %>%
        mutate(identifier = read_name,
               data_type = "Bulk Spanning/Junction Read Pair",
               y_position_offset = 0,
               curvature = -0.25) %>%
        select(identifier, data_type, chr14_position, chr11_position, shifted_chr14, shifted_chr11, category, fusion)  
    }
    
    sc_plot_df <- sc_chimeric_transcripts %>%
      mutate(category = case_when((min_start.x <= chr14_breakpoint & max_end.x > chr14_breakpoint) | (min_start.y < chr11_breakpoint & max_end.y > chr11_breakpoint) ~ "Supports General IGH Fusion with CCND1",
                                  min_start.x > chr14_breakpoint & min_start.y > chr11_breakpoint ~ "Supports IGH--CCND1 Fusion",
                                  max_end.x <= chr14_breakpoint & max_end.y <= chr11_breakpoint ~ "Supports Reciprocal CCND1--IGH Fusion",
                                  max_end.x <= chr14_breakpoint & min_start.y > chr11_breakpoint ~ "Supports General IGH Fusion with CCND1",
                                  TRUE ~ "Supports General IGH Fusion with CCND1")) %>%
      mutate(chr11_position = case_when(category == "Supports IGH--CCND1 Fusion" ~ min_start.y,
                                        category == "Supports Reciprocal CCND1--IGH Fusion" ~ max_end.y,
                                        category == "Supports General IGH Fusion with CCND1" ~ min_start.y),
             chr14_position = case_when(category == "Supports IGH--CCND1 Fusion" ~ min_start.x,
                                        category == "Supports Reciprocal CCND1--IGH Fusion" ~ max_end.x,
                                        category == "Supports General IGH Fusion with CCND1" ~ min_start.x)) %>%
      mutate(category = factor(category, levels = c("Supports IGH--CCND1 Fusion", "Supports Reciprocal CCND1--IGH Fusion", "Supports General IGH Fusion with CCND1"), ordered = TRUE)) %>%
      mutate(shifted_chr11 = chr11_position + chr11_shift,
             shifted_chr14 = chr14_position + chr14_shift) %>%
      mutate(identifier = str_c(cell_barcode, ":", molecular_barcode),
             data_type = "Single Cell Chimeric Transcript",
             y_position_offset = -0.5,
             curvature = 0.25,
             fusion = NA) %>%
      ungroup() %>%
      select(identifier, data_type, chr14_position, chr11_position, shifted_chr14, shifted_chr11, category, fusion)
    
    bind_rows(bulk_plot_df, sc_plot_df) %>% 
      mutate(color_column = category) %>%
      return()
  }
  plot_bulk_sc_77570 <- function(bulk_sc_plot_df, genes = gene_spans, dir, id, use_fusion_names = TRUE, bulk_color = TRUE, color_indicates = "Fusion Name"){
    reported_translocation_breakpoints <- tribble(~trans, ~chrom,      ~pos,
                                                  "t(11;14)",     11,  69404964,
                                                  "t(11;14)",     14, 105741943,
                                                  "t(11;14)",     14, 105859858,
                                                  "t(11;14)",     14, 106269143)
    chr11_breakpoint <- 69404964
    chr14_breakpoint1 <- 105741943
    chr14_breakpoint2 <- 105859858
    chr14_breakpoint3 <- 106269143
    between_genes <- 0.4e6
    chr11_offset_y <- 0.5
    
    bulk_sc_plot_df <- bulk_sc_plot_df %>% 
      mutate(color_column = case_when(!bulk_color ~ "black",
                                      use_fusion_names ~ str_remove_all(str_remove_all(fusion, "CCND1"), "--"),
                                      TRUE ~ as.character(str_detect(fusion, "IGH"))),
             reciprocal = str_detect(fusion, "^CCND1"))
    
    overall_chr14_min <- bulk_sc_plot_df %>% pull(chr14_position) %>% min()
    overall_chr14_max <- bulk_sc_plot_df %>% pull(chr14_position) %>% max() 
    overall_chr11_min <- bulk_sc_plot_df %>% pull(chr11_position) %>% min()
    overall_chr11_max <- bulk_sc_plot_df %>% pull(chr11_position) %>% max()
    
    chr11_shift = 0 - overall_chr11_min + overall_chr14_max - overall_chr14_min + between_genes
    chr14_shift = 0 - overall_chr14_min
    
    shifted_chr11_breakpoint = chr11_breakpoint + chr11_shift
    shifted_chr14_breakpoint1 = chr14_breakpoint1 + chr14_shift
    shifted_chr14_breakpoint2 = chr14_breakpoint2 + chr14_shift
    shifted_chr14_breakpoint3 = chr14_breakpoint3 + chr14_shift
    
    chr14_gene_spans <- genes %>% 
      filter(chromosome == "chr14") %>% 
      filter( (start <= overall_chr14_min & end >= overall_chr14_min) | 
                (start >= overall_chr14_min & end <= overall_chr14_max) | 
                (start <= overall_chr14_max & end >= overall_chr14_max) ) %>% 
      filter(str_detect(gene_name, "IGH")) %>%
      filter(!str_detect(gene_name, "@")) %>%
      mutate(shifted_chr14_start = start + chr14_shift,
             shifted_chr14_end = end + chr14_shift)
    
    chr11_gene_spans <- genes %>% 
      filter(chromosome == "chr11") %>% 
      filter( gene_name %in% c("CCND1") |
                (start <= overall_chr11_min & end >= overall_chr11_min) | 
                (start >= overall_chr11_min & end <= overall_chr11_max) | 
                (start <= overall_chr11_max & end >= overall_chr11_max) ) %>% 
      filter(strand == "+", type == "protein_coding") %>%
      mutate(shifted_chr11_start = start + chr11_shift,
             shifted_chr11_end = end + chr11_shift)
    
    shifted_gene_boundaries <- c(chr14_gene_spans %>% pull(shifted_chr14_start) %>% min(),
                                 chr14_gene_spans %>% pull(shifted_chr14_end) %>% max(),
                                 chr11_gene_spans %>% pull(shifted_chr11_start) %>% min(),
                                 chr11_gene_spans %>% pull(shifted_chr11_start) %>% max(),
                                 chr11_gene_spans %>% pull(shifted_chr11_end) %>% max())
    
    gene_boundaries <- c(chr14_gene_spans %>% pull(start) %>% min(),
                         chr14_gene_spans %>% pull(end) %>% max(),
                         chr11_gene_spans %>% pull(start) %>% min(),
                         chr11_gene_spans %>% pull(start) %>% max(),
                         chr11_gene_spans %>% pull(end) %>% max())
    
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
      #chr11 genes
      geom_rect(data = chr11_gene_spans,
                aes(xmin = shifted_chr11_start, 
                    xmax = shifted_chr11_end,
                    ymin = 0 - chr11_offset_y, 
                    ymax = 1 - chr11_offset_y,
                    fill = "chr11"), 
                #fill = "#feb24c",
                alpha = 1)
    
    if (bulk_sc_plot_df %>% 
        filter(data_type == "Bulk Spanning/Junction Read Pair") %>% 
        nrow() > 0) {
      # draw bulk read curves
      p <- p + geom_curve(data = bulk_sc_plot_df %>% 
                            filter(data_type == "Bulk Spanning/Junction Read Pair"),
                          aes(x = shifted_chr14, 
                              xend = shifted_chr11,
                              color = color_column,
                              linetype = reciprocal,
                              y = 1,
                              yend = 1 - chr11_offset_y),
                          curvature = -0.25,
                          ncp = 10,
                          lineend = "round",
                          alpha = 0.25)
    } 
    # draw sc read curves
    p <- p + geom_curve(data = bulk_sc_plot_df %>% 
                          filter(data_type == "Single Cell Chimeric Transcript"),
                        aes(x = shifted_chr14, 
                            xend = shifted_chr11,
                            y = 0,
                            yend = 0 - chr11_offset_y,
                            color = "black"),
                        curvature = 0.25,
                        ncp = 10,
                        lineend = "round",
                        show.legend = FALSE,
                        alpha = 0.25) +
      # label chr11 genes
      geom_text(data = chr11_gene_spans,
                aes(x = (shifted_chr11_end + shifted_chr11_start)/2,
                    y = 1 - chr11_offset_y - 0.5,
                    label = gene_name),
                fontface = "italic") +
      # translocation breakpoints 
      geom_vline(xintercept = shifted_chr11_breakpoint,
                 linetype = 2) +
      geom_vline(xintercept = shifted_chr14_breakpoint1,
                 linetype = 2) +
      geom_vline(xintercept = shifted_chr14_breakpoint2,
                 linetype = 2) +
      geom_vline(xintercept = shifted_chr14_breakpoint3,
                 linetype = 2) +
      annotate(geom = "text", x = shifted_chr11_breakpoint + between_genes*0.02, y = -1.5, label = chr11_breakpoint, hjust = 0, vjust = 1) +
      annotate(geom = "text", x = shifted_chr14_breakpoint1 - between_genes*0.02, y = -1.5, label = chr14_breakpoint1, hjust = 1, vjust = 1) +
      annotate(geom = "text", x = shifted_chr14_breakpoint2 - between_genes*0.02, y = -1.5, label = chr14_breakpoint2, hjust = 1, vjust = 1) +
      annotate(geom = "text", x = shifted_chr14_breakpoint3 - between_genes*0.02, y = -1.5, label = chr14_breakpoint3, hjust = 1, vjust = 1) +
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
    
    ggsave(str_c(dir, id, ".discordant_read_positions.pdf"),
           p,
           width = 14.5, height = 12,
           useDingbats = FALSE)
    
    ggsave(str_c(dir, id, ".discordant_read_positions.no_legend.pdf"),
           p + guides(color = FALSE, linetype = FALSE, fill = FALSE),
           width = 7.25, height = 12,
           useDingbats = FALSE)
  }
  
  # 77570
  if (TRUE) {
    
    bulk_fusion_reads_77570 <- get_bulk_fusion_reads_77570(bulk_reads = bulk_reads_77570, star_fusion = NULL, use_SF_only = FALSE)
    plot_sc_read_length_distribution(dis_reads = dis_reads_77570_discover, cutoff = 1024, id = "77570", dir = str_c(paper_supp, "77570/"))
    sc_chimeric_transcripts_77570 <- get_sc_chimeric_transcripts_77570(dis_reads_77570_discover)
    bulk_sc_plot_df_77570 <- get_bulk_sc_plot_df_77570(bulk_fusion_reads_77570, sc_chimeric_transcripts_77570)
    
    plot_bulk_sc_77570(bulk_sc_plot_df_77570, id = "77570", bulk_color = FALSE, dir = str_c(paper_supp, "77570/"))
    
    plot_cell_chimeric_transcripts(bulk_sc = bulk_sc_plot_df_77570, 
                                   tsne_umap = get_tsne_umap(cell_types = cell_types_77570, seurat_object = seurat_object_77570), 
                                   id = "77570", 
                                   reduction = "UMAP", 
                                   dir = str_c(paper_main, "77570/"))
    plot_cell_chimeric_transcripts(bulk_sc = bulk_sc_plot_df_77570, 
                                   tsne_umap = get_tsne_umap(cell_types = cell_types_77570, seurat_object = seurat_object_77570), 
                                   id = "77570", 
                                   reduction = "t-SNE", 
                                   dir = str_c(paper_supp, "77570/"))
    plot_gene_expression_violin(bulk_sc = bulk_sc_plot_df_77570, 
                                seurat_object = seurat_object_77570, 
                                tsne_umap = get_tsne_umap(cell_types = cell_types_77570, 
                                                          seurat_object = seurat_object_77570), 
                                ensg = "ENSG00000110092", 
                                id = "77570", 
                                gene = "CCND1", 
                                dir = str_c(paper_supp, "77570/"),
                                max = 5)
  }
}

# ==============================================================================
# Work with 81012
# ==============================================================================
if (TRUE) {
  get_bulk_fusion_reads_81012 <- function(bulk_reads, star_fusion = NULL, use_SF_only = FALSE){
    
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
      filter(chromosome_donor %in% c("chr11", "chr14"), 
             chromosome_acceptor %in% c("chr11", "chr14")) %>%
      select(read_name, 
             chromosome_donor, first_base_donor, 
             chromosome_acceptor, first_base_acceptor) %>% 
      unique() %>% 
      mutate(chr14_position = case_when(chromosome_donor == "chr14" ~ first_base_donor,
                                        TRUE ~ first_base_acceptor),
             chr11_position = case_when(chromosome_donor == "chr11" ~ first_base_donor,
                                        TRUE ~ first_base_acceptor)) %>%
      filter(chr14_position > 105500000,
             chr14_position < 107000000) %>%
      filter(chr11_position > 69640000,
             chr11_position < 69660000) %>%
      select(read_name, chr14_position, chr11_position) %>%
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
  get_sc_chimeric_transcripts_81012 <- function(dis_reads){
    
    sc_chr14_min_max <- dis_reads %>% 
      filter(chrom == 14) %>% 
      filter(end - start < 1024) %>%
      group_by(cell_barcode, molecular_barcode) %>% 
      summarize(min_start = min(start), max_start = max(start), 
                min_end = min(end), max_end = max(end), 
                n_reads = n())
    
    sc_chr11_min_max <- dis_reads %>% 
      filter(chrom == 11) %>% 
      filter(end - start < 1024) %>% 
      group_by(cell_barcode, molecular_barcode) %>% 
      summarize(min_start = min(start), max_start = max(start), 
                min_end = min(end), max_end = max(end), 
                n_reads = n())
    
    sc_chr14_min_max %>% 
      full_join(sc_chr11_min_max, 
                by = c("cell_barcode", "molecular_barcode")) %>%
      filter(!any(is.na(min_start.x), is.na(min_start.y))) %>%
      return()
  }
  get_bulk_sc_plot_df_81012 <- function(bulk_fusion_reads, sc_chimeric_transcripts){
    
    chr11_breakpoint <- NULL
    chr14_breakpoint <- NULL
    between_genes <- 0.2e6
    
    #if (bulk_fusion_reads %>% nrow() == 0) {
    #  bulk_fusion_reads <- NULL
    #}
    
    if (!identical(bulk_fusion_reads, NULL)) {
      bulk_chr14_min <- bulk_fusion_reads %>% pull(chr14_position) %>% min()
      bulk_chr14_max <- bulk_fusion_reads %>% pull(chr14_position) %>% max()
      bulk_chr11_min <- bulk_fusion_reads %>% pull(chr11_position) %>% min()
      bulk_chr11_max <- bulk_fusion_reads %>% pull(chr11_position) %>% max()
    } else {
      bulk_chr11_min <- NULL
      bulk_chr11_max <- NULL
      bulk_chr14_min <- NULL
      bulk_chr14_max <- NULL
    }
    
    sc_chr14_min <- sc_chimeric_transcripts %>% pull(min_start.x) %>% min()
    sc_chr14_max <- sc_chimeric_transcripts %>% pull(min_start.x) %>% max()
    sc_chr11_min <- sc_chimeric_transcripts %>% pull(min_start.y) %>% min()
    sc_chr11_max <- sc_chimeric_transcripts %>% pull(min_start.y) %>% max()
    
    overall_chr14_min <- min(bulk_chr14_min, sc_chr14_min)
    overall_chr14_max <- max(bulk_chr14_max, sc_chr14_max)
    overall_chr11_min <- min(bulk_chr11_min, sc_chr11_min)
    overall_chr11_max <- max(bulk_chr11_max, sc_chr11_max)
    
    chr11_shift = 0 - overall_chr11_min + overall_chr14_max - overall_chr14_min + between_genes
    chr14_shift = 0 - overall_chr14_min
    
    #shifted_chr11_breakpoint = chr11_breakpoint + chr11_shift
    #shifted_chr14_breakpoint = chr14_breakpoint + chr14_shift
    shifted_chr11_breakpoint <- NULL
    shifted_chr11_breakpoint <- NULL
    
    if (identical(bulk_fusion_reads, NULL)) {
      bulk_plot_df <- NULL
    } else {
      bulk_plot_df <- bulk_fusion_reads %>%
        #mutate(category = case_when(chr11_position <= chr11_breakpoint & chr14_position <= chr14_breakpoint ~ "Supports Reciprocal CCND1--IGH Fusion",
        #                            chr11_position > chr11_breakpoint & chr14_position <= chr14_breakpoint ~ "Supports General IGH Fusion with CCND1",
        #                            chr11_position <= chr11_breakpoint & chr14_position > chr14_breakpoint ~ "Supports General IGH Fusion with CCND1",
        #                            chr11_position > chr11_breakpoint & chr14_position > chr14_breakpoint ~ "Supports IGH--CCND1 Fusion")) %>%
        #mutate(category = factor(category, levels = c("Supports IGH--CCND1 Fusion", "Supports Reciprocal CCND1--IGH Fusion", "Supports General IGH Fusion with CCND1"), ordered = TRUE)) %>%
        mutate(category = "None") %>%
        mutate(shifted_chr11 = chr11_position + chr11_shift,
               shifted_chr14 = chr14_position + chr14_shift) %>%
        mutate(identifier = read_name,
               data_type = "Bulk Spanning/Junction Read Pair",
               y_position_offset = 0,
               curvature = -0.25) %>%
        select(identifier, data_type, chr14_position, chr11_position, shifted_chr14, shifted_chr11, category, fusion)  
    }
    
    sc_plot_df <- sc_chimeric_transcripts %>%
      #mutate(category = case_when((min_start.x <= chr14_breakpoint & max_end.x > chr14_breakpoint) | (min_start.y < chr11_breakpoint & max_end.y > chr11_breakpoint) ~ "Supports General IGH Fusion with CCND1",
      #                            min_start.x > chr14_breakpoint & min_start.y > chr11_breakpoint ~ "Supports IGH--CCND1 Fusion",
      #                            max_end.x <= chr14_breakpoint & max_end.y <= chr11_breakpoint ~ "Supports Reciprocal CCND1--IGH Fusion",
      #                            max_end.x <= chr14_breakpoint & min_start.y > chr11_breakpoint ~ "Supports General IGH Fusion with CCND1",
      #                            TRUE ~ "Supports General IGH Fusion with CCND1")) %>%
      #mutate(chr11_position = case_when(category == "Supports IGH--CCND1 Fusion" ~ min_start.y,
      #                                  category == "Supports Reciprocal CCND1--IGH Fusion" ~ max_end.y,
      #                                  category == "Supports General IGH Fusion with CCND1" ~ min_start.y),
      #       chr14_position = case_when(category == "Supports IGH--CCND1 Fusion" ~ min_start.x,
      #                                  category == "Supports Reciprocal CCND1--IGH Fusion" ~ max_end.x,
      #                                  category == "Supports General IGH Fusion with CCND1" ~ min_start.x)) %>%
      #mutate(category = factor(category, levels = c("Supports IGH--CCND1 Fusion", "Supports Reciprocal CCND1--IGH Fusion", "Supports General IGH Fusion with CCND1"), ordered = TRUE)) %>%
      mutate(chr11_position = min_start.y,
             chr14_position = min_start.x,
             category = "None") %>%
      mutate(shifted_chr11 = chr11_position + chr11_shift,
             shifted_chr14 = chr14_position + chr14_shift) %>%
      mutate(identifier = str_c(cell_barcode, ":", molecular_barcode),
             data_type = "Single Cell Chimeric Transcript",
             y_position_offset = -0.5,
             curvature = 0.25,
             fusion = NA) %>%
      ungroup() %>%
      select(identifier, data_type, chr14_position, chr11_position, shifted_chr14, shifted_chr11, category, fusion)
    
    bind_rows(bulk_plot_df, sc_plot_df) %>% 
      mutate(color_column = category) %>%
      return()
  }
  plot_bulk_sc_81012 <- function(bulk_sc_plot_df, genes = gene_spans, dir, id, use_fusion_names = TRUE, bulk_color = TRUE, color_indicates = "Fusion Name"){
    #reported_translocation_breakpoints <- tribble(~trans, ~chrom,      ~pos,
    #                                              "t(11;14)",     11,  69404964,
    #                                              "t(11;14)",     14, 105741943,
    #                                              "t(11;14)",     14, 105859858,
    #                                              "t(11;14)",     14, 106269143)
    chr11_breakpoint <- NULL
    chr14_breakpoint <- NULL
    between_genes <- 0.2e6
    chr11_offset_y <- 0.5
    
    bulk_sc_plot_df <- bulk_sc_plot_df %>% 
      mutate(color_column = case_when(!bulk_color ~ "black",
                                      use_fusion_names ~ str_remove_all(str_remove_all(fusion, "CCND1"), "--"),
                                      TRUE ~ as.character(str_detect(fusion, "IGH"))),
             reciprocal = str_detect(fusion, "^CCND1"))
    
    overall_chr14_min <- bulk_sc_plot_df %>% pull(chr14_position) %>% min()
    overall_chr14_max <- bulk_sc_plot_df %>% pull(chr14_position) %>% max() 
    overall_chr11_min <- bulk_sc_plot_df %>% pull(chr11_position) %>% min()
    overall_chr11_max <- bulk_sc_plot_df %>% pull(chr11_position) %>% max()
    
    chr11_shift = 0 - overall_chr11_min + overall_chr14_max - overall_chr14_min + between_genes
    chr14_shift = 0 - overall_chr14_min
    
    #shifted_chr11_breakpoint = chr11_breakpoint + chr11_shift
    #shifted_chr14_breakpoint1 = chr14_breakpoint1 + chr14_shift
    #shifted_chr14_breakpoint2 = chr14_breakpoint2 + chr14_shift
    #shifted_chr14_breakpoint3 = chr14_breakpoint3 + chr14_shift
    
    chr14_gene_spans <- genes %>% 
      filter(chromosome == "chr14") %>% 
      filter( (start <= overall_chr14_min & end >= overall_chr14_min) | 
                (start >= overall_chr14_min & end <= overall_chr14_max) | 
                (start <= overall_chr14_max & end >= overall_chr14_max) ) %>% 
      filter(str_detect(gene_name, "IGH")) %>%
      filter(!str_detect(gene_name, "@")) %>%
      mutate(shifted_chr14_start = start + chr14_shift,
             shifted_chr14_end = end + chr14_shift)
    
    chr11_gene_spans <- genes %>% 
      filter(chromosome == "chr11") %>% 
      filter( gene_name %in% c("CCND1") |
                (start <= overall_chr11_min & end >= overall_chr11_min) | 
                (start >= overall_chr11_min & end <= overall_chr11_max) | 
                (start <= overall_chr11_max & end >= overall_chr11_max) ) %>% 
      filter(strand == "+", type == "protein_coding") %>%
      mutate(shifted_chr11_start = start + chr11_shift,
             shifted_chr11_end = end + chr11_shift)
    
    shifted_gene_boundaries <- c(chr14_gene_spans %>% pull(shifted_chr14_start) %>% min(),
                                 chr14_gene_spans %>% pull(shifted_chr14_end) %>% max(),
                                 chr11_gene_spans %>% pull(shifted_chr11_start) %>% min(),
                                 chr11_gene_spans %>% pull(shifted_chr11_start) %>% max(),
                                 chr11_gene_spans %>% pull(shifted_chr11_end) %>% max())
    
    gene_boundaries <- c(chr14_gene_spans %>% pull(start) %>% min(),
                         chr14_gene_spans %>% pull(end) %>% max(),
                         chr11_gene_spans %>% pull(start) %>% min(),
                         chr11_gene_spans %>% pull(start) %>% max(),
                         chr11_gene_spans %>% pull(end) %>% max())
    
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
      #chr11 genes
      geom_rect(data = chr11_gene_spans,
                aes(xmin = shifted_chr11_start, 
                    xmax = shifted_chr11_end,
                    ymin = 0 - chr11_offset_y, 
                    ymax = 1 - chr11_offset_y,
                    fill = "chr11"), 
                #fill = "#feb24c",
                alpha = 1)
    
    if (bulk_sc_plot_df %>% 
        filter(data_type == "Bulk Spanning/Junction Read Pair") %>% 
        nrow() > 0) {
      # draw bulk read curves
      p <- p + geom_curve(data = bulk_sc_plot_df %>% 
                            filter(data_type == "Bulk Spanning/Junction Read Pair"),
                          aes(x = shifted_chr14, 
                              xend = shifted_chr11,
                              color = color_column,
                              linetype = reciprocal,
                              y = 1,
                              yend = 1 - chr11_offset_y),
                          curvature = -0.25,
                          ncp = 10,
                          lineend = "round",
                          alpha = 0.25)
    } 
    # draw sc read curves
    p <- p + geom_curve(data = bulk_sc_plot_df %>% 
                          filter(data_type == "Single Cell Chimeric Transcript"),
                        aes(x = shifted_chr14, 
                            xend = shifted_chr11,
                            y = 0,
                            yend = 0 - chr11_offset_y,
                            color = "black"),
                        curvature = 0.25,
                        ncp = 10,
                        lineend = "round",
                        show.legend = FALSE,
                        alpha = 0.25) +
      # label chr11 genes
      geom_text(data = chr11_gene_spans,
                aes(x = (shifted_chr11_end + shifted_chr11_start)/2,
                    y = 1 - chr11_offset_y - 0.5,
                    label = gene_name),
                fontface = "italic") +
      # translocation breakpoints 
      #geom_vline(xintercept = shifted_chr11_breakpoint,
      #           linetype = 2) +
      #geom_vline(xintercept = shifted_chr14_breakpoint1,
      #           linetype = 2) +
      #geom_vline(xintercept = shifted_chr14_breakpoint2,
      #           linetype = 2) +
      #geom_vline(xintercept = shifted_chr14_breakpoint3,
      #           linetype = 2) +
      #annotate(geom = "text", x = shifted_chr11_breakpoint + between_genes*0.02, y = -1.5, label = chr11_breakpoint, hjust = 0, vjust = 1) +
      #annotate(geom = "text", x = shifted_chr14_breakpoint1 - between_genes*0.02, y = -1.5, label = chr14_breakpoint1, hjust = 1, vjust = 1) +
      #annotate(geom = "text", x = shifted_chr14_breakpoint2 - between_genes*0.02, y = -1.5, label = chr14_breakpoint2, hjust = 1, vjust = 1) +
      #annotate(geom = "text", x = shifted_chr14_breakpoint3 - between_genes*0.02, y = -1.5, label = chr14_breakpoint3, hjust = 1, vjust = 1) +
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
    
    ggsave(str_c(dir, id, ".discordant_read_positions.pdf"),
           p,
           width = 14.5, height = 12,
           useDingbats = FALSE)
    
    ggsave(str_c(dir, id, ".discordant_read_positions.no_legend.pdf"),
           p + guides(color = FALSE, linetype = FALSE, fill = FALSE),
           width = 7.25, height = 12,
           useDingbats = FALSE)
  }
  
  # 81012_1
  if (TRUE) {
    plot_sc_read_length_distribution(dis_reads = dis_reads_81012_1_discover, cutoff = 1024, id = "81012_1", dir = str_c(paper_supp, "81012_1/"))
    sc_chimeric_transcripts_81012_1 <- get_sc_chimeric_transcripts_81012(dis_reads_81012_1_discover)
    bulk_sc_plot_df_81012_1 <- get_bulk_sc_plot_df_81012(bulk_fusion_reads = NULL, sc_chimeric_transcripts_81012_1)
    
    plot_bulk_sc_81012(bulk_sc_plot_df_81012_1, id = "81012_1", bulk_color = FALSE, dir = str_c(paper_supp, "81012_1/"))
    
    plot_cell_chimeric_transcripts(bulk_sc = bulk_sc_plot_df_81012_1, 
                                   tsne_umap = get_tsne_umap(cell_types = cell_types_81012_1, seurat_object = seurat_object_81012_1), 
                                   id = "81012_1", 
                                   reduction = "UMAP", 
                                   dir = str_c(paper_main, "81012_1/"))
    plot_cell_chimeric_transcripts(bulk_sc = bulk_sc_plot_df_81012_1, 
                                   tsne_umap = get_tsne_umap(cell_types = cell_types_81012_1, seurat_object = seurat_object_81012_1), 
                                   id = "81012_1", 
                                   reduction = "t-SNE", 
                                   dir = str_c(paper_supp, "81012_1/"))
    plot_gene_expression_violin(bulk_sc = bulk_sc_plot_df_81012_1, 
                                seurat_object = seurat_object_81012_1, 
                                tsne_umap = get_tsne_umap(cell_types = cell_types_81012_1, 
                                                          seurat_object = seurat_object_81012_1), 
                                ensg = "ENSG00000110092", 
                                id = "81012_1", 
                                gene = "CCND1", 
                                dir = str_c(paper_supp, "81012_1/"),
                                max = 5)
  }
  
  # 81012_2
  if (TRUE) {
    plot_sc_read_length_distribution(dis_reads = dis_reads_81012_2_discover, cutoff = 1024, id = "81012_2", dir = str_c(paper_supp, "81012_2/"))
    sc_chimeric_transcripts_81012_2 <- get_sc_chimeric_transcripts_81012(dis_reads_81012_2_discover)
    bulk_sc_plot_df_81012_2 <- get_bulk_sc_plot_df_81012(bulk_fusion_reads = NULL, sc_chimeric_transcripts_81012_2)
    
    plot_bulk_sc_81012(bulk_sc_plot_df_81012_2, id = "81012_2", bulk_color = FALSE, dir = str_c(paper_supp, "81012_2/"))
    
    plot_cell_chimeric_transcripts(bulk_sc = bulk_sc_plot_df_81012_2, 
                                   tsne_umap = get_tsne_umap(cell_types = cell_types_81012_2, seurat_object = seurat_object_81012_2), 
                                   id = "81012_2", 
                                   reduction = "UMAP", 
                                   dir = str_c(paper_main, "81012_2/"))
    plot_cell_chimeric_transcripts(bulk_sc = bulk_sc_plot_df_81012_2, 
                                   tsne_umap = get_tsne_umap(cell_types = cell_types_81012_2, seurat_object = seurat_object_81012_2), 
                                   id = "81012_2", 
                                   reduction = "t-SNE", 
                                   dir = str_c(paper_supp, "81012_2/"))
    plot_gene_expression_violin(bulk_sc = bulk_sc_plot_df_81012_2, 
                                seurat_object = seurat_object_81012_2, 
                                tsne_umap = get_tsne_umap(cell_types = cell_types_81012_2, 
                                                          seurat_object = seurat_object_81012_2), 
                                ensg = "ENSG00000110092", 
                                id = "81012_2", 
                                gene = "CCND1", 
                                dir = str_c(paper_supp, "81012_2/"),
                                max = 5)
  }
}

# ==============================================================================
# Explore QC of samples
# ==============================================================================
if (TRUE) {
  bad_df <- read_tsv("data/bad_qc_reads_plot_df.tsv") %>%
    filter(chromosome %in% c(2, 14, 22)) %>%
    filter(!(chromosome == 2 & start < 88857244)) %>%
    group_by(cell_barcode, molecular_barcode, sample, test, chromosome) %>% 
    summarize(start = min(start), end = max(end)) %>%
    ungroup() %>%
    separate(test, by = "_", into = c("ig", "other_gene")) %>%
    mutate(ig_name = case_when(chromosome == 2 ~ "IGK",
                               chromosome == 14 ~ "IGH",
                               chromosome == 22 ~ "IGL"))
  
  original_df <- bind_rows(dis_reads_27522_1_discover %>% mutate(sample = "27522_1"), 
                           dis_reads_27522_4_discover %>% mutate(sample = "27522_4"), 
                           dis_reads_47499_discover %>% mutate(sample = "47499"), 
                           dis_reads_56203_1_discover %>% mutate(sample = "56203_1"), 
                           dis_reads_56203_2_discover %>% mutate(sample = "56203_2"),
                           dis_reads_77570_discover %>% mutate(sample = "77570"), 
                           dis_reads_81012_1_discover %>% mutate(sample = "81012_1"), 
                           dis_reads_81012_2_discover %>% mutate(sample = "81012_2")) %>% 
    filter(end - start < 1024) %>% 
    group_by(sample, cell_barcode, molecular_barcode, chrom) %>% 
    summarize(start = min(start), end = max(end)) %>%
    mutate(other_gene = "original") %>%
    filter(chrom %in% c(2, 14, 22)) %>%
    filter(!(chrom == 2 & start < 88857244)) %>%
    mutate(ig_name = case_when(chrom == 2 ~ "IGK",
                               chrom == 14 ~ "IGH",
                               chrom == 22 ~ "IGL"))
  
  ggplot(original_df, aes(x = start/1e6, fill = sample)) + 
    geom_histogram(binwidth = 10000/1e6) +
    #geom_vline(xintercept = 105854501) +
    #geom_vline(xintercept = 105922264) +
    facet_wrap(~ ig_name, ncol = 1, scales = "free") +
    labs(fill = "Sample ID",
         x = "Position (Mb)",
         y = "Number of 'Chimeric Transcripts'") +
    theme_bw() +
    ggsave(str_c(paper_supp, "pre_QC_discovery.pdf"), width = 10, height = 20, useDingbats = FALSE)
  
  
  ggplot(bad_df, aes(x = start/1e6, fill = other_gene)) + 
    geom_histogram(binwidth = 10000/1e6) +
    #geom_vline(xintercept = 105854501) +
    #geom_vline(xintercept = 105922264) +
    facet_wrap(~ ig_name, ncol = 1, scales = "free") +
    labs(fill = "'Partner' Gene",
         x = "Position (Mb)",
         y = "Number of 'Chimeric Transcripts'") +
    theme_bw() +
    ggsave(str_c(paper_supp, "false_positives.pdf"), width = 10, height = 20, useDingbats = FALSE)
  
}

# Draw genes
if (FALSE) {
  library(biomaRt) # 06_single_cell.R
  library(Gviz) # 06_single_cell.R
  # from Bioconductor http://bioconductor.org/packages/release/bioc/html/Gviz.html
  bm <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "hsapiens_gene_ensembl")
  
  # WHSC1 (aka NSD2)
  if (TRUE) {
    gene = "WHSC1"
    this_chrom = "chr4"
    this_start = 1.7e6
    this_end = 2.0e6
    
    pdf(str_c(paper_main, gene, ".gene_region.pdf"), width = 20, height = 10)
    axisTrack <- GenomeAxisTrack()
    biomTrack <- BiomartGeneRegionTrack(genome = "hg38",
                                        name = "ENSEMBL",
                                        chromosome = this_chrom,
                                        start = this_start,
                                        end = this_end,
                                        biomart = bm,
                                        protein_coding = "#1b9e77",
                                        utr3 = "#7570b3",
                                        utr5 = "#e7298a")
    plotTracks(list(axisTrack, biomTrack),
               col.line = NULL,
               col = NULL,
               stackHeight = 0.3,
               collapseTranscripts = "meta",
               transcriptAnnotation = "symbol")
    dev.off()
    
    pdf(str_c(paper_main, this_chrom, ".ideogram.pdf"), width = 7.5, height = 0.5)
    ideoTrack <- IdeogramTrack(genome = "hg38", chromosome = this_chrom)
    plotTracks(ideoTrack,
               from = this_start, 
               to = this_end, 
               showId = FALSE, 
               showBandId = FALSE, 
               cex.bands = 0.5)
    dev.off()
  }
  
  # IGH super locus
  if (TRUE) {
    gene = "IGH"
    this_chrom = "chr14"
    this_start = 105500000
    this_end = 107000000
    
    pdf(str_c(paper_main, gene, ".gene_region.pdf"), width = 20, height = 20)
    axisTrack <- GenomeAxisTrack()
    biomTrack <- BiomartGeneRegionTrack(genome = "hg38",
                                        name = "ENSEMBL",
                                        chromosome = this_chrom,
                                        start = this_start,
                                        end = this_end,
                                        biomart = bm)
    plotTracks(list(axisTrack, biomTrack),
               col.line = NULL,
               col = NULL,
               stackHeight = 0.2,
               collapseTranscripts = "meta",
               transcriptAnnotation = "symbol")
    dev.off()
    
    pdf(str_c(paper_main, this_chrom, ".ideogram.pdf"), width = 7.5, height = 0.5)
    ideoTrack <- IdeogramTrack(genome = "hg38", chromosome = this_chrom)
    plotTracks(ideoTrack,
               from = this_start, 
               to = this_end, 
               showId = FALSE, 
               showBandId = FALSE, 
               cex.bands = 0.5)
    dev.off()
  }

}

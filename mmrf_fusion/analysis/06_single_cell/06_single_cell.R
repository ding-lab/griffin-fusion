# ==============================================================================
# scRNA analysis
# Steven Foltz (smfoltz@wustl.edu), May 2019
# ==============================================================================

paper_main = "paper/main/06_single_cell/"
paper_supp = "paper/supplemental/06_single_cell/"

# Create directories 
dir.create(paper_main, recursive = TRUE, showWarnings = FALSE)
dir.create(paper_supp, recursive = TRUE, showWarnings = FALSE)

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
get_two_genes_expression <- function(seurat_object, tsne_umap, ensg1, engs2){
  two_genes_expression <- tibble(barcode = seurat_object@assays$RNA@data@Dimnames[[2]],
                                 expression1 = seurat_object@assays$RNA@data[ensg1,] %>% as.vector(),
                                 expression2 = seurat_object@assays$RNA@data[ensg2,] %>% as.vector())
  tsne_umap %>%
    left_join(gene_expression, by = "barcode") %>%
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
    coord_equal() +
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
  
  ggsave(str_c(dir, id, ".expression.", gene, ".pdf"),
         p, 
         width = 3.5, height = 3.5, useDingbats = FALSE)
}

# ==============================================================================
# Plot two genes expression
# ==============================================================================
plot_two_genes_expression <- function(seurat_object, tsne_umap, ensg1, ensg2, id, gene1, gene2, dir = paper_main) {
  
  plot_df <- get_two_genes_expression(seurat_object, tsne_umap, ensg1, ensg2)
  
  p <- ggplot(data = plot_df, aes(x = expression_1, y = expression_2))
  
  p <- p + labs(x = str_c(gene1, " Expression"), y = str_c(gene2, " Expression"))
  
  p <- p + geom_point(aes(color = expression),
                      shape = 16,
                      size = 1.5,
                      show.legend = FALSE,
                      alpha = 0.25) +
    theme_bw() +
    coord_fixed() +
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
plot_cell_chimeric_transcripts <- function(bulk_sc, tsne_umap, reduction = "UMAP", id, dir = paper_main){
  
  if (reduction %in% c("UMAP", "t-SNE")) {
    
    plot_df <- bulk_sc %>% 
      filter(data_type == "Single Cell Chimeric Transcript") %>% 
      separate(identifier, 
               into = c("cell_barcode", "molecular_barcode"), 
               sep = ":") %>% 
      select(cell_barcode, data_type) %>% 
      unique() %>%
      right_join(tsne_umap, 
                 by = c("cell_barcode" = "barcode")) %>%
      arrange(!is.na(data_type))
    
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
  
  p <- p + geom_point(aes(color = !is.na(data_type)),
                      shape = 16, 
                      size = 1.5, 
                      alpha = 1)
  
  p <- p + 
    labs(color = "Cell Contains\nChimeric Transcript") +
    theme_bw() +
    coord_equal() +
    scale_color_manual(values = c("#a6bddb", "#1c9099")) +
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
  
  ggsave(str_c(dir, id, ".cells_chimeric_transcripts.", reduction, ".pdf"),
         p,
         width = 3.5, height = 3.5, useDingbats = FALSE)
  
  ggsave(str_c(dir, id, ".cells_chimeric_transcripts.no_legend.", reduction, ".pdf"),
         p + guides(color = FALSE),
         width = 3.5, height = 3.5, useDingbats = FALSE)

}

# ==============================================================================
# Work with 27522 reads
# ==============================================================================
get_bulk_fusion_reads_27522 <- function(bulk_reads, star_fusion = NULL, use_SF_only = TRUE){
  
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
    mutate(shifted_chr4 = chr4_position + chr4_shift,
           shifted_chr14 = chr14_position + chr14_shift) %>%
    mutate(category = case_when(chr4_position <= chr4_breakpoint & chr14_position <= chr14_breakpoint ~ "Category 1",
                                chr4_position > chr4_breakpoint & chr14_position <= chr14_breakpoint ~ "Category 2",
                                chr4_position <= chr4_breakpoint & chr14_position > chr14_breakpoint ~ "Category 3",
                                chr4_position > chr4_breakpoint & chr14_position > chr14_breakpoint ~ "Category 4")) %>%
    mutate(identifier = read_name,
           data_type = "Bulk Spanning/Junction Read Pair",
           y_position_offset = 0,
           curvature = -0.25) %>%
    select(identifier, data_type, chr14_position, chr4_position, shifted_chr14, shifted_chr4, category, fusion) #, y_position_offset, curvature)
  
  sc_plot_df <- sc_chimeric_transcripts %>%
    mutate(chr4_position = min_start.y,
           chr14_position = min_start.x) %>%
    mutate(shifted_chr4 = chr4_position + chr4_shift,
           shifted_chr14 = chr14_position + chr14_shift) %>%
    mutate(category = case_when(chr4_position <= chr4_breakpoint & chr14_position <= chr14_breakpoint ~ "Category 1",
                                chr4_position > chr4_breakpoint & chr14_position <= chr14_breakpoint ~ "Category 2",
                                chr4_position <= chr4_breakpoint & chr14_position > chr14_breakpoint ~ "Category 3",
                                chr4_position > chr4_breakpoint & chr14_position > chr14_breakpoint ~ "Category 4")) %>%
    mutate(identifier = str_c(cell_barcode, ":", molecular_barcode),
           data_type = "Single Cell Chimeric Transcript",
           y_position_offset = -0.5,
           curvature = 0.25,
           fusion = NA) %>%
    ungroup() %>%
    select(identifier, data_type, chr14_position, chr4_position, shifted_chr14, shifted_chr4, category, fusion) #, y_position_offset, curvature)
  
  bind_rows(bulk_plot_df, sc_plot_df) %>% return()
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
    filter( (start <= overall_chr4_min & end >= overall_chr4_min) | 
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


#ANALYSIS (all)
if (FALSE) {
  

# ==============================================================================
# ANALYSIS
# ==============================================================================

# ==============================================================================
# Plot cell types of each sample (4)
# UMAP in paper_main; t-SNE in paper_supp
# ==============================================================================
if (TRUE) {
  plot_cell_types(cell_types_27522_1, seurat_object_27522_1, id = "27522_1")
  plot_cell_types(cell_types_27522_4, seurat_object_27522_4, id = "27522_4")
  plot_cell_types(cell_types_56203_1, seurat_object_56203_1, id = "56203_1")
  plot_cell_types(cell_types_56203_2, seurat_object_56203_2, id = "56203_2")
  
  plot_cell_types(cell_types_27522_1, seurat_object_27522_1, id = "27522_1", reduction = "t-SNE", dir = paper_supp)
  plot_cell_types(cell_types_27522_4, seurat_object_27522_4, id = "27522_4", reduction = "t-SNE", dir = paper_supp)
  plot_cell_types(cell_types_56203_1, seurat_object_56203_1, id = "56203_1", reduction = "t-SNE", dir = paper_supp)
  plot_cell_types(cell_types_56203_2, seurat_object_56203_2, id = "56203_2", reduction = "t-SNE", dir = paper_supp)  
}

# ==============================================================================
# Plot interesting gene expressions
# 27522: WHSC1 (ENSG00000109685) and FGFR3 (ENSG00000068078)
# 56203: MYC (ENSG00000136997) and PVT1 (ENSG00000249859)
# ==============================================================================
if (TRUE) {
  plot_gene_expression(seurat_object_27522_1, 
                       get_tsne_umap(cell_types = cell_types_27522_1, 
                                     seurat_object = seurat_object_27522_1), 
                       ensg = "ENSG00000109685", id = "27522_1", gene = "WHSC1")
  plot_gene_expression(seurat_object_27522_1, 
                       get_tsne_umap(cell_types = cell_types_27522_1, 
                                     seurat_object = seurat_object_27522_1), 
                       ensg = "ENSG00000068078", id = "27522_1", gene = "FGFR3")
  plot_gene_expression(seurat_object_27522_4, 
                       get_tsne_umap(cell_types = cell_types_27522_4, 
                                     seurat_object = seurat_object_27522_4), 
                       ensg = "ENSG00000109685", id = "27522_4", gene = "WHSC1")
  plot_gene_expression(seurat_object_27522_4, 
                       get_tsne_umap(cell_types = cell_types_27522_4, 
                                     seurat_object = seurat_object_27522_4), 
                       ensg = "ENSG00000068078", id = "27522_4", gene = "FGFR3")
  
  plot_gene_expression(seurat_object_56203_1, 
                       get_tsne_umap(cell_types = cell_types_56203_1, 
                                     seurat_object = seurat_object_56203_1), 
                       ensg = "ENSG00000136997", id = "56203_1", gene = "MYC")
  plot_gene_expression(seurat_object_56203_1, 
                       get_tsne_umap(cell_types = cell_types_56203_1, 
                                     seurat_object = seurat_object_56203_1), 
                       ensg = "ENSG00000249859", id = "56203_1", gene = "PVT1")
  plot_gene_expression(seurat_object_56203_2, 
                       get_tsne_umap(cell_types = cell_types_56203_2, 
                                     seurat_object = seurat_object_56203_2), 
                       ensg = "ENSG00000136997", id = "56203_2", gene = "MYC")
  plot_gene_expression(seurat_object_56203_2, 
                       get_tsne_umap(cell_types = cell_types_56203_2, 
                                     seurat_object = seurat_object_56203_2), 
                       ensg = "ENSG00000249859", id = "56203_2", gene = "PVT1")
}

# ==============================================================================
# Plot FGFR3 CNV in 25722_1
# Misleading because FGFR3 is upregulated in plasma cells, values dont' make sense
# ==============================================================================
plot_gene_cnv(infercnv_27522_1, 
              get_tsne_umap(cell_types = cell_types_27522_1, 
                            seurat_object = seurat_object_27522_1), 
              id = "27522_1", gene = "FGFR3")

# ==============================================================================
# Plot chr4 CNV in 25722_1
# ==============================================================================
plot_chr_cnv(infercnv_27522_1, 
             get_tsne_umap(cell_types = cell_types_27522_1, 
                           seurat_object = seurat_object_27522_1), 
             id = "27522_1", chr = "chr4")
}

# 27522
if (TRUE) {
  bulk_fusion_reads_27522_1 <- get_bulk_fusion_reads_27522(bulk_reads_27522_1, star_fusion_calls_27522_1)
  plot_sc_read_length_distribution(dis_reads = dis_reads_27522_1_discover, cutoff = 1024, id = "27522_1")
  sc_chimeric_transcripts_27522_1 <- get_sc_chimeric_transcripts_27522(dis_reads_27522_1_discover)
  bulk_sc_plot_df_27522_1 <- get_bulk_sc_plot_df_27522(bulk_fusion_reads_27522_1, sc_chimeric_transcripts_27522_1)
  plot_bulk_sc_27522(bulk_sc_plot_df_27522_1, dir = paper_main, id = "27522_1", bulk_color = FALSE)
  
  plot_cell_chimeric_transcripts(bulk_sc = bulk_sc_plot_df_27522_1, tsne_umap = get_tsne_umap(cell_types = cell_types_27522_1, seurat_object = seurat_object_27522_1), id = "27522_1")
  plot_cell_chimeric_transcripts(bulk_sc = bulk_sc_plot_df_27522_1, tsne_umap = get_tsne_umap(cell_types = cell_types_27522_1, seurat_object = seurat_object_27522_1), id = "27522_1", reduction = "t-SNE", dir = paper_supp)
  
  bulk_fusion_reads_27522_4 <- get_bulk_fusion_reads_27522(bulk_reads_27522_4)
  plot_sc_read_length_distribution(dis_reads = dis_reads_27522_4_discover, cutoff = 1024, id = "27522_4")
  sc_chimeric_transcripts_27522_4 <- get_sc_chimeric_transcripts_27522(dis_reads_27522_4_discover)
  bulk_sc_plot_df_27522_4 <- get_bulk_sc_plot_df_27522(bulk_fusion_reads_27522_4, sc_chimeric_transcripts_27522_4)
  plot_bulk_sc_27522(bulk_sc_plot_df_27522_4, dir = paper_main, id = "27522_4", bulk_color = FALSE)
  
  plot_cell_chimeric_transcripts(bulk_sc = bulk_sc_plot_df_27522_4, tsne_umap = get_tsne_umap(cell_types = cell_types_27522_4, seurat_object = seurat_object_27522_4), id = "27522_4")
  plot_cell_chimeric_transcripts(bulk_sc = bulk_sc_plot_df_27522_4, tsne_umap = get_tsne_umap(cell_types = cell_types_27522_4, seurat_object = seurat_object_27522_4), id = "27522_4", reduction = "t-SNE", dir = paper_supp)
}



  
  


  
# ============================================================================
# Cells with fusion
# ============================================================================
  
  p <- plot_df_breakpoints %>% group_by(cell_barcode) %>% 
    summarize(cats = str_c(sort(unique(category)), collapse = "; ")) %>%
    right_join(umap_coords_27522_1, by = c("cell_barcode" = "barcode")) %>%
    filter(!is.na(UMAP_1)) %>%
    arrange(!is.na(cats)) %>%
    ggplot(aes(x = UMAP_1, y = UMAP_2)) + 
    geom_point(aes(color = !is.na(cats)), shape = 16, size = 1.5) +
    labs(x = "UMAP 1", y = "UMAP 2") +
    theme_bw() +
    coord_equal() +
    #scale_color_brewer(palette = "Purples", drop = FALSE, direction = 1) +
    scale_color_manual(values = c("#a6bddb", "#1c9099")) +
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
  
  ggsave(str_c(paper_main, "chimeric_transcripts.27522_1.umap.pdf"), p,
         width = 3.5, height = 3.5, useDingbats = FALSE)
  ggsave(str_c(paper_main, "chimeric_transcripts.27522_1.umap.no_legend.pdf"), p + guides(color = FALSE),
         width = 3.5, height = 3.5, useDingbats = FALSE)
  
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

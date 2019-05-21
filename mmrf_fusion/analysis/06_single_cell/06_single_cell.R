# ==============================================================================
# scRNA analysis
# Steven Foltz (smfoltz@wustl.edu), May 2019
# ==============================================================================

paper_main = "paper/main/06_single_cell/"
paper_supp = "paper/supplemental/06_single_cell/"

# Create directories 
dir.create(paper_main, recursive = TRUE, showWarnings = FALSE)
dir.create(paper_supp, recursive = TRUE, showWarnings = FALSE)

# 27552_1_fusion
if (TRUE) {
  
  # ==============================================================================
  # Remove redundancy in discordant reads data
  # ==============================================================================
  cell_barcodes_with_fusions <- dis_reads_27522_1_fusion %>%
    mutate(fusion2 = case_when(str_detect(fusion, "NSD2") ~ "IGH--WHSC1",
                               TRUE ~ fusion)) %>%
    group_by(cell_barcode) %>% 
    summarize(fusions = str_c(unique(sort(fusion2)), collapse = ", "))
  
  # ==============================================================================
  # Smartly combine relevant data in single tibble
  # Relevant data: cell barcode, tSNE-coords, cell type,
  # WHSC1/FGFR3/CSNK1E expression, t(4;14) fusion, other fusion
  # ==============================================================================
  
  tsne_coords_27522_1 <- Embeddings(object = seurat_object_27522_1, 
                                    reduction = "tsne") %>% 
    as_tibble() %>% 
    mutate(barcode = row.names(Embeddings(object = seurat_object_27522_1, 
                                               reduction = "tsne"))) %>% 
    select(barcode, tSNE_1, tSNE_2)
  
  umap_coords_27522_1 <- Embeddings(object = seurat_object_27522_1, 
                                    reduction = "umap") %>% 
    as_tibble() %>% 
    mutate(barcode = row.names(Embeddings(object = seurat_object_27522_1, 
                                               reduction = "umap"))) %>% 
    select(barcode, UMAP_1, UMAP_2)
  
  
  key_gene_expression <- tibble(barcode = seurat_object_27522_1@assays$RNA@data@Dimnames[[2]],
                                WHSC1 = seurat_object_27522_1@assays$RNA@data["ENSG00000109685",] %>% as.vector(),
                                FGFR3 = seurat_object_27522_1@assays$RNA@data["ENSG00000068078",] %>% as.vector(),
                                CSNK1E = seurat_object_27522_1@assays$RNA@data["ENSG00000213923",] %>% as.vector())
  
  plot_df <- cell_types_27522_1 %>% 
    left_join(tsne_coords_27522_1, by = "barcode") %>%
    left_join(umap_coords_27522_1, by = "barcode") %>%
    left_join(cell_barcodes_with_fusions, 
              by = c("barcode" = "cell_barcode")) %>%
    replace_na(list(fusions = "None Detected")) %>%
    left_join(key_gene_expression, by = "barcode") %>%
    arrange(cell_type)
  
  cell_type_positions <- tribble(~cell_type,     ~cell_type_long, ~tSNE_1, ~tSNE_2, ~UMAP_1, ~UMAP_2,
                                 "B",                  "B Cells",      -5,      30,       1,      10,
                                 "CD4+T",        "CD4+\nT Cells",     -40,      15,      -6,      -9,
                                 "CD8+T",        "CD8+\nT Cells",     -42,     -20,    -0.5,    -9.5,
                                 "DC",        "Dendritic\nCells",      34,     -20,      12,       5,
                                 "CD16+Mono",      "Macrophages",      NA,      NA,      NA,      NA,
                                 "CD14+Mono",        "Monocytes",      25,      30,      13,      -4,
                                 "NK",   "Natural Killer\nCells",      NA,      NA,      NA,      NA,
                                 "Plasma",       "Plasma\nCells",      18,     -32,    -0.5,      -1)
  
  # ==============================================================================
  # Plot cell types
  # ==============================================================================
  
  ggplot(data = plot_df, aes(x = tSNE_1, y = tSNE_2)) + 
    geom_point(aes(color = cell_type), shape = 16, size = 1.5, show.legend = FALSE, alpha = 0.25) +
    geom_text(data = cell_type_positions, aes(x = tSNE_1, y = tSNE_2, 
                                              label = cell_type_long, 
                                              color = cell_type),
              size = 3.5,
              show.legend = FALSE) +
    annotate("text", label = "Cell Types", x = -Inf, y = Inf, 
             vjust = 1, hjust = 0, 
             size = 3.5) + 
    labs(x = "t-SNE 1", y = "t-SNE 2") +
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
          legend.text = element_text(size = 10)) +
    ggsave(str_c(paper_supp, "cell_types.27522_1.tsne.pdf"),
           width = 3.5, height = 3.5, useDingbats = FALSE)
  
  ggplot(data = plot_df, aes(x = UMAP_1, y = UMAP_2)) + 
    geom_point(aes(color = cell_type), shape = 16, size = 1.5, show.legend = FALSE, alpha = 0.25) +
    geom_text(data = cell_type_positions, aes(x = UMAP_1, y = UMAP_2, 
                                              label = cell_type_long, 
                                              color = cell_type),
              size = 3.5,
              show.legend = FALSE) +
    annotate("text", label = "Cell Types", x = -Inf, y = Inf, 
             vjust = 1, hjust = 0, 
             size = 3.5) + 
    labs(x = "UMAP 1", y = "UMAP 2") +
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
          legend.text = element_text(size = 10)) +
    ggsave(str_c(paper_main, "cell_types.27522_1.umap.pdf"),
           width = 3.5, height = 3.5, useDingbats = FALSE)
  
  # ==============================================================================
  # WHSC1 expression
  # ==============================================================================
  
  ggplot(data = plot_df %>% arrange(WHSC1), aes(x = UMAP_1, y = UMAP_2)) + 
    geom_point(aes(color = WHSC1), shape = 16, size = 1.5, show.legend = FALSE, alpha = 0.25) +
    annotate("text", label = "WHSC1", x = -Inf, y = Inf, 
             vjust = 1, hjust = 0,
             size = 3.5, fontface = "italic") + 
    labs(x = "UMAP 1", y = "UMAP 2") +
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
          legend.text = element_text(size = 10)) +
    ggsave(str_c(paper_supp, "expression.WHSC1.27522_1.pdf"),
           width = 3.5, height = 3.5, useDingbats = FALSE)
  
  # ==============================================================================
  # FGFR3 expression
  # ==============================================================================
  
  ggplot(data = plot_df %>% arrange(FGFR3), aes(x = UMAP_1, y = UMAP_2)) + 
    geom_point(aes(color = FGFR3), shape = 16, size = 1.5, show.legend = FALSE, alpha = 0.25) +
    annotate("text", label = "FGFR3", x = -Inf, y = Inf, 
             vjust = 1, hjust = 0, 
             size = 3.5, fontface = "italic") + 
    labs(x = "UMAP 1", y = "UMAP 2") +
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
          legend.text = element_text(size = 10)) +
    ggsave(str_c(paper_supp, "expression.FGFR3.27522_1.pdf"),
           width = 3.5, height = 3.5, useDingbats = FALSE)
  
  # ==============================================================================
  # CSNK1E expression
  # ==============================================================================
  
  ggplot(data = plot_df %>% arrange(CSNK1E), aes(x = UMAP_1, y = UMAP_2)) + 
    geom_point(aes(color = CSNK1E), shape = 16, size = 1.5, show.legend = FALSE, alpha = 0.25) +
    annotate("text", label = "CSNK1E", x = -Inf, y = Inf, 
             vjust = 1, hjust = 0,
             size = 3.5, fontface = "italic") + 
    labs(x = "UMAP 1", y = "UMAP 2") +
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
          legend.text = element_text(size = 10)) +
    ggsave(str_c(paper_supp, "expression.CSNK1E.27522_1.pdf"),
           width = 3.5, height = 3.5, useDingbats = FALSE)
  
}


# 27552_1_discover
if (TRUE) {
  
  # ============================================================================
  # Smartly combine relevant data in single tibble
  # Relevant data: cell barcode, tSNE-coords, cell type,
  # WHSC1/FGFR3/CSNK1E expression, t(4;14) fusion, other fusion
  # ============================================================================
  
  tsne_coords_27522_1 <- Embeddings(object = seurat_object_27522_1, 
                                    reduction = "tsne") %>% 
    as_tibble() %>% 
    mutate(barcode = row.names(Embeddings(object = seurat_object_27522_1, 
                                          reduction = "tsne"))) %>% 
    select(barcode, tSNE_1, tSNE_2)
  
  umap_coords_27522_1 <- Embeddings(object = seurat_object_27522_1, 
                                    reduction = "umap") %>% 
    as_tibble() %>% 
    mutate(barcode = row.names(Embeddings(object = seurat_object_27522_1, 
                                          reduction = "umap"))) %>% 
    select(barcode, UMAP_1, UMAP_2)
  
  # ============================================================================
  # Transcript breakpoints
  # ============================================================================
  
  reported_translocation_breakpoints <- tribble(~translocation, ~chrom,      ~pos,
                                                     "t(4;14)",      4,   1871964,
                                                     "t(4;14)",     14, 105858090)
  
  chr4_breakpoint <- 1871964
  chr14_breakpoint <- 105858090
  
  reported_fusion_breakpoints <- tribble(      ~fusion,   ~chrom,      ~pos, ~fusion_n,
                                         "IGHJ1--NSD2",  "chr14", 105865407,         1,
                                         "IGHJ1--NSD2",   "chr4",   1900626,         1,
                                         "IGHJ6--NSD2",  "chr14", 105863198,         2,
                                         "IGHJ6--NSD2",   "chr4",   1900626,         2,
                                          "NSD2--IGHM",   "chr4",   1871542,         3,
                                          "NSD2--IGHM",  "chr14", 105856217,         3,
                                         "IGHJ5--NSD2",  "chr14", 105863814,         4,
                                         "IGHJ5--NSD2",   "chr4",   1900626,         4,
                                         "IGHJ3--NSD2",  "chr14", 105864587,         5,
                                         "IGHJ3--NSD2",   "chr4",   1900626,         5,
                                          "IGH@--NSD2",  "chr14", 105864215,         6,
                                          "IGH@--NSD2",   "chr4",   1900626,         6,
                                          "IGH@--NSD2",  "chr14", 105865199,         7,
                                          "IGH@--NSD2",   "chr4",   1900626,         7) %>%
    mutate(label_fusion = str_c(fusion_n, fusion, sep = ": "),
           chrom = str_remove_all(chrom, "chr"))

  chr14_min_max <- dis_reads_27522_1_discover %>% 
    filter(chrom == 14) %>% 
    filter(end - start < 100) %>% 
    group_by(cell_barcode, molecular_barcode) %>% 
    summarize(min_start = min(start), max_start = max(start), 
              min_end = min(end), max_end = max(end), 
              n_reads = n())
  
  chr4_min_max <- dis_reads_27522_1_discover %>% 
    filter(chrom == 4) %>% 
    filter(end - start < 100) %>% 
    group_by(cell_barcode, molecular_barcode) %>% 
    summarize(min_start = min(start), max_start = max(start), 
              min_end = min(end), max_end = max(end), 
              n_reads = n())
  
  plot_df_breakpoints <- chr14_min_max %>% full_join(chr4_min_max, 
                              by = c("cell_barcode", "molecular_barcode")) %>%
    filter(!any(is.na(min_start.x), is.na(min_start.y))) %>%
    mutate(cell_molecular_barcode = str_c(cell_barcode, molecular_barcode, sep = ":"),
           category = case_when(min_start.y > chr4_breakpoint & min_start.x > chr14_breakpoint ~ "Canonical IGH--WHSC1",
                                max_end.y < chr4_breakpoint & max_end.x < chr14_breakpoint ~ "Reciprocal WHSC1--IGH",
                                TRUE ~ "Something else?")) %>%
    mutate(chr14_bp = case_when(category == "Canonical IGH--WHSC1" ~ min_start.x - 105600000,
                                category == "Reciprocal WHSC1--IGH" ~ max_end.x - 105600000,
                                category == "Something else?" ~ min_start.x - 105600000),
           chr4_bp = case_when(category == "Canonical IGH--WHSC1" ~ min_start.y - 1800000 + 2e5,
                               category == "Reciprocal WHSC1--IGH" ~ max_end.y - 1800000 + 2e5,
                               category == "Something else?" ~ min_start.y - 1800000 + 2e5))
  
  ggplot(data = plot_df_breakpoints,
         aes(x = chr14_bp, xend = chr4_bp, y = 14, yend = 4)) +
    scale_y_continuous(breaks = c(4,14), labels = c("chr4", "chr14")) +
    geom_segment(aes(color = category), 
                 show.legend = FALSE) + 
    geom_segment(aes(x = 1871964 - 1800000 + 2e5, xend = 1871964 - 1800000 + 2e5,
                     y = 3.5, yend = 4.5), linetype = 1) +
    geom_segment(aes(x = 105858090 - 105600000, xend = 105858090 - 105600000,
                     y = 13.5, yend = 14.5), linetype = 1) +
    annotate("text", x = 0, y = 14.5, label = "105.6 Mb") +
    annotate("text", x = 200000, y = 14.5, label = "105.8 Mb") +
    annotate("text", x = 400000, y = 14.5, label = "106.0 Mb") +
    annotate("text", x = 600000, y = 14.5, label = "106.2 Mb") +
    annotate("text", x = 200000, y = 3.5, label = "1.8 Mb") +
    annotate("text", x = 400000, y = 3.5, label = "2.0 Mb") +
    facet_wrap(~ category, ncol = 1, strip.position = "right") +
    labs(x = NULL, y = NULL) +
    theme_bw() +
    theme(panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          strip.background = element_blank()) +
    ggsave(str_c(paper_main, "reads.discovered.27522_1.pdf"),
           width = 7.25, height = 5, useDingbats = FALSE)
  
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


# 56203_1_discover
if (FALSE) {
  
  # ==============================================================================
  # Remove redundancy in discordant reads data
  # ==============================================================================
  cell_barcodes_with_fusions <- dis_reads_56203_1_discover %>%
    select(cell_barcode) %>%
    unique() %>%
    mutate(fusions = "IGH--MYC")
  
  # ==============================================================================
  # Smartly combine relevant data in single tibble
  # Relevant data: cell barcode, tSNE-coords, cell type,
  # WHSC1/FGFR3/CSNK1E expression, t(4;14) fusion, other fusion
  # ==============================================================================
  
  key_gene_expression <- tibble(barcode = seurat_object_56203_1@assays$RNA@data@Dimnames[[2]],
                                MYC = seurat_object_56203_1@assays$RNA@data["ENSG00000136997",] %>% as.vector(),
                                PVT1 = seurat_object_56203_1@assays$RNA@data["ENSG00000249859",] %>% as.vector())
  
  plot_df <- cell_types_56203_1 %>% left_join(tsne_coords_56203_1, 
                                              by = c("barcode" = "Barcode")) %>%
    left_join(cell_barcodes_with_fusions, 
              by = c("barcode" = "cell_barcode")) %>%
    replace_na(list(fusions = "None Detected")) %>%
    left_join(key_gene_expression, by = "barcode")
  
  cell_type_positions <- tribble(~cell_type,    ~cell_type_long, ~tSNE_1, ~tSNE_2,
                                 "B",                 "B Cells",       7,      39,
                                 "CD14+Mono", "CD14+ Monocytes",      25,      27,
                                 "CD4+T",        "CD4+ T Cells",     -40,      10,
                                 "CD8+T",        "CD8+ T Cells",     -40,     -15,
                                 "DC",        "Dendritic\nCells",     35,     -15,
                                 "Plasma",       "Plasma Cells",      18,     -30)
  
  # ==============================================================================
  # Plot cell types
  # ==============================================================================
  
  ggplot(data = plot_df, aes(x = tSNE_1, y = tSNE_2)) + 
    geom_point(aes(color = cell_type), shape = 16, size = 3) + #, show.legend = FALSE) +
    #geom_text(data = cell_type_positions, aes(x = tSNE_1, y = tSNE_2, 
    #                                          label = cell_type_long, 
    #                                          color = cell_type),
    #          size = 5,
    #          show.legend = FALSE) +
    labs(x = "t-SNE 1", y = "t-SNE 2") +
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
  
  # ==============================================================================
  # MYC expression
  # ==============================================================================
  
  ggplot(data = plot_df %>% arrange(MYC), aes(x = tSNE_1, y = tSNE_2)) + 
    geom_point(aes(color = MYC), shape = 16, size = 3, show.legend = FALSE) +
    annotate("text", label = "MYC", x = -30, y = 40, size = 5, fontface = "italic") + 
    labs(x = "t-SNE 1", y = "t-SNE 2") +
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
  
  # ==============================================================================
  # PVT1 expression
  # ==============================================================================
 
   ggplot(data = plot_df %>% arrange(PVT1), aes(x = tSNE_1, y = tSNE_2)) + 
    geom_point(aes(color = PVT1), shape = 16, size = 3, show.legend = FALSE) +
    annotate("text", label = "PVT1", x = -30, y = 40, size = 5, fontface = "italic") + 
    labs(x = "t-SNE 1", y = "t-SNE 2") +
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
  
  # ==============================================================================
  # Cells with fusion
  # ==============================================================================
  
  ggplot(data = plot_df %>% arrange(fusions != "None Detected"), 
         aes(x = tSNE_1, y = tSNE_2)) + 
    geom_point(aes(color = fusions != "None Detected"), 
               shape = 16, size = 3) +
    labs(x = "t-SNE 1", y = "t-SNE 2", color = "Fusion Detected") +
    theme_bw() +
    coord_fixed() +
    scale_color_manual(values = c("#9ecae1", "#3182bd")) +
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
  
  # ==============================================================================
  # Transcript breakpoints
  # ==============================================================================
  
  myc_pos <- dis_reads_56203_1_discover %>%
    mutate(diff = end - start) %>% 
    filter(diff < 1000) %>% 
    filter(chrom == 8) %>% 
    group_by(cell_barcode, molecular_barcode) %>% 
    mutate(myc_pos = min(start)) %>% 
    ungroup() %>% 
    select(cell_barcode, molecular_barcode, myc_pos) %>%
    unique()
  
  igh_pos <- dis_reads_56203_1_discover %>%
    mutate(diff = end - start) %>% 
    filter(diff < 1000) %>% 
    filter(chrom == 14) %>% 
    group_by(cell_barcode, molecular_barcode) %>% 
    mutate(igh_pos = min(start)) %>% 
    ungroup() %>% 
    select(cell_barcode, molecular_barcode, igh_pos) %>% 
    unique()
  
  myc_pos %>% 
    left_join(igh_pos, by = c("cell_barcode", "molecular_barcode")) %>% 
    ggplot(aes(x = myc_pos/1e6, y = igh_pos/1e6)) +
    geom_point(shape = 16, size = 3, alpha = 0.5, show.legend = FALSE) +
    labs(x = "Mapping Position (Mb) (MYC, chr8)", 
         y = "Mapping Position (Mb) (IGH, chr14)") +
    theme_bw() +
    coord_equal() +
    scale_x_continuous() + 
    scale_y_continuous() + 
    theme(panel.background = element_blank(),
          panel.border = element_blank(),
          plot.background = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_text(size = 8),
          axis.title = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10))
  
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

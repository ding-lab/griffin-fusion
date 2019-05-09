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
  
  key_gene_expression <- tibble(barcode = seurat_object_27522_1@assays$RNA@data@Dimnames[[2]],
                                WHSC1 = seurat_object_27522_1@assays$RNA@data["ENSG00000109685",] %>% as.vector(),
                                FGFR3 = seurat_object_27522_1@assays$RNA@data["ENSG00000068078",] %>% as.vector(),
                                CSNK1E = seurat_object_27522_1@assays$RNA@data["ENSG00000213923",] %>% as.vector())
  
  plot_df <- cell_types_27522_1 %>% left_join(tsne_coords_27522_1, 
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
                                 "DC",        "Dendritic\nCells",      35,     -15,
                                 "Plasma",       "Plasma Cells",      18,     -30)
  
  # ==============================================================================
  # Plot cell types
  # ==============================================================================
  
  ggplot(data = plot_df, aes(x = tSNE_1, y = tSNE_2)) + 
    geom_point(aes(color = cell_type), shape = 16, size = 3, show.legend = FALSE) +
    geom_text(data = cell_type_positions, aes(x = tSNE_1, y = tSNE_2, 
                                              label = cell_type_long, 
                                              color = cell_type),
              size = 5,
              show.legend = FALSE) +
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
  # WHSC1 expression
  # ==============================================================================
  
  ggplot(data = plot_df %>% arrange(WHSC1), aes(x = tSNE_1, y = tSNE_2)) + 
    geom_point(aes(color = WHSC1), shape = 16, size = 3, show.legend = FALSE) +
    annotate("text", label = "WHSC1", x = -40, y = 40, size = 5, fontface = "italic") + 
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
  # FGFR3 expression
  # ==============================================================================
  
  ggplot(data = plot_df %>% arrange(FGFR3), aes(x = tSNE_1, y = tSNE_2)) + 
    geom_point(aes(color = FGFR3), shape = 16, size = 3, show.legend = FALSE) +
    annotate("text", label = "FGFR3", x = -40, y = 40, size = 5, fontface = "italic") + 
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
  # CSNK1E expression
  # ==============================================================================
  
  ggplot(data = plot_df %>% arrange(CSNK1E), aes(x = tSNE_1, y = tSNE_2)) + 
    geom_point(aes(color = CSNK1E), shape = 16, size = 3, show.legend = FALSE) +
    annotate("text", label = "CSNK1E", x = -40, y = 40, size = 5, fontface = "italic") + 
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
  
  whsc1_pos <- dis_reads_27522_1_fusion %>% filter(str_detect(fusion, "NSD2")) %>% 
    mutate(diff = end - start) %>% 
    filter(diff < 1000) %>% 
    filter(chrom == 4) %>% 
    group_by(cell_barcode, molecular_barcode) %>% 
    mutate(whsc1_pos = min(start)) %>% 
    ungroup() %>% 
    select(cell_barcode, molecular_barcode, whsc1_pos) %>%
    unique()
  
  igh_pos <- dis_reads_27522_1_fusion %>% filter(str_detect(fusion, "NSD2")) %>% 
    mutate(diff = end - start) %>% 
    filter(diff < 1000) %>% 
    filter(chrom == 14) %>% 
    group_by(cell_barcode, molecular_barcode) %>% 
    mutate(igh_pos = min(start)) %>% 
    ungroup() %>% 
    select(cell_barcode, molecular_barcode, igh_pos) %>% 
    unique()
  
  whsc1_pos %>% 
    left_join(igh_pos, by = c("cell_barcode", "molecular_barcode")) %>% 
    ggplot(aes(x = whsc1_pos/1e6, y = igh_pos/1e6)) +
    geom_point(shape = 16, size = 3, alpha = 0.5, show.legend = FALSE) +
    labs(x = "Mapping Position (Mb) (WHSC1, chr4)", 
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
  
  
  # ==============================================================================
  # Connect discordant reads
  # ==============================================================================
  
  left_shift_igh <- igh_pos %>% pull(igh_pos) %>% min() %>% plyr::round_any(10000, floor)
  left_shift_whsc1 <- whsc1_pos %>% pull(whsc1_pos) %>% min() %>% plyr::round_any(10000, floor)

  igh_labels = 
  
  upper_limit_whsc1 <- plot_df %>% pull(new_whsc1_pos) %>% max() %>% plyr::round_any(10000, ceiling)
  whsc1_labels <- seq(whsc1_pos %>% pull(whsc1_pos) %>% min() %>% 
                        plyr::round_any(10000, floor) - 10000, whsc1_pos %>% 
                        pull(whsc1_pos) %>% max() %>% 
                        plyr::round_any(10000, ceiling), 10000)/1e6

  plot_df <- whsc1_pos %>% 
    left_join(igh_pos, by = c("cell_barcode", "molecular_barcode")) %>% 
    mutate(new_igh_pos = igh_pos - left_shift_igh,
           new_whsc1_pos = whsc1_pos - left_shift_whsc1 + 10000)
  ggplot(plot_df) +
    geom_segment(aes(x = new_igh_pos, xend = new_whsc1_pos,
                     y = 2, yend = 1), linetype = 2, color = "grey50") +
    geom_point(aes(x = new_igh_pos, y = 2), shape = 16) +
    geom_point(aes(x = new_whsc1_pos, y = 1), shape = 16) +
    scale_x_continuous(limits = c(0, upper_limit_whsc1),
                       labels = whsc1_labels) +
    theme_bw()
    
  
}


# 27552_1_discover
if (TRUE) {
  
  # ==============================================================================
  # Remove redundancy in discordant reads data
  # ==============================================================================
  cell_barcodes_with_fusions <- dis_reads_27522_1_discover %>%
    select(cell_barcode) %>%
    unique() %>%
    mutate(fusions = "IGH--WHSC1")
  
  # ==============================================================================
  # Smartly combine relevant data in single tibble
  # Relevant data: cell barcode, tSNE-coords, cell type,
  # WHSC1/FGFR3/CSNK1E expression, t(4;14) fusion, other fusion
  # ==============================================================================
  
  key_gene_expression <- tibble(barcode = seurat_object_27522_1@assays$RNA@data@Dimnames[[2]],
                                WHSC1 = seurat_object_27522_1@assays$RNA@data["ENSG00000109685",] %>% as.vector(),
                                FGFR3 = seurat_object_27522_1@assays$RNA@data["ENSG00000068078",] %>% as.vector(),
                                CSNK1E = seurat_object_27522_1@assays$RNA@data["ENSG00000213923",] %>% as.vector())
  
  plot_df <- cell_types_27522_1 %>% left_join(tsne_coords_27522_1, 
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
                                 "DC",        "Dendritic\nCells",      35,     -15,
                                 "Plasma",       "Plasma Cells",      18,     -30)
  
  # ==============================================================================
  # Plot cell types
  # ==============================================================================
  
  ggplot(data = plot_df, aes(x = tSNE_1, y = tSNE_2)) + 
    geom_point(aes(color = cell_type), shape = 16, size = 3, show.legend = FALSE) +
    geom_text(data = cell_type_positions, aes(x = tSNE_1, y = tSNE_2, 
                                              label = cell_type_long, 
                                              color = cell_type),
              size = 5,
              show.legend = FALSE) +
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
  # WHSC1 expression
  # ==============================================================================
  
  ggplot(data = plot_df %>% arrange(WHSC1), aes(x = tSNE_1, y = tSNE_2)) + 
    geom_point(aes(color = WHSC1), shape = 16, size = 3, show.legend = FALSE) +
    annotate("text", label = "WHSC1", x = -40, y = 40, size = 5, fontface = "italic") + 
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
  # FGFR3 expression
  # ==============================================================================
  
  ggplot(data = plot_df %>% arrange(FGFR3), aes(x = tSNE_1, y = tSNE_2)) + 
    geom_point(aes(color = FGFR3), shape = 16, size = 3, show.legend = FALSE) +
    annotate("text", label = "FGFR3", x = -40, y = 40, size = 5, fontface = "italic") + 
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
  # CSNK1E expression
  # ==============================================================================
  
  ggplot(data = plot_df %>% arrange(CSNK1E), aes(x = tSNE_1, y = tSNE_2)) + 
    geom_point(aes(color = CSNK1E), shape = 16, size = 3, show.legend = FALSE) +
    annotate("text", label = "CSNK1E", x = -40, y = 40, size = 5, fontface = "italic") + 
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
  
  whsc1_pos <- dis_reads_27522_1_discover %>%
    mutate(diff = end - start) %>% 
    filter(diff < 1000) %>% 
    filter(chrom == 4) %>% 
    group_by(cell_barcode, molecular_barcode) %>% 
    mutate(whsc1_pos = min(start)) %>% 
    ungroup() %>% 
    select(cell_barcode, molecular_barcode, whsc1_pos) %>%
    unique()
  
  igh_pos <- dis_reads_27522_1_discover %>%
    mutate(diff = end - start) %>% 
    filter(diff < 1000) %>% 
    filter(chrom == 14) %>% 
    group_by(cell_barcode, molecular_barcode) %>% 
    mutate(igh_pos = min(start)) %>% 
    ungroup() %>% 
    select(cell_barcode, molecular_barcode, igh_pos) %>% 
    unique()
  
  whsc1_pos %>% 
    left_join(igh_pos, by = c("cell_barcode", "molecular_barcode")) %>% 
    ggplot(aes(x = whsc1_pos/1e6, y = igh_pos/1e6)) +
    geom_point(shape = 16, size = 3, alpha = 0.5, show.legend = FALSE) +
    labs(x = "Mapping Position (Mb) (WHSC1, chr4)", 
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


# 27552_1_discover
if (TRUE) {
  
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
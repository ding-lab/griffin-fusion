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
# Smartly combine relevant data in single tibble
# Relevant data: cell barcode, tSNE-coords, cell type,
# WHSC1/FGFR3/CSNK1E expression, t(4;14) fusion, other fusion
# ==============================================================================



# ==============================================================================
# Plot cell types
# ==============================================================================

cell_types %>% View()

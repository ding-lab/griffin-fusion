# ==============================================================================
# Reads in data for analysis
# Steven Foltz (smfoltz@wustl.edu), August 2018
# ==============================================================================

# ==============================================================================
# Load all packages at beginning
# ==============================================================================
library(tidyverse)
library(readxl)
source("analysis/common_functions.R")

# Other R pacakges used in downstream scripts
library(pheatmap) # 01_overview.R
library(RColorBrewer) # 01_overview.R
library(UpSetR) # 01_overview.R

# ==============================================================================
# Sample lists (MMRF, SRR) of data related to samples in this analysis
# ==============================================================================
samples_all <- read_tsv("data/sample_list.806.txt", 
                        col_names = c("mmrf","srr"))
samples_primary <- read_tsv("data/sample_list.primary.txt",
                            col_names = c("mmrf","srr"))

# ==============================================================================
# Fusion tibbles. fusions_primary contains fusions from primary time points
# This is BEFORE removing significantly undervalidated fusions
# ==============================================================================
fusions_all <- read_tsv("data/fusion_df.txt") %>%
  mutate(fusion = str_remove(fusion, "@"),
         geneA = str_remove(geneA, "@"),
         geneB = str_remove(geneB, "@"))
fusions_primary <- fusions_all %>% filter(srr %in% samples_primary$srr)

#fusions_hard_all <- read_tsv("data/Hard_Filtered_Fusions.tsv")
#fusions_hard_primary <- fusions_hard_all %>%
#  filter(Sample %in% samples_primary$srr)

# ==============================================================================
# Flag and remove potential false positives (significantly undervalidated)
# This REMOVES significantly undervalidated fusions from fusions data frames
# ==============================================================================
wgs_discordant_read_validation_rate <- 
fusions_primary %>% 
  filter(!is.na(n_discordant), 
         Overlap != "Overlapping_regions") %>% 
  mutate(validated = n_discordant >= 1) %>%
  pull(validated) %>% mean()

significantly_under_validated_fusions <- fusions_primary %>% 
  filter(!is.na(n_discordant), 
         Overlap != "Overlapping_regions") %>% 
  mutate(validated = n_discordant >= 1) %>% 
  group_by(fusion) %>% 
  summarize(n = n(), n_validated = sum(validated)) %>%
  arrange(desc(n)) %>%
  mutate(p_value = pbinom(q = n_validated, size = n, 
                          prob = wgs_discordant_read_validation_rate)) %>% 
  filter(p_value < 0.15) %>% pull(fusion)

fusions_all <- fusions_all %>% 
  filter(!(fusion %in% significantly_under_validated_fusions))
fusions_primary <- fusions_primary %>% 
  filter(!(fusion %in% significantly_under_validated_fusions))

# ==============================================================================
# seqFISH and clinical information
# ==============================================================================
seqfish_clinical_info <- read_tsv("data/seqfish_clinical.txt")

# ==============================================================================
# Gene expression data
# ==============================================================================
# read in expression data but only keep data from samples in analysis
expression_all <- read_tsv("data/mmy_gene_expr_with_fusions.tsv") %>%
  filter(srr %in% samples_all$srr)
# keep primary samples only
expression_primary <- expression_all %>% filter(srr %in% samples_primary$srr)

# ==============================================================================
# Sample output file locations
# ==============================================================================
file_locations <- read_tsv("data/sample_list.with_file_names.txt",
                           col_names = FALSE)

# ==============================================================================
# ENSG gene names
# ==============================================================================
# read in list of ENSGs and gene names used in this study
ensg_gene_list <- read_tsv("data/ensg_gene_list.tsv")

# ==============================================================================
# TCGA Pan-cancer fusion analysis results
# ==============================================================================
pancan_fusions <- read_excel("data/tcga_pancancer_fusions.xlsx",
                             sheet = "Final fusion call set")
names(pancan_fusions) <- pancan_fusions[1,]
pancan_fusions <- pancan_fusions[-1,]

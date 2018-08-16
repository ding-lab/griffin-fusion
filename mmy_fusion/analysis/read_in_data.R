# Read in data for analysis

# Load all pacakges at beginning
library(tidyverse)

# Sample lists (MMRF, SRR) of data related to samples in this analysis
samples_all <- read_tsv("data/sample_list.806.txt", col_names = c("MMRF","SRR"))
samples_primary <- read_tsv("data/sample_list.primary.txt", col_names = c("MMRF","SRR"))

# Fusion tibbles. fusions_primary contains fusions from primary time points
fusions_all <- read_tsv("data/fusion_df.txt")
fusions_primary <- fusions_all %>% filter(srr %in% samples_primary$SRR)

# seqFISH and clinical information
seqfish_clinical_info <- read_tsv("data/seqfish_clinical.txt")

# Gene expression data
expr <- read_tsv("data/mmy_gene_expr_with_fusions.tsv")
expr <- expr %>% filter(srr %in% samples_all$SRR) # only keep data from samples in analysis
primary_expr <- expr %>% filter(srr %in% samples_primary$SRR) # filter out non-primary samples
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
library(ggrepel) # 02_expression.R
library(gridExtra) # 01_overview.R
library(pheatmap) # 01_overview.R
library(RColorBrewer) # 01_overview.R
library(survival) # 01_overview.R
library(survminer) # 01_overview.R
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

fusions_hard_all <- read_tsv("data/Hard_Filtered_Fusions.tsv")
fusions_hard_primary <- fusions_hard_all %>%
  filter(Sample %in% samples_primary$srr)

# ==============================================================================
# Flag and remove potential false positives (significantly undervalidated)
# This REMOVES significantly undervalidated fusions from fusions data frames
# ==============================================================================
wgs_discordant_read_validation_rate <- fusions_primary %>% 
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

undervalidated <- fusions_all %>% 
  filter(fusion %in% significantly_under_validated_fusions) %>%
  select(fusion, geneA, geneB) %>%
  group_by(fusion, geneA, geneB) %>%
  summarize(fusion_count = n()) %>%
  ungroup() %>%
  mutate(filter = "Undervalidated") %>%
  rename("FusionName" = "fusion")

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
# Information about kinases
# ==============================================================================
kinases <- read_tsv("data/Kinase_fusion_info.txt") %>% 
  mutate(Fusion = str_remove(Fusion, "@")) %>% 
  right_join(fusions_primary, by = c("PatientID" = "mmrf", 
                                 "SampleID" = "srr", 
                                 "Fusion" = "fusion")) %>% 
  filter(!is.na(KinaseDomain)) %>%
  mutate(kinase_group_full_name = case_when(Group == "TK" ~ "Tyrosine\nKinase",
                                            Group == "OTHER" ~ "Other",
                                            Group == "TKL" ~ "Tyrosine\nKinase-Like",
                                            TRUE ~ Group))

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

# ==============================================================================
# DEPO database
# ==============================================================================
depo <- read_tsv("data/DEPO_final_20170206.txt")

# ==============================================================================
# Soft filtering
# ==============================================================================
soft_columns <- c("FusionName",	"LeftBreakpoint",	"RightBreakpoint", "Cancer", 
                  "Sample", "JunctionReadCount", "SpanningFragCount", "FFPM", 
                  "PROT_FUSION_TYPE", "GTEx", "Callers", "CallerNumber")
efi <- read_tsv("data/Fusions_EFI.tsv", col_names = soft_columns)
efi <- efi %>% mutate(filter = "EFI")
low_count <- read_tsv("data/Fusions_low_count.tsv", col_names = soft_columns)
low_count <- low_count %>% mutate(filter = "Low Count")
many_partners <- read_tsv("data/Fusions_with_many_partners.tsv", col_names = soft_columns)
many_partners <- many_partners %>% mutate(filter = "Many Partners")
within_300kb <- read_tsv("data/Fusions_within_300kb.tsv", col_names = soft_columns)
within_300kb <- within_300kb %>% mutate(filter = "Within 300Kb")
soft_filtered <- bind_rows(efi, low_count, many_partners, within_300kb)

# ==============================================================================
# Mutation calls from Hua
# ==============================================================================
mutation_calls <- read_tsv("data/wxs_bm_data.withmutect.merged.maf.rc.caller.renamed.Bone_Marrow.tsv",
                           col_types = cols_only(Hugo_Symbol = "c",
                                                 Chromosome	= "n",
                                                 Start_Position	= "n",
                                                 End_Position	= "n",
                                                 Variant_Classification = "c",
                                                 Variant_Type = "c",
                                                 Reference_Allele	= "c",
                                                 Tumor_Seq_Allele1 = "c",
                                                 Tumor_Seq_Allele2 = "c",
                                                 Tumor_Sample_Barcode	= "c",
                                                 Matched_Norm_Sample_Barcode = "c",
                                                 HGVSc	= "c",
                                                 HGVSp	= "c",
                                                 HGVSp_Short = "c",
                                                 Transcript_ID = "c",
                                                 t_depth = "n",
                                                 t_ref_count = "n",
                                                 t_alt_count = "n"))

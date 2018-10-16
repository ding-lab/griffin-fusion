# ==============================================================================
# Sample characteristics (seqFISH heatmaps and clinical variables)
# Steven Foltz (smfoltz@wustl.edu), September 2018
# ==============================================================================

plot_dir = "analysis/fusion_summaries/sample_characteristics/"

# ==============================================================================
# Load necessary libraries
# ==============================================================================

library(pheatmap)
library(RColorBrewer)
library(UpSetR)

# ==============================================================================
# Function to create heatmap of binary seqFISH data
# ==============================================================================

plot_seqfish_heatmap <- function(input_tibble, output_filename){
  
  arranged_info <- input_tibble %>% 
    arrange(desc(seqfish_Hyperdiploidy), 
            desc(seqfish_CN_del_13q14),
            desc(seqfish_CN_del_13q34),
            desc(seqfish_CN_del_17p13), 
            desc(seqfish_CN_gain_1q21),
            desc(seqfish_Translocation_WHSC1_4_14),
            desc(seqfish_Translocation_CCND3_6_14),
            desc(seqfish_Translocation_MYC_8_14), 
            desc(seqfish_Translocation_MAFA_8_14),
            desc(seqfish_Translocation_CCND1_11_14), 
            desc(seqfish_Translocation_CCND2_12_14),
            desc(seqfish_Translocation_MAF_14_16), 
            desc(seqfish_Translocation_MAFB_14_20))
  
  pheatmap_df <- data.frame(
    arranged_info %>%
      select(seqfish_Hyperdiploidy, seqfish_CN_del_13q14, seqfish_CN_del_13q34, 
             seqfish_CN_del_17p13, seqfish_CN_gain_1q21, 
             seqfish_Translocation_WHSC1_4_14, seqfish_Translocation_CCND3_6_14, 
             seqfish_Translocation_MYC_8_14, seqfish_Translocation_MAFA_8_14, 
             seqfish_Translocation_CCND1_11_14, 
             seqfish_Translocation_CCND2_12_14, seqfish_Translocation_MAF_14_16, 
             seqfish_Translocation_MAFB_14_20) %>% 
      mutate_all(funs(replace(., is.na(.), 2))) )
  
  names(pheatmap_df) <- c("Hyperdiploidy", "CNV del(13q14)", "CNV del(13q34)", 
    "CNV_del(17p13)", "CNV gain(1q21)", "Translocation t(4;14) (WHSC1)",
    "Translocation t(6;14) (CCND3)", "Translocation t(8;14) (MYC)", 
    "Translocation t(8;14) (MAFA)", "Translocation t(11;14) (CCND1)", 
    "Translocation t(12;14) (CCND2)", "Translocation t(14;16) (MAF)", 
    "Translocation t(14;20) (MAFB)")
  
  pheatmap_mat <- t(as.matrix(pheatmap_df))
  colnames(pheatmap_mat) <- arranged_info %>% pull(mmrf)
  
  annotation_col_df <- data.frame(
    arranged_info %>%
      select(age_ge_66, Female, race, ISS_Stage) %>%
      mutate_all(factor))
  names(annotation_col_df) <- c("Age >= 66", "Sex", "Ethnicity", "ISS Stage")
  levels(annotation_col_df$Sex) <- c("Male", "Female")
  levels(annotation_col_df$Ethnicity) <- c("White", "Black", "Other")
  levels(annotation_col_df$`Age >= 66`) <- c("No", "Yes")
  rownames(annotation_col_df) <- arranged_info %>% pull(mmrf)
  
  pheatmap(pheatmap_mat, color = brewer.pal(n = 3, name = "Blues"), 
           annotation_col = annotation_col_df, gaps_row = c(5),
           cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = FALSE, 
           cellwidth = 10, cellheight = 10, height = 10, width = 55,  
           legend = FALSE, filename = output_filename)
}

# ==============================================================================
# Function to plot seqFISH as UpSetR
# ==============================================================================

plot_seqfish_upsetr <- function(plot_tibble, output_filename){
  upsetr_df <- plot_tibble %>% 
    select(seqfish_Hyperdiploidy, seqfish_CN_del_13q14, seqfish_CN_del_13q34, 
           seqfish_CN_del_17p13, seqfish_CN_gain_1q21, 
           seqfish_Translocation_WHSC1_4_14, seqfish_Translocation_CCND3_6_14, 
           seqfish_Translocation_MYC_8_14, seqfish_Translocation_MAFA_8_14, 
           seqfish_Translocation_CCND1_11_14, 
           seqfish_Translocation_CCND2_12_14, seqfish_Translocation_MAF_14_16, 
           seqfish_Translocation_MAFB_14_20) %>% 
    rowwise() %>% mutate(none = as.numeric(
      sum(seqfish_Hyperdiploidy, seqfish_CN_del_13q14, seqfish_CN_del_13q34, 
          seqfish_CN_del_17p13, seqfish_CN_gain_1q21, 
          seqfish_Translocation_WHSC1_4_14, seqfish_Translocation_CCND3_6_14, 
          seqfish_Translocation_MYC_8_14, seqfish_Translocation_MAFA_8_14, 
          seqfish_Translocation_CCND1_11_14, seqfish_Translocation_CCND2_12_14, 
          seqfish_Translocation_MAF_14_16, seqfish_Translocation_MAFB_14_20) 
      == 0)) %>% as.data.frame()
  names(upsetr_df) <- c("Hyperdiploidy", "CNV del(13q14)", "CNV del(13q34)", 
                        "CNV_del(17p13)", "CNV gain(1q21)", 
                        "Translocation t(4;14) (WHSC1)",
                        "Translocation t(6;14) (CCND3)", 
                        "Translocation t(8;14) (MYC)", 
                        "Translocation t(8;14) (MAFA)", 
                        "Translocation t(11;14) (CCND1)", 
                        "Translocation t(12;14) (CCND2)", 
                        "Translocation t(14;16) (MAF)", 
                        "Translocation t(14;20) (MAFB)", "No seqFISH events")
  
  pdf(file = output_filename, width = 40, height = 20)
  upset(upsetr_df, nsets = ncol(upsetr_df), nintersects = NA, order.by = "freq",
        text.scale = 2, point.size = 3, line.size = 1)
  dev.off()
}

# ==============================================================================
# Plot two different heatmaps and upsetrs: hyperdiploid, non-hyperdiploid
# ==============================================================================

hyperdiploid_status <- c(1, 0)
hyperdiploid_string <- c("hyperdiploid", "non-hyperdiploid")
for (i in 1:2) {
  h_status <- hyperdiploid_status[i]
  h_string <- hyperdiploid_string[i]
  plot_tibble <- seqfish_clinical_info %>% 
    filter(seqfish_Hyperdiploidy == h_status)
  plot_seqfish_heatmap(plot_tibble, 
                       str_c(plot_dir, "heatmap.", h_string, ".pdf"))
  plot_seqfish_upsetr(plot_tibble, 
                      str_c(plot_dir, "upsetr.", h_string, ".pdf")) 
  
}

# ==============================================================================
# Plot one upsetr for all samples
# ==============================================================================
plot_tibble <- seqfish_clinical_info %>% filter(!is.na(seqfish_Hyperdiploidy))
plot_seqfish_upsetr(plot_tibble, str_c(plot_dir, "upsetr.both.pdf"))

# ==============================================================================
# Number of non/hyperdiploid samples
# ==============================================================================
n_hyperdiploid_samples <- seqfish_clinical_info %>% 
  filter(seqfish_Hyperdiploidy == 1) %>% nrow()
n_nonhyperdiploid_samples <- seqfish_clinical_info %>% 
  filter(seqfish_Hyperdiploidy == 0) %>% nrow()
n_na_hyperdiploid_samples <- seqfish_clinical_info %>% 
  filter(is.na(seqfish_Hyperdiploidy)) %>% nrow()

# ==============================================================================
# Important clinical features
# ==============================================================================
age_summary <- seqfish_clinical_info %>% select(Age) %>% summary()
age_na <- seqfish_clinical_info %>% filter(is.na(Age)) %>% nrow()
sex_summary <- seqfish_clinical_info %>% 
  mutate(sex = factor(Female, labels = c("Male", "Female"))) %>% 
  select(sex) %>% summary()
sex_na <- seqfish_clinical_info %>% filter(is.na(Female)) %>% nrow()
race_summary <- seqfish_clinical_info %>% 
  mutate(race_name = factor(race, labels = c("White", "Black", "Other"))) %>% 
  select(race_name) %>% summary()
race_na <- seqfish_clinical_info %>% filter(is.na(race)) %>% nrow()
ecog_summary <- seqfish_clinical_info %>% mutate_at("ECOG", factor) %>%
  select(ECOG) %>% summary()
ecog_na <- seqfish_clinical_info %>% filter(is.na(ECOG)) %>% nrow()
plasma_summary <- seqfish_clinical_info %>% 
  select(BM_Plasma_Cell_Percent) %>% summary()
plasma_na <- seqfish_clinical_info %>% filter(is.na(BM_Plasma_Cell_Percent)) %>% 
  nrow()
stage_summary <- seqfish_clinical_info %>% mutate_at("ISS_Stage", factor) %>% 
  select(ISS_Stage) %>% summary() 
stage_na <- seqfish_clinical_info %>% filter(is.na(ISS_Stage)) %>% nrow()
ldh_summary <- seqfish_clinical_info %>% 
  select(LDH) %>% summary()
ldh_na <- seqfish_clinical_info %>% filter(is.na(LDH)) %>% nrow()
bone_summary <- seqfish_clinical_info %>% mutate_at("Bone_lesions", factor) %>% 
  select(Bone_lesions) %>% summary()
bone_na <- seqfish_clinical_info %>% filter(is.na(Bone_lesions)) %>% nrow()
plasmacytoma_summary <- seqfish_clinical_info %>% 
  mutate_at("Plasmacytoma", factor) %>% 
  select(Plasmacytoma) %>% summary()
plasmacytoma_na <- seqfish_clinical_info %>% 
  filter(is.na(Plasmacytoma)) %>% nrow()

# ==============================================================================
# Plot continuous clinical variables
# ==============================================================================
pdf(str_c(plot_dir, "clinical_variables.age.pdf"), height = 10, width = 15)
n_missing <- seqfish_clinical_info %>% filter(is.na(Age)) %>% nrow()
p <- seqfish_clinical_info %>% 
  ggplot(aes(x = Age)) + geom_histogram() + 
  labs(x = "Age at onset (years)", y = "Number of patients", 
       caption = str_c("Number missing = ", n_missing)) +
  ggplot2_standard_additions()
print(p)
dev.off()

pdf(str_c(plot_dir, "clinical_variables.bm_pct.pdf"), height = 10, width = 15)
n_missing <- seqfish_clinical_info %>% 
  filter(is.na(BM_Plasma_Cell_Percent)) %>% nrow()
p <- seqfish_clinical_info %>% 
  ggplot(aes(x = BM_Plasma_Cell_Percent)) + geom_histogram() + 
  labs(x = "Bone marrow plasma cell (%)", y = "Number of patients",
       caption = str_c("Number missing = ", n_missing)) +
  ggplot2_standard_additions()
print(p)
dev.off()

pdf(str_c(plot_dir, "clinical_variables.ldh.pdf"), height = 10, width = 15)
n_missing <- seqfish_clinical_info %>% filter(is.na(LDH)) %>% nrow()
p <- seqfish_clinical_info %>% 
  ggplot(aes(x = LDH)) + geom_histogram() +
  labs(x = "Lactate dehydrogenase (LDH) (U/L)", y = "Number of patients",
       caption = str_c("Number missing = ", n_missing)) +
  ggplot2_standard_additions()
print(p)
dev.off()

# ==============================================================================
# Data type summary
# ==============================================================================
n_samples_total <- samples_all %>% nrow()
n_samples_primary <- samples_primary %>% nrow()
n_samples_seqfish <- seqfish_clinical_info %>% 
  filter(!is.na(seqfish_Study_Visit_ID)) %>% nrow()
n_samples_seqfish_na <- seqfish_clinical_info %>% 
  filter(is.na(seqfish_Study_Visit_ID)) %>% nrow()
# X14 is the column of SRR,SRR for WGS tumor,normal
n_samples_wgs <- file_locations %>% filter(X14 != "NA,NA") %>% nrow()
n_samples_wgs_na <- file_locations %>% filter(X14 == "NA,NA") %>% nrow()

# ==============================================================================
# Summary table output
# ==============================================================================
summary_tibble <- tribble(
  ~Category, ~Subcategory, ~N, ~Missing, 
  ~Percentage, ~Min, ~`25%`, ~Median, ~Mean, ~`75%`, ~Max,
  "Hyperdiploid Status", 
  "Hyperdiploid",
  n_hyperdiploid_samples,
  ".", 
  n_hyperdiploid_samples/n_samples_primary,
  ".",  ".",    ".",     ".",   ".",    ".",
  
  "Hyperdiploid Status", 
  "Non-Hyperdiploid",
  n_nonhyperdiploid_samples,
  ".", 
  n_nonhyperdiploid_samples/n_samples_primary,
  ".",  ".",    ".",     ".",   ".",    ".",  
  
  "Hyperdiploid Status", 
  "Not available",
  n_na_hyperdiploid_samples,
  ".", 
  n_na_hyperdiploid_samples/n_samples_primary,
  ".",  ".",    ".",     ".",   ".",    ".",  
  
  "Age",
  ".",
  n_samples_primary - age_na,
  age_na,
  ".",
  as.numeric(str_split(age_summary[1,1], ":", simplify = TRUE)[1,2]),
  as.numeric(str_split(age_summary[2,1], ":", simplify = TRUE)[1,2]),
  as.numeric(str_split(age_summary[3,1], ":", simplify = TRUE)[1,2]),
  as.numeric(str_split(age_summary[4,1], ":", simplify = TRUE)[1,2]),
  as.numeric(str_split(age_summary[5,1], ":", simplify = TRUE)[1,2]),
  as.numeric(str_split(age_summary[6,1], ":", simplify = TRUE)[1,2]),
  
  "Bone marrow plasma cell (%)",
  ".",
  n_samples_primary - plasma_na,
  plasma_na,
  ".",
  as.numeric(str_split(plasma_summary[1,1], ":", simplify = TRUE)[1,2]),
  as.numeric(str_split(plasma_summary[2,1], ":", simplify = TRUE)[1,2]),
  as.numeric(str_split(plasma_summary[3,1], ":", simplify = TRUE)[1,2]),
  as.numeric(str_split(plasma_summary[4,1], ":", simplify = TRUE)[1,2]),
  as.numeric(str_split(plasma_summary[5,1], ":", simplify = TRUE)[1,2]),
  as.numeric(str_split(plasma_summary[6,1], ":", simplify = TRUE)[1,2]),
  
  "LDH",
  ".",
  n_samples_primary - ldh_na,
  ldh_na,
  ".",
  as.numeric(str_split(ldh_summary[1,1], ":", simplify = TRUE)[1,2]),
  as.numeric(str_split(ldh_summary[2,1], ":", simplify = TRUE)[1,2]),
  as.numeric(str_split(ldh_summary[3,1], ":", simplify = TRUE)[1,2]),
  as.numeric(str_split(ldh_summary[4,1], ":", simplify = TRUE)[1,2]),
  as.numeric(str_split(ldh_summary[5,1], ":", simplify = TRUE)[1,2]),
  as.numeric(str_split(ldh_summary[6,1], ":", simplify = TRUE)[1,2]),
  
  "Sex",
  "Female",
  as.numeric(str_split(sex_summary[2,1], ":", simplify = TRUE)[1,2]),
  ".",
  as.numeric(str_split(sex_summary[2,1], ":", 
                       simplify = TRUE)[1,2])/n_samples_primary,
  ".", ".", ".", ".", ".", ".",
  
  "Sex",
  "Male",
  as.numeric(str_split(sex_summary[1,1], ":", simplify = TRUE)[1,2]),
  ".",
  as.numeric(str_split(sex_summary[1,1], ":", 
                       simplify = TRUE)[1,2])/n_samples_primary,
  ".", ".", ".", ".", ".", ".",
  
  "Race",
  "White",
  as.numeric(str_split(race_summary[1,1], ":", simplify = TRUE)[1,2]),
  ".",
  as.numeric(str_split(race_summary[1,1], ":", 
                       simplify = TRUE)[1,2])/n_samples_primary,
  ".", ".", ".", ".", ".", ".",
  
  "Race",
  "Black",
  as.numeric(str_split(race_summary[2,1], ":", simplify = TRUE)[1,2]),
  ".",
  as.numeric(str_split(race_summary[2,1], ":", 
                       simplify = TRUE)[1,2])/n_samples_primary,
  ".", ".", ".", ".", ".", ".",
  
  "Race",
  "Other",
  as.numeric(str_split(race_summary[3,1], ":", simplify = TRUE)[1,2]),
  ".",
  as.numeric(str_split(race_summary[3,1], ":", 
                       simplify = TRUE)[1,2])/n_samples_primary,
  ".", ".", ".", ".", ".", ".",
  
  "ECOG",
  "0 = Fully active",
  as.numeric(str_split(ecog_summary[1,1], ":", simplify = TRUE)[1,2]),
  ".",
  as.numeric(str_split(ecog_summary[1,1], ":", 
                       simplify = TRUE)[1,2])/n_samples_primary,
  ".", ".", ".", ".", ".", ".",
  
  "ECOG",
  "1 = Restricted in physically strenuous activity",
  as.numeric(str_split(ecog_summary[2,1], ":", simplify = TRUE)[1,2]),
  ".",
  as.numeric(str_split(ecog_summary[2,1], ":", 
                       simplify = TRUE)[1,2])/n_samples_primary,
  ".", ".", ".", ".", ".", ".",
  
  "ECOG",
  "2 = Ambulatory and capable of all self-care",
  as.numeric(str_split(ecog_summary[3,1], ":", simplify = TRUE)[1,2]),
  ".",
  as.numeric(str_split(ecog_summary[3,1], ":", 
                       simplify = TRUE)[1,2])/n_samples_primary,
  ".", ".", ".", ".", ".", ".",
  
  "ECOG",
  "3 = Capable of only limited self-care",
  as.numeric(str_split(ecog_summary[4,1], ":", simplify = TRUE)[1,2]),
  ".",
  as.numeric(str_split(ecog_summary[4,1], ":", 
                       simplify = TRUE)[1,2])/n_samples_primary,
  ".", ".", ".", ".", ".", ".",
  
  "ECOG",
  "4 = Completely disabled",
  as.numeric(str_split(ecog_summary[5,1], ":", simplify = TRUE)[1,2]),
  ".",
  as.numeric(str_split(ecog_summary[5,1], ":", 
                       simplify = TRUE)[1,2])/n_samples_primary,
  ".", ".", ".", ".", ".", ".",
  
  "ECOG",
  "Not available",
  as.numeric(str_split(ecog_summary[6,1], ":", simplify = TRUE)[1,2]),
  ".",
  as.numeric(str_split(ecog_summary[6,1], ":", 
                       simplify = TRUE)[1,2])/n_samples_primary,
  ".", ".", ".", ".", ".", ".", 
  
  "ISS Stage",
  "I",
  as.numeric(str_split(stage_summary[1,1], ":", simplify = TRUE)[1,2]),
  ".",
  as.numeric(str_split(stage_summary[1,1], ":", 
                       simplify = TRUE)[1,2])/n_samples_primary,
  ".", ".", ".", ".", ".", ".", 
  
  "ISS Stage",
  "II",
  as.numeric(str_split(stage_summary[2,1], ":", simplify = TRUE)[1,2]),
  ".",
  as.numeric(str_split(stage_summary[2,1], ":", 
                       simplify = TRUE)[1,2])/n_samples_primary,
  ".", ".", ".", ".", ".", ".", 
  
  "ISS Stage",
  "III",
  as.numeric(str_split(stage_summary[3,1], ":", simplify = TRUE)[1,2]),
  ".",
  as.numeric(str_split(stage_summary[3,1], ":", 
                       simplify = TRUE)[1,2])/n_samples_primary,
  ".", ".", ".", ".", ".", ".", 
  
  "ISS Stage",
  "Not available",
  as.numeric(str_split(stage_summary[4,1], ":", simplify = TRUE)[1,2]),
  ".",
  as.numeric(str_split(stage_summary[4,1], ":", 
                       simplify = TRUE)[1,2])/n_samples_primary,
  ".", ".", ".", ".", ".", ".",
  
  "Bone lesions",
  "No",
  as.numeric(str_split(bone_summary[1,1], ":", simplify = TRUE)[1,2]),
  ".",
  as.numeric(str_split(bone_summary[1,1], ":", 
                       simplify = TRUE)[1,2])/n_samples_primary,
  ".", ".", ".", ".", ".", ".",
  
  "Bone lesions",
  "Yes",
  as.numeric(str_split(bone_summary[2,1], ":", simplify = TRUE)[1,2]),
  ".",
  as.numeric(str_split(bone_summary[2,1], ":", 
                       simplify = TRUE)[1,2])/n_samples_primary,
  ".", ".", ".", ".", ".", ".",
  
  "Plasmacytoma",
  "No",
  as.numeric(str_split(plasmacytoma_summary[1,1], ":", simplify = TRUE)[1,2]),
  ".",
  as.numeric(str_split(plasmacytoma_summary[1,1], ":", 
                       simplify = TRUE)[1,2])/n_samples_primary,
  ".", ".", ".", ".", ".", ".",
  
  "Plasmacytoma",
  "Yes",
  as.numeric(str_split(bone_summary[2,1], ":", simplify = TRUE)[1,2]),
  ".",
  as.numeric(str_split(bone_summary[2,1], ":", 
                       simplify = TRUE)[1,2])/n_samples_primary,
  ".", ".", ".", ".", ".", ".", 
  
  "Number of samples",
  "Total (includes multiple time points) (with RNA-seq)",
  n_samples_total, 
  ".",
  ".",
  ".", ".", ".", ".", ".", ".",
  
  "Number of samples",
  "Primary (with RNA-seq)",
  n_samples_primary,
  ".",
  n_samples_primary/n_samples_primary,
  ".", ".", ".", ".", ".", ".",
  
  "Number of samples",
  "Primary (with RNA-seq + seqFISH)",
  n_samples_seqfish,
  ".",
  n_samples_seqfish/n_samples_primary,
  ".", ".", ".", ".", ".", ".",
  
  "Number of samples",
  "Primary (with RNA-seq + seqFISH missing)",
  n_samples_seqfish_na,
  ".",
  n_samples_seqfish_na/n_samples_primary,
  ".", ".", ".", ".", ".", ".",
  
  "Number of samples",
  "Primary (with RNA-seq + WGS)",
  n_samples_wgs,
  ".",
  n_samples_wgs/n_samples_primary,
  ".", ".", ".", ".", ".", ".",
  
  "Number of samples",
  "Primary (with RNA-seq + WGS missing)",
  n_samples_wgs_na,
  ".",
  n_samples_wgs_na/n_samples_primary,
  ".", ".", ".", ".", ".", "."
  )

summary_tibble <- summary_tibble %>% arrange(Category, Subcategory)
write_tsv(summary_tibble, 
          "analysis/fusion_summaries/sample_characteristics/summary_table.txt", 
          na = "NA", append = FALSE, col_names = TRUE)
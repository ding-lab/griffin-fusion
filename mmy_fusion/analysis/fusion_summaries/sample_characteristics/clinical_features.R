# ==============================================================================
# Sample characteristics (seqFISH heatmaps and clinical variables)
# Steven Foltz (smfoltz@wustl.edu), September 2018
# ==============================================================================

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
    mutate(seqfish_Nonhyperdiploidy = !as.logical(seqfish_Hyperdiploidy)) %>%
    rowwise() %>% mutate(none = as.numeric(
      sum(seqfish_Hyperdiploidy, seqfish_CN_del_13q14, seqfish_CN_del_13q34, 
          seqfish_CN_del_17p13, seqfish_CN_gain_1q21, 
          seqfish_Translocation_WHSC1_4_14, seqfish_Translocation_CCND3_6_14, 
          seqfish_Translocation_MYC_8_14, seqfish_Translocation_MAFA_8_14, 
          seqfish_Translocation_CCND1_11_14, seqfish_Translocation_CCND2_12_14, 
          seqfish_Translocation_MAF_14_16, seqfish_Translocation_MAFB_14_20) 
      == 0)) %>% as.data.frame()
  names(upsetr_df) <- c("Hyperdiploidy", "CNV del(13q14)", "CNV del(13q34)", 
                        "CNV_del(17p13)", "CNV gain(1q21)", "Translocation t(4;14) (WHSC1)",
                        "Translocation t(6;14) (CCND3)", "Translocation t(8;14) (MYC)", 
                        "Translocation t(8;14) (MAFA)", "Translocation t(11;14) (CCND1)", 
                        "Translocation t(12;14) (CCND2)", "Translocation t(14;16) (MAF)", 
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
  plot_tibble <- seqfish_clinical_info %>% filter(seqfish_Hyperdiploidy == h_status)
  plot_seqfish_heatmap(plot_tibble, 
              str_c("analysis/fusion_summaries/sample_characteristics/heatmap.", 
                   h_string, ".pdf"))
  plot_seqfish_upsetr(plot_tibble, 
              str_c("analysis/fusion_summaries/sample_characteristics/upsetr.", 
                   h_string, ".pdf")) 
  
}

# ==============================================================================
# Plot one upsetr for all samples
# ==============================================================================
plot_tibble <- seqfish_clinical_info %>% filter(!is.na(seqfish_Hyperdiploidy))
plot_seqfish_upsetr(plot_tibble, 
             "analysis/fusion_summaries/sample_characteristics/upsetr.both.pdf")

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
sex_summary <- seqfish_clinical_info %>% 
  mutate(sex = factor(Female, labels = c("Male", "Female"))) %>% 
  select(sex) %>% summary()
race_summary <- seqfish_clinical_info %>% 
  mutate(race_name = factor(race, labels = c("White", "Black", "Other"))) %>% 
  select(race_name) %>% summary()
ecog_summary <- seqfish_clinical_info %>% mutate_at("ECOG", factor) %>%
  select(ECOG) %>% summary()
plasma_summary <- seqfish_clinical_info %>% 
  select(BM_Plasma_Cell_Percent) %>% summary()
stage_summary <- seqfish_clinical_info %>% mutate_at("ISS_Stage", factor) %>% 
  select(ISS_Stage) %>% summary() 
ldh_summary <- seqfish_clinical_info %>% 
  select(LDH) %>% summary()
bone_summary <- seqfish_clinical_info %>% mutate_at("Bone_lesions", factor) %>% 
  select(Bone_lesions) %>% summary()
plasmacytoma_summary <- seqfish_clinical_info %>% select(Plamacytoma) %>% 
  summary()

# ==============================================================================
# Plot continuous clinical variables
# ==============================================================================
pdf("analysis/fusion_summaries/sample_characteristics/clinical_variables.pdf", 
    height = 10, width = 15)

n_missing <- seqfish_clinical_info %>% filter(is.na(Age)) %>% nrow()
seqfish_clinical_info %>% 
  ggplot(aes(x = Age)) + geom_histogram() + 
  labs(x = "Age at onset (years)", y = "Number of patients", 
       caption = str_c("Number missing = ", n_missing)) +
  ggplot2_standard_additions()

n_missing <- seqfish_clinical_info %>% 
  filter(is.na(BM_Plasma_Cell_Percent)) %>% nrow()
seqfish_clinical_info %>% 
  ggplot(aes(x = BM_Plasma_Cell_Percent)) + geom_histogram() + 
  labs(x = "Bone marrow plasma cell (%)", y = "Number of patients",
       caption = str_c("Number missing = ", n_missing)) +
  ggplot2_standard_additions()

n_missing <- seqfish_clinical_info %>% filter(is.na(LDH)) %>% nrow()
seqfish_clinical_info %>% 
  ggplot(aes(x = LDH)) + geom_histogram() + 
  labs(x = "Lactate dehydrogenase (LDH) (U/L)", y = "Number of patients",
       caption = str_c("Number missing = ", n_missing)) +
  ggplot2_standard_additions()
dev.off()

# ==============================================================================
# Sample characteristics (seqFISH heatmaps and clinical variables)
# Steven Foltz (smfoltz@wustl.edu), September 2018
# ==============================================================================

# ==============================================================================
# Load necessary libraries
# ==============================================================================

library(pheatmap)

# ==============================================================================
# Heatmap of binary seqFISH data
# ==============================================================================

arranged_info <- seqfish_clinical_info %>% 
  filter(seqfish_Hyperdiploidy == 1) %>%
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

pheatmap_mat <- t(as.matrix(data.frame(
  arranged_info %>%
    select(seqfish_Hyperdiploidy, seqfish_CN_del_13q14, seqfish_CN_del_13q34, 
seqfish_CN_del_17p13, seqfish_CN_gain_1q21, seqfish_Translocation_WHSC1_4_14,
seqfish_Translocation_CCND3_6_14, seqfish_Translocation_MYC_8_14, 
seqfish_Translocation_MAFA_8_14, seqfish_Translocation_CCND1_11_14, 
seqfish_Translocation_CCND2_12_14, seqfish_Translocation_MAF_14_16, 
seqfish_Translocation_MAFB_14_20)) %>%
  mutate_all(funs(replace(., is.na(.), 2)))))

colnames(pheatmap_mat) <- arranged_info %>% pull(mmrf)

annotation_col_df <- data.frame(
  arranged_info %>%
    select(seqfish_Hyperdiploidy, age_ge_66, Female, race, ISS_Stage) %>%
  mutate_all(factor))
rownames(annotation_col_df) <- arranged_info %>% pull(mmrf)

pheatmap(pheatmap_mat, color = brewer.pal(n = 3, name = "Blues"), annotation_col = annotation_col_df, cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = FALSE, filename = "~/Desktop/test2.pdf", cellwidth = 10, cellheight = 10, legend = FALSE, heigh = 10)

#cellwidth = 1, cellheight = 1, legend = FALSE

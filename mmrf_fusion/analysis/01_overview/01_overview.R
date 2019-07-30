# ==============================================================================
# Overview (MMRF Fusions)
# Steven Foltz (smfoltz@wustl.edu), April 2019
# ==============================================================================

paper_main = "paper/main/01_overview/"
paper_supp = "paper/supplemental/01_overview/"

# Create directories 
dir.create(paper_main, recursive = TRUE, showWarnings = FALSE)
dir.create(paper_supp, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# Clinical features of MMRF patients -- Supplemental
# Originally written September 2018, Updated April 2019
# ==============================================================================

if (TRUE) {
  
  # ============================================================================
  # Function to create heatmap of binary seqFISH data
  # ============================================================================
  
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
        select(seqfish_Hyperdiploidy, 
               seqfish_CN_del_13q14, 
               seqfish_CN_del_13q34, 
               seqfish_CN_del_17p13, 
               seqfish_CN_gain_1q21, 
               seqfish_Translocation_WHSC1_4_14,
               seqfish_Translocation_CCND3_6_14, 
               seqfish_Translocation_MYC_8_14, 
               seqfish_Translocation_MAFA_8_14, 
               seqfish_Translocation_CCND1_11_14, 
               seqfish_Translocation_CCND2_12_14, 
               seqfish_Translocation_MAF_14_16, 
               seqfish_Translocation_MAFB_14_20) %>% 
        mutate_all(list(~replace(., is.na(.), 2))) )
    
    names(pheatmap_df) <- c("Hyperdiploidy", 
                            "CNV del(13q14)", 
                            "CNV del(13q34)", 
                            "CNV_del(17p13)", 
                            "CNV gain(1q21)", 
                            "Translocation t(4;14) (WHSC1)",
                            "Translocation t(6;14) (CCND3)", 
                            "Translocation t(8;14) (MYC)", 
                            "Translocation t(8;14) (MAFA)", 
                            "Translocation t(11;14) (CCND1)", 
                            "Translocation t(12;14) (CCND2)", 
                            "Translocation t(14;16) (MAF)", 
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
  
  # ============================================================================
  # Function to plot seqFISH as UpSetR
  # ==============================================================================
  
  plot_seqfish_upsetr <- function(plot_tibble, output_filename){
    upsetr_df <- plot_tibble %>% 
      select(seqfish_Hyperdiploidy, 
             seqfish_CN_del_13q14, 
             seqfish_CN_del_13q34, 
             seqfish_CN_del_17p13, 
             seqfish_CN_gain_1q21, 
             seqfish_Translocation_WHSC1_4_14, 
             seqfish_Translocation_CCND3_6_14, 
             seqfish_Translocation_MYC_8_14, 
             seqfish_Translocation_MAFA_8_14, 
             seqfish_Translocation_CCND1_11_14, 
             seqfish_Translocation_CCND2_12_14,
             seqfish_Translocation_MAF_14_16, 
             seqfish_Translocation_MAFB_14_20) %>% 
      rowwise() %>% mutate(none = as.numeric(
        sum(seqfish_Hyperdiploidy, 
            seqfish_CN_del_13q14, 
            seqfish_CN_del_13q34, 
            seqfish_CN_del_17p13, 
            seqfish_CN_gain_1q21, 
            seqfish_Translocation_WHSC1_4_14,
            seqfish_Translocation_CCND3_6_14, 
            seqfish_Translocation_MYC_8_14, 
            seqfish_Translocation_MAFA_8_14, 
            seqfish_Translocation_CCND1_11_14,
            seqfish_Translocation_CCND2_12_14, 
            seqfish_Translocation_MAF_14_16, 
            seqfish_Translocation_MAFB_14_20) 
        == 0)) %>% as.data.frame()
    names(upsetr_df) <- c("Hyperdiploidy", 
                          "CNV del(13q14)", 
                          "CNV del(13q34)", 
                          "CNV_del(17p13)", 
                          "CNV gain(1q21)", 
                          "Translocation t(4;14) (WHSC1)",
                          "Translocation t(6;14) (CCND3)", 
                          "Translocation t(8;14) (MYC)", 
                          "Translocation t(8;14) (MAFA)", 
                          "Translocation t(11;14) (CCND1)", 
                          "Translocation t(12;14) (CCND2)", 
                          "Translocation t(14;16) (MAF)", 
                          "Translocation t(14;20) (MAFB)", 
                          "No seqFISH events")
    
    pdf(file = output_filename, width = 40, height = 20, useDingbats = FALSE)
    upset(upsetr_df, 
          nsets = ncol(upsetr_df), 
          nintersects = NA, 
          order.by = "freq",
          #text.scale = c(intersection size title, 
          #               intersection size tick labels, 
          #               set size title, 
          #               set size tick labels, 
          #               set names, 
          #               numbers above bars),
          text.scale = c(4, 
                         3, 
                         4, 
                         3, 
                         3, 
                         3),
          point.size = 4, 
          line.size = 2)
    dev.off()
  }
  
  # ============================================================================
  # Plot two different heatmaps and upsetrs: hyperdiploid, non-hyperdiploid
  # ============================================================================
  
  hyperdiploid_status <- c(1, 0)
  hyperdiploid_string <- c("hyperdiploid", "non-hyperdiploid")
  for (i in 1:2) {
    h_status <- hyperdiploid_status[i]
    h_string <- hyperdiploid_string[i]
    plot_tibble <- seqfish_clinical_info %>% 
      filter(seqfish_Hyperdiploidy == h_status)
    plot_seqfish_heatmap(plot_tibble, 
                         str_c(paper_supp, "heatmap.", h_string, ".pdf"))
    plot_seqfish_upsetr(plot_tibble, 
                        str_c(paper_supp, "upsetr.", h_string, ".pdf")) 
    
  }
  
  # ============================================================================
  # Plot one upsetr for all samples
  # ============================================================================
  
  plot_tibble <- seqfish_clinical_info %>% filter(!is.na(seqfish_Hyperdiploidy))
  plot_seqfish_upsetr(plot_tibble, str_c(paper_supp, "upsetr.both.pdf"))
  
  # ============================================================================
  # Number of non/hyperdiploid samples
  # ============================================================================
  
  n_hyperdiploid_samples <- seqfish_clinical_info %>% 
    filter(seqfish_Hyperdiploidy == 1) %>% nrow()
  n_nonhyperdiploid_samples <- seqfish_clinical_info %>% 
    filter(seqfish_Hyperdiploidy == 0) %>% nrow()
  n_na_hyperdiploid_samples <- seqfish_clinical_info %>% 
    filter(is.na(seqfish_Hyperdiploidy)) %>% nrow()
  
  # ============================================================================
  # Important clinical features
  # ============================================================================
  
  age_summary <- seqfish_clinical_info %>% 
    select(Age) %>% summary()
  age_na <- seqfish_clinical_info %>% 
    filter(is.na(Age)) %>% nrow()
  sex_summary <- seqfish_clinical_info %>% 
    mutate(sex = factor(Female, labels = c("Male", "Female"))) %>% 
    select(sex) %>% summary()
  sex_na <- seqfish_clinical_info %>% 
    filter(is.na(Female)) %>% nrow()
  race_summary <- seqfish_clinical_info %>% 
    mutate(race_name = factor(race, labels = c("White", "Black", "Other"))) %>% 
    select(race_name) %>% summary()
  race_na <- seqfish_clinical_info %>%
    filter(is.na(race)) %>% nrow()
  ecog_summary <- seqfish_clinical_info %>% 
    mutate_at("ECOG", factor) %>% select(ECOG) %>% summary()
  ecog_na <- seqfish_clinical_info %>% 
    filter(is.na(ECOG)) %>% nrow()
  plasma_summary <- seqfish_clinical_info %>% 
    select(BM_Plasma_Cell_Percent) %>% summary()
  plasma_na <- seqfish_clinical_info %>% 
    filter(is.na(BM_Plasma_Cell_Percent)) %>% nrow()
  stage_summary <- seqfish_clinical_info %>%
    mutate_at("ISS_Stage", factor) %>% select(ISS_Stage) %>% summary() 
  stage_na <- seqfish_clinical_info %>% 
    filter(is.na(ISS_Stage)) %>% nrow()
  ldh_summary <- seqfish_clinical_info %>% 
    select(LDH) %>% summary()
  ldh_na <- seqfish_clinical_info %>% 
    filter(is.na(LDH)) %>% nrow()
  bone_summary <- seqfish_clinical_info %>% 
    mutate_at("Bone_lesions", factor) %>% select(Bone_lesions) %>% summary()
  bone_na <- seqfish_clinical_info %>%
    filter(is.na(Bone_lesions)) %>% nrow()
  plasmacytoma_summary <- seqfish_clinical_info %>% 
    mutate_at("Plasmacytoma", factor) %>% select(Plasmacytoma) %>% summary()
  plasmacytoma_na <- seqfish_clinical_info %>% 
    filter(is.na(Plasmacytoma)) %>% nrow()
  
  # ============================================================================
  # Plot continuous clinical variables
  # ============================================================================
  
  pdf(str_c(paper_supp, "clinical_variables.age.pdf"), height = 10, width = 15,
      useDingbats = FALSE)
  n_missing <- seqfish_clinical_info %>% filter(is.na(Age)) %>% nrow()
  p <- seqfish_clinical_info %>% 
    ggplot(aes(x = Age)) + geom_histogram() + 
    labs(x = "Age at onset (years)", y = "Number of patients", 
         caption = str_c("Number missing = ", n_missing)) +
    theme_bw(base_size = 20)
  print(p)
  dev.off()
  
  pdf(str_c(paper_supp, "clinical_variables.bm_pct.pdf"), height = 10, 
      width = 15, useDingbats = FALSE)
  n_missing <- seqfish_clinical_info %>% 
    filter(is.na(BM_Plasma_Cell_Percent)) %>% nrow()
  p <- seqfish_clinical_info %>% 
    ggplot(aes(x = BM_Plasma_Cell_Percent)) + geom_histogram() + 
    labs(x = "Bone marrow plasma cell (%)", y = "Number of patients",
         caption = str_c("Number missing = ", n_missing)) +
    theme_bw(base_size = 20)
  print(p)
  dev.off()
  
  pdf(str_c(paper_supp, "clinical_variables.ldh.pdf"), height = 10, width = 15,
      useDingbats = FALSE)
  n_missing <- seqfish_clinical_info %>% filter(is.na(LDH)) %>% nrow()
  p <- seqfish_clinical_info %>% 
    ggplot(aes(x = LDH)) + geom_histogram() +
    labs(x = "Lactate dehydrogenase (LDH) (U/L)", y = "Number of patients",
         caption = str_c("Number missing = ", n_missing)) +
    theme_bw(base_size = 20)
  print(p)
  dev.off()
  
  # ============================================================================
  # Data type summary
  # ============================================================================
  n_samples_total <- samples_all %>% nrow()
  n_samples_primary <- samples_primary %>% nrow()
  n_samples_seqfish <- seqfish_clinical_info %>% 
    filter(!is.na(seqfish_Study_Visit_ID)) %>% nrow()
  n_samples_seqfish_na <- seqfish_clinical_info %>% 
    filter(is.na(seqfish_Study_Visit_ID)) %>% nrow()
  # X14 is the column of SRR,SRR for WGS tumor,normal
  n_samples_wgs <- file_locations %>% filter(X14 != "NA,NA") %>% nrow()
  n_samples_wgs_na <- file_locations %>% filter(X14 == "NA,NA") %>% nrow()
  
  print(str_c("Number of patients: ", n_samples_primary))
  print(str_c("Number with WGS: ", n_samples_wgs, "/", n_samples_primary, " = ", round(100*n_samples_wgs/n_samples_primary, digits = 1), "%"))
  print(str_c("Number with additional samples: ", samples_all %>% 
                group_by(mmrf) %>% 
                summarize(count = n()) %>% 
                filter(count > 1) %>% 
                nrow()))
  print(str_c("Number of total RNA samples: ", n_samples_total))
  
  # ============================================================================
  # Summary table output
  # ============================================================================
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
    as.numeric(str_split(plasmacytoma_summary[2,1], ":", simplify = TRUE)[1,2]),
    ".",
    as.numeric(str_split(plasmacytoma_summary[2,1], ":", 
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
  write_tsv(summary_tibble, str_c(paper_supp, "summary_table.txt"), 
            na = "NA", append = FALSE, col_names = TRUE)
}

# ==============================================================================
# Basic survival plots of MMRF patients -- Supplemental
# Originally written September 2018, Updated April 2019
# ==============================================================================

if (TRUE) {
  
  # ============================================================================
  # Create EFS survival object
  # ============================================================================
  
  # Use EFS_censor == 0 because that is TRUE for death, FALSE for censored
  # Stratify by Stage
  EFS_tibble <- seqfish_clinical_info %>% 
    filter(!is.na(ISS_Stage), !is.na(EFS_censor))
  EFS_fit <- survfit(Surv(EFS, EFS_censor == 0) ~ ISS_Stage, 
                     data = EFS_tibble)
  
  # Plot survival curve stratified by Stage
  pdf(str_c(paper_supp, "event_free_survival.pdf"), width = 7.25, height = 5,
      useDingbats = FALSE)
  print(ggsurvplot(EFS_fit, data = EFS_tibble,  conf.int = TRUE,
             surv.median.line = "hv", pval = TRUE,
             legend.labs = c("ISS Stage I", "ISS Stage II", "ISS Stage III"),
             xlab = "Time (days)",
             ggtheme = theme_survminer()))
  dev.off()
  
  # Some survival stats
  # Number of samples necessary data
  seqfish_clinical_info %>% 
    filter(is.na(ISS_Stage) | is.na(EFS_censor)) %>% nrow()
  # Number of censored samples
  summary(EFS_fit)$table[,"n.start"] - summary(EFS_fit)$table[,"events"]
  # Number of samples with event (progression, death)
  summary(EFS_fit)$table[,"events"]
  # Median years of event-free survival
  summary(EFS_fit)$table[,"median"]
  
  # ============================================================================
  # Create Death survival object
  # ============================================================================
  
  # Use EFS_censor == 0 because that is TRUE for death, FALSE for censored
  # Stratify by Stage
  death_tibble <- seqfish_clinical_info %>% 
    filter(!is.na(ISS_Stage), !is.na(D_PT_lstalive)) %>% rowwise() %>% 
    mutate( time_on_trial = max(D_PT_deathdy, D_PT_lstalive, na.rm = TRUE), 
            death = as.numeric(!is.na(D_PT_deathdy)))
  death_fit <- survfit(Surv(time_on_trial, death) ~ ISS_Stage, 
                       data = death_tibble)
  
  # Plot survival curve stratified by Stage
  pdf(str_c(paper_supp, "overall_survival.pdf"), width = 7.25, height = 5)
  print(ggsurvplot(death_fit, data = death_tibble,  conf.int = TRUE,
             surv.median.line = "hv", pval = TRUE, 
             legend.labs = c("ISS Stage I", "ISS Stage II", "ISS Stage III"),
             xlab = "Time (days)",
             ggtheme = theme_survminer()))
  dev.off()
  
  # Some survival stats
  seqfish_clinical_info %>% 
    filter(is.na(ISS_Stage) | is.na(D_PT_lstalive)) %>% nrow()
  # Number of censored samples
  summary(death_fit)$table[,"n.start"] - summary(death_fit)$table[,"events"]
  # Number of samples with event (progression, death)
  summary(death_fit)$table[,"events"]
  # Median years of event-free survival
  summary(death_fit)$table[,"median"]
  
  # ============================================================================
  # Survival table 
  # ============================================================================
  
  summary_tibble <- tribble(
    ~`Category`, ~`ISS Stage`, ~`N Samples`, ~`N Events`, ~`N Censored`, ~`Median Survival (Days)`, ~`95% Confidence Interval (Days)`,
    "Event Free", "Stage I", 
    summary(EFS_fit)$table["ISS_Stage=1","records"],
    summary(EFS_fit)$table["ISS_Stage=1","events"],
    summary(EFS_fit)$table["ISS_Stage=1","records"] - 
      summary(EFS_fit)$table["ISS_Stage=1","events"],
  summary(EFS_fit)$table["ISS_Stage=1","median"],
  str_c(str_replace_na(summary(EFS_fit)$table["ISS_Stage=1","0.95LCL"]),
        str_replace_na(summary(EFS_fit)$table["ISS_Stage=1","0.95UCL"]),
        sep = " - "),
  
  "Event Free", "Stage II", 
  summary(EFS_fit)$table["ISS_Stage=2","records"],
  summary(EFS_fit)$table["ISS_Stage=2","events"],
  summary(EFS_fit)$table["ISS_Stage=2","records"] - 
    summary(EFS_fit)$table["ISS_Stage=2","events"],
  summary(EFS_fit)$table["ISS_Stage=2","median"],
  str_c(str_c(str_replace_na(summary(EFS_fit)$table["ISS_Stage=2","0.95LCL"]), 
              str_replace_na(summary(EFS_fit)$table["ISS_Stage=2","0.95UCL"]), 
              sep = " - ")),
  
  "Event Free", "Stage III", 
  summary(EFS_fit)$table["ISS_Stage=3","records"],
  summary(EFS_fit)$table["ISS_Stage=3","events"],
  summary(EFS_fit)$table["ISS_Stage=3","records"] - 
    summary(EFS_fit)$table["ISS_Stage=3","events"],
  summary(EFS_fit)$table["ISS_Stage=3","median"],
  str_c(str_c(str_replace_na(summary(EFS_fit)$table["ISS_Stage=3","0.95LCL"]), 
              str_replace_na(summary(EFS_fit)$table["ISS_Stage=3","0.95UCL"]), 
              sep = " - ")),
        
  
  "Overall", "Stage I", 
  summary(death_fit)$table["ISS_Stage=1","records"],
  summary(death_fit)$table["ISS_Stage=1","events"],
  summary(death_fit)$table["ISS_Stage=1","records"] - 
    summary(death_fit)$table["ISS_Stage=1","events"],
  summary(death_fit)$table["ISS_Stage=1","median"],
  str_c(str_replace_na(summary(death_fit)$table["ISS_Stage=1","0.95LCL"]),
        str_replace_na(summary(death_fit)$table["ISS_Stage=1","0.95UCL"]),
        sep = " - "),
  
  "Overall", "Stage II", 
  summary(death_fit)$table["ISS_Stage=2","records"],
  summary(death_fit)$table["ISS_Stage=2","events"],
  summary(death_fit)$table["ISS_Stage=2","records"] - 
    summary(death_fit)$table["ISS_Stage=2","events"],
  summary(death_fit)$table["ISS_Stage=2","median"],
  str_c(str_c(str_replace_na(summary(death_fit)$table["ISS_Stage=2","0.95LCL"]), 
              str_replace_na(summary(death_fit)$table["ISS_Stage=2","0.95UCL"]), 
              sep = " - ")),
  
  "Overall", "Stage III", 
  summary(death_fit)$table["ISS_Stage=3","records"],
  summary(death_fit)$table["ISS_Stage=3","events"],
  summary(death_fit)$table["ISS_Stage=3","records"] - 
    summary(death_fit)$table["ISS_Stage=3","events"],
  summary(death_fit)$table["ISS_Stage=3","median"],
  str_c(str_c(str_replace_na(summary(death_fit)$table["ISS_Stage=3","0.95LCL"]), 
              str_replace_na(summary(death_fit)$table["ISS_Stage=3","0.95UCL"]), 
              sep = " - "))
  
  )

  write_tsv(summary_tibble, str_c(paper_supp, "survival_table.txt"), 
            na = "NA", append = FALSE, col_names = TRUE)

}

# ==============================================================================
# Number of fusions per sample
# Originally written August 2018, Updated April 2019
# ==============================================================================

if (TRUE) {
  
  # Create data frame for plotting
  
  n_hdp <- seqfish_clinical_info %>% 
    filter(seqfish_Hyperdiploidy == 1) %>% nrow()
  n_not_hdp <- seqfish_clinical_info %>% 
    filter(seqfish_Hyperdiploidy == 0) %>% nrow()
  n_na_hdp <- seqfish_clinical_info %>% 
    filter(is.na(seqfish_Hyperdiploidy)) %>% nrow()
  hpd_key <- tribble(~seqfish_Hyperdiploidy, ~hyperdiploid_categories, ~count,
                     0, str_c("Non-Hyperdiploid (", n_not_hdp, ")"), n_not_hdp,
                     1, str_c("Hyperdiploid (", n_hdp, ")"), n_hdp,
                     NA, str_c("Not Available (", n_na_hdp, ")"), n_na_hdp
  )
  
  plot_df <- seqfish_clinical_info %>% select(mmrf, seqfish_Hyperdiploidy) %>% 
    left_join(fusions_primary, by = "mmrf") %>%
    group_by(mmrf, seqfish_Hyperdiploidy, srr) %>% 
    summarize(n = n(), 
              n_ig = sum(geneA %in% c("IGH", "IGL", "IGK") | 
                           geneB %in% c("IGH", "IGL", "IGK"))) %>%
    ungroup() %>%
    mutate(n_fusions = n - is.na(srr)) %>%
    left_join(hpd_key, by = "seqfish_Hyperdiploidy")
  
  # t-test to compare number of fusions between hyperdiploidy and not
  
  n_fusions_with_hyperdiploid_info <- seqfish_clinical_info %>% 
    left_join(fusions_primary, by = "mmrf") %>% 
    filter(!is.na(seqfish_Hyperdiploidy)) %>% 
    group_by(mmrf, srr, seqfish_Hyperdiploidy)  %>% 
    summarize(n = n(),
              n_ig = sum(geneA %in% c("IGH", "IGL", "IGK") | 
                           geneB %in% c("IGH", "IGL", "IGK"))) %>% 
    ungroup() %>% 
    mutate(n_corrected = n - is.na(srr))
  
  n_fusions_hyperdiploid <- n_fusions_with_hyperdiploid_info %>% 
    filter(seqfish_Hyperdiploidy == 1) %>% pull(n_corrected)
  n_fusions_nonhyperdiploid <- n_fusions_with_hyperdiploid_info %>% 
    filter(seqfish_Hyperdiploidy == 0) %>% pull(n_corrected)
  n_fusions_hyperdiploid_ig <- n_fusions_with_hyperdiploid_info %>% 
    filter(seqfish_Hyperdiploidy == 1) %>% pull(n_ig)
  n_fusions_nonhyperdiploid_ig <- n_fusions_with_hyperdiploid_info %>%
    filter(seqfish_Hyperdiploidy == 0) %>% pull(n_ig)
  if ( t.test(n_fusions_hyperdiploid, n_fusions_nonhyperdiploid)$p.value < 0.05) {
    print("Number of fusions significantly different between non- and hyperdiploid.")
  } else {
    print("Number of fusions not significantly different between non- and hyperdiploid.")
  }
  if ( t.test(n_fusions_hyperdiploid_ig, n_fusions_nonhyperdiploid_ig)$p.value < 0.05) {
    print("Number of IG fusions significantly different between non- and hyperdiploid.")
  } else {
    print("Number of IG fusions not significantly different between non- and hyperdiploid.")
  }
  
  t.test(n_fusions_with_hyperdiploid_info %>% 
    filter(seqfish_Hyperdiploidy == 1) %>% pull(n_ig),
  n_fusions_with_hyperdiploid_info %>% 
    filter(seqfish_Hyperdiploidy == 0) %>% pull(n_ig))
  
  # Plot number of fusions per sample
  
  plot_df %>% 
    ggplot(aes(x = n_fusions)) + 
    geom_histogram(binwidth = 1, center = 0) + 
    facet_wrap(~ fct_reorder(hyperdiploid_categories, -count), ncol = 1) +
    labs(x = "Number of Fusions Detected", y = "Number of Samples") +
    theme_bw() +
    scale_x_continuous(expand = c(0,0)) +
    theme(panel.background = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.border = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title = element_text(size = 12),
          axis.text.x = element_text(size = 8),
          strip.background = element_blank(),
          strip.text = element_text(size = 12)) +
    ggsave(str_c(paper_supp, "histogram_n_fusions_per_sample.pdf"), 
           device = "pdf", width = 7.25, height = 5, useDingbats = FALSE)
  
  # Plot frequency of number of fusions per sample
  
  n_fusion_tibble <- plot_df %>% group_by(hyperdiploid_categories) %>% 
    summarize(`Median` = median(n_fusions), 
              `Mean` = round(mean(n_fusions),1), 
              `Max` = max(n_fusions),
              `IG` = round(mean(n_ig),1))
  
  n_total_samples <- plot_df %>% nrow()
  
  overall_n_fusion_tibble <- plot_df %>% 
    summarize(hyperdiploid_categories = str_c("Overall (", n_total_samples, ")"),
              `Median` = median(n_fusions),
              `Mean` = round(mean(n_fusions),1),
              `Max` = max(n_fusions),
              `IG` = round(mean(n_ig),1))
  
  n_fusion_all <- n_fusion_tibble %>% bind_rows(overall_n_fusion_tibble)
  
  plot_df %>% ggplot(aes(x = n_fusions, y = ..density..)) + 
    geom_vline(xintercept = c(0, 3), color = "grey70") + 
    geom_freqpoly(aes(color = fct_reorder(hyperdiploid_categories, count)), 
                  binwidth = 1, center = 0, size = 2, show.legend = FALSE) +
    geom_freqpoly(binwidth = 1, center = 0, size = 2, show.legend = FALSE) +
    labs(x = "Number of Fusions Detected (per Sample)", 
         y = "Sample Density",
         color = "Hyperdiploid Category") +
    scale_color_brewer(palette = "Paired", direction = -1) +
    theme_bw() +
    scale_x_continuous(limits = c(plot_df %>% pull(n_fusions) %>% min(),
                                  plot_df %>% pull(n_fusions) %>% max()),
                       breaks = c(0, 3, 20, 40, 60),
                       labels = c(0, 3, 20, 40, 60)) +
    scale_y_continuous() +
    annotation_custom(tableGrob(n_fusion_all, rows = NULL, 
                                cols = c("HRD Status", "Median", "Mean", "Max", "IG"),
                                theme = ttheme_default(core = list(fg_params = list(col = matrix(c("#a6cee3", "#1f78b4", "#b2df8a", rep("#000000", 17)), nrow = 4, byrow = FALSE),
                                                                                    fontface = matrix(rep(c("bold", "plain", "plain", "plain", "plain"), 4), nrow = 4, byrow = TRUE),
                                                                                    fontsize = 8)), 
                                                       colhead = list(fg_params = list( fontsize = 10)))),
                      xmax = 80, ymax = 0.25) + 
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title = element_text(size = 12),
          axis.text.x = element_text(size = 8)) +
    ggsave(str_c(paper_main, "freqpoly_n_fusions_per_sample.pdf"), 
           device = "pdf", width = 5, height = 5/1.618, useDingbats = FALSE)
  
} # change median IG to mean

# ==============================================================================
# Top recurrent fusions (with validation)
# April 2019
# ==============================================================================

if (TRUE) {
  keep_fusions <- fusions_primary %>% group_by(fusion) %>% 
    summarize(count = n(), 
              n_not_na = sum(!is.na(n_discordant)), 
              n_validated = sum(!is.na(n_discordant) & n_discordant >= 3)) %>% 
    filter(n_not_na > 1) %>%
    mutate(validation_pct = 100*n_validated/n_not_na) %>%
    filter(n_validated >= 1) %>% pull(fusion)
  
  total_each_fusion <- fusions_primary %>% filter(fusion %in% keep_fusions) %>%
    mutate(fusion = case_when(geneB %in% c("IGH", "IGK", "IGL") ~ str_c(fusion, "*"), # mark as reciprocal
                              TRUE ~ fusion)) %>%
    group_by(fusion) %>% summarize(total = n())
  
  total_by_status <- fusions_primary %>% filter(fusion %in% keep_fusions) %>% 
    mutate(fusion = case_when(geneB %in% c("IGH", "IGK", "IGL") ~ str_c(fusion, "*"), # mark as reciprocal
                              TRUE ~ fusion)) %>%
    select(fusion, n_discordant) %>% 
    mutate(validation_status = case_when(is.na(n_discordant) ~ "Not Available", 
                                         n_discordant >= 3 ~ "WGS Validated", 
                                         TRUE ~ "Not Validated" )) %>% 
    group_by(fusion, validation_status) %>% 
    summarize(count = n())
  
  plot_df <- total_each_fusion %>% left_join(total_by_status, by = "fusion")
  
  p <- ggplot(data = plot_df, aes(x = fct_reorder(fusion, total), 
                             y = count, 
                             fill = validation_status)) + 
    geom_col() +
    coord_flip() +
    theme_bw() +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0),
                       breaks = seq(0, 100, 20),
                       labels = seq(0, 100, 20),
                       position = "right") +
    scale_fill_manual(values = c("#cbc9e2", "#9e9ac8", "#6a51a3"),
                      breaks = c("WGS Validated", "Not Validated", "Not Available")) +
    theme(panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = 8),
          axis.title = element_text(size = 12),
          axis.text.y = element_text(face = "italic", size = 10),
          legend.position = "bottom") +
    labs(x = NULL, fill = NULL, y = "Top Recurrent Validated Fusions (Number of Samples)")
    
  ggsave(str_c(paper_main, "top_recurrent_validated_fusions.pdf"), p,
           device = "pdf", width = 7.25, height = 7.25/1.618, useDingbats = FALSE)
  ggsave(str_c(paper_main, "top_recurrent_validated_fusions.no_legend.pdf"), 
         p + guides(fill = FALSE), 
         device = "pdf", width = 7.25, height = 7.25/1.618, useDingbats = FALSE)
  
}

# ==============================================================================
# Create a plot showing the overlap/agreement between tools
# Originally written September 2018, Updated April 2019
# ==============================================================================

if (TRUE) {
  
  # Plot tool overlap (upsetr)
  
  upsetr_df <- data.frame(fusions_primary %>% select(starts_with("called_by")))
  names(upsetr_df) <- c("EricScript", "FusionCatcher", "INTEGRATE", 
                        "PRADA", "STAR-Fusion")
  pdf(file = str_c(paper_supp, "tool_overlap.upsetr.pdf"), 
     width = 7.25, height = 5, useDingbats = FALSE)
  upset(upsetr_df, nsets = ncol(upsetr_df), nintersects = NA, order.by = "freq",
        set_size.angles = 90,
        text.scale = 1.5, point.size = 3, line.size = 1)
  dev.off()
}

# ==============================================================================
# Create a summary table of soft filtering
# Written April 2019
# ==============================================================================
soft_filtered %>% 
  group_by(FusionName, filter) %>% 
  summarize(fusion_count = n()) %>% 
  separate(col = "FusionName", into = c("geneA", "geneB"), sep = "--", remove = FALSE) %>%
  select(FusionName, geneA,	geneB, fusion_count, filter) %>%
  bind_rows(undervalidated) %>%
  rename("Fusion Name" = "FusionName", 
         "Fusion Count" = "fusion_count", 
         "Filter" = "filter") %>%
  arrange(desc(`Fusion Count`), `Filter`) %>%
  write_tsv(str_c(paper_supp, "soft_filtering.tsv"))

# ==============================================================================
# Fusion overview paragraph output
# ==============================================================================
n_igh_whsc1 <- fusions_primary %>% 
  filter(fusion %in% c("IGH--WHSC1", "WHSC1--IGH")) %>% 
  pull(srr) %>% unique() %>% length()
n_igh_whsc1_with_wgs <- fusions_primary %>% 
  filter(fusion %in% c("IGH--WHSC1")) %>% 
  filter(!is.na(n_discordant)) %>% 
  pull(srr) %>% unique() %>% length()
n_igh_whsc1_validated <- fusions_primary %>% 
  filter(fusion %in% c("IGH--WHSC1")) %>% 
  filter(!is.na(n_discordant), n_discordant >= 3) %>% 
  pull(srr) %>% unique() %>% length()
n_fusions_with_wgs <- fusions_primary %>% filter(!is.na(n_discordant)) %>% nrow()
n_fusions_validated <- fusions_primary %>% filter(!is.na(n_discordant), n_discordant >= 3) %>% nrow()
post_filtering_validation_rate <- n_fusions_validated/n_fusions_with_wgs
n_igh_fusions <- fusions_primary %>% 
  filter(geneA %in% c("IGH", "IGK", "IGL") | 
           geneB %in% c("IGH", "IGK", "IGL")) %>% nrow()
n_fusions_total <- fusions_primary %>% nrow()
prop_igh_fusions <- n_igh_fusions/n_fusions_total
myc_pvt1_ig_fusions <- fusions_primary %>% 
  filter(geneA %in% c("MYC", "PVT1") & geneB %in% c("IGH", "IGK", "IGL")) %>% 
  pull(fusion) %>% table() %>% sort()
print(str_c("IGH WHSC1 fusions: ", n_igh_whsc1, "/", n_samples_primary, " = ", round(100*n_igh_whsc1/n_samples_primary, 1)))
print(str_c("IGH WHSC1 validated: ", n_igh_whsc1_validated, "/", n_igh_whsc1_with_wgs, " = ", round(100*n_igh_whsc1_validated/n_igh_whsc1_with_wgs, 1)))
print(str_c("Overall fusion validation rate: ", round(100*post_filtering_validation_rate, 1), "%"))
print(str_c("Fusions involving IGH/IGK/IGL: ", n_igh_fusions, "/", n_fusions_total, " = ", round(100*prop_igh_fusions, 1), "%"))
print("MYC or PVT1 and IG fusions: ")
print(myc_pvt1_ig_fusions)
print(str_c("Number of fusion tools: (out of ", n_fusions_total, ")"))
fusions_primary %>% pull(CallerN) %>% table()
fusions_primary %>% pull(CallerN) %>% table()/n_fusions_total
print(str_c("Number of significantly undervalidated fusions: ", significantly_under_validated_fusions %>% length()))

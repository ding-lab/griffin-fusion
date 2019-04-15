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
        mutate_all(funs(replace(., is.na(.), 2))) )
    
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
    ggplot2_standard_additions()
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
    ggplot2_standard_additions()
  print(p)
  dev.off()
  
  pdf(str_c(paper_supp, "clinical_variables.ldh.pdf"), height = 10, width = 15,
      useDingbats = FALSE)
  n_missing <- seqfish_clinical_info %>% filter(is.na(LDH)) %>% nrow()
  p <- seqfish_clinical_info %>% 
    ggplot(aes(x = LDH)) + geom_histogram() +
    labs(x = "Lactate dehydrogenase (LDH) (U/L)", y = "Number of patients",
         caption = str_c("Number missing = ", n_missing)) +
    ggplot2_standard_additions()
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
  pdf(str_c(paper_supp, "event_free_survival.pdf"), width = 20, height = 15,
      useDingbats = FALSE)
  ggsurvplot(EFS_fit, data = EFS_tibble,  conf.int = TRUE,
             surv.median.line = "hv", pval = TRUE,
             legend.labs = c("ISS Stage I", "ISS Stage II", "ISS Stage III"),
             xlab = "Time (days)",
             ggtheme = theme_bw(base_size = 20))
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
  pdf(str_c(paper_supp, "overall_survival.pdf"), width = 20, height = 15)
  ggsurvplot(death_fit, data = death_tibble,  conf.int = TRUE,
             surv.median.line = "hv", pval = TRUE, 
             legend.labs = c("ISS Stage I", "ISS Stage II", "ISS Stage III"),
             xlab = "Time (days)",
             ggtheme = theme_bw(base_size = 20))
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
# Explore the landscape of fusions in samples with multiple timepoints
# Originally written September 2018, Updated April 2019
# ==============================================================================

if (TRUE) {
  
  # Samples with multiple timepoints
  
  samples_with_multiple_timepoints <- fusions_all %>% 
    filter(has_secondary == 1) %>% pull(mmrf) %>% unique()
  
  # Fusions involving important genes
  
  fusions_with_important_genes <- fusions_all %>%
    filter( fusion == "IGH--WHSC1" | 
              !(geneA %in% c("IGH", "IGK", "IGL") |
                  geneB %in% c("IGH", "IGK", "IGL"))) %>%
    filter( geneA_oncogene | geneA_tsg | geneA_kinase | 
              geneA_mmy_known | geneA_driver |
              geneB_oncogene | geneB_tsg | geneB_kinase |
              geneB_mmy_known | geneB_driver |
              drug_fusion | drug_geneA | drug_geneB) %>% 
    select(mmrf, fusion) %>% unique() %>% group_by(fusion) %>% 
    summarize(n()) %>% #filter(`n()` > 1) %>% 
    pull(fusion)
  
  # plot it 
  keep_these_mmrfs <- fusions_all %>% filter(fusion %in% fusions_with_important_genes, has_secondary) %>% pull(mmrf) %>% unique()
  
  keep_these_srrs <- samples_all %>% group_by(mmrf) %>% summarize(n()) %>% 
    filter(`n()` > 1) %>% 
    left_join(samples_all, by = "mmrf") %>% 
    mutate(mmrf_srr = str_c(mmrf, srr, sep = ": ")) %>% group_by(mmrf) %>% 
    mutate(ticker = row_number()) %>% select(mmrf, srr, mmrf_srr, ticker) %>%
    ungroup() %>% filter(mmrf %in% keep_these_mmrfs)
  
  keep_these_srrs %>% left_join(fusions_all, by = c("mmrf", "srr")) %>%
    filter(fusion %in% fusions_with_important_genes) %>%
    right_join(keep_these_srrs, by = c("mmrf", "srr")) %>% 
    mutate(mmrf_number_only = str_remove_all(mmrf, "MMRF_")) %>%
    replace_na(list(fusion = "Zero detected")) %>%
    ggplot(aes(y = factor(srr), x = fusion, fill = log10(FFPM + 1))) + 
    geom_tile(color = "black") +
    geom_text(aes(label = sample_number), size = 1.75, vjust = 0.5) +
    theme_bw() +
    facet_wrap(~ mmrf_number_only , ncol = 1, strip.position = "right", dir = "h", scales = "free_y") +
    theme(panel.grid.major = element_line(size = 0.1),
          panel.border = element_rect(size = 0.1),
          axis.ticks = element_line(size = 0.1),
          strip.background = element_blank(),
          axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1, size = 5, face = "italic"),
          axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 5),
          strip.text.y = element_text(angle = 0),
          legend.position = "bottom",
          panel.spacing.y = unit(0.025, units = "inches")) +
    scale_fill_gradient(low = "#ffedbc", high = "#ed4264",
                        limits = c(0,1.05), 
                        breaks = seq(0, 1, 0.25),
                        labels = seq(0, 1, 0.25)) + 
    scale_x_discrete(drop = FALSE) +
    labs(y = "Sample Number", x = NULL, fill = "Scaled FFPM") +
    ggsave(str_c(paper_supp, "multiple_timepoints.pdf"),
         device = "pdf", width = 8.5, height = 11)
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
    summarize(n = n()) %>% ungroup() %>%
    mutate(n_fusions = n - is.na(srr)) %>%
    left_join(hpd_key, by = "seqfish_Hyperdiploidy")
  
  # Plot number of fusions per sample
  
  plot_df %>% 
    ggplot(aes(x = n_fusions)) + 
    geom_histogram(binwidth = 1, center = 0) + 
    facet_wrap(~ fct_reorder(hyperdiploid_categories, -count), ncol = 1) +
    labs(x = "Number of Fusions Detected", y = "Number of Samples") +
    ggplot2_standard_additions() +
    ggsave(str_c(paper_supp, "histogram_n_fusions_per_sample.pdf"), 
           device = "pdf", width = 10, height = 5)
  
  # Plot frequency of number of fusions per sample
  
  n_fusion_tibble <- plot_df %>% group_by(hyperdiploid_categories) %>% 
    summarize(`Median` = median(n_fusions), 
              `Mean` = round(mean(n_fusions),1), 
              `Max` = max(n_fusions))
  
  overall_n_fusion_tibble <- plot_df %>% 
    summarize(hyperdiploid_categories = "Overall",
              `Median` = median(n_fusions),
              `Mean` = round(mean(n_fusions),1),
              `Max` = max(n_fusions))
  
  n_fusion_all <- n_fusion_tibble %>% bind_rows(overall_n_fusion_tibble)
  
  plot_df %>% ggplot(aes(x = n_fusions, y = ..density..)) + 
    geom_freqpoly(aes(color = fct_reorder(hyperdiploid_categories, -count)), 
                  binwidth = 1, center = 0, size = 2, show.legend = FALSE) + 
    labs(x = "Number of Fusions Detected", 
         y = "Sample Proportion",
         color = "Hyperdiploid Category") +
    scale_color_brewer(palette = "Set2") +
    ggplot2_standard_additions() +
    xlim(0,70) +
    theme(legend.position = "bottom",
          legend.direction = "vertical") +
    scale_y_continuous() +
    annotation_custom(tableGrob(n_fusion_all, rows = NULL, 
                                cols = c("", "Median", "Mean", "Max"),
                                theme = ttheme_default(core = list(fg_params = list(col = matrix(c("#66c2a5", "#fc8d62", "#8da0cb" ,rep("#000000", 13)), nrow = 4, byrow = FALSE))))),
                      xmax = 80, ymax = 0.25) + 
    ggsave(str_c(paper_main, "freqpoly_n_fusions_per_sample.pdf"), 
           device = "pdf", width = 6, height = 6)
  
}

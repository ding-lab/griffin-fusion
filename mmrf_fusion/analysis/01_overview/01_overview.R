# ==============================================================================
# Overview (MMRF Fusions)
# Steven Foltz (github: envest)
# ==============================================================================

paper_main = "paper/main/01_overview/"
paper_supp = "paper/supplementary/01_overview/"

# Create directories 
dir.create(paper_main, recursive = TRUE, showWarnings = FALSE)
dir.create(paper_supp, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# Clinical features of MMRF patients
# ==============================================================================

if (TRUE) {
  
  # ============================================================================
  # Number of non/hyperdiploid samples
  # ============================================================================
  
  n_hyperdiploid_samples <- seqfish_clinical_info %>% 
    filter(SeqWGS_Cp_Hyperdiploid_Call == 1) %>%
    nrow()
  n_nonhyperdiploid_samples <- seqfish_clinical_info %>% 
    filter(SeqWGS_Cp_Hyperdiploid_Call == 0) %>% 
    nrow()
  n_na_hyperdiploid_samples <- seqfish_clinical_info %>% 
    filter(is.na(SeqWGS_Cp_Hyperdiploid_Call)) %>% 
    nrow()
  
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
  
  # What drug treatments patients received
  treatment_summary <- seqfish_clinical_info %>% 
    mutate_at("D_PT_therclass", factor) %>% 
    select(D_PT_therclass) %>% 
    summary(maxsum = length(unique(seqfish_clinical_info$D_PT_therclass)))
  treatment_na <- seqfish_clinical_info %>% 
    filter(is.na(D_PT_therclass)) %>% 
    nrow()
  
  # How many with bone marrow transplant
  bmt_summary <- seqfish_clinical_info %>% 
    mutate_at("BMT", factor) %>% select(BMT) %>% summary()
  bmt_na <- seqfish_clinical_info %>% 
    filter(is.na(BMT)) %>% nrow()
  
  # ============================================================================
  # Data type summary
  # ============================================================================
  
  n_samples_total <- samples_all %>% nrow()
  n_samples_primary <- samples_primary %>% nrow()
  n_samples_seqfish <- seqfish_clinical_info %>% 
    filter(!is.na(seqfish_Study_Visit_ID)) %>% nrow()
  n_samples_seqfish_na <- seqfish_clinical_info %>% 
    filter(is.na(seqfish_Study_Visit_ID)) %>% nrow()
  
  n_samples_t_seqfish <- seqfish_clinical_info %>% 
    select(updated_seqfish_t_IGH_WHSC1) %>% 
    filter(!is.na(updated_seqfish_t_IGH_WHSC1)) %>% 
    nrow()
  n_samples_t_seqfish_na <- seqfish_clinical_info %>% 
    select(updated_seqfish_t_IGH_WHSC1) %>% 
    filter(is.na(updated_seqfish_t_IGH_WHSC1)) %>% 
    nrow()
  
  n_samples_cnv_seqfish <- seqfish_clinical_info %>% 
    select(SeqWGS_Cp_12p13_20percent) %>% 
    filter(!is.na(SeqWGS_Cp_12p13_20percent)) %>% 
    nrow()
  n_samples_cnv_seqfish_na <- seqfish_clinical_info %>% 
    select(SeqWGS_Cp_12p13_20percent) %>% 
    filter(is.na(SeqWGS_Cp_12p13_20percent)) %>% 
    nrow()
  
  # X14 is the column of SRR,SRR for WGS tumor,normal
  n_samples_wgs <- file_locations %>%
    filter(X14 != "NA,NA") %>%
    nrow()
  n_samples_wgs_na <- file_locations %>% 
    filter(X14 == "NA,NA") %>% 
    nrow()
  
  # ============================================================================
  # Summary table output
  # ============================================================================
  
  if (TRUE) {
    summary_tibble <- tribble(
      ~Category, ~Subcategory, ~N, ~Missing, 
      ~Percentage, ~Min, ~`25%`, ~Median, ~Mean, ~`75%`, ~Max,
      "Hyperdiploid Status", 
      "Hyperdiploid",
      n_hyperdiploid_samples,
      ".", 
      n_hyperdiploid_samples/(n_samples_primary - n_na_hyperdiploid_samples),
      ".",  ".",    ".",     ".",   ".",    ".",
      
      "Hyperdiploid Status", 
      "Non-Hyperdiploid",
      n_nonhyperdiploid_samples,
      ".", 
      n_nonhyperdiploid_samples/(n_samples_primary - n_na_hyperdiploid_samples),
      ".",  ".",    ".",     ".",   ".",    ".",  
      
      "Hyperdiploid Status", 
      NA,
      n_na_hyperdiploid_samples,
      ".", 
      ".",
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
                           simplify = TRUE)[1,2])/(n_samples_primary - sex_na),
      ".", ".", ".", ".", ".", ".",
      
      "Sex",
      "Male",
      as.numeric(str_split(sex_summary[1,1], ":", simplify = TRUE)[1,2]),
      ".",
      as.numeric(str_split(sex_summary[1,1], ":", 
                           simplify = TRUE)[1,2])/(n_samples_primary - sex_na),
      ".", ".", ".", ".", ".", ".",
      
      "Sex",
      NA,
      sex_na,
      ".",
      ".",
      ".", ".", ".", ".", ".", ".",
      
      "Race",
      "White",
      as.numeric(str_split(race_summary[1,1], ":", simplify = TRUE)[1,2]),
      ".",
      as.numeric(str_split(race_summary[1,1], ":", 
                           simplify = TRUE)[1,2])/(n_samples_primary - race_na),
      ".", ".", ".", ".", ".", ".",
      
      "Race",
      "Black",
      as.numeric(str_split(race_summary[2,1], ":", simplify = TRUE)[1,2]),
      ".",
      as.numeric(str_split(race_summary[2,1], ":", 
                           simplify = TRUE)[1,2])/(n_samples_primary - race_na),
      ".", ".", ".", ".", ".", ".",
      
      "Race",
      "Other",
      as.numeric(str_split(race_summary[3,1], ":", simplify = TRUE)[1,2]),
      ".",
      as.numeric(str_split(race_summary[3,1], ":", 
                           simplify = TRUE)[1,2])/(n_samples_primary - race_na),
      ".", ".", ".", ".", ".", ".",
      
      "Race",
      NA,
      race_na,
      ".",
      ".",
      ".", ".", ".", ".", ".", ".",
      
      "Bone Marrow Transplant",
      "No",
      as.numeric(str_split(bmt_summary[1,1], ":", simplify = TRUE)[1,2]),
      ".",
      as.numeric(str_split(bmt_summary[1,1], ":", simplify = TRUE)[1,2])/(n_samples_primary - bmt_na),
      ".", ".", ".", ".", ".", ".",
      
      "Bone Marrow Transplant",
      "Yes",
      as.numeric(str_split(bmt_summary[2,1], ":", simplify = TRUE)[1,2]),
      ".",
      as.numeric(str_split(bmt_summary[2,1], ":", simplify = TRUE)[1,2])/(n_samples_primary - bmt_na),
      ".", ".", ".", ".", ".", ".",
      
      "Bone Marrow Transplant",
      NA,
      bmt_na,
      ".",
      ".",
      ".", ".", ".", ".", ".", ".",
      
      "Treatment",
      "Bortezomib-based",
      as.numeric(str_split(treatment_summary[1,1], ":", simplify = TRUE)[1,2]),
      ".",
      as.numeric(str_split(treatment_summary[1,1], ":", 
                           simplify = TRUE)[1,2])/(n_samples_primary - treatment_na),
      ".", ".", ".", ".", ".", ".",
      
      "Treatment",
      "Carfilzomib-based",
      as.numeric(str_split(treatment_summary[2,1], ":", simplify = TRUE)[1,2]),
      ".",
      as.numeric(str_split(treatment_summary[2,1], ":", 
                           simplify = TRUE)[1,2])/(n_samples_primary - treatment_na),
      ".", ".", ".", ".", ".", ".",
      
      "Treatment",
      "combined bortezomib/carfilzomib-based",
      as.numeric(str_split(treatment_summary[3,1], ":", simplify = TRUE)[1,2]),
      ".",
      as.numeric(str_split(treatment_summary[3,1], ":", 
                           simplify = TRUE)[1,2])/(n_samples_primary - treatment_na),
      ".", ".", ".", ".", ".", ".",
      
      "Treatment",
      "combined bortezomib/IMIDs-based",
      as.numeric(str_split(treatment_summary[4,1], ":", simplify = TRUE)[1,2]),
      ".",
      as.numeric(str_split(treatment_summary[4,1], ":", 
                           simplify = TRUE)[1,2])/(n_samples_primary - treatment_na),
      ".", ".", ".", ".", ".", ".",
      
      "Treatment",
      "combined bortezomib/IMIDs/carfilzomib-based",
      as.numeric(str_split(treatment_summary[5,1], ":", simplify = TRUE)[1,2]),
      ".",
      as.numeric(str_split(treatment_summary[5,1], ":", 
                           simplify = TRUE)[1,2])/(n_samples_primary - treatment_na),
      ".", ".", ".", ".", ".", ".",
      
      "Treatment",
      "combined IMIDs/carfilzomib-based",
      as.numeric(str_split(treatment_summary[6,1], ":", simplify = TRUE)[1,2]),
      ".",
      as.numeric(str_split(treatment_summary[6,1], ":", 
                           simplify = TRUE)[1,2])/(n_samples_primary - treatment_na),
      ".", ".", ".", ".", ".", ".",
      
      "Treatment",
      "IMIDs-based",
      as.numeric(str_split(treatment_summary[7,1], ":", simplify = TRUE)[1,2]),
      ".",
      as.numeric(str_split(treatment_summary[7,1], ":", 
                           simplify = TRUE)[1,2])/(n_samples_primary - treatment_na),
      ".", ".", ".", ".", ".", ".",
      
      "Treatment",
      NA,
      treatment_na,
      ".",
      ".",
      ".", ".", ".", ".", ".", ".",
      
      
      "ECOG",
      "0 = Fully active",
      as.numeric(str_split(ecog_summary[1,1], ":", simplify = TRUE)[1,2]),
      ".",
      as.numeric(str_split(ecog_summary[1,1], ":", 
                           simplify = TRUE)[1,2])/(n_samples_primary - as.numeric(str_split(ecog_summary[6,1], ":", simplify = TRUE)[1,2])),
      ".", ".", ".", ".", ".", ".",
      
      "ECOG",
      "1 = Restricted in physically strenuous activity",
      as.numeric(str_split(ecog_summary[2,1], ":", simplify = TRUE)[1,2]),
      ".",
      as.numeric(str_split(ecog_summary[2,1], ":", 
                           simplify = TRUE)[1,2])/(n_samples_primary - as.numeric(str_split(ecog_summary[6,1], ":", simplify = TRUE)[1,2])),
      ".", ".", ".", ".", ".", ".",
      
      "ECOG",
      "2 = Ambulatory and capable of all self-care",
      as.numeric(str_split(ecog_summary[3,1], ":", simplify = TRUE)[1,2]),
      ".",
      as.numeric(str_split(ecog_summary[3,1], ":", 
                           simplify = TRUE)[1,2])/(n_samples_primary - as.numeric(str_split(ecog_summary[6,1], ":", simplify = TRUE)[1,2])),
      ".", ".", ".", ".", ".", ".",
      
      "ECOG",
      "3 = Capable of only limited self-care",
      as.numeric(str_split(ecog_summary[4,1], ":", simplify = TRUE)[1,2]),
      ".",
      as.numeric(str_split(ecog_summary[4,1], ":", 
                           simplify = TRUE)[1,2])/(n_samples_primary - as.numeric(str_split(ecog_summary[6,1], ":", simplify = TRUE)[1,2])),
      ".", ".", ".", ".", ".", ".",
      
      "ECOG",
      "4 = Completely disabled",
      as.numeric(str_split(ecog_summary[5,1], ":", simplify = TRUE)[1,2]),
      ".",
      as.numeric(str_split(ecog_summary[5,1], ":", 
                           simplify = TRUE)[1,2])/(n_samples_primary - as.numeric(str_split(ecog_summary[6,1], ":", simplify = TRUE)[1,2])),
      ".", ".", ".", ".", ".", ".",
      
      "ECOG",
      NA,
      as.numeric(str_split(ecog_summary[6,1], ":", simplify = TRUE)[1,2]),
      ".",
      ".",
      ".", ".", ".", ".", ".", ".", 
      
      "ISS Stage",
      "I",
      as.numeric(str_split(stage_summary[1,1], ":", simplify = TRUE)[1,2]),
      ".",
      as.numeric(str_split(stage_summary[1,1], ":", 
                           simplify = TRUE)[1,2])/(n_samples_primary - as.numeric(str_split(stage_summary[4,1], ":", simplify = TRUE)[1,2])),
      ".", ".", ".", ".", ".", ".", 
      
      "ISS Stage",
      "II",
      as.numeric(str_split(stage_summary[2,1], ":", simplify = TRUE)[1,2]),
      ".",
      as.numeric(str_split(stage_summary[2,1], ":", 
                           simplify = TRUE)[1,2])/(n_samples_primary - as.numeric(str_split(stage_summary[4,1], ":", simplify = TRUE)[1,2])),
      ".", ".", ".", ".", ".", ".", 
      
      "ISS Stage",
      "III",
      as.numeric(str_split(stage_summary[3,1], ":", simplify = TRUE)[1,2]),
      ".",
      as.numeric(str_split(stage_summary[3,1], ":", 
                           simplify = TRUE)[1,2])/(n_samples_primary - as.numeric(str_split(stage_summary[4,1], ":", simplify = TRUE)[1,2])),
      ".", ".", ".", ".", ".", ".", 
      
      "ISS Stage",
      NA,
      as.numeric(str_split(stage_summary[4,1], ":", simplify = TRUE)[1,2]),
      ".",
      ".",
      ".", ".", ".", ".", ".", ".",
      
      "Bone lesions",
      "No",
      as.numeric(str_split(bone_summary[1,1], ":", simplify = TRUE)[1,2]),
      ".",
      as.numeric(str_split(bone_summary[1,1], ":", 
                           simplify = TRUE)[1,2])/(n_samples_primary - bone_na),
      ".", ".", ".", ".", ".", ".",
      
      "Bone lesions",
      "Yes",
      as.numeric(str_split(bone_summary[2,1], ":", simplify = TRUE)[1,2]),
      ".",
      as.numeric(str_split(bone_summary[2,1], ":", 
                           simplify = TRUE)[1,2])/(n_samples_primary - bone_na),
      ".", ".", ".", ".", ".", ".",
      
      "Bone lesions",
      NA,
      bone_na,
      ".",
      ".",
      ".", ".", ".", ".", ".", ".",
      
      "Plasmacytoma",
      "No",
      as.numeric(str_split(plasmacytoma_summary[1,1], ":", simplify = TRUE)[1,2]),
      ".",
      as.numeric(str_split(plasmacytoma_summary[1,1], ":", 
                           simplify = TRUE)[1,2])/(n_samples_primary - plasmacytoma_na),
      ".", ".", ".", ".", ".", ".",
      
      "Plasmacytoma",
      "Yes",
      as.numeric(str_split(plasmacytoma_summary[2,1], ":", simplify = TRUE)[1,2]),
      ".",
      as.numeric(str_split(plasmacytoma_summary[2,1], ":", 
                           simplify = TRUE)[1,2])/(n_samples_primary - plasmacytoma_na),
      ".", ".", ".", ".", ".", ".", 
      
      "Plasmacytoma",
      NA,
      plasmacytoma_na,
      ".",
      ".",
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
      "Primary (with RNA-seq + seqFISH translocations)",
      n_samples_t_seqfish,
      n_samples_t_seqfish_na,
      n_samples_t_seqfish/n_samples_primary,
      ".", ".", ".", ".", ".", ".",
      
      "Number of samples",
      "Primary (with RNA-seq + seqFISH CNV)",
      n_samples_cnv_seqfish,
      n_samples_cnv_seqfish_na,
      n_samples_cnv_seqfish/n_samples_primary,
      ".", ".", ".", ".", ".", ".",
      
      "Number of samples",
      "Primary (with RNA-seq + WGS)",
      n_samples_wgs,
      n_samples_wgs_na,
      n_samples_wgs/n_samples_primary,
      ".", ".", ".", ".", ".", "."
    )
    
    summary_tibble <- summary_tibble %>% arrange(Category, Subcategory)
    write_tsv(summary_tibble, str_c(paper_supp, "summary_table.txt"), 
              na = "NA", append = FALSE, col_names = TRUE)  
  }
  
}

# ==============================================================================
# Basic survival plots of MMRF patients
# ==============================================================================

if (TRUE) {
  
  # ============================================================================
  # Create EFS survival object
  # ============================================================================
  
  # Use EFS_censor == 0 because that is TRUE for death, FALSE for censored
  # Stratify by Stage
  EFS_tibble <- seqfish_clinical_info %>% 
    filter(!is.na(ISS_Stage), !is.na(EFS))
  EFS_fit <- survfit(Surv(EFS, EFS_censor == 0) ~ ISS_Stage, 
                     data = EFS_tibble)
  
  # Plot survival curve stratified by Stage
  pdf(str_c(paper_supp, "event_free_survival.pdf"), width = 3.5, height = 3.5,
      useDingbats = FALSE)
  
  print(ggsurvplot(EFS_fit, data = EFS_tibble, conf.int = TRUE,
                   surv.median.line = "hv", pval = TRUE,
                   legend.labs = c("ISS Stage I", "ISS Stage II", "ISS Stage III"),
                   legend = "bottom",
                   xlab = "Time (days)", 
                   ylab = "Progression-Free Survival Probability",
                   ggtheme = theme_survminer(base_size = 12,
                                             base_family = "",
                                             font.main = c(12, "plain", "black"),
                                             font.submain = c(12, "plain", "black"),
                                             font.x = c(12, "plain", "black"),
                                             font.y = c(12, "plain", "black"),
                                             font.caption = c(12, "plain", "black"),
                                             font.tickslab = c(8, "plain", "black"),
                                             legend = c("top", "bottom", "left", "right", "none"),
                                             font.legend = c(8, "plain", "black")),
                   conf.int.alpha = 0.1))
  
  dev.off()
  
  pdf(str_c(paper_supp, "event_free_survival.no_legend.pdf"), width = 3.5, height = 3.5,
      useDingbats = FALSE)
  
  print(ggsurvplot(EFS_fit, data = EFS_tibble, conf.int = TRUE,
                   surv.median.line = "hv", pval = TRUE,
                   legend.labs = c("ISS Stage I", "ISS Stage II", "ISS Stage III"),
                   legend = "none",
                   xlab = "Time (days)", 
                   ylab = "Progression-Free Survival Probability",
                   ggtheme = theme_survminer(base_size = 12,
                                             base_family = "",
                                             font.main = c(12, "plain", "black"),
                                             font.submain = c(12, "plain", "black"),
                                             font.x = c(12, "plain", "black"),
                                             font.y = c(12, "plain", "black"),
                                             font.caption = c(12, "plain", "black"),
                                             font.tickslab = c(8, "plain", "black"),
                                             legend = c("top", "bottom", "left", "right", "none"),
                                             font.legend = c(8, "plain", "black")),
                   conf.int.alpha = 0.1))
  
  dev.off()
  
  # ============================================================================
  # Create Death survival object
  # ============================================================================
  
  # Use OS_censor == 0 because that is TRUE for death, FALSE for censored
  # Stratify by Stage
  death_tibble <- seqfish_clinical_info %>% 
    filter(!is.na(ISS_Stage), !is.na(OS))
  death_fit <- survfit(Surv(OS, OS_censor == 0) ~ ISS_Stage, data = death_tibble)
  
  # Plot survival curve stratified by Stage
  pdf(str_c(paper_supp, "overall_survival.pdf"), width = 3.5, height = 3.5,
      useDingbats = FALSE)
  
  print(ggsurvplot(death_fit, data = death_tibble, conf.int = TRUE,
                   surv.median.line = "hv", pval = TRUE,
                   legend.labs = c("ISS Stage I", "ISS Stage II", "ISS Stage III"),
                   legend = "bottom",
                   xlab = "Time (days)", 
                   ylab = "Overall Survival Probability",
                   ggtheme = theme_survminer(base_size = 12,
                                             base_family = "",
                                             font.main = c(12, "plain", "black"),
                                             font.submain = c(12, "plain", "black"),
                                             font.x = c(12, "plain", "black"),
                                             font.y = c(12, "plain", "black"),
                                             font.caption = c(12, "plain", "black"),
                                             font.tickslab = c(8, "plain", "black"),
                                             legend = c("top", "bottom", "left", "right", "none"),
                                             font.legend = c(8, "plain", "black")),
                   conf.int.alpha = 0.1))
  
  dev.off()
  
  pdf(str_c(paper_supp, "overall_survival.no_legend.pdf"), width = 3.5, height = 3.5,
      useDingbats = FALSE)
  
  print(ggsurvplot(death_fit, data = death_tibble, conf.int = TRUE,
                   surv.median.line = "hv", pval = TRUE,
                   legend.labs = c("ISS Stage I", "ISS Stage II", "ISS Stage III"),
                   legend = "none",
                   xlab = "Time (days)", 
                   ylab = "Overall Survival Probability",
                   ggtheme = theme_survminer(base_size = 12,
                                             base_family = "",
                                             font.main = c(12, "plain", "black"),
                                             font.submain = c(12, "plain", "black"),
                                             font.x = c(12, "plain", "black"),
                                             font.y = c(12, "plain", "black"),
                                             font.caption = c(12, "plain", "black"),
                                             font.tickslab = c(8, "plain", "black"),
                                             legend = c("top", "bottom", "left", "right", "none"),
                                             font.legend = c(8, "plain", "black")),
                   conf.int.alpha = 0.1))
  
  dev.off()
  
  # ============================================================================
  # Survival table 
  # ============================================================================
  
  if (TRUE) {
    survival_tibble <- tribble(
      ~`Category`, ~`ISS Stage`, ~`N Samples`, ~`N Events`, ~`N Censored`, 
      ~`Median Survival (Days)`, ~`95% Confidence Interval (Days)`,
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
    
    write_tsv(survival_tibble, str_c(paper_supp, "survival_table.txt"), 
              na = "NA", append = FALSE, col_names = TRUE)
  }
  
  # ============================================================================
  # Some overall survival figures
  # ============================================================================
  
  n_patients_progressed <- seqfish_clinical_info %>% 
    filter(!is.na(EFS), !is.na(EFS_censor)) %>% 
    filter(EFS_censor == 0) %>% 
    nrow()
  n_patients_with_pfs_data <- seqfish_clinical_info %>% 
    filter(!is.na(EFS), !is.na(EFS_censor)) %>% 
    nrow()
  n_patients_died <- seqfish_clinical_info %>% 
    filter(!is.na(OS), !is.na(OS_censor)) %>% 
    filter(OS_censor == 0) %>% 
    nrow()
  n_patients_with_os_data <- seqfish_clinical_info %>% 
    filter(!is.na(OS), !is.na(OS_censor)) %>% 
    nrow()
  summary_os_followup <- seqfish_clinical_info %>% 
    filter(!is.na(OS), !is.na(OS_censor)) %>% 
    pull(OS) %>% 
    summary()
  summary_efs_followup <- seqfish_clinical_info %>% 
    filter(!is.na(EFS), !is.na(EFS_censor)) %>%
    pull(EFS) %>%
    summary()
}

# ==============================================================================
# Number of fusions per sample
# ==============================================================================

if (TRUE) {
  
  # Create data frame for plotting
  
  n_hdp <- seqfish_clinical_info %>% 
    filter(SeqWGS_Cp_Hyperdiploid_Call == 1) %>% nrow()
  n_not_hdp <- seqfish_clinical_info %>% 
    filter(SeqWGS_Cp_Hyperdiploid_Call == 0) %>% nrow()
  n_na_hdp <- seqfish_clinical_info %>% 
    filter(is.na(SeqWGS_Cp_Hyperdiploid_Call)) %>% nrow()
  hpd_key <- tribble(~SeqWGS_Cp_Hyperdiploid_Call, ~hyperdiploid_categories, ~count,
                     0, str_c("Non-Hyperdiploid (", n_not_hdp, ")"), n_not_hdp,
                     1, str_c("Hyperdiploid (", n_hdp, ")"), n_hdp,
                     NA, str_c("Not Available (", n_na_hdp, ")"), n_na_hdp)
  
  # MWU to compare number of fusions between hyperdiploidy and not
  
  n_fusions_with_hyperdiploid_info <- seqfish_clinical_info %>% 
    left_join(fusions_primary, by = "mmrf") %>% 
    filter(!is.na(SeqWGS_Cp_Hyperdiploid_Call)) %>% 
    group_by(mmrf, srr, SeqWGS_Cp_Hyperdiploid_Call)  %>% 
    summarize(n = n(),
              n_ig = sum(geneA %in% c("IGH", "IGL", "IGK") | 
                           geneB %in% c("IGH", "IGL", "IGK"))) %>% 
    ungroup() %>% 
    mutate(n_corrected = n - is.na(srr)) # 0 for those without any calls
  
  n_fusions_hyperdiploid <- n_fusions_with_hyperdiploid_info %>% 
    filter(SeqWGS_Cp_Hyperdiploid_Call == 1) %>% 
    pull(n_corrected)
  n_fusions_nonhyperdiploid <- n_fusions_with_hyperdiploid_info %>% 
    filter(SeqWGS_Cp_Hyperdiploid_Call == 0) %>% 
    pull(n_corrected)
  n_fusions_hyperdiploid_ig <- n_fusions_with_hyperdiploid_info %>% 
    filter(SeqWGS_Cp_Hyperdiploid_Call == 1) %>% 
    pull(n_ig)
  n_fusions_nonhyperdiploid_ig <- n_fusions_with_hyperdiploid_info %>%
    filter(SeqWGS_Cp_Hyperdiploid_Call == 0) %>% 
    pull(n_ig)
  
  wilcox_n_fusions <- wilcox.test(n_fusions_hyperdiploid, n_fusions_nonhyperdiploid)
  
  wilcox_n_ig_fusions <- wilcox.test(n_fusions_hyperdiploid_ig, n_fusions_nonhyperdiploid_ig)
  
  # Plot number of fusions per sample

  plot_df <- seqfish_clinical_info %>% 
    select(mmrf, SeqWGS_Cp_Hyperdiploid_Call) %>% 
    left_join(fusions_primary, by = "mmrf") %>%
    group_by(mmrf, SeqWGS_Cp_Hyperdiploid_Call, srr) %>% 
    summarize(n = n(), 
              n_ig = sum(geneA %in% c("IGH", "IGL", "IGK") | 
                           geneB %in% c("IGH", "IGL", "IGK"))) %>%
    ungroup() %>%
    mutate(n_fusions = n - is.na(srr)) %>%
    left_join(hpd_key, by = "SeqWGS_Cp_Hyperdiploid_Call")
    
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
           device = "pdf", width = 7.25, height = 4, useDingbats = FALSE)
  
  # Plot frequency of number of fusions per sample
  
  n_fusion_tibble <- plot_df %>% 
    group_by(hyperdiploid_categories) %>% 
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
  
  plot_df %>% 
    ggplot(aes(x = n_fusions, y = ..density..)) + 
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
                                cols = c("Hyperdiploid Status", "Median", "Mean", "Max", "IG"),
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
  
}

# ==============================================================================
# Top recurrent fusions (with WGS support)
# ==============================================================================

if (TRUE) {
  keep_fusions <- fusions_primary %>% 
    group_by(fusion) %>% 
    summarize(count = n(), 
              n_not_na = sum(!is.na(n_discordant)), 
              n_validated = sum(!is.na(n_discordant) & n_discordant >= 3)) %>% 
    filter(n_not_na > 1) %>%
    mutate(validation_pct = 100*n_validated/n_not_na) %>%
    filter(n_validated >= 1) %>% pull(fusion)
  
  total_each_fusion <- fusions_primary %>% 
    filter(fusion %in% keep_fusions) %>%
    mutate(fusion = case_when(geneB %in% c("IGH", "IGK", "IGL") ~ str_c(fusion, "*"), # mark as reciprocal
                              TRUE ~ fusion)) %>%
    group_by(fusion) %>% summarize(total = n())
  
  total_by_status <- fusions_primary %>% filter(fusion %in% keep_fusions) %>% 
    mutate(fusion = case_when(geneB %in% c("IGH", "IGK", "IGL") ~ str_c(fusion, "*"), # mark as reciprocal
                              TRUE ~ fusion)) %>%
    select(fusion, n_discordant) %>% 
    mutate(validation_status = case_when(is.na(n_discordant) ~ "Not Available", 
                                         n_discordant >= 3 ~ "WGS Supported", 
                                         TRUE ~ "Not Supported" )) %>% 
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
                      breaks = c("WGS Supported", "Not Supported", "Not Available")) +
    theme(panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = 8),
          axis.title = element_text(size = 12),
          axis.text.y = element_text(face = "italic", size = 10),
          legend.position = "bottom") +
    labs(x = NULL, fill = NULL, y = "Top Recurrent Fusions Supported by WGS (Number of Samples)")
    
  ggsave(str_c(paper_main, "top_recurrent_validated_fusions.pdf"), p,
           device = "pdf", width = 7.25, height = 7.25/1.618, useDingbats = FALSE)
  ggsave(str_c(paper_main, "top_recurrent_validated_fusions.no_legend.pdf"), 
         p + guides(fill = FALSE), 
         device = "pdf", width = 7.25, height = 7.25/1.618, useDingbats = FALSE)
  
}

# ==============================================================================
# Create a plot showing the overlap/agreement between tools
# ==============================================================================

if (TRUE) {
  
  # Plot tool overlap (upsetr)
  
  upsetr_df <- data.frame(fusions_primary %>% select(starts_with("called_by")))
  
  names(upsetr_df) <- c("EricScript", "FusionCatcher", "INTEGRATE", 
                        "PRADA", "STAR-Fusion")
  
  pdf(file = str_c(paper_supp, "tool_overlap.upsetr.pdf"), 
     width = 7.25, height = 4, useDingbats = FALSE)
  
  upset(upsetr_df, 
        nsets = ncol(upsetr_df),
        nintersects = NA,
        order.by = "freq",
        set_size.angles = 90,
        text.scale = 1.5,
        point.size = 3,
        line.size = 1)
  
  dev.off()
  
}

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
n_fusions_with_wgs <- fusions_primary %>% 
  filter(!is.na(n_discordant)) %>% nrow()
n_fusions_validated <- fusions_primary %>% 
  filter(!is.na(n_discordant), n_discordant >= 3) %>% nrow()
post_filtering_validation_rate <- n_fusions_validated/n_fusions_with_wgs
n_igh_fusions <- fusions_primary %>% 
  filter(geneA %in% c("IGH", "IGK", "IGL") | 
           geneB %in% c("IGH", "IGK", "IGL")) %>% nrow()
n_fusions_total <- fusions_primary %>% nrow()
prop_igh_fusions <- n_igh_fusions/n_fusions_total
myc_pvt1_ig_fusions <- fusions_primary %>% 
  filter((geneA %in% c("MYC", "PVT1") & geneB %in% c("IGH", "IGK", "IGL")) |
           (geneB %in% c("MYC", "PVT1") & geneA %in% c("IGH", "IGK", "IGL"))) %>%
  pull(fusion) %>% table() %>% sort()

print(str_c("Number of patients: ", n_samples_primary))
print("Tissue Sources:")
print(samples_primary %>%
        left_join(samples_all, by = "srr") %>%
        pull(tissue_source) %>%
        table())
print(str_c("Primary = pre-treatment: ", 
            length(mmrf_primary_pretreatment), "/", n_samples_primary, " = ", 
            round(100*length(mmrf_primary_pretreatment)/n_samples_primary, 2), "%"))
print(str_c("Number with additional samples: ", samples_all %>% 
              group_by(mmrf) %>% 
              summarize(count = n()) %>% 
              filter(count > 1) %>% 
              nrow()))
print(str_c("Number of total RNA samples: ", n_samples_total))

print("Patient Age:")
print(summary_tibble %>% filter(Category == "Age"))
print("Patient Stage:")
print(summary_tibble %>% filter(Category == "ISS Stage"))
print("Summary of EFS followup:")
print(summary_efs_followup)
print(summary_efs_followup/365.25)
print(str_c("Number of patients progressed: ", 
            n_patients_progressed, "/", n_patients_with_pfs_data, " = ",
            round(100*n_patients_progressed/n_patients_with_pfs_data, 3), "%"))
print("Summary of OS followup:")
print(summary_os_followup)
print(summary_os_followup/365.25)
print(str_c("Number of patients died: ",
            n_patients_died, "/", n_patients_with_os_data, " = ", 
            round(100*n_patients_died/n_patients_with_os_data, 3), "%"))
print("Survival by ISS Stage:")
print(survival_tibble %>% filter(Category == "Event Free")) %>% 
  pull(`Median Survival (Days)`)/365.25
print("Hyperdiploid Status:")
print(summary_tibble %>% filter(Category == "Hyperdiploid Status"))
print("Patient Ancestry:")
print(summary_tibble %>% filter(Category == "Race"))
print("Treatment Regimens:")
print(summary_tibble %>% filter(Category == "Treatment"))
print("Bone Marrow Transplants:")
print(summary_tibble %>% filter(Category == "Bone Marrow Transplant"))

print(str_c("IGH WHSC1 fusions: ", 
            n_igh_whsc1, "/", n_samples_primary, " = ",
            round(100*n_igh_whsc1/n_samples_primary, 1), "%"))
print(str_c("IGH WHSC1 validated: ",
            n_igh_whsc1_validated, "/", n_igh_whsc1_with_wgs, " = ", 
            round(100*n_igh_whsc1_validated/n_igh_whsc1_with_wgs, 1), "%"))
print(str_c("Fusions involving IGH/IGK/IGL: ", 
            n_igh_fusions, "/", n_fusions_total, " = ", 
            round(100*prop_igh_fusions, 1), "%"))
print("MYC or PVT1 and IG fusions: ")
print(myc_pvt1_ig_fusions)

print(str_c("Number of fusions HRD: ", 
            round(mean(n_fusions_hyperdiploid), 1)))
print(str_c("Number of fusions non-HRD: ", 
            round(mean(n_fusions_nonhyperdiploid), 1)))
print(wilcox_n_fusions)

print(str_c("Number of Ig fusions HRD: ", 
            round(mean(n_fusions_hyperdiploid_ig), 1)))
print(str_c("Number of Ig fusions non-HRD: ", 
            round(mean(n_fusions_nonhyperdiploid_ig), 1)))
print(wilcox_n_ig_fusions)

print("Number of fusions overall:")
print(seqfish_clinical_info %>% pull(total_fusions) %>% summary())

print(str_c("Number of significantly undervalidated fusions: ", 
            significantly_under_validated_fusions %>% length()))
print(str_c("Overall fusion validation rate: ", 
            round(100*post_filtering_validation_rate, 1), "%"))
print(str_c("TCGA Pancancer fusion validation rate: ", 
            round(100*tcga_validation_rate, 1), "%"))
print(str_c("Number of fusion tools: (out of ", n_fusions_total, ")"))
print(fusions_primary %>% pull(CallerN) %>% table())
print(fusions_primary %>% pull(CallerN) %>% table()/n_fusions_total)

# Make supplementary tables

supp_tables <- "paper/tables/supplementary/"
dir.create(supp_tables, recursive = TRUE, showWarnings = FALSE)

# Table S1 patient characteristics and data availability

# summary_tibble from 01_overview.R
write_tsv(summary_tibble, 
          str_c(supp_tables, "Table_S1.Patient_Info.tsv"), 
          na = "NA", append = FALSE, col_names = TRUE)

# Table S2 fusions reported from all samples
write_tsv(fusions_all %>% 
            select(-wgs_bam) %>%
            select(-discordant_reads),
          str_c(supp_tables, "Table_S2.Fusion_Calls.tsv"), 
          na = "NA", append = FALSE, col_names = TRUE)

# Table S3 overexpressed fusions
overexpressed_table <- testing_tbl_pvalue_adjusted %>% 
  filter(fdr < 0.05 | median_value > 0.9, 
         event_type %in% c("Fusion Expression", 
                           "Fusion Expression Outlier")) %>% 
  select(-c("event2", "n_samples_with_tested", 
            "n_samples_without_tested", "n_samples_na_tested"))
write_tsv(overexpressed_table, 
          str_c(supp_tables, "Table_S3.Overexpressed.tsv"), 
          na = "NA", append = FALSE, col_names = TRUE)

# Table S4 kinase information
kinase_table <- read_tsv("data/Kinase_fusion_info.txt") %>% 
  mutate(Fusion = str_remove(Fusion, "@")) %>% 
  mutate(kinase_group_full_name = case_when(Group == "TK" ~ "Tyrosine Kinase",
                                            Group == "OTHER" ~ "Other",
                                            Group == "TKL" ~ "Tyrosine Kinase-Like",
                                            TRUE ~ Group))
write_tsv(kinase_table, 
          str_c(supp_tables, "Table_S4.Kinases.tsv"), 
          na = "NA", append = FALSE, col_names = TRUE)

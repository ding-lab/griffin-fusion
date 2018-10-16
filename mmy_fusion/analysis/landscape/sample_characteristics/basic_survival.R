# ==============================================================================
# Sample survival
# Steven Foltz (smfoltz@wustl.edu), September 2018
# ==============================================================================

plot_dir = "analysis/landscape/sample_characteristics/"

# ==============================================================================
# Load necessary libraries
# ==============================================================================

library(survival)
library(survminer)

# ==============================================================================
# Create EFS survival object
# ==============================================================================

# Use EFS_censor == 0 because that is TRUE for death, FALSE for censored
# Stratify by Stage
EFS_tibble <- seqfish_clinical_info %>% 
  filter(!is.na(ISS_Stage), !is.na(EFS_censor))
EFS_fit <- survfit(Surv(EFS, EFS_censor == 0) ~ ISS_Stage, 
                   data = EFS_tibble)

# Plot survival curve stratified by Stage
pdf(str_c(plot_dir, "event_free_survival.pdf"), width = 20, height = 15)
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

# ==============================================================================
# Create Death survival object
# ==============================================================================

# Use EFS_censor == 0 because that is TRUE for death, FALSE for censored
# Stratify by Stage
death_tibble <- seqfish_clinical_info %>% 
  filter(!is.na(ISS_Stage), !is.na(D_PT_lstalive)) %>% rowwise() %>% 
  mutate( time_on_trial = max(D_PT_deathdy, D_PT_lstalive, na.rm = TRUE), 
          death = as.numeric(!is.na(D_PT_deathdy)))
death_fit <- survfit(Surv(time_on_trial, death) ~ ISS_Stage, 
                   data = death_tibble)

# Plot survival curve stratified by Stage
pdf(str_c(plot_dir, "overall_survival.pdf"), width = 20, height = 15)
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
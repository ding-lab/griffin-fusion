# ==============================================================================
# Fusions involving oncogenes and tumor suppressors
# Steven Foltz (smfoltz@wustl.edu), November 2018
# ==============================================================================

fused_oncogenes_geneA <- fusions_primary %>% 
  filter(geneA_oncogene == 1) %>% select(geneA) %>% unique()

fused_oncogenes_geneB <- fusions_primary %>% 
  filter(geneB_oncogene == 1) %>% select(geneB) %>% unique()

fused_tsg_geneA <- fusions_primary %>% 
  filter(geneA_tsg == 1) %>% select(geneA) %>% unique()

fused_tsg_geneB <- fusions_primary %>% 
  filter(geneB_tsg == 1) %>% select(geneB) %>% unique()

fused_kinase_geneA <- fusions_primary %>% 
  filter(geneA_kinase == 1) %>% select(geneA) %>% unique()

fused_kinase_geneB <- fusions_primary %>% 
  filter(geneB_kinase == 1) %>% select(geneB) %>% unique()
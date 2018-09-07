# ==============================================================================
# Perform GO enrichment on genes involved in fusions
# Steven Foltz (smfoltz@wustl.edu), September 2018
# ==============================================================================

# ==============================================================================
# Gene list
# ==============================================================================

# There is a conflict with select() if AnnotationDbi library is loaded
fusion_genes_gt2 <- fusions_primary %>% 
  gather(geneA, geneB, key = "geneAB", value = "fusion_gene") %>%
  dplyr::select(mmrf, srr, fusion_gene) %>% distinct() %>% 
  group_by(fusion_gene) %>% summarize(count = n()) %>% filter(count >= 3) %>% 
  arrange(desc(count))

# ==============================================================================
# Load relevant packages
# ==============================================================================
library(biomaRt)
library(org.Hs.eg.db)
library(topGO)

# ==============================================================================
# Set up biomaRt for gene list conversion
# ==============================================================================

# Obtain the mart
mart <- useDataset("hsapiens_gene_ensembl", mart = useMart("ensembl"))

# Get Ensembl gene IDs
all_ensembl_gene_ids <- getBM(attributes = "ensembl_gene_id", 
                             values = "*", mart = mart)
names(all_ensembl_gene_ids) <- c("ensg")

# Give interesting genes value 1, rest get 0
interesting_ensg_genes <- ensg_gene_list %>% 
  left_join(all_ensembl_gene_ids, by = "ensg") %>%
  mutate(interesting = as.numeric(gene %in% fusion_genes_gt2$fusion_gene))

# Give interesting genes value 1. The rest get value 0.
all <- interesting_ensg_genes %>% pull(interesting) %>% factor()
names(all) <- interesting_ensg_genes %>% pull(ensg)
  
# GO enrichment using genes with value 1
GOdata <- new("topGOdata", ontology = "CC", allGenes = all, 
              geneSel = function(p) p == 1,
              description = "Genes in MMRF Fusions",
              annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "Ensembl")



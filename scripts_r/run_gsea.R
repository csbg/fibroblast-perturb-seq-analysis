library(fgsea)
library(tidyverse)
library(fs)
source("scripts_r/utils.R")



# Load data ---------------------------------------------------------------

dge <- 
  dir_ls("data_generated/", regex = "DEG_") %>% 
  read_csv(id = "file") %>% 
  separate_wider_regex(file, c(".+DEG_", comparison = ".*", "\\.csv")) %>% 
  select(comparison, group, gene = names,
         logFC = logfoldchanges, p_adj = pvals_adj)

gene_map <-
  read_tsv("metadata/biomart_human_mouse_20210727.tsv") %>% 
  select(human_gene = `Human gene name`, mouse_gene = `Gene name`) %>% 
  distinct()

enrichr_genesets <- read_rds("data_generated/enrichr_genesets_mouse.rds")



# Perform enrichment ------------------------------------------------------

# returns dataframe, comprising columns "db", "comparison", "group",
# as well as all columns in the result of fgsea().
run_gsea <- function(comparison, group, db) {
  ranked_genes <-
    dge %>%
    filter(comparison == {{comparison}}, group == {{group}}) %>%
    arrange(desc(logFC)) %>% 
    select(gene, logFC) %>%
    deframe()
  ranked_genes <- ranked_genes[!is.na(ranked_genes)]
  
  fgsea(
    enrichr_genesets[[db]],
    ranked_genes,
    eps = 0
  ) %>%
    as_tibble() %>%
    mutate(
      db = {{db}},
      comparison = {{comparison}},
      group = {{group}},
      .before = 1
    )
}


gsea_results <-
  expand_grid(
    comparison = unique(dge$comparison),
    group = unique(dge$group),
    db = names(enrichr_genesets)
  ) %>% 
  pmap(run_gsea, .progress = TRUE) %>% 
  list_rbind()



# Save results ------------------------------------------------------------

write_rds(dge, "data_generated/dge.rds")
write_rds(gsea_results, "data_generated/gsea.rds")

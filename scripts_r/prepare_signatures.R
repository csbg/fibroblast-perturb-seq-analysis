library(tidyverse)
library(readxl)
source("scripts_r/utils.R")



# Fibroblast signatures ---------------------------------------------------

cluster_names_forte <- c(
  "0"  = "Homeostatic Epicardial Derived",
  "1"  = "Late Resolution",
  "2"  = "Progenitor Like State",
  "3"  = "Myofibroblasts",
  "4"  = "Phagocytic",
  "5"  = "Injury Response",
  "6"  = "Matrifibrocyte",
  "7"  = "Endocardial Derived",
  "8"  = "Proliferative",
  "9"  = "Epicardium",
  "10" = "Interferon Response",
  "11" = "Dendritic Cells"
)

gene_map <-
  read_tsv("metadata/biomart_human_mouse_20210727.tsv") %>% 
  select(gene_human = `Human gene name`, gene_mouse = `Gene name`) %>% 
  distinct()

markers_human <-
  bind_rows(
    .id = "ref",
    
    Amrute =
      read_xlsx("data_raw/signatures/signatures_Amrute.xlsx", sheet = 1, skip = 1) %>% 
      select(
        cluster,
        gene_human = ...1,
        logFC = avg_log2FC,
        p_adj = p_val_adj
      ),
    
    Chaffin =
      read_xlsx("data_raw/signatures/signatures_Chaffin.xlsx", sheet = 1) %>% 
      select(
        cluster = `Sub-Cluster`,
        gene_human = Gene,
        logFC = `limma-voom: logFC`,
        p_adj = `limma-voom: Adjusted P-Value`
      ),
    
    Fu =
      read_xlsx("data_raw/signatures/signatures_Fu.xlsx", sheet = 1) %>% 
      select(
        cluster,
        gene_human = gene,
        logFC = avg_logFC,
        p_adj = p_val_adj
      ),
    
    Koenig = 
      read_xlsx("data_raw/signatures/markers_Koenig.xlsx", sheet = 3) %>% 
      select(
        cluster,
        gene_human = gene,
        logFC = avg_log2FC,
        p_adj = p_val_adj
      ),
    
    Kuppe =
      c("Fib1", "Fib2", "Fib3", "Fib4") %>% 
      set_names() %>% 
      map(\(s) read_xlsx("data_raw/signatures/signatures_Kuppe.xlsx", sheet = s)) %>% 
      list_rbind(names_to = "cluster") %>% 
      select(
        cluster, 
        gene_human = gene, 
        logFC = avg_log2FC,
        p_adj = p_val_adj
      )
  ) %>% 
  left_join(
    gene_map,
    by = "gene_human",
    relationship = "many-to-many"
  )

markers_mouse <-
  bind_rows(
    .id = "ref",
    
    # the Buechler dataset can be simply loaded
    Buechler = 
      map(
        1:10,
        \(s) read_xlsx("data_raw/signatures/markers_Buechler.xlsx", sheet = s)
      ) %>% 
      list_rbind() %>% 
      select(
        cluster,
        gene_mouse = Gene,
        logFC = avg_logFC,
        p_adj = p_val_adj
      ),
    
    # the Forte dataset contains numbered clusters,
    # to which we assign meaningful names
    Forte =
      map(
        1:11,
        \(s) read_xlsx("data_raw/signatures/markers_Forte.xlsx", sheet = s)
      ) %>% 
      list_rbind() %>%
      select(
        cluster,
        gene_mouse = gene,
        logFC = avg_logFC,
        p_adj = p_val_adj
      ) %>%
      mutate(cluster = recode(cluster, !!!cluster_names_forte))
  ) %>% 
  left_join(
    gene_map,
    by = "gene_mouse",
    relationship = "many-to-many"
  )

markers_own <- 
  read_csv("data_raw/signatures/markers_murine.csv", comment = "#") %>% 
  left_join(
    gene_map,
    by = "gene_mouse",
    relationship = "many-to-many"
  )


# table with six columns: dataset, cluster (cell type), marker gene (human + mouse),
# log fold change, and adjusted p value
# filtering: only keep significant genes (padj < 0.01)
# and the 100 genes with highest logFC per cluster
markers <- 
  bind_rows(markers_human, markers_mouse) %>% 
  filter(p_adj < 0.01) %>%
  slice_max(logFC, n = 100, by = c(ref, cluster)) %>%
  bind_rows(markers_own) %>% 
  relocate(gene_mouse, .after = gene_human)

write_rds(markers, "data_generated/fibroblast_markers.rds")



# Gene sets for GSEA ------------------------------------------------------

enrichr_databases <- c(
  # pathways
  "KEGG_2019_Mouse",
  "MSigDB_Hallmark_2020",
  "NCI-Nature_2016",
  "Reactome_2022",
  "WikiPathways_2019_Mouse",
  
  # transcription factor targets
  "ChEA_2022",
  "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X",
  "ENCODE_TF_ChIP-seq_2015",
  "TRANSFAC_and_JASPAR_PWMs",
  "TRRUST_Transcription_Factors_2019"
)

# download Enrichr databases in a format that can be used by fgsea:
# assemble a named list (names = databases)
#   of named lists (names = enrichment terms)
#   of character vectors (all genes associated with the respective term).
enrichr_genesets_human <-
  enrichr_databases %>% 
  set_names() %>% 
  map(\(db) {
    info("Downloading {db}")
    url <- paste0(
      "https://maayanlab.cloud/Enrichr/geneSetLibrary",
      "?mode=text&libraryName=",
      db
    )
    db <- read_lines(url)
    m <- str_match(db, "(.+?)\\t\\t(.+)")
    terms <- m[, 2]
    genes <- m[, 3] %>% str_split("\\t")
    genes %>% 
      map(stringi::stri_remove_empty) %>% 
      set_names(terms)
  })

write_rds(enrichr_genesets_human, "data_generated/enrichr_genesets_human.rds")

# convert human to mouse genes
enrichr_genesets_mouse <- 
  enrichr_genesets_human %>% 
  modify(\(db) {
    modify(db, \(gene_set) {
      gene_map %>% 
        filter(gene_human %in% gene_set) %>% 
        pull(gene_mouse) %>% 
        unique()   
    })
  })

write_rds(enrichr_genesets_human, "data_generated/enrichr_genesets_mouse.rds")

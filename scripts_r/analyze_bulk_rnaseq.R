library(ComplexHeatmap)
library(edgeR)
library(tidyverse)
source("scripts_r/utils.R")



# Human -------------------------------------------------------------------

## Load data ----

counts_human <- read_tsv("data_raw/rna-seq/DataRaw_human.txt")
samples_human <-
  read_csv("data_raw/rna-seq/sample_data.csv", comment = "#") %>% 
  filter(organism == "human")


rna_data_unfiltered <- DGEList(
  counts =
    counts_human %>%
    select(gene_name, Kat5_pos_1_2015_S15:Mock_pos_2017_S17),
  samples = 
    samples_human %>%
    select(!sample),
  group = samples_human$condition,
  genes = 
    counts_human %>% 
    select(gene_id, gene_type, level:strand)
)


## Perform DGE ----

genes_to_keep <- filterByExpr(rna_data_unfiltered)
rna_data <- rna_data_unfiltered[genes_to_keep,, keep.lib.sizes = FALSE]

rbind(
  cpm(rna_data_unfiltered, log = TRUE) %>% 
    magrittr::set_rownames(rep("original", nrow(.))),
  cpm(rna_data, log = TRUE) %>% 
    magrittr::set_rownames(rep("filtered", nrow(.)))
) %>% 
  as_tibble(rownames = "dataset") %>% 
  mutate(dataset = as_factor(dataset)) %>% 
  pivot_longer(!dataset, names_to = "sample", values_to = "log_cpm") %>% 
  ggplot(aes(log_cpm, color = sample)) + 
  geom_density() +
  facet_wrap(vars(dataset)) +
  theme_pub()



rna_data <- calcNormFactors(rna_data)

plotMDS(cpm(rna_data, log = TRUE))

design <- model.matrix(~condition, data = rna_data$samples)
colnames(design) <- c("(Intercept)", "Kat5_inh", "TGFb")
design

efit <-
  voom(rna_data, design, plot = TRUE) %>% 
  lmFit(design) %>% 
  eBayes()

plotSA(efit)

dge <- bind_rows(
  .id = "comparison",
  TGFb =
    topTable(efit, coef = "TGFb", number = Inf) %>%
    filter(gene_type == "protein_coding") %>% 
    as_tibble() %>% 
    select(gene = gene_name, logFC, p = P.Value, p_adj = adj.P.Val),
  Kat5_inh =
    topTable(efit, coef = "Kat5_inh", number = Inf) %>%
    filter(gene_type == "protein_coding") %>% 
    as_tibble() %>% 
    select(gene = gene_name, logFC, p = P.Value, p_adj = adj.P.Val)
)
  
dge



## Analyze results ----

dge %>% 
  filter(p_adj <= 0.5)

ggplot(dge, aes(logFC, -log10(p))) +
  geom_point() +
  facet_wrap(vars(comparison))


top_genes <- c(
  dge %>%
    filter(comparison == "TGFb") %>%
    slice_max(logFC, n = 50, with_ties = FALSE) %>%
    pull(gene),
  dge %>%
    filter(comparison == "TGFb") %>%
    slice_min(logFC, n = 50, with_ties = FALSE) %>%
    pull(gene)
)

mat_human <- 
  cpm(rna_data, log = TRUE) %>%
  magrittr::set_rownames(rna_data$genes$gene_name) %>% 
  t() %>% 
  scale() %>% 
  t()

selected_samples <- 
  samples_human %>% 
  filter(condition != "Kat5_inhibitor1") %>% 
  pull(sample)

Heatmap(
  mat_human[top_genes, selected_samples],
  cluster_rows = FALSE
)


## GSEA ----

run_gsea <- function(comparison, db) {
  ranked_genes <-
    dge %>%
    filter(comparison == {{comparison}}) %>%
    arrange(desc(logFC)) %>% 
    select(gene, logFC) %>%
    deframe()
  ranked_genes <- ranked_genes[!is.na(ranked_genes)]
  
  # info("GSEA of comparison {comparison}, group {group}, ",
  #      "db {db} ({length(ranked_genes)} genes)")
  
  fgsea(
    enrichr_genesets[[db]],
    ranked_genes,
    eps = 0
  ) %>%
    as_tibble() %>%
    mutate(
      db = {{db}},
      comparison = {{comparison}},
      .before = 1
    )
}

gsea_results <-
  expand_grid(
    comparison = unique(dge$comparison),
    db = c(
      "MSigDB_Hallmark_2020",
      "WikiPathways_2019_Mouse",
      "Reactome_2022",
      "KEGG_2019_Mouse"
    )
  ) %>% 
  pmap(run_gsea, .progress = TRUE) %>% 
  list_rbind()

ggplot(gsea_results, aes(NES, -log10(padj))) +
  geom_point(alpha = 0.25)

gsea_results %>% 
  filter(padj <= 0.05)





# Mouse -------------------------------------------------------------------

counts_mouse <- read_tsv("data_raw/rna-seq/DataRaw_mouse.txt")
samples_mouse <-
  read_csv("data_raw/rna-seq/sample_data.csv", comment = "#") %>% 
  filter(organism == "mouse") %>% 
  mutate(condition = factor(condition) %>% fct_relevel("PBS"))


rna_data_unfiltered <- DGEList(
  counts =
    counts_mouse %>%
    select(gene_name, IK22_1A_S4:IK33_3D_S3),
  samples = 
    samples_mouse %>%
    select(!sample),
  group = samples_mouse$condition,
  genes = 
    counts_mouse %>% 
    select(gene_id, gene_type, level:strand)
)

genes_to_keep <- filterByExpr(rna_data_unfiltered)
rna_data <- rna_data_unfiltered[genes_to_keep,, keep.lib.sizes = FALSE]

rbind(
  cpm(rna_data_unfiltered, log = TRUE) %>% 
    magrittr::set_rownames(rep("original", nrow(.))),
  cpm(rna_data, log = TRUE) %>% 
    magrittr::set_rownames(rep("filtered", nrow(.)))
) %>% 
  as_tibble(rownames = "dataset") %>% 
  mutate(dataset = as_factor(dataset)) %>% 
  pivot_longer(!dataset, names_to = "sample", values_to = "log_cpm") %>% 
  ggplot(aes(log_cpm, color = sample)) + 
  geom_density() +
  facet_wrap(vars(dataset)) +
  theme_pub()



rna_data <- calcNormFactors(rna_data)

plotMDS(cpm(rna_data, log = TRUE))

design <- model.matrix(~condition, data = rna_data$samples)
colnames(design) <- c("(Intercept)", "Kat5_inh1", "Kat5_inh2")
design

efit <-
  voom(rna_data, design, plot = TRUE) %>% 
  lmFit(design) %>% 
  eBayes()

plotSA(efit)

dge <-
  topTable(efit, coef = "Kat5_inh2", number = Inf) %>%
  as_tibble()

dge %>% 
  filter(adj.P.Val <= 0.5)

ggplot(dge, aes(logFC, -log10(P.Value))) +
  geom_point()

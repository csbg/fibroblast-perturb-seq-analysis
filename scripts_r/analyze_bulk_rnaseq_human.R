library(ComplexHeatmap)
library(edgeR)
library(fgsea)
library(tidyverse)
source("scripts_r/utils.R")

ht_opt(
  simple_anno_size = unit(1.5, "mm"),
  COLUMN_ANNO_PADDING = unit(1, "pt"),
  DENDROGRAM_PADDING = unit(1, "pt"),
  HEATMAP_LEGEND_PADDING = unit(1, "mm"),
  ROW_ANNO_PADDING = unit(1, "pt"),
  TITLE_PADDING = unit(2, "mm"),
  heatmap_row_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
  heatmap_row_names_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
  heatmap_column_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
  heatmap_column_names_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
  legend_labels_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
  legend_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
  legend_border = FALSE
)



# Load data ---------------------------------------------------------------

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



# Perform DGE -------------------------------------------------------------

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
  geom_density(key_glyph = "point") +
  facet_wrap(vars(dataset)) +
  theme_pub() +
  theme(
    legend.key.height = unit(2, "mm"),
    legend.key.width = unit(2, "mm")
  )
ggsave_default("rnaseq/human_logcpm_dist", width = 100, height = 50)


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
dge %>% save_table("rnaseq_human_dge")



# Analyze results ---------------------------------------------------------

dge %>% 
  filter(p_adj <= 0.5)

ggplot(dge, aes(logFC, -log10(p))) +
  geom_point() +
  facet_wrap(vars(comparison))
ggsave_default("rnaseq/human_volcano")

# top 30 genes with pos/neg logFC
top_genes <- c(
  dge %>%
    filter(comparison == "TGFb") %>%
    slice_max(logFC, n = 30, with_ties = FALSE) %>%
    pull(gene),
  dge %>%
    filter(comparison == "TGFb") %>%
    slice_min(logFC, n = 30, with_ties = FALSE) %>%
    pull(gene)
)

selected_samples <- 
  samples_human %>% 
  filter(condition != "Kat5_inhibitor1") %>% 
  pull(sample)

mat_human <- 
  cpm(rna_data, log = TRUE) %>%
  magrittr::set_rownames(rna_data$genes$gene_name) %>% 
  t() %>% 
  scale() %>% 
  t() %>% 
  magrittr::extract(top_genes, selected_samples)


(p <- Heatmap(
  mat_human,
  cluster_rows = FALSE,
  width = ncol(mat_human) * unit(2, "mm"),
  height = nrow(mat_human) * unit(2, "mm"),
))
ggsave_default("rnaseq/human_dge_heatmap", plot = p, width = 100)



# Perform GSEA ------------------------------------------------------------

# enrichr_genesets from run_gsea.R, saved as RDS file

enrichr_genesets <- read_rds("data_generated/enrichr_genesets_human.rds")

run_gsea <- function(comparison, db) {
  ranked_genes <-
    dge %>%
    filter(comparison == {{comparison}}) %>%
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
ggsave_default("rnaseq/human_gsea_nes")

gsea_results %>% filter(padj <= 0.05)
gsea_results %>% save_table("rnaseq_human_gsea")



# Plot GSEA ---------------------------------------------------------------

selected_terms <- tribble(
  ~db, ~pathway, ~term_display,
  "MSigDB_Hallmark_2020",
  "TNF-alpha Signaling via NF-kB",
  "TNF-alpha Signaling via NF-kB (MSigDB)",
  "MSigDB_Hallmark_2020",
  "Inflammatory Response",
  "Inflammatory Response (MSigDB)",
  "KEGG_2019_Mouse",
  "JAK-STAT signaling pathway",
  "JAK-STAT signaling pathway (KEGG)",
  "KEGG_2019_Mouse",
  "Hedgehog signaling pathway",
  "Hedgehog signaling pathway (KEGG)",
  "WikiPathways_2019_Mouse",
  "Striated Muscle Contraction WP216",
  "Striated Muscle Contraction (WikiPathways)",
  "MSigDB_Hallmark_2020",
  "Oxidative Phosphorylation",
  "Oxidative Phosphorylation (MSigDB)",
  "Reactome_2022",
  "Mitochondrial Translation R-HSA-5368287",
  "Mitochondrial Translation (Reactome)",
  "WikiPathways_2019_Mouse",
  "Retinol metabolism WP1259",
  "Retinol metabolism (WikiPathways)",
  "MSigDB_Hallmark_2020",
  "G2-M Checkpoint",
  "G2-M Checkpoint (MSigDB)"
)

plot_terms <- function(terms) {
  terms <- 
    terms %>% 
    mutate(term_display = as_factor(term_display))
  
  # prepare data
  data_vis <- 
    gsea_results %>% 
    left_join(terms, by = join_by(db, pathway)) %>% 
    filter(!is.na(term_display))
  
  # set limits for colorbar
  color_limit <- max(abs(data_vis$NES), na.rm = TRUE)
  
  # make plot  
  ggplot(data_vis, aes(comparison, term_display, size = -log10(padj))) +
    geom_point(aes(color = NES)) +
    xlab("Condition") +
    ylab("Term (database)") +
    scale_color_distiller(
      name = "Normalized\nenrichment\nscore",
      palette = "RdBu",
      direction = -1,
      limits = c(-color_limit, color_limit),
      breaks = c(-color_limit, 0, color_limit),
      labels = function(x) round(x, 2),
      guide = guide_colorbar(
        barheight = unit(15, "mm"),
        barwidth = unit(2, "mm"),
        ticks = FALSE
      )
    ) +
    scale_size_area(
      name = "-log10 padj",
      max_size = 3,
      limits = c(0, 10),
      breaks = c(0, 5, 10),
      labels = c("0", "5", "10 or higher"),
      oob = scales::oob_squish
    )  +
    coord_fixed() +
    theme_pub(TRUE) +
    theme(
      legend.key.height = unit(3, "mm"),
      legend.key.width = unit(3, "mm"),
      legend.margin = margin(),
      panel.grid = element_blank(),
      panel.spacing = unit(-.5, "pt")
    )
}

plot_terms(selected_terms)
ggsave_default("rnaseq/human_gsea", type = "pdf", width = 65)

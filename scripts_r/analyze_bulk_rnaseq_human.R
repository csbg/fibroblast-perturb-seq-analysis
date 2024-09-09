library(ComplexHeatmap)
library(edgeR)
library(fgsea)
library(tidyverse)
library(ggrepel)
library(khroma)
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

counts_raw <- read_tsv("data_raw/rna-seq/DataRaw_human.txt")
samples <-
  read_csv("data_raw/rna-seq/sample_data.csv", comment = "#") %>% 
  filter(organism == "human")


rna_data_unfiltered <- DGEList(
  counts =
    counts_raw %>%
    select(gene_name, Kat5_pos_1_2015_S15:Mock_pos_2017_S17),
  samples = 
    samples %>%
    select(!sample),
  group = samples$condition,
  genes = 
    counts_raw %>% 
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

mds <- plotMDS(cpm(rna_data, log = TRUE), plot = FALSE)
tibble(
  mds_1 = mds$x,
  mds_2 = mds$y,
  label = rownames(mds$distance.matrix.squared)
) %>% 
  ggplot(aes(x = mds_1, y = mds_2)) +
  geom_point() +
  geom_text_repel(aes(label = label), size = BASE_TEXT_SIZE_MM) +
  xlab(
    str_glue("Leading logFC dim 1 ({round(mds$var.explained[1] * 100, 1)} %)")
  ) +
  ylab(
    str_glue("Leading logFC dim 2 ({round(mds$var.explained[2] * 100, 1)} %)")
  ) +
  coord_fixed() +
  theme_pub()
ggsave_default("rnaseq/human_mds", width = 80)

design <- model.matrix(~0 + condition + patient, data = rna_data$samples)
colnames(design) <- str_replace(colnames(design), "condition", "")
design

contrasts <- makeContrasts(
  TGFb_vs_healthy = TGFb - healthy,
  Kat5i_vs_TGFb = Kat5i - TGFb,
  Kat5i_vs_healthy = Kat5i - healthy,
  levels = colnames(design)
)

efit <-
  voom(rna_data, design, plot = TRUE) %>% 
  lmFit(design) %>% 
  contrasts.fit(contrasts) %>% 
  eBayes()

plotSA(efit)

dge <- 
  colnames(contrasts) %>% 
  set_names() %>% 
  map(\(c)
    topTable(efit, coef = c, number = Inf) %>%
      filter(gene_type == "protein_coding") %>% 
      as_tibble() %>% 
      select(gene = gene_name, logFC, p = P.Value, p_adj = adj.P.Val)
  ) %>% 
  list_rbind(names_to = "comparison") %>% 
  mutate(comparison = factor(comparison) %>% fct_relevel(colnames(contrasts)))
  
dge
dge %>% save_table("rnaseq_human_dge")



# Analyze results ---------------------------------------------------------

dge %>% filter(p_adj <= 0.5)
dge %>% filter(p_adj <= 0.5) %>% count(comparison)

ggplot(dge, aes(logFC, -log10(p))) +
  geom_point(alpha = .25, size = 0.1) +
  facet_wrap(vars(comparison)) +
  theme_pub()
ggsave_default("rnaseq/human_volcano", width = 120, height = 40)

# top 20 genes with pos/neg logFC
top_genes <- c(
  dge %>%
    filter(comparison == "TGFb_vs_healthy") %>%
    slice_max(logFC, n = 20, with_ties = FALSE) %>%
    pull(gene),
  dge %>%
    filter(comparison == "TGFb_vs_healthy") %>%
    slice_min(logFC, n = 20, with_ties = FALSE) %>%
    pull(gene)
)

selected_samples <- 
  samples %>% 
  filter(condition != "Kat5_inhibitor1") %>% 
  pull(sample)

mat <- 
  cpm(rna_data, log = TRUE) %>%
  magrittr::set_rownames(rna_data$genes$gene_name) %>% 
  t() %>% 
  scale() %>% 
  t() %>% 
  magrittr::extract(top_genes, selected_samples)


(p <- Heatmap(
  mat,
  cluster_rows = FALSE,
  width = ncol(mat) * unit(2, "mm"),
  height = nrow(mat) * unit(2, "mm"),
  row_split = rep(c("up", "down"), each = 20)
))
ggsave_default("rnaseq/human_dge_heatmap", plot = p, width = 100)



# Perform GSEA ------------------------------------------------------------

enrichr_genesets <- read_rds("data_generated/enrichr_genesets_human.rds")
enrichr_genesets$fibroblast_markers <-
  read_rds("data_generated/fibroblast_markers.rds") %>% 
  summarise(.by = c(ref, cluster), genes = list(gene_human)) %>% 
  unite(ref, cluster, col = db) %>% 
  deframe()

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
      "KEGG_2019_Mouse",
      "fibroblast_markers"
    )
  ) %>% 
  pmap(run_gsea, .progress = TRUE) %>% 
  list_rbind() %>%
  mutate(comparison = factor(comparison) %>% fct_relevel(colnames(contrasts)))

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
) %>% 
  bind_rows(
    tibble(
      db = "fibroblast_markers",
      pathway = names(enrichr_genesets$fibroblast_markers),
      term_display = str_c(pathway, " (markers)")
    )
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
ggsave_default("rnaseq/human_gsea", type = "pdf", width = 70)



# Correlation heatmap -----------------------------------------------------

corr_mat <-
  limma::removeBatchEffect(
    cpm(rna_data$counts, log = TRUE),
    batch = rna_data$samples$patient,
    group = rna_data$samples$group
  ) %>%
  cor()

distance <- as.dist(1 - corr_mat)

(p <- Heatmap(
  corr_mat,
  col = circlize::colorRamp2(
    seq(min(corr_mat), max(corr_mat), length.out = 9),
    color("davos", reverse = TRUE)(9),
  ),
  
  name = "correlation of\nbatch-corrected counts",
  heatmap_legend_param = list(
    at = round(c(min(corr_mat), max(corr_mat)), 2),
    border = FALSE,
    grid_width = unit(2, "mm"),
    legend_height = unit(15, "mm")
  ),

  clustering_distance_rows = distance,
  clustering_distance_columns = distance,
  row_dend_gp = gpar(lwd = 0.5),
  row_title = "samples",
  row_title_side = "right",

  width = unit(20, "mm"),
  height = unit(20, "mm"),

  show_column_dend = FALSE,
  show_column_names = FALSE
))
ggsave_default("rnaseq/human_cor_heatmap", plot = p)

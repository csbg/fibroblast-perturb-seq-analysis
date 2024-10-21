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

gene_map <-
  read_tsv("metadata/biomart_human_mouse_20210727.tsv") %>% 
  select(gene_human = `Human gene name`, gene_mouse = `Gene name`) %>% 
  distinct()

counts_raw <- read_tsv("data_raw/rna-seq/DataRaw_human.txt")
samples <-
  read_csv("data_raw/rna-seq/sample_data.csv", comment = "#") %>% 
  filter(organism == "human")


rna_data_unfiltered <- DGEList(
  counts =
    counts_raw %>%
    select(gene_name, all_of(samples$sample)),
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
  left_join(samples, by = join_by(label == sample)) %>% 
  mutate(sample_id = str_sub(label, -3)) %>% 
  ggplot(aes(x = mds_1, y = mds_2)) +
  geom_point(aes(color = condition, shape = patient)) +
  geom_text_repel(
    aes(label = sample_id),
    size = BASE_TEXT_SIZE_MM,
    segment.size = BASE_LINEWIDTH,
    seed = 1
  ) +
  xlab(
    str_glue("Leading logFC dim 1 ({round(mds$var.explained[1] * 100, 1)} %)")
  ) +
  ylab(
    str_glue("Leading logFC dim 2 ({round(mds$var.explained[2] * 100, 1)} %)")
  ) +
  coord_fixed() +
  theme_pub() +
  theme(
    legend.key.height = unit(2, "mm"),
    legend.key.width = unit(2, "mm")
  )
ggsave_default("rnaseq/human_mds", type = "pdf", width = 80)

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

dge %>% filter(p_adj <= 0.1)
dge %>% filter(p_adj <= 0.1) %>% count(comparison)

ggplot(dge, aes(logFC, -log10(p))) +
  geom_point(alpha = .25, size = 0.1) +
  facet_wrap(vars(comparison)) +
  theme_pub()
ggsave_default("rnaseq/human_volcano", type = "pdf", width = 120, height = 40)

c("IL11", "COMP", "NPR3", "NOX4", "COL7A1", "LRRC15", "TSPAN2", "GAL",
  "ACKR2", "CILP", "ELN", "PDE4B", "DEPP1", "C7", "SYNE3", "APOL6",
  "ANKRD33", "SLC14A1", "SCARA5", "CAMK2B", "CHI3L1", "MYOZ2", "LACC1",
  "ABCA9", "CFB", "ZBTB7C")

plot_volcano <- function() {
  selected_genes <- tribble(
    ~gene_human, ~type,
    "ACTA2", "fibrotic",
    "ABCA9", "antifibrotic",
    "ACKR3", "fibrotic",
    "ANKRD33", "antifibrotic",
    "APOL6", "antifibrotic",
    "C7", "antifibrotic",
    "CAMK2B", "antifibrotic",
    "CFB", "antifibrotic",
    "CHI3L1", "antifibrotic",
    "CILP", "fibrotic",
    "COL1A1", "fibrotic",
    "COL7A1", "fibrotic",
    "COMP", "fibrotic",
    "CTHRC1", "fibrotic",
    "DEPP1", "antifibrotic",
    "ELN", "fibrotic",
    "GAL", "fibrotic",
    "IL11", "fibrotic",
    "LACC1", "antifibrotic",
    "LOX", "fibrotic",
    "LRRC15", "fibrotic",
    "MEOX1", "fibrotic",
    "MYOZ2", "antifibrotic",
    "NOX4", "fibrotic",
    "NPR3", "fibrotic",
    "PDE4B", "fibrotic",
    "PDGFRA", "antifibrotic",
    "POSTN", "fibrotic",
    "SCARA5", "antifibrotic",
    "SLC14A1", "antifibrotic",
    "SYNE3", "antifibrotic",
    "TSPAN2", "fibrotic",
    "ZBTB7C", "antifibrotic",
  )

  plot_data <- 
    dge %>% 
    filter(comparison != "Kat5i_vs_healthy")
  
  plot_data_highlight <- 
    plot_data %>%
    inner_join(selected_genes, by = join_by(gene == gene_human))
  
  ggplot(plot_data, aes(logFC, -log10(p))) +
    geom_point(
      size = 0.1,
      shape = 20,
      color = "gray80"
    ) +
    geom_point(
      data = plot_data_highlight,
      aes(color = type),
      size = 0.5,
      shape = 20,
      show.legend = FALSE
    ) +
    geom_text_repel(
      data = plot_data_highlight,
      aes(label = gene),
      size = BASE_TEXT_SIZE_MM,
      segment.size = BASE_LINEWIDTH,
      min.segment.length = 0.1,
      seed = 1,
      max.overlaps = 15
    ) +
    xlab("log-fold change") +
    ylab("-log10 p-value") +
    scale_color_manual(
      values = c(
        fibrotic = "#363594",
        antifibrotic = "#218b43",
        quiescent = "#218b43",
        x = "black"
      )
    ) +
    facet_wrap(vars(comparison)) +
    theme_pub() +
    theme(panel.grid = element_blank())
}

plot_volcano()
ggsave_default("rnaseq/6e_human_volcano_pub", type = "pdf", width = 80, height = 40)
ggsave_default("rnaseq/6e_human_volcano_pub_large", type = "pdf", width = 160, height = 80)


plot_expression_heatmap <- function(genes, samples = NULL) {
  samples <- samples %||% rownames(rna_data_unfiltered$samples)
  
  mat <- 
    cpm(rna_data_unfiltered, log = TRUE) %>%
    limma::removeBatchEffect(
      batch = rna_data_unfiltered$samples$patient,
      group = rna_data_unfiltered$samples$group
    ) %>%
    magrittr::set_rownames(rna_data_unfiltered$genes$gene_name) %>% 
    t() %>% 
    scale() %>% 
    t() %>% 
    magrittr::extract(genes, samples)
  
  Heatmap(
    mat,
    cluster_rows = FALSE,
    # cluster_columns = FALSE,
    width = ncol(mat) * unit(2, "mm"),
    height = nrow(mat) * unit(2, "mm")
    # row_split = rep(c("up", "down"), each = 20)
  )
}

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

(p <- plot_expression_heatmap(top_genes, selected_samples))
ggsave_default("rnaseq/human_dge_heatmap", plot = p, type = "pdf", width = 100)


# manually selected genes
# manual_genes <- c(
#   # fibrotic
#   "LOX", "CTHRC1", "COL1A1", "ACTA2", "POSTN", "MEOX1",
#   
#   # antifibrotic/quiescent
#   "PDGFRA", "MGP", "SPON2"
# )

gene_map <-
  read_tsv("metadata/biomart_human_mouse_20210727.tsv") %>% 
  select(gene_human = `Human gene name`, gene_mouse = `Gene name`) %>% 
  distinct()

manual_genes <- c("Spon2", "Col6a1", "Cygb", "Vcan",
               # "Igfp4",
               "Sox9", "Cebpb",
               "Sfrp1",
               # "Lgasl3",
               "Prl2c3", "Mmp3", "Prl2c2", "Meox1",
               "Bhlhe41", "Myl9", "Egr2", "Acta2", "Tagln", "Tpm1", "Ltbp2",
               "Ccn2", "Hbegf", "Serpine1", "Pdlim1", "Lox", "Ccn4", "Postn",
               "Tnc", "Loxl2", "Tgfb1", "Pmepa1", "Runx1", "Cthrc1", "Ncam1",
               "Adam12", "Clu", "Mfap4", "Prg4", "Col3a1", "Col1a1", "Col1a2",
               "Mmp14", "Adamts10", "Serpinh1", "Col4a2", "Loxl3", "Mmp23",
               "Junb", "Ddah1", "Col5a1", "Myh9", "Foxc2", "Actn1", "Myl6",
               "Jun", "Prss23", "Adam5", "Palld", "Hspg2")

manual_genes <-
  tibble(gene_mouse = manual_genes) %>% 
  left_join(gene_map) %>% 
  filter(!is.na(gene_human)) %>% 
  pull(gene_human)

(p <- plot_expression_heatmap(manual_genes))
ggsave_default("rnaseq/human_selected_genes_heatmap", plot = p)




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
ggsave_default("rnaseq/human_gsea_nes", type = "pdf")

gsea_results %>% filter(padj <= 0.05)
gsea_results %>% save_table("rnaseq_human_gsea")



# Plot GSEA ---------------------------------------------------------------

selected_terms <- bind_rows(
  tibble(
    db = "fibroblast_markers",
    pathway = names(enrichr_genesets$fibroblast_markers)
  ) %>% 
    mutate(term_display = str_replace(pathway, "(.+?)_(.+)", "\\2 (\\1)")),
  tribble(
    ~db, ~pathway, ~term_display,
    "KEGG_2019_Mouse",
      "JAK-STAT signaling pathway",
      "JAK-STAT signaling pathway (KEGG)",
    "KEGG_2019_Mouse",
      "Hedgehog signaling pathway",
      "Hedgehog signaling pathway (KEGG)",
    "MSigDB_Hallmark_2020",
      "TNF-alpha Signaling via NF-kB",
      "TNF-alpha Signaling via NF-kB (MSigDB)",
    "MSigDB_Hallmark_2020",
      "Inflammatory Response",
      "Inflammatory Response (MSigDB)",
    "MSigDB_Hallmark_2020",
      "Oxidative Phosphorylation",
      "Oxidative Phosphorylation (MSigDB)",
    "MSigDB_Hallmark_2020",
      "G2-M Checkpoint",
      "G2-M Checkpoint (MSigDB)",
    "Reactome_2022",
      "Mitochondrial Translation R-HSA-5368287",
      "Mitochondrial Translation (Reactome)",
    "WikiPathways_2019_Mouse",
      "Retinol metabolism WP1259",
      "Retinol metabolism (WikiPathways)",
    "WikiPathways_2019_Mouse",
      "Striated Muscle Contraction WP216",
      "Striated Muscle Contraction (WikiPathways)",
  )
)


plot_terms <- function(terms) {
  terms <- 
    terms %>% 
    left_join(selected_terms, by = join_by(pathway)) %>% 
    mutate(
      term_display =
        factor(term_display) %>%
        fct_inorder() %>%
        fct_rev()
    )
  
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

terms_for_main_fig <- tribble(
  ~pathway,
  "Amrute_Fib1",
  "Amrute_Fib3",
  "Amrute_Fib5",
  "Chaffin_FB-ZBTB7C",
  "Chaffin_FB-CNTNAP2",
  "Chaffin_Activated fibroblast",
  "Fu_FB0",
  "Fu_FB2",
  "Fu_FB3",
  "Koenig_Fb1 - Baseline",
  "Koenig_Fb6 - TNC",
  "Koenig_Fb8 - THBS4",
  "Koenig_Fb9 - SERPINE1",
  "Kuppe_Fib2",
  "Kuppe_Fib3",
  "Kuppe_Fib4",
)
plot_terms(terms_for_main_fig)
ggsave_default("rnaseq/6f_human_gsea_pub_main_fig", type = "pdf", height = 65)


terms_for_supp_fig <- tribble(
  ~pathway,
  "Amrute_Fib2",
  "Amrute_Fib4",
  "Amrute_Fib6",
  "Amrute_Fib7",
  "Amrute_Fib8",
  "Chaffin_FB-TLL2",
  "Chaffin_FB-X1",
  "Chaffin_FB-PTCHD4",
  "Fu_FB1",
  "Fu_FB4",
  "Fu_FB5",
  "Fu_FB6",
  "Koenig_Epicardium",
  "Koenig_Fb2 - PCOLCE2",
  "Koenig_Fb3 - GPX3",
  "Koenig_Fb4 - PLA2G2A",
  "Koenig_Fb5 - ELN",
  "Koenig_Fb7 - CCL2",
  "Kuppe_Fib1",
  # "JAK-STAT signaling pathway",
  # "Hedgehog signaling pathway",
  # "TNF-alpha Signaling via NF-kB",
  # "Inflammatory Response",
)
plot_terms(terms_for_supp_fig)
ggsave_default("rnaseq/S10g_human_gsea_pub_supp_fig", type = "pdf", width = 51)



# Correlation heatmap -----------------------------------------------------

## Counts ----

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
ggsave_default("rnaseq/human_cor_heatmap", plot = p, width = 150, type = "pdf")


## log-fold changes ----

corr_mat <-
  dge %>% 
  distinct(comparison, gene, .keep_all = TRUE) %>% 
  pivot_wider(names_from = comparison, values_from = logFC, id_cols = gene) %>% 
  mutate(TGFb_vs_Kat5i = -Kat5i_vs_TGFb) %>% 
  select(!c(gene, Kat5i_vs_TGFb)) %>%
  cor()

distance <- as.dist(1 - corr_mat)

(p <- Heatmap(
  corr_mat,
  col = circlize::colorRamp2(
    seq(min(corr_mat), max(corr_mat), length.out = 9),
    color("davos", reverse = TRUE)(9),
  ),
  
  name = "correlation of\nlog-fold changes",
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
ggsave_default("rnaseq/human_cor_heatmap_logfc", plot = p, width = 150, type = "pdf")

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

enrichr_genesets <- read_rds("data_generated/enrichr_genesets_human.rds")
enrichr_genesets$fibroblast_markers <-
  read_rds("data_generated/fibroblast_markers.rds") %>% 
  summarise(.by = c(ref, cluster), genes = list(gene_human)) %>% 
  unite(ref, cluster, col = db) %>% 
  deframe()

counts_raw <- read_tsv("data_raw/rna-seq/raw_counts.txt")
samples <- read_csv("metadata/samples_bulk_rnaseq.csv", comment = "#") 



# Normalize data ----------------------------------------------------------

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


rna_data <- calcNormFactors(rna_data)

mds <- plotMDS(cpm(rna_data, log = TRUE), plot = FALSE)
tibble(
  mds_1 = mds$x,
  mds_2 = mds$y,
  label = rownames(mds$distance.matrix.squared)
) %>%
  left_join(samples, by = join_by(label == sample)) %>% 
  mutate(
    patient = recode(patient, patient1 = "Patient 1", patient2 = "Patient 2")
  ) %>% 
  ggplot(aes(x = mds_1, y = mds_2)) +
  geom_point(aes(color = condition, shape = patient)) +
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
ggsave_default("S10c_bulkrnaseq_mds", type = "pdf", width = 80)



# Perform DGE -------------------------------------------------------------

design <- model.matrix(~0 + condition + patient, data = rna_data$samples)
colnames(design) <- str_replace(colnames(design), "condition", "")

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
  
# dge %>% save_table("rnaseq_human_dge")



# Plot DGE ----------------------------------------------------------------

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
ggsave_default("6d_bulkrnaseq_volcano", type = "pdf", width = 160, height = 80)



# Perform GSEA ------------------------------------------------------------

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

# gsea_results %>% save_table("rnaseq_human_gsea")



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
ggsave_default("6e_bulkrnaseq_signatures", type = "pdf", height = 65)


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
)
plot_terms(terms_for_supp_fig)
ggsave_default("S10e_bulkrnaseq_signatures", type = "pdf", width = 51)



# Correlation heatmap -----------------------------------------------------

corr_mat <-
  limma::removeBatchEffect(
    cpm(rna_data$counts, log = TRUE),
    batch = rna_data$samples$patient,
    group = rna_data$samples$group
  ) %>%
  cor()

distance <- as.dist(1 - corr_mat)


figure_colnames <- c(
    Kat5_pos_1_2015_S15 = "iKat5 Patient 1",
    Kat5_pos_2017_S18 = "iKat5 Patient 2",
    Mock_neg_1_2015_S13 = "Healthy Patient 1",
    Mock_neg_2017_S16 = "Healthy Patient 2",
    Mock_pos_1_2015_S14 = "TGFb Patient 1",
    Mock_pos_2017_S17 = "TGFb Patient 2"
  )
colnames(corr_mat) <- colnames(corr_mat) %>% recode(!!!figure_colnames)
rownames(corr_mat) <- rownames(corr_mat) %>% recode(!!!figure_colnames)
  
  
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
ggsave_default("S10d_bulkrnaseq_expression_correlation",
               plot = p, width = 150, type = "pdf")

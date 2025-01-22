library(SingleCellExperiment)
library(DropletUtils)
library(tidyverse)
library(scuttle)
source("scripts_r/utils.R")



# Load data ---------------------------------------------------------------

col_metadata <-
  read_rds("data_generated/rna.rds") %>% 
  colData() %>%
  as_tibble()

samples <-
  read_csv("metadata/samples_sc_rnaseq.csv") %>% 
  mutate(path = str_glue("data_raw/rna/{geo_id}_{name_orig}_filtered_feature_bc_matrix.h5"))


sce <-
  samples %>% 
  select(name_python, path) %>% 
  deframe() %>% 
  read10xCounts(intersect.genes = TRUE)

rownames(sce) <- rowData(sce)$Symbol
colnames(sce) <- str_c(
  str_sub(colData(sce)$Barcode, end = -2),
  colData(sce)$Sample
)


common_barcodes <- intersect(col_metadata$barcode, colnames(sce))

sce <-
  sce[, common_barcodes] %>% 
  logNormCounts()

col_metadata <-
  col_metadata %>%
  filter(barcode %in% common_barcodes)



# Plot gene ---------------------------------------------------------------

add_count_cols <- function(df, genes) {
  count_cols <- 
    logcounts(sce)[genes, , drop = FALSE] %>%
    t() %>%
    as.matrix() %>%
    as_tibble() %>% 
    rename_with(\(n) paste0("expression_", n))
  
  bind_cols(df, count_cols)
}


plot_data <- 
  col_metadata %>% 
  add_count_cols("Fap")


plot_gene <- function(col) {
  ggplot(plot_data, aes({{col}}, expression_Fap)) +
    geom_violin(
      color = "#74a9cf",
      fill = "#74a9cf",
      scale = "width",
      width = 0.8
    ) +
    stat_summary(geom = "point", fun = mean, size = .2) +
    theme_pub(rotate_x_labels = TRUE) + 
    theme(panel.grid = element_blank())
}


plot_gene(condition)
ggsave_default("R15a_expression_Fap_condition",
               type = "pdf", width = 20, height = 50)

plot_gene(cell_type)
ggsave_default("R15b_expression_Fap_celltype",
               type = "pdf", width = 50, height = 80)

plot_gene(target)
ggsave_default("R15c_expression_Fap_target",
               type = "pdf", width = 100, height = 80)



# Plot expression percentage ----------------------------------------------

plot_data_pct <- 
  plot_data %>% 
  summarise(
    .by = c(cell_type, condition),
    n_fap = sum(expression_Fap > 0),
    n = n()
  ) %>%
  mutate(
    pct_fap = n_fap / n,
    high_fap = pct_fap > 0.05
  )

ggplot(plot_data_pct, aes(cell_type, condition, size = n)) + 
  geom_point(aes(color = pct_fap)) +
  geom_point(
    data = plot_data_pct %>% filter(high_fap),
    shape = 1,
    show.legend = FALSE
  ) +
  xlab("fibroblast subtype") +
  scale_color_gradient(
    "% expressing FAP",
    low = "lightgrey", 
    high = "blue", 
    limits = c(0, NA),
    guide = guide_colorbar(
      barheight = unit(13, "mm"),
      barwidth = unit(2, "mm"),
      ticks = FALSE
    )
  ) + 
  scale_size_area(
    "number of cells",
    max_size = 3
  ) +
  theme_pub(rotate_x_labels = TRUE) +
  theme(
    legend.key.height = unit(3, "mm"),
    legend.key.width = unit(3, "mm"),
    legend.margin = margin(),
    legend.box = "horizontal",
    legend.box.just = "bottom"
  )
ggsave_default("pct_FAP", type = "pdf", width = 85, height = 45)

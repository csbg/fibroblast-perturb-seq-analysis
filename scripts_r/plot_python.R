library(tidyverse)
library(RColorBrewer)
library(SingleCellExperiment)
library(conflicted)
library(ComplexHeatmap)
library(viridisLite)
library(ggrepel)
source("scripts_r/utils.R")

conflicts_prefer(
  dplyr::count,
  dplyr::filter,
  dplyr::rename
)

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

col_data <- 
  read_rds("data_generated/rna.rds") %>% 
  colData() %>% 
  as_tibble()

dge <- read_rds("data_generated/dge.rds")



# 1c ----------------------------------------------------------------------

plot_data <- 
  read_csv("data_generated/python_figure_data/1c_data.csv") %>% 
  rename(guide = ...1) %>% 
    mutate(viability_score = Cas_vs_Wt_lfc, fibrotic_score = Act_vs_Quies_lfc) %>% 
  {.}
  
ggplot(plot_data, aes(viability_score, fibrotic_score)) +
  geom_point(
    data = plot_data %>% filter(Significance == "Not significant"),
    shape = 19,
    color = "gray80"
  ) +
  geom_point(
    data = plot_data %>% filter(Significance != "Not significant"),
    aes(color = Significance),
    shape = 19
  ) +
  geom_text_repel(
    data = plot_data %>% filter(Significance != "Not significant"),
    aes(label = guide),
    size = BASE_TEXT_SIZE_MM
  ) +
  theme_pub()
ggsave_default("python/1c", type = "pdf", width = 150, height = 120)

plot_data %>% 
  select(guide, viability_score, fibrotic_score, significance = Significance) %>% 
  save_table("source_data/1c")



# 1d ----------------------------------------------------------------------

plot_data <- read_csv("data_generated/python_figure_data/1d_data.csv")

plot_data %>% 
  mutate(Gene = as_factor(Gene) %>% fct_rev()) %>%
  pivot_longer(!Gene, names_to = "type", values_to = "value") %>% 
  ggplot(aes(type, Gene, color = value)) +
  geom_point(
    shape = 19,
    size = 2
  ) +
  scale_color_gradientn(
    colors = brewer.pal(9, "Spectral"),
    limits = c(-2.5, 2.5),
    oob = scales::oob_squish
  ) +
  theme_pub() +
  theme(panel.grid = element_blank())
ggsave_default("python/1d", type = "pdf", width = 50)

plot_data %>% 
  rename(fibrotic_score = Act_vs_Quies_lfc, viability_score = Cas_vs_Wt_lfc) %>% 
  save_table("source_data/1d_S2a")



# 2b, 2c, S2b -------------------------------------------------------------

ggplot(col_data, aes(x = UMAP_X, y = UMAP_Y)) +
  geom_point(
    aes(color = cell_type), # 2b: cell_type; 2c: condition; S2b: sample
    shape = 19,
    size = 0.01
  ) +
  scale_color_brewer(
    palette = "Paired",
    guide = guide_legend(override.aes = list(size = 2))
  ) +
  coord_fixed() +
  theme_pub() +
  theme(panel.grid = element_blank())
ggsave_default("python/2b", width = 150)

col_data %>% 
  select(barcode, UMAP_X, UMAP_Y, cell_type, condition, sample) %>% 
  save_table("source_data/2b_2c_S2b")



# 2d ----------------------------------------------------------------------

col_data %>% 
  ggplot(aes(condition, fill = cell_type)) +
  geom_bar(position =  "fill") +
  scale_fill_brewer(palette = "Paired")
  coord_flip() +
  theme_pub()
 
col_data %>% 
  count(condition, cell_type) %>% 
  mutate(
    .by = condition,
    frac = n / sum(n)
  ) %>% 
  save_table("source_data/2d")



# 2e ---------------------------------------------------------------------

markers <-
  read_csv("data_generated/python_figure_data/2e_markers.txt") %>% 
  pull(gene)

cell_types <- c(
  "Quiescent",
  "Transitory",
  "Stress-fiber fibroblasts",
  "Myofibroblasts",
  "Deactivated myofibroblasts",
  "Inflammatory (injury response)",
  "Inflammatory & fibrotic",
  "Tissue-repair",
  "Lamp1-fibroblasts",
  "Proliferative",
  "Proliferative myofibroblast"
)

plot_data <- 
  left_join(
    read_csv("data_generated/python_figure_data/2e_dot_color.csv") %>% 
      pivot_longer(!...1, names_to = "gene", values_to = "mean_exp"),
    read_csv("data_generated/python_figure_data/2e_dot_size.csv") %>% 
      pivot_longer(!...1, names_to = "gene", values_to = "frac"),
    by = c("...1", "gene")
  ) %>% 
  rename(cell_type = ...1) %>% 
  mutate(
    gene = fct_relevel(gene, !!!markers) %>% fct_rev(),
    cell_type = fct_relevel(cell_type, !!!cell_types) %>% fct_rev()
  )

ggplot(plot_data, aes(cell_type, gene)) +
  geom_point(aes(color = mean_exp, size = frac)) +
  # scale_color_distiller(palette = "YlGnBu", direction = 1) +
  scale_color_gradientn(colors = brewer.pal(9, "Reds")) +
  scale_radius(
    name = "-log10 padj",
    # max_size = 2,
    limits = c(0, 1),
    breaks = c(0, .5, 1),
    labels = c("0", ".5", "1"),
    range = c(0, 2)
    # oob = scales::oob_squish
  ) +
  coord_fixed() +
  theme_pub(rotate_x_labels = TRUE)

ggsave_default("python/2e", type = "pdf", width = 60)

plot_data %>% 
  save_table("source_data/2e")



# 3a ----------------------------------------------------------------------

selected_genes <-
  read_csv("data_generated/python_figure_data/3a_markers.txt") %>% 
  pull(gene)

selected_guides <- c(
  "Tgfb vs unstimulated",
  "Tgfbr1",
  "Smad3",
  "Wdr82",
  "Hcfc1",
  "Kat8",
  "Srcap",
  "Dmap1",
  "Kat5",
  "Paxip1"
)

plot_data <- 
  bind_rows(
    dge %>% 
      filter(
        comparison == "Tgfb1_Target_vs_NTC",
        group %in% selected_guides,
      ),
    dge %>%
      filter(
        comparison == "NTC_Stim_vs_Resting",
        group == "Tgfb1"
      ) %>% 
      mutate(group = "Tgfb vs unstimulated")
  ) %>% 
  filter(gene %in% selected_genes) %>% 
  mutate(
    group = fct_relevel(group, selected_guides) %>% fct_rev(),
    gene = fct_relevel(gene, selected_genes),
  ) %>% 
  select(!comparison)

plot_data %>% 
  mutate(p_adj = pmax(p_adj, 10^-10)) %>% 
  ggplot(aes(gene, group)) +
  geom_point(aes(color = logFC, size = -log10(p_adj))) +
  scale_color_gradientn(
    colors = rev(brewer.pal(9, "RdBu")),
    limits = c(-1.5, 1.5),
    oob = scales::oob_squish
  ) +
  scale_radius(
    name = "-log10 padj",
    # max_size = 2,
    limits = c(0, 10),
    range = c(0, 2.5),
    # oob = scales::oob_squish
  ) +
  coord_fixed() +
  theme_pub(rotate_x_labels = TRUE) +
  theme(panel.grid = element_blank())
ggsave_default("python/3a", type = "pdf", width = 180)

plot_data %>% 
  save_table("source_data/3a")



# S2c ---------------------------------------------------------------------

markers <-
  read_csv("data_generated/python_figure_data/S2c_markers.txt") %>% 
  pull(gene)

cell_types <- c(
  "Quiescent",
  "Transitory",
  "Stress-fiber fibroblasts",
  "Myofibroblasts",
  "Deactivated myofibroblasts",
  "Inflammatory (injury response)",
  "Inflammatory & fibrotic",
  "Tissue-repair",
  "Lamp1-fibroblasts",
  "Proliferative",
  "Proliferative myofibroblast"
)

plot_data <- 
  left_join(
    read_csv("data_generated/python_figure_data/S2c_dot_color.csv") %>% 
      pivot_longer(!...1, names_to = "gene", values_to = "mean_exp"),
    read_csv("data_generated/python_figure_data/S2c_dot_size.csv") %>% 
      pivot_longer(!...1, names_to = "gene", values_to = "frac"),
    by = c("...1", "gene")
  ) %>% 
  rename(cell_type = ...1) %>% 
  mutate(
    gene = fct_relevel(gene, !!!markers) %>% fct_rev(),
    cell_type = fct_relevel(cell_type, !!!cell_types)
  )

ggplot(plot_data, aes(cell_type, gene)) +
  geom_point(aes(color = mean_exp, size = frac)) +
  # scale_color_distiller(palette = "YlGnBu", direction = 1) +
  scale_color_gradientn(colors = brewer.pal(9, "Reds")) +
  scale_radius(
    name = "-log10 padj",
    # max_size = 2,
    limits = c(0, 1),
    breaks = c(0, .5, 1),
    labels = c("0", ".5", "1"),
    range = c(0, 2)
    # oob = scales::oob_squish
  ) +
  coord_fixed() +
  theme_pub(rotate_x_labels = TRUE)

ggsave_default("python/S2c", type = "pdf")

plot_data %>% 
  save_table("source_data/S2c")



# S3 ----------------------------------------------------------------------

read_csv("data_generated/python_figure_data/S3_data.csv") %>% 
  rename(cell_type = ...1) %>% 
  save_table("source_data/S3")



# S4a ---------------------------------------------------------------------

genotypes <- 
  dge %>% 
  filter(comparison != "NTC_Stim_vs_Resting") %>% 
  pull(group) %>%
  unique()

plot_data <- 
  dge %>% 
  filter(
    comparison == "Il1b_Target_vs_NTC",
    gene %in% genotypes,
  ) %>% 
  select(!comparison)

plot_data %>% 
  filter(p_adj < .1) %>% 
  mutate(p_adj = pmax(p_adj, 10^-10)) %>% 
  ggplot(aes(gene, group)) +
  geom_point(aes(color = logFC, size = -log10(p_adj))) +
  scale_color_gradientn(
    colors = rev(brewer.pal(9, "RdBu")),
    limits = c(-1.5, 1.5),
    oob = scales::oob_squish
  ) +
  scale_radius(
    name = "-log10 padj",
    # max_size = 2,
    limits = c(0, 10),
    range = c(0, 2.5),
    # oob = scales::oob_squish
  ) +
  coord_fixed() +
  theme_pub(rotate_x_labels = TRUE) +
  theme(panel.grid = element_blank())
ggsave_default("python/S4a", type = "pdf", width = 120)

plot_data %>% 
  save_table("source_data/S4a")



# S4b ---------------------------------------------------------------------

guides <- 
  read_csv(
    "data_generated/python_figure_data/S4b_labels.txt",
    col_names = "guide"
  ) %>% 
  pull(guide)

corr_mat_df <- 
  read_csv(
    "data_generated/python_figure_data/S4b_matrix.txt",
    col_names = guides
  ) %>% 
  mutate(
    .before = 1,
    guide = guides
  )

corr_mat <- 
  column_to_rownames("guide") %>% 
  as.matrix()

distance <- as.dist(1 - corr_mat)

Heatmap(
  corr_mat,
  col = circlize::colorRamp2(
    seq(min(corr_mat), max(corr_mat), length.out = 9),
    magma(9, direction = -1),
  ),
  
  name = "corr",
  heatmap_legend_param = list(
    at = round(c(min(corr_mat), max(corr_mat)), 2),
    border = FALSE,
    grid_width = unit(2, "mm"),
    legend_height = unit(15, "mm")
  ),
  
  clustering_distance_rows = distance,
  clustering_distance_columns = distance,
  row_dend_gp = gpar(lwd = 0.5),
  row_title = "gRNAs",
  row_title_side = "right",
  
  width = unit(2, "mm") * nrow(corr_mat),
  height = unit(2, "mm") * nrow(corr_mat),
  
  show_column_dend = FALSE,
  show_column_names = FALSE
)


corr_mat_df %>% 
  save_table("source_data/S4b")

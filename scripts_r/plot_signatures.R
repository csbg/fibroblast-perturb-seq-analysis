library(SingleCellExperiment)
library(readxl)
library(tidyverse)
library(ComplexHeatmap)
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

sce <- read_rds("data_generated/rna.rds")

markers <-
  read_rds("data_generated/fibroblast_markers.rds") %>% 
  select(!gene_human) %>% 
  rename(gene = gene_mouse)



# Calculate signatures ----------------------------------------------------

calculate_signatures <- function(markers) {
  # remove genes that are absent from the count matrix
  markers <- 
    markers %>% 
    filter(gene %in% rownames(sce))
  
  # matrix containing only marker genes; values are normalized between 0 and 1
  mat <- logcounts(sce)[unique(markers$gene), ]
  mat <- mat - rowMins(mat)
  mat <- mat / rowMaxs(mat)
  
  
  markers %>% 
    group_by(ref, cluster) %>% 
    group_map(\(df, keys) {
      info("Markers for {keys$ref}, {keys$cluster}")
      mat[df$gene, ] %>%
        colMeans() %>% 
        enframe("barcode", str_c("signature_", keys$ref, "_", keys$cluster))
    }) %>% 
    purrr::reduce(left_join, by = "barcode")
}



# selected_terms <- tribble(
#   ~ref, ~cluster,
#   "Amrute", "Fib1", 
#   "Amrute", "Fib3", 
#   "Amrute", "Fib7", 
#   "Chaffin", "Activated fibroblast", 
#   "Fu", "FB0", 
#   "Fu", "FB5", 
#   "Koenig", "Fb5 - ELN", 
#   "Koenig", "Fb7 - CCL2", 
#   "Koenig", "Fb8 - THBS4", 
#   "Kuppe", "Fib1",
#   "Kuppe", "Fib2",
#   "exvivo", "resting",
#   "exvivo", "fibrotic",
#   "exvivo", "inflammatory"
# )

ignored_terms <- tribble(
  ~ref, ~cluster,
  "Amrute", "Fib4",
)


# signatures <- calculate_signatures(markers)
# signatures <- calculate_signatures(markers %>% semi_join(selected_terms))
signatures <- calculate_signatures(markers %>% anti_join(ignored_terms))


# Plots -------------------------------------------------------------------

## UMAP ----

colData(sce) %>%
  as_tibble() %>% 
  ggplot(aes(UMAP_X, UMAP_Y)) +
  geom_point(aes(color = cell_type), size = 0.01) +
  coord_fixed() +
  theme_pub() +
  theme(
    legend.key.height = unit(2, "mm"),
    legend.key.width = unit(2, "mm"),
    panel.grid = element_blank()
  )
# ggsave_default("umap_cell_types", height = 120, width = 120)



## Signatures ----

plot_data_sig <- 
  colData(sce) %>% 
  as_tibble() %>% 
  left_join(signatures, by = "barcode") %>% 
  mutate(across(starts_with("signature"), ~ scale(.x)[,1]))



plot_signature <- function(ref, signatures = NULL, ncol = NULL) {
  col_prefix <- str_c("signature_", ref, "_")
  
  plot_data <- 
    plot_data_sig %>%
    select(barcode, UMAP_X, UMAP_Y, starts_with(col_prefix)) %>% 
    pivot_longer(
      starts_with(col_prefix),
      names_to = "signature",
      values_to = "value"
    ) %>%
    mutate(signature = str_remove(signature, col_prefix))

  if (!is.null(signatures)) {
    plot_data <- 
      plot_data %>% 
      filter(signature %in% {{signatures}})
  }
  
  ggplot(plot_data, aes(UMAP_X, UMAP_Y)) +
    stat_summary_hex(aes(z = pmin(3, value)), fun = mean, bins = 100) +
    scale_fill_gradient2(
      name = "Mean signature",
      low = "blue",
      midpoint = 0,
      high = "red",
      guide = guide_colorbar(
        barheight = unit(2, "mm"),
        barwidth = unit(15, "mm"),
        label.position = "top",
        title.vjust = 0.1,
        ticks = FALSE
      ),
    ) +
    # coord_fixed() +
    facet_wrap(vars(signature), ncol = ncol) +
    theme_pub() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      # legend.direction = "vertical",
      # legend.key.height = unit(3, "mm"),
      # legend.key.width = unit(3, "mm"),
      legend.position = "bottom",
      legend.margin = margin(),
      panel.grid = element_blank()
    )
}

plot_signature(
  "Forte",
  signatures = c("Homeostatic Epicardial Derived", "Myofibroblasts"),
  ncol = 1
)
ggsave_default("3g_signature_umaps", type = "pdf", width = 40, height = 70)


plot_signature("Buechler")
# ggsave_default("umap_signature_Buechler", width = 150, type = "pdf")

plot_signature("Koenig")
# ggsave_default("umap_signature_Koenig", width = 150, type = "pdf")

plot_signature("Forte")
# ggsave_default("umap_signature_Forte", width = 150, type = "pdf")



## Signature summary ----

cell_type_order <- c(
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

plot_signature_heatmap <- function(data,
                                   ref = NULL,
                                   color_limit = NULL,
                                   row_order = NULL) {
  col_prefix <- str_c("signature", ref, "", sep = "_")
  
  mat <- 
    data %>% 
    summarise(
      .by = cell_type,
      across(starts_with(col_prefix), ~mean(.x, na.rm = TRUE))
    ) %>% 
    column_to_rownames("cell_type") %>% 
    as.matrix()
  
  colnames(mat) <-
    colnames(mat) %>%
    str_remove(col_prefix)
  
  if (!is.null(row_order))
    mat <- mat[row_order, ]
    
  if (is.null(color_limit))
    color_limit <- max(range(mat))
  
  Heatmap(
    mat,
    col = circlize::colorRamp2(
      c(-color_limit, 0, color_limit),
      c("blue", "white", "red"),
    ),
    name = "Mean signature",
    heatmap_legend_param = list(
      at = round(c(-color_limit, 0, color_limit), 2),
      direction = "horizontal",
      grid_height = unit(2, "mm"),
      legend_width = unit(10, "mm")
    ),
    
    row_names_side = "left",
    row_title = "cell type",
    row_title_side = "left",
    cluster_rows = FALSE,

    column_title = glue::glue("Population signatures ({ref})",
                              .null = "all references"),
    column_title_side = "bottom",

    width = ncol(mat) * unit(2, "mm"),
    height = nrow(mat) * unit(2, "mm"),
  )
}

(p <- plot_signature_heatmap(plot_data_sig, "Forte",
                             color_limit = 2, row_order = cell_type_order))
ggsave_default("3f_signatures_Forte", plot = p)

(p <- plot_signature_heatmap(plot_data_sig, "Buechler"))
ggsave_default("S3_signatures_Buechler", plot = p)

(p <- plot_signature_heatmap(plot_data_sig, "Koenig"))
ggsave_default("S3_signatures_Koenig", plot = p)


(p <-
    plot_data_sig %>%
    select(!contains(c("Buechler", "Forte"))) %>%
    plot_signature_heatmap(color_limit = 2, row_order = cell_type_order))
ggsave_default("XX_signatures_all", plot = p)

walk(
  unique(markers$ref),
  \(ref) {
    p <- plot_signature_heatmap(plot_data_sig, ref = ref, color_limit = 2,
                                row_order = cell_type_order)
    ggsave_default(str_glue("XX_signatures_{ref}"), plot = p)
  }
)



## Selected signatures ----

selected_signatures <- tribble(
  ~signature,                     ~type,
  "Amrute_Fib7",                  "heart failure",
  "Koenig_Fb7 - CCL2",            "heart failure",
  "Chaffin_Activated fibroblast", "heart failure",
  "Amrute_Fib3",                  "heart failure",
  "Koenig_Fb6 - TNC",             "heart failure",
  "Fu_FB5",                       "heart failure",
  "Fu_FB0",                       "heart failure",
  "Koenig_Fb8 - THBS4",           "heart failure",
  "Buechler_Lrrc15+",             "heart failure",
  
  "Chaffin_FB-CNTNAP2",           "healthy",
  "Amrute_Fib1",                  "healthy",
  "Koenig_Fb3 - GPX3",            "healthy",
  "Chaffin_FB-ZBTB7C",            "healthy",
  
  "Koenig_Fb1 - Baseline",        "steady state",
  "Buechler_Col15a1+",            "steady state",
  "Buechler_Ccl19+",              "steady state",
) %>%
  mutate(type = fct_inorder(type))


plot_selected_signatures <- function(data,
                                     selected_signatures,
                                     color_limit = NULL,
                                     row_order = NULL) {
  data <-
    data %>%
    select(cell_type, contains(selected_signatures$signature))
  
  col_prefix <- "signature_"
  
  mat <-
    data %>%
    summarise(
      .by = cell_type,
      across(starts_with(col_prefix), ~mean(.x, na.rm = TRUE))
    ) %>%
    column_to_rownames("cell_type") %>%
    as.matrix()
  
  colnames(mat) <-
    colnames(mat) %>%
    str_remove(col_prefix)
  mat <- mat[, selected_signatures$signature]
  
  if (!is.null(row_order))
    mat <- mat[row_order, ]
  
  if (is.null(color_limit))
    color_limit <- max(range(mat))
  
  Heatmap(
    mat,
    col = circlize::colorRamp2(
      c(-color_limit, 0, color_limit),
      c("blue", "white", "red"),
    ),
    name = "Mean signature",
    heatmap_legend_param = list(
      at = round(c(-color_limit, 0, color_limit), 2),
      direction = "horizontal",
      grid_height = unit(2, "mm"),
      legend_width = unit(10, "mm")
    ),
    
    row_names_side = "left",
    row_title = "cell type",
    row_title_side = "left",
    
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    cluster_column_slices = FALSE,
    column_split = selected_signatures$type,
    show_parent_dend_line = FALSE,
    
    width = ncol(mat) * unit(2, "mm") + unit(2, "mm"),
    height = nrow(mat) * unit(2, "mm"),
    
    top_annotation = HeatmapAnnotation(
      type = selected_signatures$type,
      col = list(type = c(
        "heart failure" = "red",
        "healthy" = "green",
        "steady state" = "gray"
      )),
      show_annotation_name = FALSE,
      show_legend = FALSE,
      annotation_legend_param = list(
        type = list(
          title = "type",
          grid_width = unit(2, "mm"),
          labels_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
          title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT)
        )
      )
    )
  ) %>%
    draw(
      column_title = "Population signatures",
      column_title_side = "bottom",
      column_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT)
    )
}


(p <- plot_selected_signatures(plot_data_sig,
                               selected_signatures,
                               # color_limit = 2,
                               row_order = cell_type_order))
ggsave_default("XX_signatures_selected_no_clustering", plot = p,
               type = "pdf", width = 150)

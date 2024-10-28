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


ignored_terms <- tribble(
  ~ref, ~cluster,
  "Amrute", "Fib4",
)

signatures <-
  markers %>%
  anti_join(ignored_terms, by = join_by(ref, cluster)) %>% 
  calculate_signatures()



# Plots -------------------------------------------------------------------

plot_data_sig <- 
  colData(sce) %>% 
  as_tibble() %>% 
  left_join(signatures, by = "barcode") %>% 
  mutate(across(starts_with("signature"), ~ scale(.x)[,1]))


## Signature hexplots ----

plot_signature_hexplot <- function(ref, signatures = NULL, ncol = NULL) {
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
      legend.position = "bottom",
      legend.margin = margin(),
      panel.grid = element_blank()
    )
}

# plot_signature_hexplot(
#   "Forte",
#   signatures = c("Homeostatic Epicardial Derived", "Myofibroblasts"),
#   ncol = 1
# )
# 
# plot_signature_hexplot("Forte")
# plot_signature_hexplot("Buechler")


## Murine signature heatmaps ----

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

plot_murine_signatures <- function(data,
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
      direction = "vertical",
      grid_height = unit(2, "mm"),
      legend_width = unit(10, "mm")
    ),
    
    row_names_side = "right",
    row_title = "cell type",
    row_title_side = "right",
    cluster_rows = TRUE,

    column_title = glue::glue("Population signatures ({ref})",
                              .null = "all references"),
    column_title_side = "bottom",

    width = ncol(mat) * unit(2, "mm"),
    height = nrow(mat) * unit(2, "mm"),
  )
}

(p <- plot_murine_signatures(plot_data_sig, "Forte",
                             color_limit = 2, row_order = cell_type_order))
ggsave_default("S2c_signatures_Forte", plot = p)

(p <- plot_murine_signatures(plot_data_sig, "Buechler"))
ggsave_default("S2d_signatures_Buechler", plot = p)


## Human signature heatmaps ----

selected_signatures <- tribble(
  ~signature,                     ~type,           ~display,
  "Fu_FB5",                       "heart failure", "Fu-FB5-ACTA2+",
  "Chaffin_Activated fibroblast", "heart failure", "Chaffin-Activated-Fib",
  "Kuppe_Fib2",                   "heart failure", "Kuppe-Fib2-POSTN+TNC+",
  "Amrute_Fib3",                  "heart failure", "Amrute-Fib3-POSTN+",
  "Koenig_Fb6 - TNC",             "heart failure", "Koenig-Fb6-TNC+",
  "Koenig_Fb9 - SERPINE1",        "heart failure", "Koenig-Fb9-SERPINE1+",
  
  "Fu_FB1",                       "healthy",       "Fu-FB1-PLA2G2A+",
  "Fu_FB3",                       "healthy",       "Fu-FB3-PTGDS+",
  "Chaffin_FB-ZBTB7C",            "healthy",       "Chaffin-FB-ZBTB7C+",
  "Chaffin_FB-CNTNAP2",           "healthy",       "Chaffin-FB-CNTNAP2+",
  "Kuppe_Fib3",                   "healthy",       "Kuppe-Fib3",
  "Amrute_Fib1",                  "healthy",       "Amrute-Fib1-SCN7A+",
  "Koenig_Fb1 - Baseline",        "healthy",       "Koenig-Fb1-Basline",
  "Koenig_Fb3 - GPX3",            "healthy",       "Koenig-Fb3-GPX3+",
  
  "Amrute_Fib5",                  "inflammatory",  "Amrute-Fib5",
) %>%
  mutate(
    type = fct_inorder(type),
    display = coalesce(display, signature)
  )


plot_human_signatures <- function(data,
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
  
  colnames(mat) <- 
    selected_signatures %>% 
    select(signature, display) %>% 
    deframe() %>% 
    magrittr::extract(colnames(mat))
  
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
    cluster_columns = TRUE,
    cluster_column_slices = TRUE,
    column_split = selected_signatures$type,
    show_parent_dend_line = FALSE,
    
    width = ncol(mat) * unit(2, "mm") + unit(2, "mm"),
    height = nrow(mat) * unit(2, "mm"),
    
    top_annotation = HeatmapAnnotation(
      type = selected_signatures$type,
      col = list(type = c(
        "heart failure" = "red",
        "healthy" = "green",
        "steady state" = "gray",
        "inflammatory" = "blue"
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


(p <- plot_human_signatures(plot_data_sig,
                            selected_signatures,
                            row_order = cell_type_order))
ggsave_default("6b_signatures_human", plot = p, type = "pdf", width = 150)

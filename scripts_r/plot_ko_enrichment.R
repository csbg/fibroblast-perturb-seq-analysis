library(SingleCellExperiment)
library(tidyverse)
library(RColorBrewer)
library(patchwork)
source("scripts_r/utils.R")



# Load data ---------------------------------------------------------------

col_data <- 
  read_rds("data_generated/rna.rds") %>% 
  colData() %>% 
  as_tibble()



# Perform Fisher's exact test ---------------------------------------------

# all combinations of features for the test
comparisons <- expand_grid(
  condition = unique(col_data$condition),
  cell_type = unique(col_data$cell_type),
  target = unique(col_data$target) %>% setdiff("non-targeting"),
  ntc = c("NTC_0005", "NTC_0006"),
)

# how to perform Fisher's exact test:
# 1. filter cells from one condition, wich a given KO or NTC
# 2. calculate contingency table by counting which of these cells are in a
#    given cluster (or not) and whether they contain the NTC (or not)
run_test <- function(condition, cell_type, target, ntc) {
  mat <- 
    col_data %>%  
    filter(
      condition == {{condition}},
      sgRNA == {{ntc}} | target == {{target}}
    ) %>%
    mutate(
      is_cluster = cell_type == {{cell_type}},
      is_target = target == {{target}}
    ) %>%
    with(table(is_cluster, is_target))
  
  if (all(dim(mat) == 2))
    res <- fisher.test(mat)
  else
    res <- list(estimate = 0, p.value = 1)
  
  tibble(
    condition = {{condition}},
    cell_type = {{cell_type}},
    target = {{target}},
    ntc = {{ntc}},
    odds_ratio = res$estimate,
    p_val = res$p.value
  )
}

# for testing
# run_test("Tgfb1", "Quiescent", "Smad3", "NTC_0005")

#  detailed test results
test_results <- 
  comparisons %>%
  pmap(run_test, .progress = TRUE) %>% 
  list_rbind() %>% 
  mutate(p_adj = p.adjust(p_val, method = "BH"))


# summarise across NTCs:
# - mean log odds ratio
# - fraction of significant adjusted p values
# - log odds ratio must have the same sign for NTCs
test_results_summary <-
  test_results %>% 
  mutate(
    log2_or = log2(odds_ratio),
  ) %>% 
  summarise(
    .by = c(condition, cell_type, target),
    mean_lor = mean(log2_or),
    frac_signif = sum(p_adj < 0.01) / n(),
    same_sign = n_distinct(sign(log2_or[p_adj < 0.01])) <= 1,
  ) %>%
  mutate(
    mean_lor = if_else(same_sign, mean_lor, NA),
    frac_signif = if_else(same_sign, frac_signif, NA),
  )



# Plot --------------------------------------------------------------------

## Overview ----

test_results %>% 
  ggplot(aes(condition, ntc)) +
  geom_point(
    aes(color = log2(odds_ratio), size = -log10(p_adj)),
    shape = 16
  ) +
  scale_color_distiller(
    name = "log2(OR)",
    palette = "RdYlBu",
    limits = c(-5, 5),
    breaks = c(-5, 0, 5),
    labels = c("-5 or lower", "0", "5 or higher"),
    oob = scales::oob_squish_any
  ) +
  scale_size_area(
    name = "-log10 padj",
    limits = c(0, 5),
    oob = scales::oob_squish_any,
    max_size = 3
  ) +
  facet_grid(vars(target), vars(cell_type), space = "free", scales = "free") +
  theme_pub(rotate_x_labels = TRUE) +
  theme(
    legend.key.height = unit(2, "mm"),
    legend.key.width = unit(2, "mm"),
    panel.grid = element_blank()
  )



## Summary, selected clusters ----

selected_targets <- c(
  "Tgfbr1", "Smad3", "Smad4", "Wdr82", "Hcfc1", "Srcap", "Znhit1",
  "Yeats4", "Dmap1","Kat5", "Ino80", "Tfpt", "Chd4", "Smarca4",
  "Smarcd1", "Brd9", "Brd7", "Arid2", "Pbrm1", "Kat8", "Kansl1",
  "Paxip1", "Setd1b", "Egr2", "Rest", "Yy1", "Rnf40"
)

clusters_resting <- c(
  "Quiescent",
  "Transitory",
  "Stress-fiber fibroblasts",
  "Tissue-repair",
  "Lamp1-fibroblasts",
  "Proliferative"
) %>% rev()

clusters_tgfb1 <- c(
  "Quiescent",
  "Transitory",
  "Myofibroblasts",
  "Deactivated myofibroblasts",
  "Tissue-repair",
  "Lamp1-fibroblasts"
) %>% rev()

clusters_il1b <- c(
  "Quiescent",
  "Transitory",
  "Inflammatory & fibrotic",
  "Inflammatory (injury response)",
  "Tissue-repair",
  "Lamp1-fibroblasts"
) %>% rev()


plot_summary <- function(condition,
                         targets = NULL,
                         cell_types = NULL,
                         plot_x_axis = TRUE) {
  plot_data <- 
    test_results_summary %>% 
    filter(condition == {{condition}})
  
  # filter and reorder data
  if (!is.null(targets))
    plot_data <- 
      plot_data %>% 
      filter(target %in% {{targets}}) %>% 
      mutate(target = fct_relevel(target, !!!targets))
  
  if (!is.null(cell_types))
    plot_data <- 
      plot_data %>% 
      filter(cell_type %in% {{cell_types}}) %>% 
      mutate(cell_type = fct_relevel(cell_type, !!!cell_types))
  
  p <-
    ggplot(plot_data, aes(target, cell_type)) +
    geom_point(
      aes(color = mean_lor, size = frac_signif),
      shape = 16
    ) +
    facet_grid(vars(condition), NULL) +
    scale_x_discrete(name = NULL, position = "top") +
    scale_y_discrete(name = NULL) +
    scale_color_distiller(
      name = "mean log2(OR)",
      palette = "RdYlBu",
      limits = c(-5, 5),
      breaks = c(-5, 0, 5),
      labels = c("-5 or lower", "0", "5 or higher"),
      oob = scales::oob_squish_any
    ) +
    scale_size_area(
      name = "percentage\nsignificance",
      limits = c(0, 1),
      breaks = c(0, .5, 1),
      max_size = 3
    ) +
    coord_fixed() +
    theme_pub(rotate_x_labels = TRUE) +
    theme(
      axis.text.x = element_text(hjust = 0, vjust = 0.5, face = "italic"),
      legend.key.height = unit(2, "mm"),
      legend.key.width = unit(2, "mm"),
      panel.grid = element_blank(),
      strip.text.y = element_text(face = "bold", angle = 90, vjust = 1)
    )
  
  if (!plot_x_axis)
    p <- 
    p +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  
  p
}

wrap_plots(
  plot_summary(
    "Tgfb1",
    targets = selected_targets,
    cell_types = clusters_tgfb1
  ),
  plot_summary(
    "Il1b",
    targets = selected_targets,
    cell_types = clusters_il1b,
    plot_x_axis = FALSE
  ),
  plot_summary(
    "Resting",
    targets = selected_targets,
    cell_types = clusters_resting,
    plot_x_axis = FALSE
  ),
  ncol = 1,
  guides = "collect"
)
ggsave_default("2g_ko_enrichment_summary", type = "pdf", width = 130)


## Summary, all clusters ----

target_order <- c(
  "Tgfbr1", "Smad2", "Smad3", "Smad4", "Chd4", "Kat5", "Dmap1", "Yeats4",
  "Srcap", "Znhit1", "Ino80", "Tfpt", "Smarca4", "Smarcd1", "Brd9", "Brd7",
  "Arid2", "Pbrm1", "Kat8", "Kansl1", "Hcfc1", "Wdr82", "Paxip1", "Kmt2a",
  "Setd1b", "Setdb1", "Egr2", "Rest", "Yy1", "Rnf40"
)

cluster_order <- c(
  "Quiescent",
  "Transitory",
  "Stress-fiber fibroblasts",
  "Myofibroblasts",
  "Deactivated myofibroblasts",
  "Inflammatory & fibrotic",
  "Inflammatory (injury response)",
  "Tissue-repair",
  "Lamp1-fibroblasts",
  "Proliferative",
  "Proliferative myofibroblast"
) %>% rev()


plot_summary_all <- function(target_order = NULL, cluster_order = NULL) {
  plot_data <-
    test_results_summary %>%
    mutate(condition = fct_relevel(condition, "Resting", "Tgfb1", "Il1b"))

  # filter and reorder data
  if (!is.null(target_order))
    plot_data <-
      plot_data %>%
      mutate(target = fct_relevel(target, !!!target_order))

  if (!is.null(cluster_order))
    plot_data <-
      plot_data %>%
      mutate(cell_type = fct_relevel(cell_type, !!!cluster_order))

  ggplot(plot_data, aes(target, cell_type)) +
    geom_point(
      aes(color = mean_lor, size = frac_signif),
      shape = 16
    ) +
    ylab("cell type") +
    facet_grid(vars(condition), NULL) +
    scale_color_distiller(
      name = "mean log2(OR)",
      palette = "RdYlBu",
      limits = c(-5, 5),
      breaks = c(-5, 0, 5),
      labels = c("-5 or lower", "0", "5 or higher"),
      oob = scales::oob_squish_any
    ) +
    scale_size_area(
      name = "percentage\nsignificance",
      limits = c(0, 1),
      breaks = c(0, .5, 1),
      max_size = 3
    ) +
    coord_fixed() +
    theme_pub(rotate_x_labels = TRUE) +
    theme(
      legend.key.height = unit(2, "mm"),
      legend.key.width = unit(2, "mm"),
      panel.grid = element_blank()
    )
}

plot_summary_all(
  target_order = target_order,
  cluster_order = cluster_order
)
ggsave_default("S6a_ko_enrichment_summary_all", type = "pdf", width = 130)

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



# Plot data ---------------------------------------------------------------

plot_ko_umap <- function(target, condition = NULL, n_bins = 50) {
  if (!is.null(condition))
    col_data <- 
      col_data %>% 
      filter(condition == {{condition}})
  
  ggplot(col_data, aes(x = UMAP_X, y = UMAP_Y)) +
    stat_binhex(
      data = col_data %>% filter(target == "non-targeting"),
      mapping = aes(fill = "x"),
      bins = n_bins,
      fill = "lightgrey"
    ) +
    stat_binhex(
      data = col_data %>% filter(target == {{target}}),
      mapping = aes(fill = log10(after_stat(count))),
      bins = n_bins
    ) +
    scale_fill_distiller(
      name = "log10 cell count",
      palette = "YlOrRd",
      direction = 1,
      limits = c(0, 1.2),
      breaks = c(0, 1.2)
    ) +
    coord_fixed() +
    ggtitle(target) +
    theme_pub() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      legend.key.height = unit(2, "mm"),
      legend.key.width = unit(2, "mm"),
      panel.grid = element_blank(),
      panel.border = element_blank()
    )
}

# for testing
# plot_ko_umap("Smad3")
# plot_ko_umap("non-targeting")

all_targets <- c(
  "non-targeting",
  unique(col_data$target) %>% 
    str_sort() %>% 
    setdiff("non-targeting")
)

targets_main_figure <- 
  c(B = "non-targeting", C = "Tgfbr1", D = "Smad3", E = "Wdr82", F = "Kat5") 

wrap_plots(
  map(targets_main_figure, plot_ko_umap),
  design = "#B\nCD\nEF",
  guides = "collect"
) & theme(legend.position = "bottom")
ggsave_default("2f_umap_ko_distribution", type = "pdf", width = 60)


wrap_plots(
  map(all_targets, plot_ko_umap),
  nrow = 6,
  guides = "collect"
)
ggsave_default("S5_umap_ko_distribution", type = "pdf")

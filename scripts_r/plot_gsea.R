library(tidyverse)
source("scripts_r/utils.R")



# Load data ---------------------------------------------------------------

gsea_results <- read_rds("data_generated/gsea.rds")

selected_comparisons <- c(
  "Resting" = "Resting_Target_vs_NTC",
  "Tgfb1" = "Tgfb1_Target_vs_NTC",
  "Il1b" = "Il1b_Target_vs_NTC"
)



# Plot data ---------------------------------------------------------------

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

selected_targets <- c(
  "Tgfbr1", "Smad2", "Smad3", "Chd4", "Kat5",
  "Dmap1", "Srcap", "Smarca4", "Brd9", "Brd7",
  "Arid2", "Pbrm1", "Kat8", "Hcfc1", "Wdr82",
  "Paxip1", "Kmt2a", "Egr2", "Rest", "Rnf40"
)

plot_terms <- function(terms, targets) {
  terms <- 
    terms %>% 
    mutate(term_display = as_factor(term_display))
  
  # prepare data
  data_vis <- 
    gsea_results %>% 
    filter(
      comparison %in% selected_comparisons,
      group %in% targets
    ) %>% 
    left_join(terms, by = join_by(db, pathway)) %>% 
    filter(!is.na(term_display)) %>% 
    mutate(
      group = fct_relevel(group, !!!targets),
      comparison = 
        comparison %>% 
        fct_recode(!!!selected_comparisons) %>% 
        fct_relevel(!!!names(selected_comparisons))
    )
  
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
    ggtitle("Target vs NTC") +
    facet_wrap(vars(group), nrow = 1) +
    theme_pub(TRUE) +
    theme(
      legend.key.height = unit(3, "mm"),
      legend.key.width = unit(3, "mm"),
      legend.margin = margin(),
      panel.grid = element_blank(),
      panel.spacing = unit(-.5, "pt")
    )
}

plot_terms(selected_terms, selected_targets)
ggsave_default("S6b_gsea", type = "pdf", width = 210)


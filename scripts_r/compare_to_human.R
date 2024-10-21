library(tidyverse)
library(fs)
library(readxl)
library(fgsea)
source("scripts_r/utils.R")



# Load data ---------------------------------------------------------------

## Differentially expressed genes ----

# 3 comparisons * 30 targets * 5288 genes = 475920 rows
dge <- 
  dir_ls("data_generated/", regex = "DEG_") %>% 
  read_csv(id = "file") %>% 
  separate_wider_regex(file, c(".+DEG_", comparison = ".*", "\\.csv")) %>% 
  select(comparison, group, gene = names,
         logFC = logfoldchanges, p_adj = pvals_adj) %>% 
  filter(comparison != "NTC_Stim_vs_Resting")

gene_map <-
  read_tsv("metadata/biomart_human_mouse_20210727.tsv") %>% 
  select(human_gene = `Human gene name`, mouse_gene = `Gene name`) %>% 
  distinct()


## Signatures ----

markers <-
  bind_rows(
    .id = "ref",
    
    Amrute =
      read_xlsx("data_raw/signatures/signatures_Amrute.xlsx", sheet = 1, skip = 1) %>% 
      select(
        cluster,
        human_gene = ...1,
        logFC = avg_log2FC,
        p_adj = p_val_adj
      ),
    
    Chaffin =
      read_xlsx("data_raw/signatures/signatures_Chaffin.xlsx", sheet = 1) %>% 
      select(
        cluster = `Sub-Cluster`,
        human_gene = Gene,
        logFC = `limma-voom: logFC`,
        p_adj = `limma-voom: Adjusted P-Value`
      ),
    
    Fu =
      read_xlsx("data_raw/signatures/signatures_Fu.xlsx", sheet = 1) %>% 
      select(
        cluster,
        human_gene = gene,
        logFC = avg_logFC,
        p_adj = p_val_adj
      ),
    
    Koenig = 
      read_xlsx("data_raw/signatures/markers_Koenig.xlsx", sheet = 3) %>% 
      select(
        cluster,
        human_gene = gene,
        logFC = avg_log2FC,
        p_adj = p_val_adj
      ),
    
    Kuppe =
      c("Fib1", "Fib2", "Fib3", "Fib4") %>% 
      set_names() %>% 
      map(\(s) read_xlsx("data_raw/signatures/signatures_Kuppe.xlsx", sheet = s)) %>% 
      list_rbind(names_to = "cluster") %>% 
      select(
        cluster, 
        human_gene = gene, 
        logFC = avg_log2FC,
        p_adj = p_val_adj
      )
  ) %>% 
  left_join(
    gene_map,
    by = "human_gene",
    relationship = "many-to-many"
  ) %>% 
  select(ref, cluster, gene = mouse_gene, logFC, p_adj) %>% 
  filter(
    p_adj < 0.01,
    !is.na(gene)
  ) %>% 
  bind_rows(read_csv("data_raw/signatures/markers_murine.csv", comment = "#"))



# FGSEA -------------------------------------------------------------------

run_gsea <- function(comparison, group, markers) {
  gene_sets <- 
    markers %>%
    unite(ref, cluster, col = "signature") %>% 
    select(signature, gene) %>% 
    chop(gene) %>% 
    deframe() %>%
    as.list()
  
  ranked_genes <-
    dge %>%
    filter(comparison == {{comparison}}, group == {{group}}) %>%
    arrange(desc(logFC)) %>% 
    select(gene, logFC) %>%
    deframe()
  ranked_genes <- ranked_genes[!is.na(ranked_genes)]
  
  fgsea(
    gene_sets,
    ranked_genes,
    eps = 0
  ) %>%
    as_tibble() %>%
    mutate(
      comparison = {{comparison}},
      group = {{group}},
      .before = 1
    )
}

gsea_results <-
  expand_grid(
    c = unique(dge$comparison),
    g = unique(dge$group)
  ) %>% 
  pmap(\(c, g) run_gsea(c, g, markers), .progress = TRUE) %>% 
  list_rbind()



# Plots -------------------------------------------------------------------

plot_signatures <- function(enriched_terms) {
  selected_comparisons <- c(
    "Resting" = "Resting_Target_vs_NTC",
    "Tgfb1" = "Tgfb1_Target_vs_NTC",
    "Il1b" = "Il1b_Target_vs_NTC"
  )
  
  selected_targets <- c("Tgfbr1", "Smad3", "Smad4", "Kat8", "Kansl1",
                        "Kat5", "Dmap1", "Srcap", "Wdr82", "Hcfc1", "Paxip1")
  
  selected_terms <-
    c("Amrute", "Chaffin", "Fu", "Koenig", "Kuppe") %>% 
    str_flatten("|")
  
  # prepare data
  data_vis <- 
    enriched_terms %>% 
    filter(
      comparison %in% selected_comparisons,
      group %in% selected_targets,
      str_starts(pathway, selected_terms),
      pathway != "Fu_FB7"
    ) %>% 
    mutate(
      group = fct_relevel(group, !!!selected_targets),
      comparison = 
        comparison %>% 
        fct_recode(!!!selected_comparisons) %>% 
        fct_relevel(!!!names(selected_comparisons))
    )
  
  # set limits for colorbar
  color_limit <- max(abs(data_vis$NES), na.rm = TRUE)
  
  # make plot  
  ggplot(data_vis, aes(comparison, pathway, size = -log10(padj))) +
    geom_point(aes(color = NES)) +
    xlab("Condition") +
    ylab("Fibroblast signature") +
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
      limits = c(0, 4),
      breaks = c(0, 2, 4),
      labels = c("0", "2", "4 or higher"),
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

gsea_results %>% 
  mutate(
    pathway =
      as_factor(pathway) %>% 
      fct_rev() %>% 
      fct_relevel(
        "exvivo_resting",
        "exvivo_inflammatory",
        "exvivo_fibrotic"
      )
  ) %>% 
  plot_signatures()
ggsave_default("S10a_fibroblast_signatures_selected_targets", 
               type = "pdf", width = 130)

gsea_results %>% 
  select(comparison, target = group, signature = pathway, NES, padj) %>% 
  {split(., .$comparison)} %>% 
  save_table("fib_signatures")

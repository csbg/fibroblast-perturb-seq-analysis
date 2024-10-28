library(Matrix)
library(SingleCellExperiment)
library(DropletUtils)
library(tidyverse)
library(fs)



# Load data ---------------------------------------------------------------

## Raw data ----

samples <-
  read_csv("metadata/samples_sc_rnaseq.csv") %>% 
  mutate(path = str_glue("data_raw/rna/{name_orig}/outs/raw_feature_bc_matrix"))

sce <-
  samples %>% 
  select(name_python, path) %>% 
  deframe() %>% 
  read10xCounts(intersect.genes = TRUE)


## Data from AnnData object ----

coldata_anndata <-
  read_csv(
    "data_generated/final_adata_obs.csv",
    name_repair = "universal"
  ) %>%
  mutate(UMAP_Y = -UMAP_Y) %>% 
  rename(sample = chip_lane) %>% 
  separate_wider_delim(
    ...1,
    delim = "-",
    names = "barcode",
    too_many = "drop"
  )

logcounts_anndata <-
  read_csv("data_generated/final_adata_X.csv") %>% 
  column_to_rownames("...1") %>%
  as.matrix() %>% 
  t() %>% 
  as("dgCMatrix")



# Subsetting --------------------------------------------------------------

## Columns ----

# AnnData barcodes, sample name, and order
barcodes_anndata <- 
  coldata_anndata %>% 
  select(barcode, sample) %>% 
  mutate(order_anndata = row_number())

# raw data barcodes and sample name,
# ordered by corresponding barcode in AnnData
barcodes_filtered <- 
  colData(sce) %>% 
  as_tibble() %>% 
  separate_wider_delim(
    Barcode,
    delim = "-",
    names = "barcode",
    too_many = "drop"
  ) %>% 
  rename(sample = Sample) %>% 
  mutate(order = row_number()) %>% 
  left_join(barcodes_anndata, by = join_by(sample, barcode)) %>% 
  filter(!is.na(order_anndata)) %>% 
  arrange(order_anndata)

# subset SCE object, ensuring that the order of columns corresponds to AnnData
sce <- sce[, barcodes_filtered$order]


## Rows ----

# use (unique) gene symbols instead of Ensembl IDs
rowData(sce)$Symbol <- make.unique(rowData(sce)$Symbol, sep = "-")
rownames(sce) <- rowData(sce)$Symbol

# create mapping between raw data and AnnData gene names
# (use uppercase names since AnnData names have spurious cases)
genes <-
  rowData(sce) %>% 
  as_tibble() %>%
  mutate(symbol_upper = str_to_upper(Symbol))

genes_py <-
  tibble(symbol_py = rownames(logcounts_anndata)) %>% 
  mutate(symbol_upper = str_to_upper(symbol_py))

genes_data <- left_join(genes_py, genes, by = "symbol_upper")

# subset SCE object, ensuring that the order of rows corresponds to AnnData
sce <- sce[genes_data$Symbol, ]



# Combine data ------------------------------------------------------------

assay(sce, "logcounts", withDimnames = FALSE) <- logcounts_anndata
sce

colData(sce) <- 
  coldata_anndata %>%
  unite(barcode, sample, col = "barcode", sep = "-", remove = FALSE) %>% 
  as("DataFrame")

# add unique colnames (barcode sequence-sample name)
colnames(sce) <- colData(sce)$barcode
sce


# Save SCE object ---------------------------------------------------------

sce %>% write_rds("data_generated/rna.rds")


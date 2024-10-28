# fibroblast-perturb-seq-analysis

This code supplements the upcoming publication by Aguado Alvaro, Garitano, Esser-Skala et al.


## Folders

(Not all of these folders appear in the git repository.)

- `data_generated`: output files generated by the scripts in this repository
- `data_raw`: raw input data
- `doc`: project documentation
- `geo_upload_scripts`: scripts for uploading raw data to Gene Expression Omnibus
- `metadata`: additional required data
- `plots`: generated plots
- `renv`: R environment data
- `scripts_python`: Python scripts
- `scripts_r`: R scripts
- `tables`: tables exported from scripts



## Download data

Create a folder `data_raw` that will contain raw data in the following subfolders:

- `rna`: Download `GSE261783_RAW.tar` from GEO Series GSE261783 and extract all files.
- `rna-seq`: Download XXX from GEO Series GSExxx and extract all files.
- `signatures`:
  - `Amrute.xlsx`: Table S10 from Amrute et al (https://doi.org/10.1038/s44161-023-00260-8)
  - `Buechler.xlsx`: Table S5 from Buechler et al (https://doi.org/10.1038/s41586-021-03549-5)
  - `Chaffin.xlsx`: Table S11 from Chaffin et al (https://doi.org/10.1038/s41586-022-04817-8)
  - `Forte.xlsx`: Table S3 from Forte et al (https://doi.org/10.1016/j.celrep.2020.02.008)
  - `Fu.xlsx`: Additional File 6 from Fu et al (https://doi.org/10.1186/s12916-023-03232-8)
  - `Koenig.xlsx`: Table S27 from Koenig et al (https://doi.org/10.1038/s44161-022-00028-6)
  - `Kuppe.xlsx`: Table S13 from Kuppe et al (https://doi.org/10.1038/s41586-022-05060-x)

Optionally, obtain intermediary data: Extract the contents of `data_generated.tgz` from Zenodo repository https://doi.org/XXX to folder `data_generated`.



## Analysis workflow

### Python

Run Python scripts in order.


### R

Run the following scripts in the folder `scripts_r` in order to run the R analysis pipeline.
`utils.R` contains auxiliary functions and definitions required by several other scripts.

1. `create_sce.R`: interface to Python scripts; creates a SingleCellObject with the same data and metadata as the AnnData object
2. `prepare_signatures.R`: assemble fibroblast signatures and gene sets for enrichment analysis
3. `run_gsea.R`: perform gene set enrichment analysis
4. `plot_signature_scores.R`: fibroblast signature scores as heatmaps (figure 6b, S2c, S2d)
5. `plot_signature_enrichments.R`: fibroblast enrichments as dotplots (S10a)
6. `plot_ko_distribution.R`: UMAPs with knockout distribution (2f, S5)
7. `plot_ko_enrichment.R`: summary of knockout enrichment (2g, S6a)
8. `plot_gsea.R`: plot GSEA results (S6b)
9. `plot_egr2ko_genes.R`: Volcano plot of Egr2-KO vs NTC (S8h)
10. `plot_bulk_rnaseq.R`: bulk RNA-seq analysis of patient-derived fibroblasts (6d, e, S10c, d, e)
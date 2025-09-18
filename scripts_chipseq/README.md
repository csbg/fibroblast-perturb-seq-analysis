# ChIP-seq analysis

## ChIP Pipeline (`ChIP_pipeline.sbs`)

### Overview
The ChIP pipeline is designed for comprehensive ChIP-seq analysis, including preprocessing, mapping, peak calling, and differential binding analysis. This pipeline integrates various bioinformatics tools and custom scripts to handle a sequence of tasks from raw data processing to the identification of significantly enriched regions.

### Prerequisites
- Linux or macOS environment
- SLURM Workload Manager for job scheduling
- Conda package manager
- Required software and libraries:
  - Bowtie2
  - SAMtools
  - BEDTools
  - deepTools
  - cutadapt
  - subread (featureCounts)
  - Miniconda3
  - MACS2
  - HOMER
  - DESeq2 (for differential analysis)
- Access to computational resources with at least 40G of memory and 8 CPUs

### Installation
1. Install Miniconda3 from [Miniconda](https://docs.conda.io/en/latest/miniconda.html) if not already installed.
2. Create a new Conda environment using the provided `chip.yaml` file:

> conda env create -f scripts/chip.yaml --name ChIP

3. Activate the Conda environment:

> conda activate ChIP


### Usage
Submit the job to SLURM using the following command:

> sbatch ChIP_pipeline.sh

The pipeline requires a sample information file (`sample.info.txt`) located in the `data/` directory, detailing sample IDs, reference genomes, and paired-end FASTQ file paths.

### Configuration
- Modify `sample_info` variable to point to your sample information file.
- Adjust memory (`--mem`) and CPUs (`--cpus-per-task`) based on your computational resources and dataset size.

### Pipeline Steps
1. **Preprocessing**: Trimming adapters from raw FASTQ files using cutadapt.
2. **Mapping**: Aligning reads to the reference genome sequences with Bowtie2.
3. **Peak Calling**: Identifying enriched regions using MACS2.
4. **Post-processing**: Removing blacklisted regions and generating BigWig files for visualization.
5. **Differential Binding Analysis**: Using DESeq2 to compare binding across conditions.

### Output
- Trimmed FASTQ files in `data/fastq/`
- BAM files and sorted, deduplicated BAM files in `data/bam/`
- BigWig files for visualization in `data/bigwig_files/`
- Peak files (MACS2 output) in `data/peak_files/`
- Differential analysis results and annotated peak lists in `data/peak_files/Consensus/`

### Note
This pipeline is designed for flexibility and can be adapted to various ChIP-seq analysis needs. Adjust parameters and paths as necessary to fit your dataset and computational environment.


## Differential Binding Analysis in ChIP-Seq Data (`ChIP_pipeline.R`)

This R script outlines the process of performing differential binding analysis on ChIP-seq data using the DESeq2 package, and visualizing the results with a scatterplot highlighting differentially bound peaks based on specific thresholds for log fold change and adjusted p-value.

### Requirements

- R environment
- Packages: DESeq2, ggplot2

### Steps

1. **Load Required Libraries**: Make sure DESeq2 and ggplot2 are installed and loaded.

2. **Project Name Definition**: Specify a project name to label output files.

3. **Data Preparation**:
   - **Binding Intensity Data**: Load a table containing binding intensities per peak across samples.
   - **Conditions Specification**: Define the experimental conditions for each sample.

4. **Differential Expression Analysis**:
   - Initialize a `DESeqDataSet` object.
   - Run the DESeq analysis to perform normalization and differential analysis.

5. **Result Extraction**:
   - Apply thresholds to identify significantly differentially bound peaks. The thresholds are |log2FC| > 0.75 and adjusted p-value < 0.05.

6. **Visualization**:
   - **Scatterplot Creation**: Plot average normalized counts for two conditions, highlighting peaks with significant differential binding. Peaks with a log fold change greater than 0.75 are shown in red (upregulated), those with a log fold change less than -0.75 in blue (downregulated), and non-significant peaks in grey.

### Usage

To run the analysis, execute the R script in an environment where the necessary libraries are installed. Adjust the file path and conditions according to your dataset.

### Output

- A scatterplot visualizing the differential binding of peaks between two conditions, highlighting the significantly differentially bound peaks.

This workflow facilitates the identification and visualization of genomic regions with significant changes in binding affinity, providing insights into the regulatory mechanisms under investigation.

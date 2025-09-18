# Differential Binding Analysis in ChIP-Seq Data

This R script outlines the process of performing differential binding analysis on ChIP-seq data using the DESeq2 package, and visualizing the results with a scatterplot highlighting differentially bound peaks based on specific thresholds for log fold change and adjusted p-value.

## Requirements

- R environment
- Packages: DESeq2, ggplot2

## Steps

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

## Usage

To run the analysis, execute the R script in an environment where the necessary libraries are installed. Adjust the file path and conditions according to your dataset.

## Output

- A scatterplot visualizing the differential binding of peaks between two conditions, highlighting the significantly differentially bound peaks.

This workflow facilitates the identification and visualization of genomic regions with significant changes in binding affinity, providing insights into the regulatory mechanisms under investigation.

# Load Required Libraries
library(DESeq2)
library(ggplot2)

# Define the project name. This variable is used to name output files accordingly.
project <- "ChIP-Seq.pf"

# Load Binding Intensity Data
# Specify the path to the file containing binding intensities per peak.
file_path <- "data/total.peaks_counts.txt"
# Read the table. Assumes first column as rownames (peaks) and remaining columns as samples.
peaks_counts <- read.table(file_path, header = TRUE, sep = "\t", row.names = 1)

# Specify Conditions for the Samples
# Adjust the condition assignments based on your experimental design.
conditions <- factor(c(rep("Normal", 2), rep("TGFb", 2)))

# Differential Expression Analysis with DESeq2
# Initialize a DESeqDataSet. The design formula specifies the conditions under which the experiments were conducted.
dds <- DESeqDataSetFromMatrix(countData = peaks_countdata,
                              colData = DataFrame(condition = conditions),
                              design = ~ condition)

# Run the DESeq analysis pipeline
dds <- DESeq(dds)

# Extracting Results with Log Fold Change Threshold
# Apply log2 fold change threshold (|log2FC| > 0.75) and adjusted p-value (< 0.05) to identify significantly differentially bound peaks.
resultsData <- results(dds, contrast = c("condition", "TGFb", "Normal"))
sigResults <- resultsData[which(abs(resultsData$log2FoldChange) > 0.75 & resultsData$padj < 0.05), ]

# Visualization: Scatterplot of Differentially Bound Peaks

# Filter for significant results based on adjusted p-value < 0.05 and |log2FoldChange| > 0.75
sigResults <- subset(res, padj < 0.05 & abs(log2FoldChange) > 0.75)

# Prepare normalized counts data for plotting
normCounts <- as.data.frame(counts(dds, normalized=TRUE))
normCounts$Normal_avg <- rowMeans(normCounts[,1:2])
normCounts$TGFb_avg <- rowMeans(normCounts[,3:4])

# Assign colors based on differential expression status
# Peaks upregulated (logFC > 0.75) are colored red, downregulated (logFC < -0.75) blue,
# and all others grey to denote non-significant differential expression
normCounts$color <- ifelse(rownames(normCounts) %in% rownames(sigResults[log2FoldChange > 0.75, ]), "red2",
                           ifelse(rownames(normCounts) %in% rownames(sigResults[log2FoldChange < -0.75, ]), "blue3", "grey"))

# Create the scatterplot
ggplot(normCounts, aes(x = Normal_avg, y = TGFb_avg)) +
  geom_point(aes(color = color), alpha=0.6) +
  scale_color_identity() +  # Use actual color values specified in 'color' column
  theme_minimal() +
  labs(title = "Differentially Bound Peaks Scatterplot",
       subtitle = "Red: Upregulated (logFC > 0.75), Blue: Downregulated (logFC < -0.75), Grey: Non-significant",
       x = "Average Normalized Counts in Normal Condition",
       y = "Average Normalized Counts in TGFb Condition",
       color = "Differential Expression Status") +
  theme(legend.position = "none")  # Optionally remove the legend if colors are self-explanatory

# Note: This plot visually distinguishes differentially bound peaks that meet the specific
# criteria of adjusted p-value < 0.05 and |log2FoldChange| > 0.75, according to the defined color scheme.

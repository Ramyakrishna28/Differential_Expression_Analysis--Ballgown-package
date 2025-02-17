# Load necessary libraries
library(ballgown)
library(tidyverse)
library(genefilter)
library(dplyr)
library(ggfortify)

# Set path to the parent directory where the StringTie quantification results are located
parent_dir <- "PATH_TO_YOUR_STRINGTIE_QUANTIFICATION_FOLDER"

# List files in the parent directory
list.files(parent_dir)

# Load all samples from the parent folder
# Modify samplePattern to match the subdirectory structure
bg <- ballgown(dataDir = parent_dir, samplePattern = "SAMPLE_PATTERN", meas = 'all')

# Add phenotype data (e.g., 3 controls and 3 treatments)
pheno_data <- data.frame(
  sample = sampleNames(bg),
  condition = c("condition1", "condition1", "condition2", "condition2") # Modify based on your conditions
)
pData(bg) <- pheno_data

# Filter low-abundance genes
bg_filtered <- subset(bg, "rowMeans(texpr(bg)) > 1", genomesubset = TRUE)

# Perform differential expression analysis at the transcript level
results_transcript <- stattest(bg_filtered, feature = 'transcript', meas = 'FPKM', covariate = 'condition')

# Adjust p-values for multiple testing
results_transcript$qval <- p.adjust(results_transcript$pval, method = "BH")

# Extract transcript expression data
expr_data <- texpr(bg_filtered)  # Transcript expression data

# Get the condition information from the phenotype data
condition <- pData(bg_filtered)$condition

# Calculate mean expression per condition (e.g., "case" and "control")
mean_expr_case <- rowMeans(expr_data[, condition == 'case'])
mean_expr_control <- rowMeans(expr_data[, condition == 'control'])

# Calculate fold change (case/control)
FC <- mean_expr_case / mean_expr_control

# Add fold change and log2 fold change to results
results_transcript$FC <- FC
results_transcript$log2FC <- log2(FC)  # Log2 fold change

# Save the updated results with log2FC to the desired directory
write.csv(results_transcript, file = "PATH_TO_SAVE_RESULTS/transcript_level_differential_expression_with_logFC.csv", row.names = FALSE)

# Define thresholds
log2FC_threshold <- 1.5  # Adjust fold change threshold (e.g., ≈2.8-fold)
pval_threshold <- 0.01   # More strict p-value threshold

# Filter differentially expressed transcripts
filtered_results <- results_transcript[abs(results_transcript$log2FC) > log2FC_threshold & results_transcript$pval < pval_threshold, ]

# Save filtered results
write.csv(filtered_results, file = "PATH_TO_SAVE_RESULTS/up-down_DE_results.csv", row.names = FALSE)

# Merge with additional gene/transcript information
filtered_data <- read.csv("PATH_TO_YOUR_FILTERED_DE_RESULTS/Filtered_DE_results.csv")
t_data <- read.csv("PATH_TO_YOUR_TRANSCRIPT_DATA/t_data_cleaned.ctab", sep = "\t", header = TRUE)

# Ensure correct column names
colnames(filtered_data)
colnames(t_data)

# Merge using the correct columns (adjust column names as needed)
merged_data <- merge(
  filtered_data, 
  t_data[, c("t_id", "t_name", "gene_id", "gene_name")], 
  by.x = "id", 
  by.y = "t_id", 
  all.x = TRUE
)

# Save the merged dataset
write.csv(merged_data, "PATH_TO_SAVE_RESULTS/Trans_significant_genes_mapped.csv", row.names = FALSE)

# Clean up gene and transcript IDs (remove version numbers)
merged_data$t_name <- sub("\\..*", "", merged_data$t_name)
merged_data$gene_id <- sub("\\..*", "", merged_data$gene_id)

# Save the cleaned dataset
write.csv(merged_data, "PATH_TO_SAVE_RESULTS/significant_trans_mapped.csv", row.names = FALSE)

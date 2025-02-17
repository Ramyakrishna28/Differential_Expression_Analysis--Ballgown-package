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

# Perform differential expression analysis
results <- stattest(bg_filtered, feature = 'gene', meas = 'FPKM', covariate = 'condition')

# Adjust p-values for multiple testing
results$qval <- p.adjust(results$pval, method = "BH")

# Save results
write.csv(results, file = "PATH_TO_SAVE_RESULTS/EdgeR_differential_expression_results.csv", row.names = TRUE)

# Perform PCA on the gene expression data
pca_data <- prcomp(texpr(bg_filtered))

# Plot the PCA results
autoplot(pca_data, data = pData(bg_filtered), colour = 'condition', 
         main = "PCA Plot of Samples", frame = TRUE)

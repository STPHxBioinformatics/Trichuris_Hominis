# Load necessary library
library(ggplot2)
library(dplyr)
# Read the eigenvec file
pca_data <- read.table("2runsnpeff_missense_added_pheno_filt_added_sex_pca.eigenvec", header = FALSE, stringsAsFactors = FALSE)

# Remove duplicate first column (PLINK writes the sample ID twice)
pca_data <- pca_data[, -1]

# Rename columns
colnames(pca_data) <- c("SampleID", paste0("PC", 1:(ncol(pca_data) - 1)))

# Display first few rows
head(pca_data)


# Plot PC1 vs PC2
ggplot(pca_data, aes(x = PC1, y = PC2, label = SampleID)) +
  geom_point(alpha = 0.7, color = "blue") +
  labs(title = "PCA Plot (PC1 vs PC2)", x = "PC1", y = "PC2") +
  theme_minimal()




pheno_data <- read.table("Main_Phenotype_File_sex.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(pheno_data) <- c("SampleID", "DuplicateID", "Sex")
pheno_data <- pheno_data[, c("SampleID", "Sex")]
pheno_data$Sex <- factor(pheno_data$Sex, levels = c(1, 2), labels = c("Male", "Female"))
pca_merged <- merge(pca_data, pheno_data, by = "SampleID")
head(pca_merged)



ggplot(pca_merged, aes(x = PC1, y = PC2, color = Sex)) +
  geom_point(alpha = 0.7, size = 3) +
  labs(title = "PCA Plot Colored by Sex", x = "PC1", y = "PC2") +
  theme_minimal() +
  scale_color_manual(values = c("Male" = "blue", "Female" = "red"))  # Custom colors





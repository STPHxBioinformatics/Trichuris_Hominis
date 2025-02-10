
library(qqman)
library(tidyverse)

layout_matrix <- matrix(c(1, 2, 
                          3, 4, 
                          5, 6), 
                        nrow = 3, byrow = TRUE)

layout(layout_matrix, widths = c(2, 1))  # 2/3 for the first column, 1/3 for the second

par(mar = c(5, 6, 2, 3))


gwas_results <- read.table("GWAS_assoc_both_sex.assoc", header = TRUE, stringsAsFactors = FALSE)

gwas_results$CHR <- as.integer(sub(".*scaffold([0-9]+)", "\\1", gwas_results$CHR))
gwas_results <- gwas_results[ , c("SNP", "CHR", "BP", "P")]
gwas_results <- gwas_results %>%
  filter(!is.na(BP) & is.finite(BP) & !is.na(P) & is.finite(P))
gwas_results$Gene <- sub(".*_g(\\d+)$", "g\\1", gwas_results$SNP)
threshold=8.89e-7
highlight <- gwas_results$SNP[gwas_results$P < threshold]
highlight <- sub(".*_g(\\d+)$", "g\\1", highlight)
highlight
gene_list <- highlight
filtered_gwas_results <- gwas_results[gwas_results$Gene %in% highlight, ]
filtered_gwas_results$SNP
highlight <- filtered_gwas_results$SNP

manhattan(gwas_results,
          bp = "BP",               
          p = "P",                
          snp = "SNP",             
          xlab = "",
          ylim = c(0, 7), 
          highlight = highlight,
          suggestiveline = FALSE, 
          genomewideline = -log10(8.89e-7),
)

abline(v = 21250000, col = "grey15", lwd = 2, lty = 2)
abline(v = 53380673, col = "grey15", lwd = 2, lty = 2)
abline(v = 75600000, col = "grey15", lwd = 2, lty = 2)



qq(gwas_results$P)



gwas_results <- read.table("GWAS_assoc_female_sex.assoc", header = TRUE, stringsAsFactors = FALSE)

gwas_results$CHR <- as.integer(sub(".*scaffold([0-9]+)", "\\1", gwas_results$CHR))
gwas_results <- gwas_results[ , c("SNP", "CHR", "BP", "P")]
gwas_results <- gwas_results %>%
  filter(!is.na(BP) & is.finite(BP) & !is.na(P) & is.finite(P))
gwas_results$Gene <- sub(".*_g(\\d+)$", "g\\1", gwas_results$SNP)
threshold=8.89e-7
highlight <- gwas_results$SNP[gwas_results$P < threshold]
highlight <- sub(".*_g(\\d+)$", "g\\1", highlight)
gene_list <- c(gene_list, highlight)
filtered_gwas_results <- gwas_results[gwas_results$Gene %in% highlight, ]
filtered_gwas_results$SNP
highlight <- filtered_gwas_results$SNP

manhattan(gwas_results,
          bp = "BP",               
          p = "P",                
          snp = "SNP",             
          xlab = "",
          ylim = c(0, 7), 
          highlight = highlight,
          suggestiveline = FALSE, 
          genomewideline = -log10(8.89e-7),
)

abline(v = 21250000, col = "grey15", lwd = 2, lty = 2)
abline(v = 53380673, col = "grey15", lwd = 2, lty = 2)
abline(v = 75600000, col = "grey15", lwd = 2, lty = 2)



qq(gwas_results$P)



gwas_results <- read.table("GWAS_assoc_male_sex.assoc", header = TRUE, stringsAsFactors = FALSE)

gwas_results$CHR <- as.integer(sub(".*scaffold([0-9]+)", "\\1", gwas_results$CHR))
gwas_results <- gwas_results[ , c("SNP", "CHR", "BP", "P")]
gwas_results <- gwas_results %>%
  filter(!is.na(BP) & is.finite(BP) & !is.na(P) & is.finite(P))
gwas_results$Gene <- sub(".*_g(\\d+)$", "g\\1", gwas_results$SNP)
threshold=8.89e-7
highlight <- gwas_results$SNP[gwas_results$P < threshold]
highlight <- sub(".*_g(\\d+)$", "g\\1", highlight)
gene_list <- c(gene_list, highlight)
filtered_gwas_results <- gwas_results[gwas_results$Gene %in% highlight, ]
filtered_gwas_results$SNP
highlight <- filtered_gwas_results$SNP

manhattan(gwas_results,
          bp = "BP",               
          p = "P",                
          snp = "SNP",             
          xlab = "",
          ylim = c(0, 7), 
          highlight = highlight,
          suggestiveline = FALSE, 
          genomewideline = -log10(8.89e-7),
)

abline(v = 21250000, col = "grey15", lwd = 2, lty = 2)
abline(v = 53380673, col = "grey15", lwd = 2, lty = 2)
abline(v = 75600000, col = "grey15", lwd = 2, lty = 2)



qq(gwas_results$P)


gene_list
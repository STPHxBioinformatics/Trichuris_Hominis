if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Rsamtools")
BiocManager::install("VariantAnnotation")
library(Rsamtools)
install.packages("zoo")
library(zoo)

bamFile <- "SQK-NBD114-24_barcode09_filtered.bam"
bamCoverage <- pileup(bamFile)
coverage_df <- data.frame(pos = bamCoverage$pos, coverage = bamCoverage$count)


library(VariantAnnotation)
vcf <- readVcf("variants_btub.vcf", "hg19")
# Extract variant positions
variant_pos <- start(rowRanges(vcf))

# Extract the DP4 information from the INFO column
info_data <- info(vcf)
dp4_values <- info_data$DP4

allele_freq <- sapply(dp4_values, function(dp4) {
  if (length(dp4) == 4) {
    ref_count <- dp4[1] + dp4[2]  # Sum of reference forward and reverse reads
    alt_count <- dp4[3] + dp4[4]  # Sum of alternate forward and reverse reads
    alt_count / (ref_count + alt_count)  # Allele frequency
  } else {
    NA  # Handle cases where DP4 is not properly formatted
  }
})
variant_df <- data.frame(pos = variant_pos, allele_frequency = allele_freq)

# Display the data frame
print(variant_df)


library(ggplot2)

# Plot coverage

coverage_df$smoothed_coverage <- rollmean(coverage_df$coverage, k = 100, fill = NA)
vertical_bands <- data.frame(
  xmin = c(847, 854, 835, 771, 949),
  xmax = c(849, 856, 837, 773, 1009),
  label = c("Glu198", "Phe200", "Ser194", "Phe168", "Truncation")
)
vertical_bands$center <- (vertical_bands$xmin + vertical_bands$xmax) / 2


smoothed_plot <- ggplot(coverage_df, aes(x = pos, y = smoothed_coverage)) +
  geom_line(color = "blue") +
  geom_vline(data = variant_df, aes(xintercept = pos), color = "red", linetype = "dashed", alpha = 0.7) +  # Vertical lines for variants
  geom_point(data = variant_df, aes(x = pos, y = allele_frequency * 100), color = "purple", size = 2) +  # Points for allele frequencies (scaled)
  geom_rect(data = vertical_bands, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf), 
            fill = "darkgray", alpha = 0.5, inherit.aes = FALSE) +  # Darker rectangles
  # Add labels to the center of each vertical band
  geom_text(data = vertical_bands, aes(x = center, y = 52.5, label = label), color = "black", size = 3, angle = 90, vjust = 0.5, inherit.aes = FALSE) +
  scale_y_continuous(
    limits = c(0, 105),  # Set the limits for the primary y-axis (coverage)
    sec.axis = sec_axis(~ . / 100, name = "Allele Frequency")  # Add a secondary y-axis
  ) +
  labs(
    x = "Position on beta-tubulin gene",
    y = "Coverage",
    title = "Nanopore long read data on beta-tubulin region of interest containing known mutation positions and predicted truncation"
  ) +
  theme_minimal()


smoothed_plot

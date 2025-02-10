install.packages("esquisse")

library(tidyverse)
library(esquisse)
library(ggsignif)
library(patchwork)
library(gridExtra)
# Read the tab-separated text file
merged_coverage <- read.table("merged_coverage.txt", header = TRUE, sep = "\t")
merged_coverage2 <- read.table("merged_coverage_ref_2.txt", header = TRUE, sep = "\t")
View(merged_coverage2)
merged_coverage2 <- merged_coverage2 %>%
  mutate(species = ifelse(grepl("^GFB", Sample), "T. hominis", "others")) %>%
  mutate(
    coverage_btub = as.numeric(coverage_btub),
    coverage_autosome = as.numeric(coverage_autosome),
    ratio = as.numeric(ratio)
  )

# View the imported data
head(merged_coverage)
head(merged_coverage2)
tail(merged_coverage2)

# Read the CSV file without headers
main_phenotype <- read.csv("Main Phenotype File.csv", header = FALSE, sep = ";")

# Assign column names based on the structure
# Assuming the first column contains the sample IDs and the other columns contain relevant data
colnames(main_phenotype) <- c("Sample", "Sample2", "Col3", "Col4", "Col5", "Col6", "Treatment", "Sex")

# View the imported data
head(main_phenotype)

# Merge the two data frames by the "Sample" column
merged_df <- merge(merged_coverage, main_phenotype, by = "Sample")

# View the merged data
head(merged_df)



view(merged_df)

sample_label <- merged_df %>% 
  filter(Sample == "GFB-7577_HWNGYDSX5_3_TRICHURIS777Tiagba_S147_L003")

# Violin plot for Sex
p1 <- merged_df %>%
  filter(!(Sex %in% "Unknown")) %>%
  ggplot() +
  aes(x = Sex, y = ratio) +
  geom_violin(trim = FALSE, fill = "lightblue") +  # Violin plot
  geom_jitter(aes(color = coverage_autosome), width = 0.2, size = 1) +  # Jitter to show samples
  stat_summary(fun = mean, geom = "point", shape = 23, size = 4, fill = "red") +  # Mean point
  geom_signif(comparisons = list(c("Male", "Female")),
              map_signif_level = TRUE, test = "t.test") +  # Add p-value using t-test
  labs(title = "Coverage Ratio by Sex", x = "Sex", y = "Ratio", color = "Average Overall\nCoverage") +  # Labels
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5)) +
  scale_color_gradient(low = "blue", high = "red")

# Violin plot for Treatment
p2 <- ggplot(merged_df, aes(x = Treatment, y = ratio)) +
  geom_violin(trim = FALSE, fill = "lightblue") +  # Violin plot
  geom_jitter(aes(color = coverage_autosome), width = 0.2, size = 1) +  # Jitter to show samples
  stat_summary(fun = mean, geom = "point", shape = 23, size = 4, fill = "red") +  # Mean point
  geom_signif(comparisons = list(c("Treatment 1", "Treatment 2")),
              map_signif_level = TRUE, test = "t.test") +  # Add p-value using t-test
  labs(title = "Coverage Ratio by Treatment", x = "Treatment", y = "Ratio", color = "Average Overall\nCoverage") +  # Labels
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5)) +
  scale_color_gradient(low = "blue", high = "red")

merged_coverage2$species <- factor(merged_coverage2$species, levels = c("T. hominis", "others"))

# Violin plot for species
p3 <- ggplot(merged_coverage2, aes(x = species, y = ratio)) +
  geom_violin(trim = FALSE, fill = "lightblue") +  # Violin plot
  geom_jitter(aes(color = coverage_autosome), width = 0.2, size = 1) +  # Jitter to show samples
  stat_summary(fun = mean, geom = "point", shape = 23, size = 4, fill = "red") +  # Mean point
  geom_signif(comparisons = list(c("T. hominis", "others")),
              map_signif_level = TRUE, test = "t.test") +  # Add p-value using t-test
  labs(title = "Coverage Ratio by Species", x = "Species", y = "Ratio", color = "Average Overall\nCoverage") +  # Labels
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "none") +
  scale_color_gradient(low = "blue", high = "red") +
  scale_x_discrete(labels = c("T. hominis" = expression(italic("T. incognita")), "others" = "others"))

p3

combined_plot <- p1 + p2 + p3 + plot_layout(guides = 'collect')
combined_plot
write.csv(merged_df, "merged_output.csv", row.names = FALSE)
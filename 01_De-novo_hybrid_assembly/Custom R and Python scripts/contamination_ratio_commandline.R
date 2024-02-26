library(tidyverse)

# Retrieve the command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if an argument is provided
if (length(args) == 0) {
  stop("Please provide the path to the input file as a command line argument.")
}

# Read the input file
input_file <- args[1]
df <- read.delim(input_file, header = FALSE)

# Create the ggplot and save the plot
ggplot(df, aes(x = V2, y = V3, color = cut(V4, breaks = c(0, 0.01, 0.5, Inf), labels = c("<1% contamination", "1-50% contamination", ">50% contamination")))) +
  geom_point() +
  scale_color_manual(values = c("<1% contamination" = "blue", "1-50% contamination" = "green", ">50% contamination" = "orange")) +
  labs(x = "Contig length [bp]", y = "Contamination length [bp]", color = "Contamination Ratio") +
  theme_minimal()
ggsave("Contamination_ratio.jpg")


# Filter the data
filtered_df <- df %>%
  filter(V4 > 0.50)

# Save the filtered_df as a tab-separated file
output_file <- "filtered_data.tsv"
#write.table(filtered_df[1], file = output_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = F)

#write.table(paste("@", filtered_df$V1), file = output_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
lines <- paste0("@", filtered_df$V1)
writeLines(lines, con = output_file)
# Print the first row of the filtered data
#print(filtered_df[1])
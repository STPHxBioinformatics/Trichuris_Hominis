library(tidyverse)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
gff_file1 <- args[1]
gff_file2 <- args[2]
contig_names_file <- args[3]

getwd()

#load df and ref df
df <- read.table(file = gff_file1)
df_ref <- read.table(file = gff_file2)
df_contig_names <- read.table(file = contig_names_file, header = FALSE)
df_contig_names <- rename(df_contig_names, V1.x = V1)
df_contig_names$V1.x <- substr(df_contig_names[,1], 2, 100)
df <- df  %>%
  filter(V3 == "gene")

#create new column with gene ID's
df$V11 <- substr(df[,9], 15, 22)
df_ref$V11 <- substr(df_ref[,9], 15, 22)

#merge df to get LG
df_merged <- dplyr::inner_join(df, df_ref, by="V11")

#create image of contig assignment based on genes
ggplot(df_merged) +
 aes(x = V1.x, fill = V1.y) +
 geom_bar() +
 scale_fill_hue(direction = 1) +
 labs(x = "Contigs", 
 y = "Gene count per contig", title = "Mapped Scaffols to Linkage Groups (LG)", fill = "Linkage groups", subtitle = "LG1 = X-Chromosome, LG2 = Autosome, LG3 = Autosome") +
 theme_minimal() +
 theme(plot.title = element_text(face = "bold", hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave("linkage_group_assignment.jpg")

#data cleaning
df_unique_contig <- as.data.frame(df_merged %>%
  group_by(V1.x) %>%
  count(V1.y) %>%
  mutate(n/sum(n)) %>%
  group_by(V1.x) %>%
  top_n(1,))
df_unique_contig <- df_unique_contig %>%
  arrange(V1.y, desc(n))
df_unique_contig <- tibble::rowid_to_column(df_unique_contig, "ID")
df_unique_contig$New_name <- str_c("TTRE_chr", substr(df_unique_contig[,3], 8,8), "_scaffold", df_unique_contig[,1])
df_new_contig_names <- dplyr::full_join(df_contig_names, df_unique_contig, by="V1.x")
df_new_contig_names <- df_new_contig_names %>%
  arrange(ID)
df_new_contig_names <- df_new_contig_names %>%
  mutate(
    New_name = ifelse(is.na(ID), paste0("TTRE_unplaced_scaffold", seq(nrow(df_unique_contig)+1, nrow(df_unique_contig)+1 + sum(is.na(ID)) - 1)), New_name),
    ID = ifelse(is.na(ID), seq(nrow(df_unique_contig)+1, nrow(df_unique_contig)+1 + sum(is.na(ID)) - 1), ID),
    )
df_new_contig_names <- df_new_contig_names %>%
  arrange(ID)
df_final <-df_new_contig_names %>%
  select(6,1)
names(df_final)<-NULL
names(df_final[,2])<-NULL


write.table(df_final, file = "New_Old_names3", sep = "\t", row.names = F, col.names = F, quote = FALSE)
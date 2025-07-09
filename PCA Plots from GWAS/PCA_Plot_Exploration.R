# Load necessary library


library(ggplot2)
library(dplyr)
library(stringr)
library(ggrepel)

# Read the eigenvec file
pca_data <- read.table("2runsnpeff_missense_added_pheno_filt_added_sex_pca.eigenvec", header = FALSE, stringsAsFactors = FALSE)

# Remove duplicate first column (PLINK writes the sample ID twice)
pca_data <- pca_data[, -1]



# Rename columns
colnames(pca_data) <- c("SampleID", paste0("PC", 1:(ncol(pca_data) - 1)))


pca_data <- pca_data %>%
  mutate(village = str_extract(SampleID, "TRICHURIS\\d+([A-Za-z]+)") %>% str_replace("TRICHURIS\\d+", ""))

# Display first few rows
head(pca_data)


#add sex
pheno_data <- read.table("Main_Phenotype_File_sex.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(pheno_data) <- c("SampleID", "DuplicateID", "Sex")
pheno_data <- pheno_data[, c("SampleID", "Sex")]
pheno_data$Sex <- factor(pheno_data$Sex, levels = c(1, 2), labels = c("Male", "Female"))
pca_merged <- merge(pca_data, pheno_data, by = "SampleID")
head(pca_merged)

#plot sex pca
ggplot(pca_merged, aes(x = PC1, y = PC2, color = Sex, label = SampleID)) +
  geom_point(alpha = 0.7, size = 3) +
  labs(title = "PCA plot colored by Sex", x = "PC1", y = "PC2") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave("plot1.jpeg", width = 6, height = 4, dpi = 600)

#add treatment
treatment_data <- read.table("Main_Phenotype_File.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(treatment_data) <- c("SampleID", "DuplicateID", "Treatment")
treatment_data <- treatment_data[, c("SampleID", "Treatment")]
treatment_data$Treatment <- factor(treatment_data$Treatment, levels = c(1, 2), labels = c("T1", "T2"))
pca_merged <- merge(pca_data, treatment_data, by = "SampleID")
head(pca_merged)

#plot treatment
ggplot(pca_merged, aes(x = PC1, y = PC2, color = Treatment, label = SampleID)) +
  geom_point(alpha = 0.7, size = 3) +
  labs(title = "PCA plot colored by Treatment", x = "PC1", y = "PC2") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave("plot2.jpeg", width = 6, height = 4, dpi = 600 )


#plot village
ggplot(pca_data, aes(x = PC1, y = PC2, color = village, label = SampleID)) +
  geom_point(alpha = 0.7) +
  labs(title = "PCA plot colored by Village", x = "PC1", y = "PC2") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave("plot3.jpeg", width = 6, height = 4, dpi = 600 )









# Read the eigenvec file
pca_data <- read.table("2runsnpeff_missense_added_sex_filt_pca_male.eigenvec", header = FALSE, stringsAsFactors = FALSE)

# Remove duplicate first column (PLINK writes the sample ID twice)
pca_data <- pca_data[, -1]



# Rename columns
colnames(pca_data) <- c("SampleID", paste0("PC", 1:(ncol(pca_data) - 1)))


pca_data <- pca_data %>%
  mutate(village = str_extract(SampleID, "TRICHURIS\\d+([A-Za-z]+)") %>% str_replace("TRICHURIS\\d+", ""))

# Display first few rows
head(pca_data)


#add treatment
treatment_data <- read.table("Main_Phenotype_File.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(treatment_data) <- c("SampleID", "DuplicateID", "Treatment")
treatment_data <- treatment_data[, c("SampleID", "Treatment")]
treatment_data$Treatment <- factor(treatment_data$Treatment, levels = c(1, 2), labels = c("T1", "T2"))
pca_merged <- merge(pca_data, treatment_data, by = "SampleID")
head(pca_merged)

manova_result <- manova(cbind(PC1, PC2) ~ Treatment, data = pca_merged)
summary(manova_result, test = "Pillai")

ggplot(pca_merged, aes(x = PC1, y = PC2, color = Treatment)) +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(type = "t", level = 0.95) +  # 99.9% confidence ellipse per group
  labs(title = "PCA Clustering by Treatment Males", x = "PC1", y = "PC2") +
  theme_minimal() +
  ylim(-0.2, 0.2) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))


ggsave("plot4.jpeg", width = 6, height = 4, dpi = 600 )


anova_pc1 <- aov(PC1 ~ Treatment, data = pca_merged)
summary(anova_pc1)
p_val <- summary(anova_pc1)[[1]][["Pr(>F)"]][1]
p_val

ggplot(pca_merged, aes(x = Treatment, y = PC1, fill = Treatment)) +
  geom_boxplot() +
  labs(title = "PC1 by Treatment", x = "Treatment", y = "PC1") +
  annotate("text", x = 1.5, y = max(pca_merged$PC1, na.rm = TRUE), 
           label = paste0("ANOVA p = ", signif(p_val, 3)), size = 4) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave("plot5.jpeg", width = 6, height = 4, dpi = 600 )


anova_pc2 <- aov(PC2 ~ Treatment, data = pca_merged)
summary(anova_pc2)
p_val <- summary(anova_pc2)[[1]][["Pr(>F)"]][1]
p_val

ggplot(pca_merged, aes(x = Treatment, y = PC2, fill = Treatment)) +
  geom_boxplot() +
  labs(title = "PC2 by Treatment", x = "Treatment", y = "PC2") +
  annotate("text", x = 1.5, y = max(pca_merged$PC2, na.rm = TRUE), 
           label = paste0("ANOVA p = ", signif(p_val, 3)), size = 4) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("plot6.jpeg", width = 6, height = 4, dpi = 600 )

manova_result <- manova(cbind(PC1, PC2) ~ village, data = pca_merged)
summary(manova_result, test = "Pillai")

ggplot(pca_merged, aes(x = PC1, y = PC2, color = village)) +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(type = "t", level = 0.95) +  # 95% confidence ellipse per group
  labs(title = "PCA Clustering by Village Males", x = "PC1", y = "PC2") +
  ylim(-0.2, 0.2) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("plot7.jpeg", width = 6, height = 4, dpi = 600 )

anova_pc1 <- aov(PC1 ~ village, data = pca_merged)
summary(anova_pc1)

tukey_result <- TukeyHSD(anova_pc1)
p_val_tiagba_akakro <- tukey_result$village["Tiagba-Akakro", "p adj"]


ggplot(pca_merged, aes(x = village, y = PC1, fill = village)) +
  geom_boxplot() +
  labs(title = "PC1 by Village", x = "Village", y = "PC1") +
  annotate("text", x = 3.5, y = max(pca_merged$PC1, na.rm = TRUE), 
           label = paste0("p(Tiagba-Akakro)=", signif(p_val_tiagba_akakro, 3)), size = 4) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave("plot8.jpeg", width = 6, height = 4, dpi = 600 )


anova_pc2 <- aov(PC2 ~ village, data = pca_merged)
summary(anova_pc2)

tukey_result <- TukeyHSD(anova_pc2)
p_val_tiagba_akakro <- tukey_result$village["Tiagba-Akakro", "p adj"]


ggplot(pca_merged, aes(x = village, y = PC2, fill = village)) +
  geom_boxplot() +
  labs(title = "PC2 by Village", x = "Village", y = "PC2") +
  annotate("text", x = 3.5, y = max(pca_merged$PC1, na.rm = TRUE), 
           label = paste0("p(Tiagba-Akakro)=", signif(p_val_tiagba_akakro, 3)), size = 4) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("plot9.jpeg", width = 6, height = 4, dpi = 600 )















# Read the eigenvec file
pca_data <- read.table("2runsnpeff_missense_added_sex_filt_pca.eigenvec", header = FALSE, stringsAsFactors = FALSE)

# Remove duplicate first column (PLINK writes the sample ID twice)
pca_data <- pca_data[, -1]



# Rename columns
colnames(pca_data) <- c("SampleID", paste0("PC", 1:(ncol(pca_data) - 1)))


pca_data <- pca_data %>%
  mutate(village = str_extract(SampleID, "TRICHURIS\\d+([A-Za-z]+)") %>% str_replace("TRICHURIS\\d+", ""))

# Display first few rows
head(pca_data)


#add treatment
treatment_data <- read.table("Main_Phenotype_File.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(treatment_data) <- c("SampleID", "DuplicateID", "Treatment")
treatment_data <- treatment_data[, c("SampleID", "Treatment")]
treatment_data$Treatment <- factor(treatment_data$Treatment, levels = c(1, 2), labels = c("T1", "T2"))
pca_merged <- merge(pca_data, treatment_data, by = "SampleID")
head(pca_merged)

manova_result <- manova(cbind(PC1, PC2) ~ Treatment, data = pca_merged)
summary(manova_result, test = "Pillai")

ggplot(pca_merged, aes(x = PC1, y = PC2, color = Treatment)) +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(type = "t", level = 0.95) +  # 99.9% confidence ellipse per group
  labs(title = "PCA Clustering by Treatment Females", x = "PC1", y = "PC2") +
  theme_minimal() +
  xlim(-0.02, 0.02) +
  ylim(-0.02, 0.02) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("plot10.jpeg", width = 6, height = 4, dpi = 600 )



anova_pc1 <- aov(PC1 ~ Treatment, data = pca_merged)
summary(anova_pc1)
p_val <- summary(anova_pc1)[[1]][["Pr(>F)"]][1]
p_val

ggplot(pca_merged, aes(x = Treatment, y = PC1, fill = Treatment)) +
  geom_boxplot() +
  labs(title = "PC1 by Treatment", x = "Treatment", y = "PC1") +
  annotate("text", x = 1.5, y = 0.01, 
           label = paste0("ANOVA p = ", signif(p_val, 3)), size = 4) +
  ylim(-0.01, 0.01) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("plot11.jpeg", width = 6, height = 4, dpi = 600 )


anova_pc2 <- aov(PC2 ~ Treatment, data = pca_merged)
summary(anova_pc2)
p_val <- summary(anova_pc2)[[1]][["Pr(>F)"]][1]
p_val

ggplot(pca_merged, aes(x = Treatment, y = PC2, fill = Treatment)) +
  geom_boxplot() +
  labs(title = "PC2 by Treatment", x = "Treatment", y = "PC2") +
  annotate("text", x = 1.5, y = 0.01, 
           label = paste0("ANOVA p = ", signif(p_val, 3)), size = 4) +
  ylim(-0.01, 0.01) +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) 

ggsave("plot12.jpeg", width = 6, height = 4, dpi = 600 )

manova_result <- manova(cbind(PC1, PC2) ~ village, data = pca_merged)
summary(manova_result, test = "Pillai")

ggplot(pca_merged, aes(x = PC1, y = PC2, color = village)) +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(type = "t", level = 0.95) +  # 95% confidence ellipse per group
  labs(title = "PCA Clustering by Village Females", x = "PC1", y = "PC2") +
  xlim(-0.02, 0.02) +
  ylim(-0.02, 0.02) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) 

ggsave("plot13.jpeg", width = 6, height = 4, dpi = 600 )

anova_pc1 <- aov(PC1 ~ village, data = pca_merged)
summary(anova_pc1)

tukey_result <- TukeyHSD(anova_pc1)
tukey_result
p_val_tiagba_akakro <- tukey_result$village["Tiagba-Akakro", "p adj"]


ggplot(pca_merged, aes(x = village, y = PC1, fill = village)) +
  geom_boxplot() +
  labs(title = "PC1 by Village", x = "Village", y = "PC1") +
  annotate("text", x = 3.5, y = 0.01, 
           label = paste0("p(Tiagba-Akakro)=", signif(p_val_tiagba_akakro, 3)), size = 4) +
  theme_minimal() +
  ylim(-0.01, 0.01) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave("plot14.jpeg", width = 6, height = 4, dpi = 600 )


anova_pc2 <- aov(PC2 ~ village, data = pca_merged)
summary(anova_pc2)

tukey_result <- TukeyHSD(anova_pc2)
tukey_result
p_val_tiagba_akakro <- tukey_result$village["Tiagba-Akakro", "p adj"]


ggplot(pca_merged, aes(x = village, y = PC2, fill = village)) +
  geom_boxplot() +
  labs(title = "PC2 by Village", x = "Village", y = "PC2") +
  annotate("text", x = 3.5, y = 0.01, 
           label = paste0("p(Tiagba-Akakro)=", signif(p_val_tiagba_akakro, 3)), size = 4) +
  theme_minimal() +
  ylim(-0.01, 0.01) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))


ggsave("plot15.jpeg", width = 6, height = 4, dpi = 600 )































# 
# 
# 
# #plot village
# ggplot(pca_data, aes(x = PC1, y = PC2, color = village, label = SampleID)) +
#   geom_point(alpha = 0.7) +
#   labs(title = "PCA plot colored by Village", x = "PC1", y = "PC2") +
#   theme_minimal() +
#   theme(
#     plot.title = element_text(hjust = 0.5, face = "bold")
#   )
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # Read the eigenvec file
# pca_data <- read.table("2runsnpeff_missense_added_sex_filt_pca.eigenvec", header = FALSE, stringsAsFactors = FALSE)
# 
# # Remove duplicate first column (PLINK writes the sample ID twice)
# pca_data <- pca_data[, -1]
# 
# 
# 
# # Rename columns
# colnames(pca_data) <- c("SampleID", paste0("PC", 1:(ncol(pca_data) - 1)))
# 
# 
# pca_data <- pca_data %>%
#   mutate(village = str_extract(SampleID, "TRICHURIS\\d+([A-Za-z]+)") %>% str_replace("TRICHURIS\\d+", ""))
# 
# # Display first few rows
# head(pca_data)
# 
# 
# #add sex
# pheno_data <- read.table("Main_Phenotype_File_sex.txt", header = FALSE, stringsAsFactors = FALSE)
# colnames(pheno_data) <- c("SampleID", "DuplicateID", "Sex")
# pheno_data <- pheno_data[, c("SampleID", "Sex")]
# pheno_data$Sex <- factor(pheno_data$Sex, levels = c(1, 2), labels = c("Male", "Female"))
# pca_merged <- merge(pca_data, pheno_data, by = "SampleID")
# head(pca_merged)
# 
# #plot sex pca
# ggplot(pca_merged, aes(x = PC1, y = PC2, color = Sex, label = SampleID)) +
#   geom_point(alpha = 0.7, size = 3) +
#   labs(title = "PCA plot colored by Sex", x = "PC1", y = "PC2") +
#   theme_minimal() +
#   scale_color_manual(values = c("Male" = "blue", "Female" = "red")) +
#   theme(
#     plot.title = element_text(hjust = 0.5, face = "bold")
#   )
# 
# 
# 
# 
# 
# 
# #add treatment
# treatment_data <- read.table("Main_Phenotype_File.txt", header = FALSE, stringsAsFactors = FALSE)
# colnames(treatment_data) <- c("SampleID", "DuplicateID", "Treatment")
# treatment_data <- treatment_data[, c("SampleID", "Treatment")]
# treatment_data$Treatment <- factor(treatment_data$Treatment, levels = c(1, 2), labels = c("T1", "T2"))
# pca_merged <- merge(pca_data, treatment_data, by = "SampleID")
# head(pca_merged)
# 
# #plot treatment
# ggplot(pca_merged, aes(x = PC1, y = PC2, color = Treatment, label = SampleID)) +
#   geom_point(alpha = 0.7, size = 3) +
#   labs(title = "PCA plot colored by treatment", x = "PC1", y = "PC2") +
#   theme_minimal() +
#   scale_color_manual(values = c("T1" = "blue", "T2" = "red")) +
#   theme(
#     plot.title = element_text(hjust = 0.5, face = "bold")
#   )
# 
# 
# manova_result <- manova(cbind(PC1, PC2) ~ Treatment, data = pca_merged)
# summary(manova_result, test = "Pillai")
# 
# #plot village
# ggplot(pca_data, aes(x = PC1, y = PC2, color = village, label = SampleID)) +
#   geom_point(alpha = 0.7) +
#   labs(title = "PCA plot colored by Village", x = "PC1", y = "PC2") +
#   theme_minimal() +
#   theme(
#     plot.title = element_text(hjust = 0.5, face = "bold")
#   )
# 
# 
# 
# 
# 
# 
# 
# #Significance in 2D space?
# manova_result <- manova(cbind(PC1, PC2) ~ village, data = pca_merged)
# summary(manova_result, test = "Pillai")
# 
# 
# ggplot(pca_merged, aes(x = PC1, y = PC2, color = Treatment)) +
#   geom_point(size = 3, alpha = 0.7) +
#   stat_ellipse(type = "t", level = 0.999) +  # 95% confidence ellipse per group
#   labs(title = "PCA Clustering by Village", x = "PC1", y = "PC2") +
#   theme_minimal() +
#   xlim(-0.0575, 0.05) +
#   ylim(-0.0575, 0.1) +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# #Significance in 1D (each pc individually) space
# 
# anova_pc1 <- aov(PC1 ~ Treatment, data = pca_merged)
# summary(anova_pc1)
# 
# anova_pc2 <- aov(PC2 ~ village, data = pca_data)
# summary(anova_pc2)
# 
# TukeyHSD(anova_pc1)
# TukeyHSD(anova_pc2)
# 
# 
# ggplot(pca_merged, aes(x = Treatment, y = PC1, fill = Treatment)) +
#   geom_boxplot() +
#   labs(title = "PC1 by Treatment", x = "Village", y = "PC1") +
#   theme_minimal()
# 
# ggplot(pca_data, aes(x = village, y = PC2, fill = village)) +
#   geom_boxplot() +
#   labs(title = "PC2 by Village", x = "Village", y = "PC2") +
#   theme_minimal()
# 
# 
# pheno_data <- read.table("Main_Phenotype_File_sex.txt", header = FALSE, stringsAsFactors = FALSE)
# colnames(pheno_data) <- c("SampleID", "DuplicateID", "Sex")
# pheno_data <- pheno_data[, c("SampleID", "Sex")]
# pheno_data$Sex <- factor(pheno_data$Sex, levels = c(1, 2), labels = c("Male", "Female"))
# pca_merged <- merge(pca_data, pheno_data, by = "SampleID")
# head(pca_merged)
# 
# treatment_data <- read.table("Main_Phenotype_File.txt", header = FALSE, stringsAsFactors = FALSE)
# colnames(treatment_data) <- c("SampleID", "DuplicateID", "Treatment")
# treatment_data <- treatment_data[, c("SampleID", "Treatment")]
# treatment_data$Treatment <- factor(treatment_data$Treatment, levels = c(1, 2), labels = c("T1", "T2"))
# pca_merged <- merge(pca_data, treatment_data, by = "SampleID")
# head(pca_merged)
# 
# 
# ggplot(pca_merged, aes(x = PC1, y = PC2, color = Treatment, label = SampleID)) +
#   geom_point(alpha = 0.7, size = 3) +
# #  geom_text_repel(size = 3, max.overlaps = 30) +
# 
#   labs(title = "PCA plot colored by treatment", x = "PC1", y = "PC2") +
#   theme_minimal() +
#   scale_color_manual(values = c("T1" = "blue", "T2" = "red")) +
#   theme(
#     plot.title = element_text(hjust = 0.5, face = "bold")
#   )
# 
# ggplot(pca_merged, aes(x = PC1, y = PC2, color = Sex, label = SampleID)) +
#   geom_point(alpha = 0.7, size = 3) +
#   labs(title = "PCA plot colored by Sex", x = "PC1", y = "PC2") +
#   theme_minimal() +
#   scale_color_manual(values = c("Male" = "blue", "Female" = "red")) +
#   theme(
#     plot.title = element_text(hjust = 0.5, face = "bold")
#   )
# 
# 
# # PCA Plot for Male Worms
# ggplot(pca_merged %>% filter(Sex == "Male"), aes(x = PC1, y = PC2)) +
#   geom_point(color = "blue", alpha = 0.7, size = 3) +
#   labs(title = "PCA Plot for Male Worms", x = "PC1", y = "PC2") +
#   xlim(-0.0575, -0.05)+
#   theme_minimal()
# 
# # PCA Plot for Female Worms
# ggplot(pca_merged %>% filter(Sex == "Female"), aes(x = PC1, y = PC2)) +
#   geom_point(color = "red", alpha = 0.7, size = 3) +
#   labs(title = "PCA Plot for Female Worms", x = "PC1", y = "PC2") +
#   xlim(0.025, 0.035)+
#   theme_minimal()

library(data.table)
library(VennDiagram)
library(tidyverse)


x <- fread("Orthogroups.GeneCount.tsv")

View(x)
head(x)


#transcripts:

trichinella_spiralis_PRJNA12603 <- x[trichinella_spiralis_PRJNA12603 > 0]
trichuris_cote_divoire <- x[trichuris_cote_divoire > 0]
trichuris_muris_PRJEB126 <- x[trichuris_muris_PRJEB126 > 0]
trichuris_suis_PRJNA208415 <- x[trichuris_suis_PRJNA208415 > 0]
trichuris_trichiura_doyle <- x[trichuris_trichiura_doyle > 0]





trichuris_cote_divoire_orthogroups <- trichuris_cote_divoire$Orthogroup
trichuris_muris_PRJEB126_orthogroups <- trichuris_muris_PRJEB126$Orthogroup
trichuris_suis_PRJNA208415_orthogroups <- trichuris_suis_PRJNA208415$Orthogroup
trichuris_trichiura_doyle_orthogroups <- trichuris_trichiura_doyle$Orthogroup

View(trichuris_suis_PRJNA208415_orthogroups)
View(trichuris_muris_PRJEB126_orthogroups)

colors <- c("#440154ff", '#21908dff', '#fde725ff', '#6b7fff')

list(trichuris_cote_divoire_orthogroups, trichuris_muris_PRJEB126_orthogroups, trichuris_suis_PRJNA208415_orthogroups, 
     trichuris_suis_PRJNA208416_orthogroups, trichuris_trichiura_PRJEB535_orthogroups, trichuris_trichiura_doyle_orthogroups)

paste0(length(trichuris_cote_divoire_orthogroups))
paste0(length(trichuris_muris_PRJEB126_orthogroups))
paste0(length(trichuris_suis_PRJNA208415_orthogroups))
paste0(length(trichuris_trichiura_doyle_orthogroups))

venn.diagram(x = list(trichuris_cote_divoire_orthogroups, trichuris_muris_PRJEB126_orthogroups, trichuris_suis_PRJNA208415_orthogroups, trichuris_trichiura_doyle_orthogroups),
             category.names = c(expression(atop(italic("T. Côte d'Ivoire"), plain("8,001"))), expression(atop(italic("T. muris"), plain("6,145"))),
                                expression(atop(italic("T. suis"), plain("9,479"))), expression(atop(italic("T. trichiura"),plain("9,651")))),
             filename = "venn_diagram_with_trichuris_relatives.png",
             output=TRUE,
             imagetype="png" ,
             compression = "lzw",
             lwd = 1,
             col = c("#440154ff", '#21908dff', '#F6BE00', '#6b7fff'),
             fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#F6BE00',0.3),  alpha('#6b7fff',0.3)),
             cex = 0.6,
             fontfamily = "sans",
             cat.cex = 0.6,
             cat.default.pos = "outer",
             cat.fontfamily = "sans",
             cat.col = c("#440154ff", '#21908dff', '#F6BE00', '#6b7fff'),
             margin = 0.20,
             euler.d = FALSE)



datatable <- as.data.frame(x)
head(datatable[, c(1,2, 4, 6,8)])

# Assuming you have already selected columns 3, 5, and 6 and stored them in selected_columns
datatable[, c(2, 4, 6,8)] <- lapply(datatable[, c(2, 4, 6,8)], as.numeric)


#data filtering for only trichiura and cote divoire:
filtered_datatable <- datatable %>%
  filter(Trichuris_cote_divoire_freeze_genes > 0 &
           trichuris_trichiura_doyle_genes > 0 &
           trichuris_suis_PRJNA208415_genes == 0)

orthogroup_list <- filtered_datatable$Orthogroup
orthogroup_list <- paste0(orthogroup_list, ".fa")
orthogroup_list
writeLines(as.character(orthogroup_list), "orthogroup_list_new.txt")


#data filtering for only trichiura and cote divoire with muris:
filtered_datatable2 <- datatable %>%
  filter(trichuris_cote_divoire > 0 &
           trichuris_trichiura_doyle > 0 &
           trichinella_spiralis_PRJNA12603 == 0 &
           trichuris_muris_PRJEB126 > 0 &
           trichuris_suis_PRJNA208416 == 0)


orthogroup_list2 <- filtered_datatable2$Orthogroup
orthogroup_list2 <- paste0(orthogroup_list2, ".fa")
orthogroup_list2
writeLines(as.character(orthogroup_list2), "orthogroup_list2_muris.txt")

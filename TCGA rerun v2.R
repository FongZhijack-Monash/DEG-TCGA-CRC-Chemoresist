install.packages("tidyverse")
install.packages("BiocManager")
BiocManager::install("TCGAbiolinks")
BiocManager::install("edgeR")
BiocManager::install("EDASeq")
BiocManager::install("TCGAutils")
BiocManager::install("maftools")
install.packages("tibble")
install.packages("GeoTcgaData")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")

library(tidyverse)
library(TCGAbiolinks)
library(edgeR)
library(SummarizedExperiment)
library(EDASeq)
library(TCGAutils)
library(maftools)
library(tibble)
library(GeoTcgaData)
library(clusterProfiler)
library(org.Hs.eg.db)


#DEG analysis for TCGA colon adenocarcinoma

# Data download
COAD_chemo <- GDCquery(project = "TCGA-COAD", 
                       data.category = "Clinical",
                       data.type = "Clinical Supplement",
                       data.format = "BCR Biotab",
                       file.type = "drug")
GDCdownload(COAD_chemo)

# Downloading HTseq counts file from TCGA

COAD_HT_query <- GDCquery(project = "TCGA-COAD",
                          data.category = "Transcriptome Profiling",
                          data.type = "Gene Expression Quantification",
                          experimental.strategy = "RNA-Seq",
                          workflow.type = "HTSeq - Counts")

GDCdownload(COAD_HT_query)

# Preparation of HTseq counts file for downstream DGE analysis

COADprepare <- GDCprepare(COAD_HT_query)

COADmatrix <- assays(COADprepare, "raw_count", withDimnames = TRUE)

COADpreprocess <- TCGAanalyze_Preprocessing(object = COADprepare)

COADnorm <- TCGAanalyze_Normalization(tabDF = COADpreprocess, geneInfo = geneInfoHT)

# Filtering of patients' TCGA barcode (manual filtering of patient's barcode by clinical drug response in Excel)

COAD_sens <- read.csv(file = "COAD_sensitive_barcode.csv")
COAD_Res <- read.csv(file = "COAD_resistant_barcode.csv")

COAD_sens_barcode <- COAD_sens[,1]

COAD_res_barcode <- COAD_Res[,1]

COADnorm_tib <- as_tibble(COADnorm)

COAD_res_select_barcode <- dplyr::select(COADnorm_tib,
                                         starts_with(COAD_res_barcode))

COAD_Resist <- colnames(COAD_res_select_barcode)

COAD_sens_select_barcode <- dplyr::select(COADnorm_tib,
                                          starts_with(COAD_sens_barcode))

COAD_sensitive <- colnames(COAD_sens_select_barcode)

# DGE analysis

COAD_Degs <- TCGAanalyze_DEA(mat1 = COADnorm[, COAD_sensitive],
                             mat2 = COADnorm[, COAD_Resist],
                             Cond1type = "Sensitive",
                             Cond2type = "Resistant",
                             fdr.cut = 0.05,
                             logFC.cut = 1)

COAD_DEGs_tab <- TCGAanalyze_LevelTab(COAD_Degs, "Sensitive", "Resistant",
                                      COADnorm[,COAD_sensitive],
                                      COADnorm[,COAD_Resist])

write.csv(COAD_DEGs_tab, "G:/My Drive/Year 2 (2021)/Bioinformatics/Rstudio/Rerun TCGA v2/TCGA Rerun V2/COAD_DEGs.csv")




# DEG code for Rectum adenocarcinoma

# Clinical Data download

READ_chemo <- GDCquery(project = "TCGA-READ", 
                       data.category = "Clinical",
                       data.type = "Clinical Supplement",
                       data.format = "BCR Biotab",
                       file.type = "drug")
GDCdownload(READ_chemo)

# Downloading HTseq counts file from TCGA

READ_HT_query <- GDCquery(project = "TCGA-READ",
                          data.category = "Transcriptome Profiling",
                          data.type = "Gene Expression Quantification",
                          experimental.strategy = "RNA-Seq",
                          workflow.type = "HTSeq - Counts")

GDCdownload(READ_HT_query)

# Preparation of HTseq counts file for downstream DGE analysis

READprepare <- GDCprepare(READ_HT_query)

READmatrix <- assays(READprepare, "raw_count", withDimnames = TRUE)

READpreprocess <- TCGAanalyze_Preprocessing(object = READprepare)

READnorm <- TCGAanalyze_Normalization(tabDF = READpreprocess, geneInfo = geneInfoHT)

# Filtering of patients' TCGA barcode (manual filtering of patient's barcode by clinical drug response in Excel)

READ_sensitive <- read.csv(file = "READ_sensitive_barcode.csv")
READ_resist <- read.csv(file = "READ_resistant_barcode.csv")

READ_sens_bar <- READ_sensitive[,1]

READ_resist_bar <- READ_resist[,1]

READnormtib <- as_tibble(READnorm)

Select_READ_Resist_bar <- dplyr::select(READnormtib,
                                        starts_with(READ_resist_bar))

Resistant_READ <- colnames(Select_READ_Resist_bar)

Select_READ_sensitive_bar <- dplyr::select(READnormtib,
                                           starts_with(READ_sens_bar))

Sensitive_READ <- colnames(Select_READ_sensitive_bar)

# DGE analysis

READ_raw_DEGs <- TCGAanalyze_DEA(mat1 = READnorm[, Sensitive_READ],
                                 mat2 = READnorm[, Resistant_READ],
                                 Cond1type = "Sensitive",
                                 Cond2type = "Resistant",
                                 fdr.cut = 0.05,
                                 logFC.cut = 1)

READ_DEGs_tab <- TCGAanalyze_LevelTab(READ_raw_DEGs, "Sensitive", "Resistant",
                                      READnorm[,Sensitive_READ],
                                      READnorm[,Resistant_READ])

write.csv(READ_DEGs_tab, "G:/My Drive/Year 2 (2021)/Bioinformatics/Rstudio/Rerun TCGA v2/TCGA Rerun V2/READ_DEGs.csv")

########## Try filtered barcodes (COAD)

COAD_sens_filt <- read.csv(file = "COAD_sensitive_filt.csv")
COAD_Res_filt <- read.csv(file = "COAD_resistant_filt.csv")

COAD_sens_barcode_filt <- COAD_sens_filt[,1]

COAD_res_barcode_filt <- COAD_Res_filt[,1]

COAD_res_select_barcode_filt <- dplyr::select(COADnorm_tib,
                                         starts_with(COAD_res_barcode_filt))

COAD_Resist_filt <- colnames(COAD_res_select_barcode_filt)

COAD_sens_select_barcode_filt <- dplyr::select(COADnorm_tib,
                                          starts_with(COAD_sens_barcode_filt))

COAD_sensitive_filt <- colnames(COAD_sens_select_barcode_filt)

# DGE analysis

COAD_Degs_filt <- TCGAanalyze_DEA(mat1 = COADnorm[, COAD_sensitive_filt],
                             mat2 = COADnorm[, COAD_Resist_filt],
                             Cond1type = "Sensitive",
                             Cond2type = "Resistant",
                             fdr.cut = 0.05,
                             logFC.cut = 1)

COAD_DEGs_tab_filt <- TCGAanalyze_LevelTab(COAD_Degs_filt, "Sensitive", "Resistant",
                                      COADnorm[,COAD_sensitive_filt],
                                      COADnorm[,COAD_Resist_filt])

write.csv(COAD_DEGs_tab_filt, "G:/My Drive/Year 2 (2021)/Bioinformatics/Rstudio/Rerun TCGA v2/TCGA Rerun V2/COAD_DEGs_filt.csv")


########## Try filtered barcodes (READ)

READ_sensitive_filt <- read.csv(file = "READ_sensitive_filtered.csv")
READ_resist_filt <- read.csv(file = "READ_resistant_filtered.csv")

READ_sens_bar_filt <- READ_sensitive_filt[,1]

READ_resist_bar_filt <- READ_resist_filt[,1]

Select_READ_Resist_bar_filt <- dplyr::select(READnormtib,
                                        starts_with(READ_resist_bar_filt))

Resistant_READ_filt <- colnames(Select_READ_Resist_bar_filt)

Select_READ_sensitive_bar_filt <- dplyr::select(READnormtib,
                                           starts_with(READ_sens_bar_filt))

Sensitive_READ_filt <- colnames(Select_READ_sensitive_bar_filt)

# DGE analysis

READ_raw_DEGs_filt <- TCGAanalyze_DEA(mat1 = READnorm[, Sensitive_READ_filt],
                                 mat2 = READnorm[, Resistant_READ_filt],
                                 Cond1type = "Sensitive",
                                 Cond2type = "Resistant",
                                 fdr.cut = 0.05,
                                 logFC.cut = 1)

READ_DEGs_tab_filt <- TCGAanalyze_LevelTab(READ_raw_DEGs_filt, "Sensitive", "Resistant",
                                      READnorm[,Sensitive_READ],
                                      READnorm[,Resistant_READ])

write.csv(READ_DEGs_tab_filt, "G:/My Drive/Year 2 (2021)/Bioinformatics/Rstudio/Rerun TCGA v2/TCGA Rerun V2/READ_DEGs_filt.csv")



# Export normalized counts file

write.csv(COADnorm, "COADnorm.csv")


#### Try subclustering (COAD no AA)


COAD_sensnoaa <- read.csv(file = "COAD_sen_noaa.csv")
COAD_Resnoaa <- read.csv(file = "COAD_res_noaa.csv")

COAD_sensnoaa_barcode <- COAD_sensnoaa[,1]

COAD_resnoaa_barcode <- COAD_Resnoaa[,1]

COAD_resnoaa_select_barcode <- dplyr::select(COADnorm_tib,
                                         starts_with(COAD_resnoaa_barcode))

COAD_Resistnoaa <- colnames(COAD_resnoaa_select_barcode)

COAD_sensnoaa_select_barcode <- dplyr::select(COADnorm_tib,
                                          starts_with(COAD_sensnoaa_barcode))

COAD_sensitivenoaa <- colnames(COAD_sensnoaa_select_barcode)

# DGE analysis

COAD_Degsnoaa <- TCGAanalyze_DEA(mat1 = COADnorm[, COAD_sensitivenoaa],
                             mat2 = COADnorm[, COAD_Resistnoaa],
                             Cond1type = "Sensitive",
                             Cond2type = "Resistant",
                             fdr.cut = 0.05,
                             logFC.cut = 1)

COAD_DEGs_tabnoaa <- TCGAanalyze_LevelTab(COAD_Degsnoaa, "Sensitive", "Resistant",
                                      COADnorm[,COAD_sensitivenoaa],
                                      COADnorm[,COAD_Resistnoaa])

write.csv(COAD_DEGs_tabnoaa, "G:/My Drive/Year 2 (2021)/Bioinformatics/Rstudio/Rerun TCGA v2/TCGA Rerun V2/COAD_DEGsnoaa.csv")


#### Try subclustering (AA only)


COAD_sensaa <- read.csv(file = "COAD_sen_aa.csv")
COAD_Resaa <- read.csv(file = "COAD_res_aa.csv")

COAD_sensaa_barcode <- COAD_sensaa[,1]

COAD_resaa_barcode <- COAD_Resaa[,1]

COAD_resaa_select_barcode <- dplyr::select(COADnorm_tib,
                                             starts_with(COAD_resaa_barcode))

COAD_Resistaa <- colnames(COAD_resaa_select_barcode)

COAD_sensaa_select_barcode <- dplyr::select(COADnorm_tib,
                                              starts_with(COAD_sensaa_barcode))

COAD_sensitiveaa <- colnames(COAD_sensaa_select_barcode)

# DGE analysis

COAD_Degsaa <- TCGAanalyze_DEA(mat1 = COADnorm[, COAD_sensitiveaa],
                                 mat2 = COADnorm[, COAD_Resistaa],
                                 Cond1type = "Sensitive",
                                 Cond2type = "Resistant",
                                 fdr.cut = 0.05,
                                 logFC.cut = 1)

COAD_DEGs_tabaa <- TCGAanalyze_LevelTab(COAD_Degsaa, "Sensitive", "Resistant",
                                          COADnorm[,COAD_sensitiveaa],
                                          COADnorm[,COAD_Resistaa])

write.csv(COAD_DEGs_tabaa, "G:/My Drive/Year 2 (2021)/Bioinformatics/Rstudio/Rerun TCGA v2/TCGA Rerun V2/COAD_DEGsaa.csv")

####### Try heatmap

install.packages("ggplot2")
BiocManager::install("ComplexHeatmap")

library(org.Hs.eg.db)
library(ComplexHeatmap)
library(tibble)

# Convert ENSEMBL to GeneID

coad.df <- as.data.frame(COAD_Degs)
coad.df$symbol <- mapIds(org.Hs.eg.db, keys = rownames(coad.df), keytype = "ENSEMBL", column = "SYMBOL")
write.csv(coad.df, "G:/My Drive/Year 2 (2021)/Bioinformatics/Rstudio/Rerun TCGA v2/TCGA Rerun V2/COAD_DEG_GeneID.csv")

coad_geneid <- as.data.frame(COAD_DEGs_tab)
coad_geneid$symbol <- mapIds(org.Hs.eg.db, keys = rownames(coad_geneid), keytype = "ENSEMBL", column = "SYMBOL")
write.csv(coad_geneid, "G:/My Drive/Year 2 (2021)/Bioinformatics/Rstudio/Rerun TCGA v2/TCGA Rerun V2/COAD_DEG_GeneID.csv")

read_geneid <- as.data.frame(READ_DEGs_tab)
read_geneid$symbol <- mapIds(org.Hs.eg.db, keys = rownames(read_geneid), keytype = "ENSEMBL", column = "SYMBOL")
write.csv(read_geneid, "G:/My Drive/Year 2 (2021)/Bioinformatics/Rstudio/Rerun TCGA v2/TCGA Rerun V2/READ_DEG_GeneID.csv")

coad_aa_geneid <- as.data.frame(COAD_DEGs_tabaa)
coad_aa_geneid$symbol <- mapIds(org.Hs.eg.db, keys = rownames(coad_aa_geneid), keytype = "ENSEMBL", column = "SYMBOL")
write.csv(coad_aa_geneid, "G:/My Drive/Year 2 (2021)/Bioinformatics/Rstudio/Rerun TCGA v2/TCGA Rerun V2/COAD_DEGsaa_GeneID.csv")

coad_noaa_geneid <- as.data.frame(COAD_DEGs_tabnoaa)
coad_noaa_geneid$symbol <- mapIds(org.Hs.eg.db, keys = rownames(coad_noaa_geneid), keytype = "ENSEMBL", column = "SYMBOL")
write.csv(coad_noaa_geneid, "G:/My Drive/Year 2 (2021)/Bioinformatics/Rstudio/Rerun TCGA v2/TCGA Rerun V2/COAD_DEGsnoaa_GeneID.csv")


######### Generating heatmap (COAD)

install.packages("pheatmap")
install.packages("RColorBrewer")

library("pheatmap")
library("RColorBrewer")
library("dplyr")

# Extract normalized expression for significant genes

coaddegbarcode <- c(COAD_Resist, COAD_sensitive)
write.csv(coaddegbarcode, "G:/My Drive/Year 2 (2021)/Bioinformatics/Rstudio/Rerun TCGA v2/TCGA Rerun V2/Heatmap/coaddegbarcode.csv")
coaddegbarcode.fixed <- read.csv("G:/My Drive/Year 2 (2021)/Bioinformatics/Rstudio/Rerun TCGA v2/TCGA Rerun V2/Heatmap/coaddegbarcode_fixed.csv")[,1]
norm.coad.a <- COADnorm[rownames(COAD_Degs), coaddegbarcode.fixed]


# Try heatmap

heat.colors <- brewer.pal(9, "YlOrRd")

pheatmap(norm.coad.a, color = heat.colors, cluster_rows = T, show_rownames = T, border_color = NA, fontsize = 3, scale = "row", fontsize_row = 3, height = 50)

# Try annotation
coad.sens.annot <- c("Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens", "Sens")
coad.res.annot <- c("Resist","Resist","Resist","Resist","Resist","Resist","Resist","Resist","Resist","Resist","Resist","Resist","Resist","Resist","Resist","Resist","Resist","Resist","Resist","Resist","Resist","Resist","Resist","Resist","Resist","Resist","Resist","Resist","Resist","Resist")
annot.coad <- data.frame(Group1 = c(coad.res.annot, coad.sens.annot))
row.names(annot.coad) <- coaddegbarcode.fixed


pheatmap(norm.coad.a, color = heat.colors, annotation = annot.coad, cluster_rows = T, cluster_cols = F, show_rownames = F, border_color = NA, fontsize = 10, scale = "row", fontsize_col = 5, height = 20)


# Heatmap with top 100 genes (Ranked by FDR) only

coad.100.deg <- read.csv("G:/My Drive/Year 2 (2021)/Bioinformatics/Rstudio/Rerun TCGA v2/TCGA Rerun V2/Heatmap/COAD_DEG_100.csv")[,1]
norm.coad.100 <- COADnorm[coad.100.deg, coaddegbarcode.fixed]
pheatmap(norm.coad.100, color = heat.colors, annotation = annot.coad, cluster_rows = T, cluster_cols = F, show_rownames = T, border_color = NA, fontsize = 5, scale = "row", fontsize_row = 3, height = 20, clustering_distance_rows = "maximum", clustering_distance_cols = "correlation")


##### Generating heatmap (READ)

# Extract normalized expression for significant genes

readdegbarcode <- c(Resistant_READ, Sensitive_READ)
write.csv(readdegbarcode, "G:/My Drive/Year 2 (2021)/Bioinformatics/Rstudio/Rerun TCGA v2/TCGA Rerun V2/Heatmap/readdegbarcode.csv")
readdegbarcode.fixed <- read.csv("G:/My Drive/Year 2 (2021)/Bioinformatics/Rstudio/Rerun TCGA v2/TCGA Rerun V2/Heatmap/readdegbarcode_fixed.csv")[,1]
norm.read.a <- READnorm[rownames(READ_raw_DEGs), readdegbarcode.fixed]

# Annotation

read.sens.annot <- c("Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens","Sens")
read.res.annot <- c("Res","Res","Res","Res","Res","Res","Res","Res","Res","Res")
annot.read <- data.frame(Group1 = c(read.res.annot, read.sens.annot))
rownames(annot.read) <- readdegbarcode.fixed

# Heatmap

heat.colors.read <- brewer.pal(9, "YlOrRd")
pheatmap(norm.read.a, color = heat.colors.read, cluster_rows = T, cluster_cols = T, show_rownames = T, border_color = NA, fontsize = 5, scale = "row", fontsize_row = 3, height = 20, clustering_distance_cols = "correlation", clustering_distance_rows = "maximum")

# Heatmap with annotation

pheatmap(norm.read.a, color = heat.colors.read, annotation = annot.read, cluster_rows = T, cluster_cols = T, show_rownames = T, border_color = NA, fontsize = 5, scale = "row", fontsize_row = 3, height = 20, clustering_distance_cols = "correlation", clustering_distance_rows = "maximum")

# Heatmap with annotation, without column clustering

pheatmap(norm.read.a, color = heat.colors.read, annotation = annot.read, cluster_rows = T, cluster_cols = F, show_rownames = F, border_color = NA, fontsize = 10, scale = "row", fontsize_col = 5, height = 20, clustering_distance_cols = "correlation", clustering_distance_rows = "maximum")

# Heatmap with top 100 genes (Ranked by FDR) only

read.100.deg <- read.csv("G:/My Drive/Year 2 (2021)/Bioinformatics/Rstudio/Rerun TCGA v2/TCGA Rerun V2/Heatmap/READ_DEG_100.csv")[,1]
norm.read.100 <- READnorm[read.100.deg, readdegbarcode.fixed]
pheatmap(norm.read.100, color = heat.colors.read, annotation = annot.read, cluster_rows = T, cluster_cols = F, show_rownames = T, border_color = NA, fontsize = 5, scale = "row", fontsize_row = 3, height = 20, clustering_distance_cols = "correlation", clustering_distance_rows = "maximum")




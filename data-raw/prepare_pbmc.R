# Setup ----
pacman::p_load(tidyverse, Seurat, SeuratData, this.path)
this.path::here() %>% setwd()
rm(list=ls())


# Load data ----
InstallData("pbmc3k")
pbmc.demo <- LoadData("pbmc3k")


# Add artificial batches for demonstration ----
pbmc.demo$batch <- paste0("B", 1:3)
table(pbmc.demo$batch)


# Gross quality control ----
pbmc.demo$percent.mt <- PercentageFeatureSet(pbmc.demo, pattern = "^MT-|^Mt-")

VlnPlot(pbmc.demo, c("nCount_RNA", "nFeature_RNA", "percent.mt"))

pbmc.demo <- subset(pbmc.demo,
                      nFeature_RNA > 200 &
                      nFeature_RNA < 2500 &
                      percent.mt   < 5)

VlnPlot(pbmc.demo, c("nCount_RNA", "nFeature_RNA", "percent.mt"))

pbmc.demo  # 2638 samples


# Ground truth ----
table(pbmc.demo$seurat_annotations, exclude=NULL)
#  Naive CD4 T Memory CD4 T   CD14+ Mono            B        CD8 T
#          697          483          480          344          271
# FCGR3A+ Mono           NK           DC     Platelet
#          162          155           32           14

table(pbmc.demo$seurat_annotations, pbmc.demo$batch, exclude=NULL)
#               B1  B2  B3
# Naive CD4 T  219 248 230
# Memory CD4 T 174 133 176
# CD14+ Mono   162 166 152
# B            127 111 106
# CD8 T         92  87  92
# FCGR3A+ Mono  55  56  51
# NK            44  63  48
# DC             7  12  13
# Platelet       4   4   6


# Save ----
usethis::use_data(pbmc.demo, overwrite = TRUE, compress = "xz")

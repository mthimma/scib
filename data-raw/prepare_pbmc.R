pacman::p_load(tidyverse, Seurat, SeuratData, this.path)
this.path::here() %>% setwd()

## Load data
InstallData("pbmc3k")
pbmc.demo <- LoadData("pbmc3k")
table(pbmc.demo$seurat_annotations, exclude=NULL)

pbmc.demo <- pbmc.demo[ , !is.na(pbmc.demo$seurat_annotations)]
table(pbmc.demo$seurat_annotations, exclude=NULL)
pbmc.demo  # 2638 samples

## Add artificial batches for demonstration
pbmc.demo$batch <- paste0("B", 1:3)

usethis::use_data(pbmc.demo, overwrite = TRUE, compress = "xz")

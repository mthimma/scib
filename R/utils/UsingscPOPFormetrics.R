## Using scPOP for scmetrics
# install.packages("devtools")
devtools::install_github("vinay-swamy/scPOP")  # from GitHub
devtools::install_github('satijalab/seurat-data')
pacman::p_load(scPOP, Seurat, SeuratData, tidyverse, patchwork)
AvailableData()
#InstallData("pbmc3k") this didnt work, hence installed from package as given below

#install.packages("C://Users/kumar/Downloads/pbmc3k.SeuratData_3.1.4.tar.gz", repos=NULL, type="source")
#InstallData("pbmc3k")
#LoadData("pbmc3k")
#data("pbmc3k")

pbmc <- pbmc3k.SeuratData::pbmc3k
pbmc.updated = UpdateSeuratObject(object = pbmc)
pbmc.updated@images <- list()

pbmc.updated[["percent.mt"]] <- PercentageFeatureSet(pbmc.updated, pattern = "^MT-")
pbmc1 <- pbmc.updated

VlnPlot(pbmc1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

pbmc1 <- subset(pbmc1, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

pbmc1 <- NormalizeData(pbmc1)

pbmc1 <- FindVariableFeatures(pbmc1, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(pbmc1)
pbmc1 <- ScaleData(pbmc1, features = all.genes)

pbmc1 <- RunPCA(pbmc1, features = VariableFeatures(object = pbmc1))

# use PCA embeddings as before
emb <- Embeddings(pbmc1, "pca")[, 1:20]

cell_labels <- pbmc1$seurat_annotations   # biological labels
#batch_labels <- pbmc$batch            # technical batches

# ensure they are factors or character
cell_labels <- as.factor(cell_labels)
batch_labels <- as.factor(batch_labels)

load("C://Users/kumar/Downloads/sceiad_subset_data.rda")
emb1 <- sceiad_subset_data %>% 
  select(Barcode, starts_with("scviDim")) %>% 
  column_to_rownames("Barcode")
metdata <- sceiad_subset_data %>% 
           select(Barcode, batch, cluster, subcluster, CellType, CellType_predict) %>% 
           column_to_rownames("Barcode")

metrics <- run_all_metrics(emb1, 
                           metadata = metdata,
                           batch_key = 'batch',
                           label1_key = 'CellType_predict',
                           #label2_key = 'cluster', 
                           run_name = 'example_fromscPOP')

metrics

# You’ll get a data frame including (names may look like):
#   
#   ASW_bio
# 
# ASW_batch
# 
# ARI, NMI, LISI…
# 
# These are exactly the metrics used in batch‑correction benchmarks for scRNA‑seq.
# 
# If your version uses different function names (e.g. asw_batch() / asw_bio()), the logic is the same: pass the embedding, batch labels, and biological labels.

library(Seurat)
library(scPOP)
library(cluster)

# 1. Integrate your datasets (e.g., with Seurat, Harmony, etc.)
# ... produce an integrated Seurat object, say 'obj_int'

# 2. Get embeddings, labels, and batch
emb <- Embeddings(obj_int, "pca")[, 1:30]
cell_labels <- obj_int$celltype  # or clusters
batch_labels <- obj_int$batch

# 3. ASW manually (for intuition)
dist_mat <- dist(emb)
asw_bio  <- mean(silhouette(as.integer(as.factor(cell_labels)), dist_mat)[, "sil_width"])
asw_batch <- mean(silhouette(as.integer(as.factor(batch_labels)), dist_mat)[, "sil_width"])

asw_bio
asw_batch

# 4. ASW + other metrics via scPOP
metrics <- scPOP::batch_metrics(
  embeddings = emb,
  batch = batch_labels,
  label = cell_labels
)

metrics


## Setup
setwd("C:/Manjula/scIntegration_local/datasets/MouseEmbryonicDevelopment/")

pacman::p_load(tidyverse, janitor, broom, patchwork, stringi,
               Seurat, reticulate, scib, ggplot2, aricode)
source("../../scripts/PlotUtilities.R")

rm(list=ls())

## Read data in 
load("data/MED_Intestines.rda")
#load("./data/GepLiver_NAFLD.rda")

dataset_name <- "MED_Intestines"
resultdir <- paste0("results/", dataset_name, "/")

dim(annot) # 53032  7


######################
## Find Optimal nPCs ##
######################

#gt_numclust <- length(unique(annot$Cell_type))

## Create seurat object with count and annot
obj <- CreateSeuratObject(counts = cts, 
                          meta.data = annot,
                          min.cells = 3, 
                          min.features = 200,
                          project = dataset_name)

obj <- DietSeurat(obj, assays = "RNA")
obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)

obj <- RunPCA(obj, npcs = 50)

v <- obj@reductions$pca@stdev^2
nPCs_50 <- which(cumsum( v/sum(v) ) > 0.80)[1]
opt_nPCs  <- nPCs_50
ggsave(paste0(resultdir, dataset_name, "_ElbowPlot_50pcs.png"),
       plot = ElbowPlot(obj, ndims=50) + geom_vline(xintercept = nPCs_50, lty=3))

## for first 50 pcs
df <- data.frame(coords = 1:50, vars = cumsum( v/sum(v)))

p <- pca_plot(df, dataset_name, opt_nPCs)

ggsave(paste0(resultdir, dataset_name,"_LinePlot_optPCs_First50pcs_testing.pdf"), plot = p, width = 8, height = 6)

rm(cts, annot, v, df, p) 
################################################
## Find Optimal clustering resolution opt_res ##
################################################

## Find NN and SNN
obj <- FindNeighbors(obj, dims = 1:opt_nPCs)
verbose = FALSE
## Community detection using Louvain
#res_try  <- c(0.001, 0.002, 0.003, 0.004, 0.005, 0.006,0.008, 0.009,)
res_try  <- c(0.001, 0.002, 0.003, 0.004, 0.005, 0.006,0.008, 0.009,0.01, 0.02, 0.03, 0.035, 0.04,0.045, 0.05, 0.10)
              #, 0.11, 0.12, 0.125, 0.25, 0.3, 0.35, 0.5, 0.45)

#g1 <- list()
gc()

obj$Cell_type2 <- stri_replace_all_regex(obj$Cell_type,
                                  pattern=c("Extraembryonic visceral endoderm","Intestinal enteroendocrine cells","Intestinal goblet cells", "Midgut/Hindgut epithelial cells" ),
                                  replacement=c("Extraembryonic\nvisceral\nendoderm", "Intestinal enteroendocrine\ncells", "Intestinal goblet\ncells", "Midgut/Hindgut\nepithelial\ncells" ),
                                  vectorize=FALSE)

obj$Cell_type2 <- factor(obj$Cell_type2)

r <- 0.02

for (r in res_try) {
  obj <- FindClusters(obj, resolution = r, verbose = verbose)
  obj <- RunUMAP(obj, reduction = "pca", dims = 1:opt_nPCs)
  aa <- DimPlot_highlight(obj, group.by = "Cell_type2", label.size = 4)
  cn <- DimPlot_highlight(obj, group.by = paste0("RNA_snn_res.", r))
  out_all <- aa[["all"]] | cn[["all"]]
  ggsave(file = paste0(resultdir, dataset_name, "Umap_global_test.png"), plot = out_all, height=12, width=12)
  
  author_all <- aa[["all"]]
  author_visceral <- aa[["Extraembryonic\nvisceral\nendoderm"]]
  author_enter    <- aa[["Intestinal enteroendocrine\ncells"]]
  author_goblet   <- aa[["Intestinal goblet\ncells"]]
  author_epi      <- aa[["Midgut/Hindgut\nepithelial\ncells"]]
  
  ggsave(file=paste0(resultdir, dataset_name, "Umap_global_authorcelltypes.png"), plot = author_all, height=12, width=12)
  ggsave(file=paste0(resultdir, dataset_name, "Umap_authorcelltypes_visceral.png"), plot = author_visceral, height=12, width=12)
  ggsave(file=paste0(resultdir, dataset_name, "Umap_authorcelltypes_entero.png"), plot = author_enter, height=12, width=12)
  ggsave(file=paste0(resultdir, dataset_name, "Umap_authorcelltypes_goblet.png"), plot = author_goblet, height=12, width=12)
  ggsave(file=paste0(resultdir, dataset_name, "Umap_authorcelltypes_epithelial.png"), plot = author_epi, height=12, width=12)
  
  predicted_all <- cn[["all"]]
  predicted_0   <- cn[["0"]]
  predicted_1   <- cn[["1"]]
  predicted_2   <- cn[["2"]]
  
  ggsave(file=paste0(resultdir, dataset_name, "Umap_global_predcelltypes.png"), plot = predicted_all, height=12, width=12)
  ggsave(file=paste0(resultdir, dataset_name, "Umap_predcelltypes_0.png"), plot = predicted_0, height=12, width=12)
  ggsave(file=paste0(resultdir, dataset_name, "Umap_predcelltypes_1.png"), plot = predicted_1, height=12, width=12)
  ggsave(file=paste0(resultdir, dataset_name, "Umap_predcelltypes_2.png"), plot = predicted_2, height=12, width=12)
  
  ## getting author cell type and predicted clusters global
  ## visualising merging/split of clusters
  #cluster1 <- (aa[["Intestinal enteroendocrine\ncells"]] /  aa[["Intestinal goblet\ncells"]])
  #out_1    <-   cluster1 | cn[["1"]] + plot_layout(widths = c(2,1))
  #ggsave(file = paste0(resultdir, dataset_name, "Umap_cluster1_test.png"), plot = out_1, height=12, width=12)
}

rm(aa, cn, cluster0, cluster1, out_1, out_all)

clusters <- obj@meta.data[ , paste0("RNA_snn_res.", res_try)]
clusters <- clusters %>% remove_constant()

res_ari <- apply(clusters, 2, function(cl) mclust::adjustedRandIndex(obj$Cell_type, cl))

res_nmi <- apply(clusters, 2, function(cl) aricode::NMI(obj$Cell_type, cl))

#numcells <- nrow(annot)

res_asw <- apply(clusters, 2, function(cl){
  if (nrow( obj@reductions$umap@cell.embeddings)  > 65536) {
    metrics_ASW( cl,  obj@reductions$umap@cell.embeddings[1:65536, ])
  } else {
    metrics_ASW( cl, obj@reductions$umap@cell.embeddings)
  }
})

#res_lisi_pca  <- apply(clusters, 2, function(cl) 1 - metrics_LISI(obj$Cell_type, obj@reductions$pca@cell.embeddings)$lisi_mean)

res_lisi_umap  <- apply(clusters, 2, function(cl) 1 - metrics_LISI(obj$Cell_type, obj@reductions$umap@cell.embeddings)$lisi_mean)

res_asw_pca <- apply(clusters, 2,function(cl) {
  if (nrow( obj@reductions$umap@cell.embeddings) > 65536) {
    metrics_ASW(cl, obj@reductions$pca@cell.embeddings[1:65536, ])
  } else {
    metrics_ASW(cl, obj@reductions$pca@cell.embeddings[1:2000, ])
  }
}
)

res <- data.frame( res = gsub("RNA_snn_res.", "", names(res_ari)) %>% as.numeric(), 
                   ari = res_ari,
                   nmi = res_nmi,
                   #asw = res_asw,
                   #aswpca = res_asw_pca,
                   #lisipca = res_lisi_pca,
                   lisi = res_lisi_umap)

## save metrics table for selection for opt_res
write.table( res, file=paste0(resultdir, dataset_name,"_comparemetrics_new.csv"),row.names = FALSE,
             sep=",", quote=FALSE)

##Get Contingency table for Author cell types Vs Clusers of unintegrated
table(obj$Cell_type)
# out <- aricode::sortPairs(obj$RNA_snn_res.0.02, obj$Cell_type, spMat = TRUE)
# out
cont_tbl <- table(obj$Cell_type, obj$RNA_snn_res.0.02)
tmp <- addmargins(cont_tbl)
#t#marginSums(t, 1) # grand total
#rowwise proportions in each column
proportions(t, 1)

#table(obj$Cell_type, obj$RNA_snn_res.0.01)
#sorted_data <- aricode::sortPairs(obj$Cell_type, obj$RNA_snn_res.0.01)
#sorted_matrix <- as.matrix(sorted_data$spMat) 

write.table(tmp, 
            file=paste0(resultdir, dataset_name,"_AuthorCelltypeClustersContingency_detailed2.csv"), row.names = TRUE,
            sep=',', quote = FALSE)

#g_scatter <- saveScatterPlot(res, paste0(getwd(), "results/MED_CardioMyocytes/", dataset_name))

## Get line plot of metrics used to select optimal clustering resolution
opt_res <- 0.02
g1 <- saveLinePlot(res, "ari", opt_res)
g2 <- saveLinePlot(res, "nmi", opt_res)

ggsave(paste0(resultdir, dataset_name,"_ARI_NMI_MetricsforOptRes_Lineplot_improved.pdf"), width = 10, height = 10,
        plot = g1 /g2 )

#rm(res_try, clusters, res_ari, res_asw, res, g1, g2)
#rm(list=ls())

#######################################################
## Create stratification and prepare for integration ##
#######################################################

#cts <- obj@assays$RNA@counts

#annot <- obj@meta.data %>% 
#   select(Sample, Sample_Origin, Cell_type) %>% 
#   mutate( nbatch_3  = stratified_sampling(Cell_type, 3),
#           nbatch_5  = stratified_sampling(Cell_type, 5),
#           nbatch_10 = stratified_sampling(Cell_type, 10))

#opt_nPCs <- nPCs_50 #12 for NAFLD Lung only 19L for MED_Intenstines
#opt_res  <- 0.02 # for MED majorcellcluster Cardiomyocytes

save(cts, annot, opt_nPCs, opt_res, obj, file = "data/MED_mcc_Cardiomyocytes.rda", compress = TRUE)

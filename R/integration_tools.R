BBKNNIntegration  <- function(object, features = NULL, layers = "counts", conda_env = NULL,
                              new.reduction = "bbknn", ndims = 30, ...) {

  reticulate::use_condaenv(conda_env, required = TRUE)

  bbknn   <- reticulate::import('bbknn')
  anndata <- reticulate::import("anndata")
  sc      <- reticulate::import("scanpy")
  scipy   <- reticulate::import("scipy")
  ig      <- reticulate::import("igraph")

  print(layers)

  batches <- SeuratWrappers:::.FindBatches(object, layers = layers)
  object <- JoinLayers(object = object, layers = "counts")

  ## BBKNN takes all genes as input, hence no need to subset to the HVG
  adata <- sc$AnnData(X   = scipy$sparse$csr_matrix(Matrix::t(LayerData(object, layer = "counts"))),
                      obs = batches)
  adata$var_names <- rownames(object)
  rm(object, batches); suppressMessages( gc(verbose = FALSE) )

  sc$tl$pca(adata, svd_solver="arpack")

  adata_bbknn <- bbknn$bbknn(adata, copy=T, trim=0, batch_key='batch')
  latent <- adata_bbknn$obsm[['X_pca']][, 1:ndims]
  rownames(latent) <- adata_bbknn$obs_names$to_list()
  rm(adata, adata_bbknn); suppressMessages( gc(verbose = FALSE) )

  suppressWarnings(latent.dr <- CreateDimReducObject(embeddings = latent,
                                                     key = new.reduction))
  output.list <- list(latent.dr)
  names(output.list) <- c(new.reduction)
  return(output.list)
}

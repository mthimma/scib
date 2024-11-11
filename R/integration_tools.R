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


CONOSIntegration <- function(object, resolution = NULL,  features = NULL, layers = "counts",
                             new.reduction = "conos", ndims = NULL, verbose = FALSE, ...) {

  if(is.null(resolution)) stop("resolution must be specified")
  if(is.null(ndims))      stop("ndims must be specified")

  ds_new <- list()
  batches <- gsub("counts.|data.", "", layers)  %>% sort()

  for(batch in batches){

    out <- GetAssayData(object, paste0("counts.", batch)) %>%
      CreateSeuratObject()

    out <- SetAssayData(out, layer = "data", GetAssayData(object, paste0("data.", batch)))

    out <- out %>%
      FindVariableFeatures(verbose = verbose) %>%
      ScaleData(verbose = verbose) %>%
      RunPCA(npcs = ndims, verbose = verbose) %>%
      FindNeighbors(dims = 1:ndims, verbose = verbose) %>%
      FindClusters(resolution = resolution, verbose = verbose)

    ds_new[[ batch ]] <- out
    rm(out, batch); gc()
  }

  require(conos)
  con <- Conos$new(ds_new, n.cores=4)
  con$buildGraph(k=20, k.self=30, space='PCA', ncomps=ndims, n.odgenes=2000, matching.method='mNN',
                 metric='angular', score.component.variance=TRUE, verbose=verbose)
  con$findCommunities(method=leiden.community, resolution=resolution, verbose=verbose)
  con$embedGraph(min.dist=0.01, spread=15, min.prob.lower=1e-3, target.dims = ndims)

  suppressWarnings(dr <- CreateDimReducObject(embeddings = con$embedding,
                                              key = new.reduction))
  output.list <- list(dr)
  names(output.list) <- c(new.reduction)
  return(output.list)
}

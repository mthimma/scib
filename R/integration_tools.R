IntegrateLayers2 <- function(object,
                             method,
                             orig.reduction = "pca",
                             new.reduction,
                             verbose        = FALSE,
                             resolution     = 0.80,
                             ndims          = ndims,
                             conda_env      = "/home/manjula/miniforge3/envs/sc" ) {
  # original methods
  if(method %in% c("HarmonyIntegration", "CCAIntegration", "RPCAIntegration", "FastMNNIntegration") ){
    out <- IntegrateLayers(object = object, method = method, orig.reduction = orig.reduction, new.reduction = new.reduction)
  }
  if (method %in% c("BBKNNIntegration", "SCANORAMAIntegration", "scVIIntegration2")) {
   out <-  IntegrateLayers(object = object, method = method,  orig.reduction = orig.reduction, new.reduction = new.reduction,
                    conda_env = conda_env)
  }

   if(method %in% c("CONOSIntegration", "LIGERIntegration","SCMCIntegration") ){
    out <- IntegrateLayers(object = object, method = method, orig.reduction = orig.reduction, new.reduction = new.reduction,
                           resolution = resolution, ndims = ndims)
  }

  return(out)
}


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



scVIIntegration2 <- function(object, features = NULL, layers = "counts",
                             conda_env = NULL, new.reduction = "integrated.dr", ndims = 30,
                             nlayers = 2, gene_likelihood = "nb", max_epochs = NULL, num_workers = 7, ...){

  reticulate::use_condaenv(conda_env, required = TRUE)
  sc <- reticulate::import("scanpy", convert = FALSE)
  scvi <- reticulate::import("scvi", convert = FALSE)
  scvi$settings$dl_num_workers <- as.integer(num_workers)
  scvi$settings$num_threads    <- as.integer(num_workers)

  anndata <- reticulate::import("anndata", convert = FALSE)
  scipy <- reticulate::import("scipy", convert = FALSE)
  if (is.null(max_epochs)) {
    max_epochs <- reticulate::r_to_py(max_epochs)
  } else {
    max_epochs <- as.integer(max_epochs)
  }

  batches <- SeuratWrappers:::.FindBatches(object, layers = layers)

  object <- JoinLayers(object = object, layers = "counts")
  adata <- sc$AnnData(X = scipy$sparse$csr_matrix(Matrix::t(LayerData(object, layer = "counts")[features, ])),
                      obs = batches,
                      var = object[[]][features,
                      ])
  scvi$model$SCVI$setup_anndata(adata, batch_key = "batch")
  model <- scvi$model$SCVI(adata = adata, n_latent = as.integer(x = ndims),
                           n_layers = as.integer(x = nlayers), gene_likelihood = gene_likelihood)
  model$train(max_epochs = max_epochs)
  latent <- model$get_latent_representation()
  latent <- as.matrix(latent)
  rownames(latent) <- reticulate::py_to_r(adata$obs$index$values)
  colnames(latent) <- paste0(new.reduction, "_", 1:ncol(latent))
  suppressWarnings(latent.dr <- CreateDimReducObject(embeddings = latent,
                                                     key = new.reduction))
  output.list <- list(latent.dr)
  names(output.list) <- new.reduction
  return(output.list)
}

SCMCIntegration  <- function(object, features = NULL, layers = "counts", new.reduction = "scmc",
                             ndims = ndims, resolution = resolution, ...){

  object <- SplitObject(object, split.by = "batch")

  ## Step2. Identify putative clusters for each dataset
  object.list <- identifyNeighbors(object, force.pca=FALSE, nDims.pca=ndims)
  object.list <- identifyClusters(object.list, resolution = resolution, nDims.consensus = ndims)

  ## Step3. Detect cluster-specific cells with high confident
  features.integration <- identifyIntegrationFeatures(object.list)
  object.list <- identifyConfidentCells(object.list, features.integration)

  ## Step4. Identify marker genes associated with the putative cell clusters in each dataset
  object.list <- identifyMarkers(object.list)

  ## Step5. Learn technical variation between any two datasets
  structured_mat <- learnTechnicalVariation(object.list, features.integration)

  ## Step 6. Learn a shared embedding of cells across all datasets after removing technical variation
  combined <- merge(object.list[[1]], object.list[2:length( object.list)])
  combined@meta.data <- combined@meta.data %>% select(-contains("_snn_res"))
  VariableFeatures(combined) <- features.integration

  combined <- JoinLayers(combined)
  combined <- scMC::integrateData(combined, structured_mat, nDims.scMC = ndims)
  # output.list <- list(Embeddings(combined, reduction = "scMC"))
  # names(output.list) <- new.reduction
   return(combined)
}



## LIGER
LIGERIntegration <- function(object, ndims = NULL, new.reduction = "liger", ...) {

  if(is.null(ndims)) stop("ndims must be supplied")

  obj <- object %>%
    CreateSeuratObject() %>%
    rliger::normalize(layer = "counts") %>%
    rliger::selectGenes() %>%
    rliger::scaleNotCenter()

  obj <- obj %>%
    runINMF(k = ndims) %>%
    quantileNorm()

  latent <- as.matrix(obj@reductions$inmfNorm@cell.embeddings)

  colnames(latent) <- paste0(new.reduction, "_", 1:ncol(latent))

  suppressWarnings(latent.dr <- CreateDimReducObject(embeddings = latent,
                                                     key = new.reduction))
  output.list <- list(latent.dr)
  names(output.list) <- c(new.reduction)
  return(output.list)
}

SCANORAMAIntegration  <- function( object, ndims = ndims, new.integration = "integrated.scanorama", layers="counts",
                                   features = NULL, conda_env = conda_env,...) {
  reticulate::use_condaenv(conda_env, required = TRUE)
  scr    <- reticulate::import('scanorama')
  scipy  <- reticulate:: import('scipy', convert = FALSE)
  datasets <- list()
  genes_list <- list()
  features <- VariableFeatures(object)
  datasets <- lapply(X = layers, FUN = function(layer) {
    layer.data <- LayerData(object = object, layer = layer, features = features)
    scipy$sparse$csr_matrix(Matrix::t(layer.data))
  })
  genes_list <- Features(object, layer=layers)

  out <- scr$correct(datasets, genes_list,
                     return_dimred=TRUE, return_dense=TRUE)
  combined <- t( do.call("rbind", out[[2]]) )

  dimnames(combined) <- list(out[[3]], unlist( sapply(datasets, rownames) ))

  emb_pca <- do.call(rbind, out[[1]])
  dimnames(emb_pca) <- list( unlist( sapply(datasets, rownames) ), paste0("PC_", 1:ncol(emb_pca)) )

  ## Create Seurat object
  combined <- CreateSeuratObject(counts=combined, assay="pano")

  combined[["pca"]] <- CreateDimReducObject(embeddings=emb_pca,
                                            stdev=apply(emb_pca, 2, sd),
                                            key="PC_", assay="pano")
  output.list <- list(Embeddings(combined, reduction = "pca"))
  names(output.list) <- new.reduction
  return(output.list)
}



IntegrateLayers2 <- function(object,
                             method,
                             orig.reduction = "pca",
                             new.reduction,
                             verbose        = FALSE,
                             resolution     = NULL,
                             ndims          = NULL,
                             conda_env      = NULL) {

  if(method == "scVIIntegration") method <- "scVIIntegration2"

  if(method %in% c("HarmonyIntegration", "CCAIntegration", "RPCAIntegration", "FastMNNIntegration") ){
    # methods without conda env
    out <- IntegrateLayers(object = object, method = method, orig.reduction = orig.reduction, new.reduction = new.reduction)
  }

  if (method %in% c("BBKNNIntegration", "SCANORAMAIntegration", "scVIIntegration2")) {
    # methods with conda env specified
    if( is.null(conda_env) ) stop("conda_env must be specified for ", method)
    out <-  IntegrateLayers(object = object, method = method, orig.reduction = orig.reduction, new.reduction = new.reduction,
                            conda_env = conda_env)
  }

  if(method %in% c("CONOSIntegration", "LIGERIntegration","SCMCIntegration") ){
    # methods which require number of PC and resolution parameters
    if( is.null(ndims) )      stop("ndims must be specified for ", method)
    if( is.null(resolution) ) stop("resolution must be specified for ", method)

    out <- IntegrateLayers(object = object, method = method, orig.reduction = orig.reduction, new.reduction = new.reduction,
                           resolution = resolution, ndims = ndims)

  }

  return(out)
}

IntegrateLayers_input <- function (object, method, orig.reduction = "pca", assay = NULL,
          features = NULL, layers = NULL, scale.layer = "scale.data",
          ...) {

  if (rlang::is_quosure(x = method)) {
    method <- eval(expr = quo_get_expr(quo = method), envir = quo_get_env(quo = method))
  }

  if (is.character(x = method)) {
    method <- get(x = method)
  }

  if (!is.function(x = method)) {
    abort(message = "'method' must be a function for integrating layers")
  }

  assay <- assay %||% DefaultAssay(object = object)

  if (inherits(x = object[[assay]], what = "SCTAssay")) {
    layers <- "data"
    scale.layer <- "scale.data"
    features <- features %||% SelectSCTIntegrationFeatures(object = object,
                                                           assay = assay)
  } else if (inherits(x = object[[assay]], what = "StdAssay")) {
    layers <- Layers(object = object, assay = assay, search = layers %||%
                       "data")
    scale.layer <- Layers(object = object, search = scale.layer)
    features <- features %||% VariableFeatures(object = object,
                                               assay = assay, nfeatures = 2000L)
  } else {
    abort(message = "'assay' must be a v5 or SCT assay")
  }

  if (!is.null(scale.layer)) {
    features <- intersect(x = features, y = Features(x = object,
                                                     assay = assay, layer = scale.layer))
  }

  if (!length(x = features)) {
    abort(message = "None of the features provided are found in this assay")
  }

  if (!is.null(orig.reduction)) {

    orig.reduction <- orig.reduction %||% DefaultDimReduc(object = object, assay = assay)

    if (!orig.reduction %in% Reductions(object = object)) {
      abort(message = paste(sQuote(x = orig.reduction),
                            "is not a dimensional reduction"))
    }

    obj.orig <- object[[orig.reduction]]

    if (is.null(x = DefaultAssay(object = obj.orig))) {
      DefaultAssay(object = obj.orig) <- assay
    }

  }

  out <- list(object = object[[assay]], assay = assay,
              orig = obj.orig, layers = layers, scale.layer = scale.layer,
              features = features, ...)

  return(out)
}


BBKNNIntegration  <- function(object, features = NULL, layers = "counts", conda_env,
                              new.reduction = "bbknn", ndims = 30, ...) {

  use_condaenv(conda_env, required = TRUE)

  bbknn   <- import('bbknn')
  anndata <- import("anndata")
  sc      <- import("scanpy")
  scipy   <- import("scipy")
  ig      <- import("igraph")

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
                             new.reduction = "conos", ndims = NULL, verbose = FALSE, n.cores = 1, ...) {

  if(is.null(resolution)) stop("resolution must be specified")
  if(is.null(ndims))      stop("ndims must be specified")
  require(conos)

  ds_new  <- list()
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

  con <- Conos$new(ds_new, n.cores=n.cores)
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


SCMCIntegration <- function(object, features = NULL, layers = "counts", new.reduction = "scmc",
                             ndims = ndims, resolution = resolution, ...){

  require(scMC)

  ds_new  <- list()
  batches <- gsub("counts.|data.", "", layers)  %>% sort()

  for(batch in batches){

    out <- GetAssayData(object, paste0("counts.", batch)) %>%
      CreateSeuratObject()

    out <- SetAssayData(out, layer = "data", GetAssayData(object, paste0("data.", batch)))

    out <- out %>%
      FindVariableFeatures(verbose = FALSE) %>%
      ScaleData(verbose = FALSE) %>%
      RunPCA(npcs = ndims, verbose = FALSE)

    ds_new[[ batch ]] <- out
    rm(out, batch); gc()
  }

  ## Step2. Identify putative clusters for each dataset
  object.list <- identifyNeighbors(ds_new, force.pca=FALSE, nDims.pca=ndims)
  object.list <- identifyClusters (object.list, resolution = resolution, nDims.consensus = ndims)

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

  latent <- as.matrix(Embeddings(combined, reduction = "scMC"))
  suppressWarnings(latent.dr <- CreateDimReducObject(embeddings = latent,
                                                     key = new.reduction))
  output.list <- list(latent.dr)
  names(output.list) <- new.reduction
  return(output.list)
}


LIGERIntegration <- function(object, ndims = NULL, new.reduction = "liger", ...) {

  if(is.null(ndims)) stop("ndims must be supplied")
  require(rliger)

  obj <- object %>%
    CreateSeuratObject() %>%
    normalize(layer = "counts") %>%
    selectGenes() %>%
    scaleNotCenter()

  obj <- obj %>%
    runINMF(k = ndims) %>%
    rliger::quantileNorm()

  latent <- as.matrix(obj@reductions$inmfNorm@cell.embeddings)

  colnames(latent) <- paste0(new.reduction, "_", 1:ncol(latent))

  suppressWarnings(latent.dr <- CreateDimReducObject(embeddings = latent,
                                                     key = new.reduction))
  output.list <- list(latent.dr)
  names(output.list) <- c(new.reduction)
  return(output.list)
}


SCANORAMAIntegration  <- function(object, layers, new.reduction = new.reduction,
                                  features = NULL, conda_env = conda_env, ...) {

  use_condaenv(conda_env, required = TRUE)
  scr   <- import('scanorama')
  scipy <- import('scipy')
  ad    <- import('anndata')

  adatas <- lapply(X = layers, FUN = function(layer) {

    layer.data <- LayerData(object = object, layer = layer, features = features)

    adata = ad$AnnData( X = scipy$sparse$csr_matrix(Matrix::t(layer.data)) )
    adata$obs_names <- colnames(layer.data)
    adata
  })

  scr$integrate_scanpy(adatas, dimred = 100L)

  latent <- lapply( adatas, function(adata){
    out <- adata$obsm['X_scanorama']
    rownames(out) <- adata$obs_names$to_list()
    return(out)
  })

  latent <- do.call("rbind", latent)

  colnames(latent) <- paste0(new.reduction, "_", 1:ncol(latent))

  suppressWarnings(latent.dr <- CreateDimReducObject(embeddings = latent,
                                                     key = new.reduction))
  output.list <- list(latent.dr)
  names(output.list) <- c(new.reduction)
  return(output.list)
}


scVIIntegration2 <- function(object, features = NULL, layers = "counts",
                             conda_env = NULL, new.reduction = "integrated.dr", ndims = 30,
                             nlayers = 2, gene_likelihood = "nb", max_epochs = 10, num_workers = 7, ...){

  ## modified from SeuratWrapper::scVIIIntegration to incorporate num_workers
  use_condaenv(conda_env, required = TRUE)
  sc   <- import("scanpy", convert = FALSE)
  scvi <- import("scvi", convert = FALSE)
  scvi$settings$dl_num_workers <- as.integer(num_workers)
  scvi$settings$num_threads    <- as.integer(num_workers)

  anndata <- import("anndata", convert = FALSE)
  scipy <- import("scipy", convert = FALSE)
  if (is.null(max_epochs)) {
    max_epochs <- r_to_py(max_epochs)
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
  rownames(latent) <- py_to_r(adata$obs$index$values)
  colnames(latent) <- paste0(new.reduction, "_", 1:ncol(latent))
  suppressWarnings(latent.dr <- CreateDimReducObject(embeddings = latent,
                                                     key = new.reduction))
  output.list <- list(latent.dr)
  names(output.list) <- new.reduction
  return(output.list)
}


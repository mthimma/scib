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


CCAIntegration2 <- function(
    object = NULL,
    assay = NULL,
    layers = NULL,
    orig = NULL,
    new.reduction = 'integrated.dr',
    reference = NULL,
    features = NULL,
    normalization.method = c("LogNormalize", "SCT"),
    dims = 1:30,
    k.filter = NA,
    scale.layer = 'scale.data',
    dims.to.integrate = NULL,
    k.weight = 100,
    weight.reduction = NULL,
    sd.weight = 1,
    sample.tree = NULL,
    preserve.order = FALSE,
    verbose = TRUE,
    ...
) {
  op <- options(Seurat.object.assay.version = "v3", Seurat.object.assay.calcn = FALSE)
  on.exit(expr = options(op), add = TRUE)
  normalization.method <- match.arg(arg = normalization.method)
  features <- features %||% SelectIntegrationFeatures5(object = object)
  assay <- assay %||% 'RNA'
  layers <- layers %||% Layers(object, search = 'data')
  if (normalization.method == 'SCT') {
    #create grouping variables
    groups <- CreateIntegrationGroups(object, layers = layers, scale.layer = scale.layer)
    object.sct <- CreateSeuratObject(counts = object, assay = 'SCT')
    object.sct$split <- groups[,1]
    object.list <- SplitObject(object = object.sct,split.by = 'split')
    object.list  <- PrepSCTIntegration(object.list, anchor.features = features)
  } else {
    object.list <- list()
    for (i in seq_along(along.with = layers)) {
      if (inherits(x = object[layers[i]], what = "IterableMatrix")) {
        warning("Converting BPCells matrix to dgCMatrix for integration ",
                "as on-disk CCA Integration is not currently supported", call. = FALSE, immediate. = TRUE)
        counts <- as(object = object[layers[i]][features, ],
                     Class = "dgCMatrix")
      }
      else {
        counts <- object[layers[i]][features, ]
      }
      object.list[[i]] <- CreateSeuratObject(counts = counts)
      if (inherits(x = object[scale.layer], what = "IterableMatrix")) {
        scale.data.layer <- as.matrix(object[scale.layer][features,
                                                          Cells(object.list[[i]])])
        object.list[[i]][["RNA"]]$scale.data <- scale.data.layer
      }
      else {
        object.list[[i]][["RNA"]]$scale.data <- object[scale.layer][features, Cells(object.list[[i]])]
      }
      object.list[[i]][['RNA']]$counts <- NULL
    }
  }

  anchor <- FindIntegrationAnchors(object.list = object.list,
                                   anchor.features = features,
                                   scale = FALSE,
                                   reduction = 'cca',
                                   normalization.method = normalization.method,
                                   dims = dims,
                                   k.filter = k.filter,
                                   reference = reference,
                                   verbose = verbose                        # removed the ... here
  )

  suppressWarnings({
    anchor@object.list <- lapply(anchor@object.list, function(x) {
      x <- DietSeurat(x, features = features[1:2])
      return(x)
    })
  }, classes = "dimWarning")

  object_merged <- IntegrateEmbeddings(anchorset = anchor,
                                       reductions = orig,
                                       new.reduction.name = new.reduction,
                                       dims.to.integrate = dims.to.integrate,
                                       k.weight = k.weight,
                                       weight.reduction = weight.reduction,
                                       sd.weight = sd.weight,
                                       sample.tree = sample.tree,
                                       preserve.order = preserve.order,
                                       verbose = verbose
  )
  output.list <- list(object_merged[[new.reduction]])
  names(output.list) <- c(new.reduction)
  return(output.list)
}

RPCAIntegration2 <- function(
    object = NULL,
    assay = NULL,
    layers = NULL,
    orig = NULL,
    new.reduction = 'integrated.dr',
    reference = NULL,
    features = NULL,
    normalization.method = c("LogNormalize", "SCT"),
    dims = 1:30,
    k.filter = NA,
    scale.layer = 'scale.data',
    dims.to.integrate = NULL,
    k.weight = 100,
    weight.reduction = NULL,
    sd.weight = 1,
    sample.tree = NULL,
    preserve.order = FALSE,
    verbose = TRUE,
    ...
) {
  op <- options(Seurat.object.assay.version = "v3", Seurat.object.assay.calcn = FALSE)
  on.exit(expr = options(op), add = TRUE)
  normalization.method <- match.arg(arg = normalization.method)
  features <- features %||% SelectIntegrationFeatures5(object = object)
  assay <- assay %||% 'RNA'
  layers <- layers %||% Layers(object = object, search = 'data')
  #check that there enough cells present
  ncells <- sapply(X = layers, FUN = function(x) {ncell <-  dim(object[x])[2]
  return(ncell) })
  if (min(ncells) < max(dims))  {
    abort(message = "At least one layer has fewer cells than dimensions specified, please lower 'dims' accordingly.")
  }
  if (normalization.method == 'SCT') {
    #create grouping variables
    groups <- CreateIntegrationGroups(object, layers = layers, scale.layer = scale.layer)
    object.sct <- CreateSeuratObject(counts = object, assay = 'SCT')
    object.sct$split <- groups[,1]
    object.list <- SplitObject(object = object.sct, split.by = 'split')
    object.list <- PrepSCTIntegration(object.list = object.list, anchor.features = features)
    object.list <- lapply(X = object.list, FUN = function(x) {
      x <- RunPCA(object = x, features = features, verbose = FALSE, npcs = max(dims))
      return(x)
    }
    )
  } else {
    object.list <- list()
    for (i in seq_along(along.with = layers)) {
      object.list[[i]] <- suppressMessages(suppressWarnings(
        CreateSeuratObject(counts = NULL, data = object[layers[i]][features,])
      ))
      VariableFeatures(object =  object.list[[i]]) <- features
      object.list[[i]] <- suppressWarnings(ScaleData(object = object.list[[i]], verbose = FALSE))
      object.list[[i]] <- RunPCA(object = object.list[[i]], verbose = FALSE, npcs=max(dims))
      suppressWarnings(object.list[[i]][['RNA']]$counts <- NULL)
    }
  }
  anchor <- FindIntegrationAnchors(object.list = object.list,
                                   anchor.features = features,
                                   scale = FALSE,
                                   reduction = 'rpca',
                                   normalization.method = normalization.method,
                                   dims = dims,
                                   k.filter = k.filter,
                                   reference = reference,
                                   verbose = verbose      # Removed ... here
  )
  slot(object = anchor, name = "object.list") <- lapply(
    X = slot(
      object = anchor,
      name = "object.list"),
    FUN = function(x) {
      suppressWarnings(expr = x <- DietSeurat(x, features = features[1:2]))
      return(x)
    })
  object_merged <- IntegrateEmbeddings(anchorset = anchor,
                                       reductions = orig,
                                       new.reduction.name = new.reduction,
                                       dims.to.integrate = dims.to.integrate,
                                       k.weight = k.weight,
                                       weight.reduction = weight.reduction,
                                       sd.weight = sd.weight,
                                       sample.tree = sample.tree,
                                       preserve.order = preserve.order,
                                       verbose = verbose
  )

  output.list <- list(object_merged[[new.reduction]])
  names(output.list) <- c(new.reduction)
  return(output.list)
}

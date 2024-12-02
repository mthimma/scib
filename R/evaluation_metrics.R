metrics_ASW <- function(grp, embed){

  d <- fields::rdist( embed )
  grp <- grp %>% as.factor %>% as.numeric()
  sil <- cluster::silhouette(x = grp, dmatrix = d)
  asw <- mean(sil[ , "sil_width"])

  return(asw)
}

metrics_LISI <- function(grp, embed, perplexity = 40) {

  meta <- data.frame(group = grp)
  rownames(meta) <- rownames(embed)

  lisi_out <- lisi::compute_lisi(embed, meta, "group", perplexity = perplexity )

  x <- lisi_out[ , "group"]
  lisi_mean   <- ( mean(x)  - min(x) ) / ( max(x) - min(x) )
  lisi_median <- ( median(x) - min(x) ) / ( max(x) - min(x) )

  return(list(lisi_mean = lisi_mean, lisi_median = lisi_median))
}

metrics_kBET <- function(batch, groundtruth, embed) {

  ## Find cell types that exist in all datasets/batches
  ctls   <- tapply(groundtruth, batch, unique)
  common <- Reduce(intersect, ctls)

  ## subset common celltype data
  keep  <- which(groundtruth %in% common)
  embed <- embed[ keep, ]
  batch <- batch[ keep ]

  #### iteration
  #pct <- c(seq(from = 5, to = 25, by = 5)) ## for dataset each with less than 1000 cells
  pct <- c(seq(from = 1, to = 5, by = 1)) ## for larger datasets
  kns <- floor(pct*nrow(embed)/100)

  acceptance_rate <- sapply(kns, function(kn){
    batch.estimate <- kBET::kBET(embed, batch, k0=kn, plot=F, do.pca=F, n_repeat = 10)
    return(1 - median(batch.estimate$stats$kBET.observed))
  })

  names(acceptance_rate) <- kns
  return(median(acceptance_rate))
}

run_eval_metrics <- function(object, reduction, groundtruth.cn, predicted.cn){

  embed <- Embeddings(object, reduction = reduction)
  groundtruth <- object@meta.data[ , groundtruth.cn]
  predicted   <- object@meta.data[ , predicted.cn]
  batch       <- object@meta.data[ , "batch"]

  res <- c(
    n_celltype           = length(unique(groundtruth)),
    ari_celltype         = mclust::adjustedRandIndex(groundtruth, predicted),
    nmi_celltype         = aricode::NMI(groundtruth, predicted),
    asw_celltype         = metrics_ASW (groundtruth, embed),
    lisi_mean_celltype   = metrics_LISI(groundtruth, embed)$lisi_mean,
    lisi_median_celltype = metrics_LISI(groundtruth, embed)$lisi_median,

    n_batch              = length(unique(batch)),
    ari_batch            = mclust::adjustedRandIndex(batch, predicted),
    nmi_batch            = aricode::NMI(batch, predicted),
    asw_batch            = metrics_ASW (batch, embed),
    lisi_mean_batch      = metrics_LISI(batch, embed)$lisi_mean,
    lisi_median_batch    = metrics_LISI(batch, embed)$lisi_median,
    kbet_batch           = metrics_kBET(batch, groundtruth, embed)
  )

  return(res)
}

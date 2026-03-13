# Calculate metrics under opt_res ----
run_eval_metrics <- function (object, reduction, groundtruth.cn, predicted.cn) 
{
  embed <- Embeddings(object, reduction = reduction)
  groundtruth <- object@meta.data[, groundtruth.cn]
  predicted <- object@meta.data[, predicted.cn]
  #batch <- object@meta.data[, "batch"]
  res <- c(n_celltype = length(unique(groundtruth)), 
           ari_celltype = mclust::adjustedRandIndex(groundtruth, predicted), 
           nmi_celltype = aricode::NMI(groundtruth, predicted), 
           lisi_mean_celltype = metrics_LISI(groundtruth, embed)$lisi_mean, 
           lisi_median_celltype = metrics_LISI(groundtruth, embed)$lisi_median) 
  #n_batch = length(unique(batch)), 
  #ari_batch = mclust::adjustedRandIndex(batch, predicted), 
  #nmi_batch = aricode::NMI(batch, predicted), 
  # asw_batch = metrics_ASW(batch, embed), 
  #lisi_mean_batch = metrics_LISI(batch, embed)$lisi_mean, 
  #lisi_median_batch = metrics_LISI(batch, embed)$lisi_median, 
  #kbet_batch = metrics_kBET(batch, groundtruth, embed))
  return(res)
}

stratified_sampling <- function(labs, k){
  
  ## Get the index for each cell type
  labs.split <- split( 1:length(labs), labs )
  
  ## Randomize each cell type's index to the k datasets
  out <- lapply(labs.split, function(index){
    
    dataset <- sample( rep(1:k, len=length(index)) )
    
    data.frame(original=index, dataset=dataset)
  })
  
  out <- do.call("rbind", out) %>% arrange(original) %>% pull(dataset)
  return(out)
}


## Function to plot PCA for first 50 PCs
## for first 50 pcs
## input is df <- data.frame(coords = 1:50, vars = cumsum( v/sum(v)))
pca_plot <- function(df, dataset_name, opt_nPCs) {
p <- ggplot(df, aes(x = coords, y = vars)) +
  geom_line(linewidth = 1.2, color = "#2c3e50") + 
  geom_point(size = 3, color = "#2c3e50") +      
  theme_classic(base_size = 15) +                 
  labs(
    title = paste0("Selection of Optimal number of PCs\n for ",dataset_name),
    x = "Principal Components",
    y = "Cumulative sum of PC variance"
  ) +
  # Horizontal line at y = 0.80
  geom_hline(yintercept = 0.80, linetype = "dashed", color = "red", linewidth = 1) +
  
  # Vertical line at x = nPCs_50
  geom_vline(xintercept = opt_nPCs, linetype = "dashed", color = "red", linewidth = 1) +
  
  geom_line(linewidth = 1.2)+
  theme(
    # Bold and larger titles
    plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black", face = "bold"),
    
    # Clean up the panel
    panel.grid.major.y = element_line(color = "gray90"), # Light horizontal grids only
    axis.line = element_line(linewidth = 0.8)
  ) 
return(p)
}


## Function to plot all and highlighted clusters
DimPlot_highlight <- function(obj, group.by, label.size = 8, base_size = 14){
  
  ## Extract the metadata, umap coordinates, rename the group.by for convenience
  mdf <- data.frame(
    obj@meta.data,
    obj@reductions$umap@cell.embeddings
  ) %>%
    mutate(grp = !!sym(group.by))
  
  ## Define color schemes ahead
  n    <- nlevels(mdf$grp)
  cmap <- setNames( scales::hue_pal()(n), levels(mdf$grp) )
  print(cmap)
  
  ## Calculate centroids for labelling
  centroids <- mdf %>%
    group_by(grp) %>%
    summarize(umap_1 = mean(umap_1),
              umap_2 = mean(umap_2))
  
  ## plot the unfiltered dataset
  plots <- list()
  
  plots[["all"]] <- ggplot(mdf, aes(x = umap_1, y = umap_2)) +
    geom_point(aes(col = grp)) +
    geom_text(data = centroids,
              aes(label = grp), size = label.size) +
    scale_color_manual(values = cmap) +
    theme_void(base_size = base_size) +
    theme(legend.position = "none") +
    coord_fixed()
  
  
  ## plot highlighting each cell type
  for(ct in levels(mdf$grp)){
    
    plots[[ct]] <- ggplot(mdf, aes(x = umap_1, y = umap_2)) +
      geom_point(col = "grey80") +
      geom_point(data = mdf %>% filter(grp %in% ct), aes(col = grp)) +
      geom_text(data = centroids %>% filter(grp %in% ct),
                aes(label = grp), col = "black", size = label.size) +
      scale_color_manual(values = cmap) +
      theme_void(base_size = base_size) +
      theme(legend.position = "none") +
      coord_fixed()
  }
  
  ## Return the two graps
  rm(mdf, n, cmap, centroids)
  return(plots)
}


### To make scatter plots for metrics given res as input
## scatter plot
saveScatterPlot <- function(res,filename) {
  g1 <- ggplot(res, aes(x=res, y=ari)) + geom_point()
  g2 <- ggplot(res, aes(x=res, y=asw)) + geom_point()

  #plot = g1 / g2 /g3/ g4 /g5/g6
  ggsave(paste0(filename,"MetricsforOptRes_Scatterplot_test.png"), width = 10, height = 10,
       plot = g1 / g2)
}

## To make line plots of clustering res vs metrics
## line plots
saveLinePlot <- function(obj, metricvar, interceptval) {
#make_lineplot <- function(obj,var,interceptval) {
  varname <- stringr::str_to_upper(metricvar)
  g <- ggplot(obj, aes(x = obj[["res"]], y = obj[[var]])) +
    geom_line(linewidth = 1.2, color = "#2c3e50") +  # Thicker line, professional color
    geom_point(size = 3, color = "#2c3e50") +      # Add points for clarity
    theme_classic(base_size = 15) +                 # High base size for readability
    labs(
      title = paste0("Selection of Optimal resolution using ", varname),
      x = "Clustering resolution",
      y = paste0(varname, " on clustering resolution")
    ) +
    geom_vline(xintercept = interceptval, linetype = "dashed", color = "red", linewidth = 1) +
    
    geom_line(linewidth = 1.2)+
    theme(
      # Bold and larger titles
      plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black", face = "bold"),
      
      # Clean up the panel
      panel.grid.major.y = element_line(color = "gray90"), # Light horizontal grids only
      axis.line = element_line(linewidth = 0.8)
    ) 
  return(g)
} 






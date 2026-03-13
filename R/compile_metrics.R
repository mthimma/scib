# SETUP -------------------------------------------------------------------

pacman::p_load(tidyverse, reshape2, ggplot2, gridExtra, writexl, 
               patchwork, ComplexHeatmap, ggnewscale, ggrepel)
rm(list=ls())

root <- "C://Manjula/scIntegration_local/datasets_benchmarked/"

setwd(root)

## Functions ---------------------------------------------------------------

process_project <- function(project_path) {
  
  project_name <- basename(project_path)
  nbatch_path <- file.path(project_path, "nbatch_3")
  
  # All *_LogNormalize.csv files
  files <- list.files(nbatch_path, pattern = "_LogNormalize\\.csv$", full.names = TRUE)
  
  # Read each file
  df_list <- lapply(files, function(f) {
    tool_name <- sub("_LogNormalize\\.csv$", "", basename(f))
    
    # Read, skipping the first line
    tmp <- read.csv(f, skip = 1, stringsAsFactors = FALSE, check.names = FALSE)
    colnames(tmp) <- c("metric", "value")
    
    # Create combined column name tool_metric
    tmp <- tmp %>%
      mutate( metric = paste0(tool_name, "_", metric) ) %>%
      select(metric, value)
      
    return(tmp)
  })
  
  # Combine all metrics for this dataset
  df <- bind_rows(df_list)
  
  # Spread into wide format (1 row)
  out <- df %>%
    pivot_wider(names_from = metric, values_from = value) %>% 
    mutate(project = project_name) %>% column_to_rownames("project")
  
  return(out)
}


make_heatmap <- function(obj, V, project_levels, offset = 0.001){
  
  tool_levels <- obj %>% 
    filter(metric == V) %>%
    mutate(value = ifelse(tool == "Unintegrated", value - offset, value)) %>% 
    group_by(tool) %>% 
    summarize(m = mean(value)) %>% 
    arrange(m) %>% 
    pull(tool)
  
  obj <- obj %>% 
    filter(metric == V) %>% 
    mutate(tool    = factor(tool, levels = tool_levels),
           Project = factor(Project, levels = project_levels))
  
  g <- ggplot(obj, aes(x = Project, y = tool, fill = value, label = round(value, 2)) ) +
    geom_tile() +
    geom_text() + 
    scale_fill_viridis_c(begin = 0.5, end = 1, direction = 1) +
    labs(x = NULL, y = NULL, fill = NULL, title = V) +
    theme_minimal() +
    theme(plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
          axis.text = element_text(size = 15))
  
  return(g)
}


# LOAD DATA  -----------------------------------------------------------


dirs <- list.dirs(root, recursive = FALSE)
dirs <- grep("HumanLiverHepatocyte|MouseEmbryonicDevelopment_BCells", dirs, invert = TRUE, value = TRUE)

res <- lapply(dirs, process_project) %>% bind_rows()
write_xlsx(res, paste0(root, "/BenchmarkedDatasetswithMetrics.xlsx"))


## Long format ----
res <- res %>% 
  rownames_to_column("Project") %>% 
  select(-contains("lisi_median"), -contains("_n_batch")) %>% 
  pivot_longer(cols = -Project) %>% 
  mutate(name = str_replace(name, "lisi_mean", "lisi"),
         name = str_replace(name, "time_minutes", "minutes_time")) %>% 
  separate(name, into = c("tool", "metric", "purpose"), sep = "_")


## Short name for each project ----
dmap  <- c(
  D1 = "GepLiver_NAFLD",
  D2 = "GepLiver",
  D3 = "MouseEmbryonicDevelopment_Cardiomyocytes",
  D4 = "MouseEmbryonicDevelopment_Intestines",
  D5 = "LungAdenoCarcinoma_NormalLung",
  D6 = "LungAdenoCarcinoma_NormalLymphNode",
  D7 = "LungAdenoCarcinoma_TumourLung",
  D8 = "TS_Blood"
)

res$Project <- names(dmap)[ match(res$Project, dmap) ]

rm(dirs, process_project, dmap)



# CELLTYPE -------------------------------------------------------

res_celltype <- res %>% 
  filter(purpose == "celltype") %>% 
  mutate(value  = ifelse(metric == "lisi", 1 - value, value),
         metric = ifelse(metric == "lisi", "1 - lisi", metric)) %>% 
  select(-purpose)



## Unintegrated ranking ----------------------------------------------------

res_celltype_unint <- res_celltype %>% 
  filter(tool=="Unintegrated") %>% 
  select(-tool) %>% 
  pivot_wider(names_from = metric, values_from = value) %>% 
  mutate(average = rowMeans(across(-Project), na.rm = TRUE)) %>%
  arrange(average)

res_celltype_unint %>%
  write_csv(paste0(root,"CellTypeMetrics_unintegrated.csv"))

PL <- res_celltype_unint$Project


## Plot celltype metrics ----
plots <-  lapply( c("ari", "asw", "nmi", "1 - lisi"),
                  function(VV) make_heatmap( obj = res_celltype, VV, PL) )

wrap_plots(plots)


rm(res_celltype_unint, PL, plots)



# BATCH -------------------------------------------------------------------

to_swap <- c("ari", "asw", "nmi")

res_batch <- res %>% 
  filter(purpose == "batch") %>% 
  mutate(value  = ifelse(metric %in% to_swap, 1 - value, value),
         metric = ifelse(metric %in% to_swap, paste0("1 - ", metric), metric),
         value  = ifelse(value > 1, 1, value)) %>% 
  select(-purpose)


## Unintegrated ranking ----------------------------------------------------

res_batch_unint <- res_batch %>% 
  filter(tool=="Unintegrated") %>% 
  select(-tool) %>% 
  pivot_wider(names_from = metric, values_from = value) %>% 
  mutate(average = rowMeans(across(-Project), na.rm = TRUE)) %>%
  arrange(average)

res_batch_unint %>%
  write_csv(paste0(root,"batchMetrics_unintegrated.csv"))

PL <- res_batch_unint$Project


## Plot batch metrics ----
plots <-  lapply( c("1 - ari", "1 - asw", "1 - nmi", "lisi", "kbet"),
                  function(VV) make_heatmap( obj = res_batch, VV, PL) )

g <- wrap_plots(plots) + plot_layout(ncol = 2)
ggsave(g, height = 12, width = 12, file = "lala.png")

rm(res_batch_unint, PL, plots)


# TIME -------------------------------------------------------

res_time <- res %>% 
  filter(purpose == "time") %>% 
  select(-purpose) %>% 
  mutate(value = round(value, 1))


## Unintegrated ranking ----------------------------------------------------

res_time_unint <- res_time %>% 
  filter(tool=="Unintegrated") %>% 
  select(-tool) %>% 
  pivot_wider(names_from = metric, values_from = value) %>% 
  mutate(average = rowMeans(across(-Project), na.rm = TRUE)) %>%
  arrange(desc(average))

res_time_unint %>%
  write_csv(paste0(root,"timeMetrics_unintegrated.csv"))

PL <- res_time_unint$Project


## Plot time metrics ----

##  Note:Need to reverse the y-axis and scale as high values are bad

plots <-  lapply( c("minutes"),
                  function(VV) make_heatmap( obj = res_time, VV, PL) )

wrap_plots(plots) +
  scale_y_discrete(limits = rev) + 
  scale_fill_viridis_c(begin = 0.5, end = 1, direction = -1)

rm(res_time_unint, PL, plots)



# AGGREGATE ACROSS DATASETS----------------------------------------------------------

## Celltype ----

res_celltype <- res_celltype %>% 
  filter(tool != "Unintegrated") %>% 
  group_by(tool, metric) %>% 
  summarize(avg = mean(value)) %>% 
  ungroup()

overall <- res_celltype %>% 
  group_by(tool) %>% 
  summarize(metric = "overall celltype",
            avg = mean(avg))
  
res_celltype <- bind_rows(res_celltype, overall) %>%
  mutate(type = "celltype")

rm(overall, to_swap, g)


## batch ----

res_batch <- res_batch %>% 
  filter(tool != "Unintegrated") %>% 
  group_by(tool, metric) %>% 
  summarize(avg = mean(value)) %>% 
  ungroup()

overall <- res_batch %>% 
  group_by(tool) %>% 
  summarize(metric = "overall batch",
            avg = mean(avg))

res_batch <- bind_rows(res_batch, overall) %>%
  mutate(type = "batch")

rm(overall)


## time ----

res_time <- res_time %>% 
  filter(tool != "Unintegrated") %>% 
  group_by(tool) %>% 
  summarize(avg = mean(value)) %>% 
  ungroup() %>%
  mutate(metric = "minutes", type = "time")


# COMBINED ----------------------------------------------------------------

res <- bind_rows(res_celltype, res_batch, res_time) %>% 
  mutate(type = factor(type, levels = c("celltype", "batch", "time"))) %>% 
  arrange(tool, type, metric)


rm(res_celltype, res_batch, res_time)


## Scatterplot ----

res_wide <- res %>% 
  select(-type) %>% 
  pivot_wider(id_cols = tool, names_from = metric, values_from = avg)

g_ari <- ggplot(res_wide, aes(x = `ari`, y = `1 - ari`, color = tool)) + 
  geom_point() + 
  geom_text_repel(aes(label = tool), size=5, force = 20, max.overlaps = 100) +
  theme_classic() + 
  theme(legend.position = 'none') + 
  labs(y = "Batch Correction", x = "Bioconservation", title = "ARI")

g_asw <- ggplot(res_wide, aes(x = `asw`, y = `1 - asw`, color = tool)) + 
  geom_point() + 
  geom_text_repel(aes(label = tool), size=5, force = 20, max.overlaps = 100) +
  theme_classic() + 
  theme(legend.position = 'none') + 
  labs(y = "Batch Correction", x = "Bioconservation", title = "ASW")

g_nmi <- ggplot(res_wide, aes(x = `nmi`, y = `1 - nmi`, color = tool)) + 
  geom_point() + 
  geom_text_repel(aes(label = tool), size=5, force = 20, max.overlaps = 100) +
  theme_classic() + 
  theme(legend.position = 'none') + 
  labs(y = "Batch Correction", x = "Bioconservation", title = "NMI")

g_lisi <- ggplot(res_wide, aes(x = `1 - lisi`, y = `lisi`, color = tool)) + 
  geom_point() + 
  geom_text_repel(aes(label = tool), size=5, force = 20, max.overlaps = 100) +
  theme_classic() + 
  theme(legend.position = 'none') + 
  labs(y = "Batch Correction", x = "Bioconservation", title = "lisi")

g_overall <- ggplot(res_wide, aes(x = `overall celltype`, y = `overall batch`, color = tool)) + 
  geom_point() + 
  geom_text_repel(aes(label = tool), size=5, force = 20) +
  theme_classic() + 
  theme(legend.position = 'none') + 
  labs(y = "Batch Correction", x = "Bioconservation", title = "Overall")

ggsave("scatterplot_metrics.png", height = 12, width = 13,
       plot = wrap_plots( list(g_ari, g_asw, g_nmi, g_lisi) ) )

g_overall

rm(g_ari, g_asw, g_lisi, g_nmi, g_overall)


## Dotplot ----

tool_order <- res %>% 
  filter(metric %in% c("overall celltype", "overall batch")) %>% 
  group_by(tool) %>% 
  mutate(overall_avg = sum(avg)) %>% 
  arrange(desc(overall_avg)) %>% 
  select(tool) %>% 
  unique() %>% 
  pull()

## Cell type metrics dotplot ----

metric_order <- c("ari", "asw", "nmi", "1 - lisi", "overall celltype")

res_celltype <- res %>% 
  filter(type == "celltype")

ggplot(res_celltype, aes(x = factor(metric, levels = metric_order ), 
                         y = factor(tool, levels = rev(tool_order)), 
                         fill = avg, size = avg)) +
  geom_point(shape = 21, color = "black") +
  scale_fill_gradientn(colors = c("#fff5f0", "#fb6a4a", "#a50f15"), name = "metric Avg") +
  labs(title = "Cell type Metrics", x = "", y = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 

## Batch metrics dotplot ----

metric_order <- c("1 - ari", "1 - asw", "1 - nmi", "lisi", "kbet", "overall batch")

res_batch <- res %>% 
  filter(type == "batch")
#, metric != "overall batch")

ggplot(res_batch, aes(x = factor(metric, levels = metric_order ), 
                         y = factor(tool, levels = rev(tool_order)), 
                         fill = avg, size = avg)) +
  geom_point(shape = 21, color = "black") +
  scale_fill_gradientn(colors = c("#fcfbfd", "#807dba", "#3f007d"), name = "metric Avg") +
  labs(title = "Batch Metrics", x = "", y = "") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

## Execution Time Plot ----

ggplot(res %>% filter(type == "time"), aes(x = metric, y = factor(tool, levels = rev(tool_order)), 
                      fill = avg, size = avg)) +
  geom_point(shape = 21, color = "black") +
  scale_fill_gradientn(colors = c("#ffeda0", "#feb24c", "#f03b20"), name = "Time Avg") +
  labs(title = "Time Metrics", x = "", y = "") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))



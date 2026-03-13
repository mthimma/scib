rm(list=ls())
library(tidyverse)
setwd("/data-A/scIntegration/scripts/")

##########################
## Generate the scripts ##
##########################

#today <- gsub("-", "", as.Date( Sys.Date() ))

#datasets <- c("GepLiver/GepLiver_NAFLD.rda")
datasets <- c("TS_Blood/blood50k.rda")
#datasets <-  c("MouseEmbryonicDevelopment/MED_Cardiomyocytes.rda")
#datasets <- c("MouseEmbryonicDevelopment/MED_BCells.rda")
#datasets <- c("/HumanLiverHepatocyte/HumanLiverHepatocytes.rda")
#datasets <- c("/LungAdenoCarcinoma/LAD_Norm:walLung.rda")
#datasets <- c("/LungAdenoCarcinoma/LAD_NormalLymphNode.rda")
#datasets <- c("/LungAdenoCarcinoma/LAD_tumourLung.rda")

parameters <- expand.grid(
  itool       = c("BBKNN", "CCA", "CONOS", "FastMNN", "Harmony", "LIGER", "RPCA", "SCANORAMA", "SCMC", "scVI", "Unintegrated"),
  data_file   = paste0("/data-A/scIntegration/datasets/", datasets),
  batch_name  = paste0("nbatch_", c(3,5,10)),
  norm_method = c("LogNormalize", "SCT")
  ) %>%
  mutate_all(as.character) %>%
  mutate( id = gsub(".rda", "", basename(data_file)) )

out_dir <- file.path(
  gsub("datasets/", "results/", dirname(parameters$data_file)), ## was 
  #paste0(gsub("datasets/", "results/", dirname(parameters$data_file)),"_TumourLung"),
  parameters$batch_name)

sapply( unique(out_dir), function(x) system( paste0("mkdir -p ", x) ) )

parameters$out_fn <- paste0(out_dir, "/", parameters$itool, "_", parameters$norm_method, "_testing.csv")
rm(out_dir)

## generate the commands

cmd <- with(parameters,
            paste("R --no-save < /data-A/scIntegration/scripts/run_integration.R",
                  data_file, batch_name, norm_method, itool, out_fn,
                   " > ",  gsub(".csv", ".log", out_fn)))
cmd

## write out the commands
cat(cmd,  file="run_TSBlood50k_twonormethods.sh", sep="\n")


#! /usr/bin/Rscript

################################################################################
#     ADD MIXING METRICS FOR BOTH LEVEL 03 OBJECTS                             #
################################################################################

## load libraries ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(Seurat)
  library(scDblFinder)
  library(ggpubr)
  library(clustree)
  library(patchwork)
  library(scater)
  library(scran)
  library(qs2)
})

## read in data ----
list_objects <- list()
list_objects$yap1 <- qs_read(here('output/data/Datsets/paper/D10025_yap1_dataset_Level_03_annotated.qs2'))
list_objects$meox1 <- qs_read(here('output/data/Datsets/paper/D10051_meox1_dataset_Level_03_annotated.qs2'))


## calculate mixing metric per cluster ----
list_results <- lapply(list_objects, function(seurat){
  seurat <- SetIdent(seurat, value = 'L3_celltype')
  mixing_metrics <- lapply(seurat$L3_celltype %>% unique(), function(x){
    filtered.seurat <- subset(seurat, idents = x);
    if (length(table(filtered.seurat$Genotype)) == 2){
      if (ncol(filtered.seurat) < 300){max_ks =  ncol(filtered.seurat)} else {max_ks = 300};
      MixingMetric(filtered.seurat, grouping.var = 'Genotype', max.k = max_ks)
    }
  })
  names(mixing_metrics) <- seurat$L3_celltype %>% unique()
  return(mixing_metrics)
})

## calculate mixing metric per cluster ----
filename <- here('output/analysis/D10025_yap1_dataset/MixingMetric/D10025_yap1_dataset_Level_03_mixing_metric.qs2')
qs_save(object = list_results$yap1, file = filename)

filename <- here('output/analysis/D10051_meox1_dataset/MixingMetric/D10051_meox1_dataset_Level_03_mixing_metric.qs2')
qs_save(object = list_results$meox1, file = filename)





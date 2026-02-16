#! /usr/bin/Rscript

################################################################################
#     ADD GENE SCORES FOR BOTH LEVEL 03 OBJECTS                                #
################################################################################

# will add:
# MAPK
# HIPPO
# AP-1 SCORE
# LEC SCORE
# VEC SCORE
# TRANSITION SCORE

## read in gene sets ----
source(here("code/general_scripts/yap1_meox1_important_genes.R"))

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


## add hippo score ----
list_objects <- lapply(list_objects, function(seurat){
  seurat <- AddModuleScore(object = seurat,
                           features = list(hippo_targets),
                           assay = 'RNA',
                           seed = 100,
                           name = 'hippo_', #function will add a number after, underscore keeps it cleaner
                           search = TRUE)
  return(seurat)
})

## add MAPK score ----
list_objects <- lapply(list_objects, function(seurat){
  seurat <- AddModuleScore(object = seurat,
                           features = list(mapk_targets$MAPK_genes),
                           assay = 'RNA',
                           seed = 100,
                           name = 'mapk_', #function will add a number after, underscore keeps it cleaner
                           search = TRUE)
  return(seurat)
})

## add wnt score ----
list_objects <- lapply(list_objects, function(seurat){
  seurat <- AddModuleScore(object = seurat,
                           features = list(wnt_targets$Wnt_genes),
                           assay = 'RNA',
                           seed = 100,
                           name = 'wnt_', #function will add a number after, underscore keeps it cleaner
                           search = TRUE)
  return(seurat)
})


## add AP1 score ----
list_objects <- lapply(list_objects, function(seurat){
  seurat <- AddModuleScore(object = seurat,
                           features = list(ap1_targets),
                           assay = 'RNA',
                           seed = 100,
                           name = 'ap1_', #function will add a number after, underscore keeps it cleaner
                           search = TRUE)
  return(seurat)
})

## add LEC score ----
list_objects <- lapply(list_objects, function(seurat){
  seurat <- AddModuleScore(object = seurat,
                           features = list(lec_score),
                           assay = 'RNA',
                           seed = 100,
                           name = 'EMBO_LEC_score_', #function will add a number after, underscore keeps it cleaner
                           search = TRUE)
  return(seurat)
})

## add LEC dev score ----
list_objects <- lapply(list_objects, function(seurat){
  seurat <- AddModuleScore(object = seurat,
                           features = list(lec_dev_score),
                           assay = 'RNA',
                           seed = 100,
                           name = 'LEC_dev_score_', #function will add a number after, underscore keeps it cleaner
                           search = TRUE)
  return(seurat)
})


## add VEC score ----
list_objects <- lapply(list_objects, function(seurat){
  seurat <- AddModuleScore(object = seurat,
                           features = list(vec_score),
                           assay = 'RNA',
                           seed = 100,
                           name = 'EMBO_VEC_score_', #function will add a number after, underscore keeps it cleaner
                           search = TRUE)
  return(seurat)
})

## add transition score ----
list_objects <- lapply(list_objects, function(seurat){
  seurat$EMBO_transient_score <- seurat$EMBO_LEC_score_1 - seurat$EMBO_VEC_score_1
  return(seurat)
})

## add transition score ----
list_objects <- lapply(list_objects, function(seurat){
  seurat$EMBO_dev_transient_score <- seurat$LEC_dev_score_1 - seurat$EMBO_VEC_score_1
  return(seurat)
})


## export ----
#overwrite as we're only adding scores
filename <- here('output/data/Datsets/paper/D10025_yap1_dataset_Level_03_annotated.qs2')
qs_save(object = list_objects$yap1, file = filename)

filename <- here('output/data/Datsets/paper/D10051_meox1_dataset_Level_03_annotated.qs2')
qs_save(object = list_objects$meox1, file = filename)


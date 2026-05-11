#! /usr/bin/Rscript

################################################################################
#     MEOX1 DATASET: SPLIT LEVEL 03 INTO SAMPLES                               #
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
level_03 <- qs_read(here('output/data/Datsets/paper/D10051_meox1_dataset_Level_03_annotated.qs2'))

## write function ----
run_split_samples <- function(seurat=level_03, sample){
  seurat <- SetIdent(seurat, value = "Genotype")
  seurat <- subset(seurat, idents = sample)
  seurat <- NormalizeData(seurat) %>%   #this is technically unneccessary, but also won't hurt
    FindVariableFeatures() %>%
    ScaleData(., vars.to.regress = c("percent.mt")) %>%
    RunPCA() %>%
    RunUMAP(., dims = 1:30) %>% #defaults
    FindNeighbors(., dims = 1:30) %>%
    FindClusters(., resolution = seq(0.1, 1, 0.1)) #this will overwrite old clustering (that's okay)
  filename <- sprintf('%s/D10051_meox1_dataset_Level_04_%s_annotated.qs2', here('output/data/Datsets/paper'), sample)
  qs_save(object = seurat, file = filename)
  return(seurat)

}

## run for both samples ----
both_samples <- level_03$Genotype %>% unique()

list_out <- lapply(both_samples, function(geno) run_split_samples(sample = geno))




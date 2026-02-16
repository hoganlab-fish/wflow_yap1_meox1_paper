#! /usr/bin/Rscript

#################################################
# Run Markers for selected metadata
#################################################

# Load libraries ----
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
  library(scGate)
})

# Load data ----
level_03 <- qs_read(here('output/data/Datsets/paper/D10051_meox1_dataset_Level_03_annotated.qs2'))

# Export path ----
save_dir <- here("output/analysis/D10051_meox1_dataset/Markers")


# Write function: get markers ----
# This function should be able to run with whatever metadata column
run_markers <- function(seurat,
                        group,
                        save_dir,
                        log2_t =0,
                        min_pct = 0.01){
  markers <- FindAllMarkers(object = seurat,
                            assay = "RNA",
                            group.by = group,
                            logfc.threshold = log2_t,
                            min.pct = min_pct)

  file_name <- sprintf("%s/%s_Markers_%s_fc%0.2f_minpct_%0.2f.csv", save_dir, seurat@misc$name, group, log2_t, min_pct)
  write.csv(file = file_name, x = markers, row.names = F)

  return(markers)

}


# Run for celltypes ----
output <- run_markers(seurat = level_03, group = "L3_celltype", save_dir = save_dir)




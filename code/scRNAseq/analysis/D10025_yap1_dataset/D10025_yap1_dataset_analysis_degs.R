#! /usr/bin/Rscript

#################################################
# Run DEGs between mutant and WT for all clusters
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
level_03 <- qs_read(here('output/data/Datsets/paper/D10025_yap1_dataset_Level_03_annotated.qs2'))

# Export path ----
save_dir <- here("output/analysis/D10025_yap1_dataset/DEGs")

# genes of interest ----
hippo_targets <- c(
  "yap1",
  "meox1",
  "wwtr1",
  "ccn2a",
  "ccn2b",
  "amotl2a",
  "amotl2b",
  "ccn1",
  "ccn1l2",
  "cavin1b",
  "cavin2a",
  "cavin2b",
  "bmp4")

# Write function: generic DEG ----
# This function should be able to run with whatever metadata column
run_degs_mut_wt <- function(seurat,
                            group,
                            save_dir,
                            log2_t =0,
                            min_pct = 0.01){
  #all idents for the group
  all_results <- lapply(seurat@meta.data[,group] %>% unique(), function(unit){
    seurat_subset <- subset(seurat, cells = colnames(seurat)[seurat@meta.data[,group] == unit])
    seurat_subset$Genotype <- factor(seurat_subset$Genotype, levels = c("wildtype", "yap1_mutant"))
    if (all(table(seurat_subset$Genotype) > 5)){
      Idents(seurat_subset) <- "Genotype"
      FindMarkers(object = seurat_subset,
                  group.by = "Genotype",
                  ident.1 = "yap1_mutant",
                  ident.2 = "wildtype",
                  min.pct = min_pct,
                  assay = "RNA",
                  logfc.threshold = log2_t) %>%
        rownames_to_column("Gene") %>%
        mutate(., group = unit)
    }} ) %>% bind_rows()

  file_name <- sprintf("%s/%s_DEG_mutVSwt_%s_fc%0.2f_minpct_%0.2f.csv", save_dir, seurat@misc$name, group, log2_t, min_pct)
  write.csv(file = file_name, x = all_results, row.names = F)

  return(all_results)

}

# Write function: selected features DEG ----
# This function should be able to run with whatever metadata column, but will make no log2FC etc filtering
run_degs_mut_wt_selected_genes <- function(seurat,
                            group,
                            save_dir,
                            log2_t =0,
                            min_pct = 0,
                            features_in,
                            export_name){
  #all idents for the group
  all_results <- lapply(seurat@meta.data[,group] %>% unique(), function(unit){
    seurat_subset <- subset(seurat, cells = colnames(seurat)[seurat@meta.data[,group] == unit])
    seurat_subset$Genotype <- factor(seurat_subset$Genotype, levels = c("wildtype", "yap1_mutant"))
    if (all(table(seurat_subset$Genotype) > 5)){
      Idents(seurat_subset) <- "Genotype"
      FindMarkers(object = seurat_subset,
                  group.by = "Genotype",
                  ident.1 = "yap1_mutant",
                  ident.2 = "wildtype",
                  min.pct = min_pct,
                  assay = "RNA",
                  features =features_in,
                  logfc.threshold = log2_t) %>%
        rownames_to_column("Gene") %>%
        mutate(., group = unit)
    }} ) %>% bind_rows()

  file_name <- sprintf("%s/%s_DEG_%s_mutVSwt_%s_fc%0.2f_minpct_%0.2f.csv", save_dir, seurat@misc$name, export_name, group, log2_t, min_pct)
  write.csv(file = file_name, x = all_results, row.names = F)

  return(all_results)

}


# Run for all celltypes ----
degs <- run_degs_mut_wt(seurat = level_03, group = "L3_celltype", save_dir = save_dir)

# Run for all celltypes - hippo genes only ----
degs_hippo <- run_degs_mut_wt_selected_genes(seurat = level_03,
                                       log2_t = 0,
                                       min_pct =0,
                                       features_in = hippo_targets,
                                       group = "L3_celltype",
                                       save_dir = save_dir,
                                       export_name = "hippo_targets")






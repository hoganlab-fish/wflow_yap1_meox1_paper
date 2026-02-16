
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
# input_paths <- c(here('output/data/Datsets/D10025_yap1_dataset/D10025_yap1_dataset_Level_03_kept_doublets.qs2'),
#                  here('output/data/Datsets/D10025_yap1_dataset/D10025_yap1_dataset_Level_03.qs2'))
input_paths <- c(here('output/data/Datsets/D10025_yap1_dataset/D10025_yap1_dataset_Level_03_SCT_kept_doublets.qs2'),
                 here('output/data/Datsets/D10025_yap1_dataset/D10025_yap1_dataset_Level_03_SCT.qs2'))
list_objects <- lapply(input_paths, qs_read)

names(list_objects) <- gsub("\\.qs2", "", basename(input_paths))

#update names
# list_objects$D10025_yap1_dataset_Level_03_kept_doublets@misc$name <- "D10025_yap1_dataset_Level_03_kept_doublets"
# list_objects$D10025_yap1_dataset_Level_03@misc$name <- "D10025_yap1_dataset_Level_03"
list_objects$D10025_yap1_dataset_Level_03_SCT_kept_doublets@misc$name <- "D10025_yap1_dataset_Level_03_SCT_kept_doublets"
list_objects$D10025_yap1_dataset_Level_03_SCT@misc$name <- "D10025_yap1_dataset_Level_03_SCT"

# Export path ----
save_dir <- here("output/analysis/D10025_yap1_dataset/DEGs")

# Write function ----
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
  write.csv(file = file_name, x = all_results)

  return(all_results)

}

run_degs_mut_wt_sct <- function(seurat, #don't run this for real
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
                  assay = "SCT",
                  logfc.threshold = log2_t) %>%
        rownames_to_column("Gene") %>%
        mutate(., group = unit)
    }} ) %>% bind_rows()

  return(all_results)

}


# 26.01.2025 - Run on res 0.2 ----
list_markers <- lapply(list_objects, function(seurat){
  run_degs_mut_wt(seurat = seurat, save_dir = save_dir, group = "RNA_snn_res.0.2")
} )

# 26.01.2025 - Run on res 0.3 ----
list_markers <- lapply(list_objects, function(seurat){
  run_degs_mut_wt(seurat = seurat, save_dir = save_dir, group = "RNA_snn_res.0.3")
} )

# 26.01.2025 - Run on res 0.4 (SCT) ----
list_markers <- lapply(list_objects, function(seurat){
  run_degs_mut_wt(seurat = seurat, save_dir = save_dir, group = "SCT_snn_res.0.4")
} )
#try on SCT
list_markers_sct <- lapply(list_objects, function(seurat){
  run_degs_mut_wt_sct(seurat = seurat, save_dir = save_dir, group = "SCT_snn_res.0.4")
} )


# OLD ANALYSIS - reproduce ----

old_object <- readRDS(here("data/input/old_objects/Level_03_seuratObject_LEC_VEC.RDS"))

old_object

sct_res <- lapply(old_object$Level_03_seurat_cluster_predicted_phenotype %>% unique(), function(unit){
  seurat_subset <- subset(old_object, cells = colnames(old_object)[old_object$Level_03_seurat_cluster_predicted_phenotype == unit])
  seurat_subset$Sample_Name <- factor(seurat_subset$Sample_Name, levels = c("WT", "yap1_mutant"))
  if (all(table(seurat_subset$Sample_Name) > 5)){
    Idents(seurat_subset) <- "Sample_Name"
    FindMarkers(object = seurat_subset,
                assay = "SCT",
                group.by = "Sample_Name",
                ident.1 = "yap1_mutant",
                ident.2 = "WT",
                min.pct = 0.01,
                logfc.threshold = 0) %>%
      rownames_to_column("Gene") %>%
      mutate(., group = unit)
  }} ) %>% bind_rows()


rna_res <- lapply(old_object$Level_03_seurat_cluster_predicted_phenotype %>% unique(), function(unit){
  seurat_subset <- subset(old_object, cells = colnames(old_object)[old_object$Level_03_seurat_cluster_predicted_phenotype == unit])
  seurat_subset$Sample_Name <- factor(seurat_subset$Sample_Name, levels = c("WT", "yap1_mutant"))
  if (all(table(seurat_subset$Sample_Name) > 5)){
    Idents(seurat_subset) <- "Sample_Name"
    FindMarkers(object = seurat_subset,
                assay = "RNA",
                group.by = "Sample_Name",
                ident.1 = "yap1_mutant",
                ident.2 = "WT",
                min.pct = 0.01,
                logfc.threshold = 0) %>%
      rownames_to_column("Gene") %>%
      mutate(., group = unit)
  }} ) %>% bind_rows()


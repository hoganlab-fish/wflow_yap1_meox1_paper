#! /usr/bin/Rscript

#################################################
# Run GO analysis running enrichR
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
  library(enrichR)
})

# Load data ----
markers <- read_csv(here("output/analysis/D10025_yap1_dataset/Markers/Level_03_Markers_L3_celltype_fc0.00_minpct_0.01.csv"))
degs <- read_csv(here("output/analysis/D10025_yap1_dataset/DEGs/Level_03_DEG_mutVSwt_L3_celltype_fc0.00_minpct_0.01.csv"))

# Export path ----
save_dir <- here("output/analysis/D10025_yap1_dataset/enrichR")

# Write function ----
run_enrichr <- function(markers, databases){
  setEnrichrSite("FishEnrichr") #set to fish
  #run enrichr
  enriched <- enrichr(markers, databases)
  ids <- lapply(enriched, nrow) %>% unlist() != 0 #remove the ones that had no hits at all
  enriched <- enriched[ids]
  complete_enriched_df <- purrr::map_df(enriched, ~as.data.frame(.x), .id="cat")
  go.terms <- apply(complete_enriched_df, 1,function(x) strsplit(x[2], "[(]")[[1]][2]) %>% as.data.frame()
  complete_enriched_df <- complete_enriched_df %>%
    mutate(., GOTerm = apply(go.terms, 1,function(x) strsplit(x, "[)]")[[1]][1]))
  return(complete_enriched_df)
}

# Run for: markers (top 100 upregulated only) ----
markers_results_all <- lapply(markers$cluster %>% unique(), function(celltype){
  input <- markers %>% filter(., p_val_adj < 0.05 & cluster == celltype) %>%
    top_n(., n=100, wt = -avg_log2FC) %>% pull(gene)
  result <- run_enrichr(markers = input, databases = c('GO_Biological_Process_2018', 'GO_Biological_Process_AutoRIF', 'KEGG_2019', 'WikiPathways_2018'))
  result$cluster <- celltype
  return(result)
}) %>% bind_rows()

filename <- here("output/analysis/D10025_yap1_dataset/enrichR/D10025_yap1_dataset_markers_L3_celltype_top100_enrichR.qs2")
qs_save(object = markers_results_all, file = filename)


# Run for: DEGs (top 100 upregulated only) ----
# no running IL_VEC because its only mutant
degs_up_results_all <- lapply(  c("preLEC", "mVEC",  "hmVEC", "cVEC", "LEC"), function(celltype){
  input <- degs %>% filter(., p_val_adj < 0.05 & group == celltype & avg_log2FC > 0) %>%
    top_n(., n=100, wt = -avg_log2FC) %>% pull(Gene)
  result <- run_enrichr(markers = input, databases = c('GO_Biological_Process_2018', 'GO_Biological_Process_AutoRIF', 'KEGG_2019', 'WikiPathways_2018'))
  result$cluster <- celltype
  return(result)
}) %>% bind_rows()

filename <- here("output/analysis/D10025_yap1_dataset/enrichR/D10025_yap1_dataset_DEGs_L3_celltype_mutVSwt_top100_UP_enrichR.qs2")
qs_save(object = degs_up_results_all, file = filename)

# Run for: DEGs (top 100 downregulated only) ----
degs_down_results_all <- lapply( c("preLEC", "mVEC",  "hmVEC", "cVEC", "LEC"), function(celltype){
  input <- degs %>% filter(., p_val_adj < 0.05 & group == celltype & avg_log2FC < 0) %>%
    top_n(., n=100, wt = avg_log2FC) %>% pull(Gene)
  result <- run_enrichr(markers = input, databases = c('GO_Biological_Process_2018', 'GO_Biological_Process_AutoRIF', 'KEGG_2019', 'WikiPathways_2018'))
  result$cluster <- celltype
  return(result)
}) %>% bind_rows()

filename <- here("output/analysis/D10025_yap1_dataset/enrichR/D10025_yap1_dataset_DEGs_L3_celltype_mutVSwt_top100_DOWN_enrichR.qs2")
qs_save(object = degs_down_results_all, file = filename)
# Run for: markers (top 50 upregulated only) ----
markers_results_all <- lapply(markers$cluster %>% unique(), function(celltype){
  input <- markers %>% filter(., p_val_adj < 0.05 & cluster == celltype) %>%
    top_n(., n=50, wt = -avg_log2FC) %>% pull(gene)
  result <- run_enrichr(markers = input, databases = c('GO_Biological_Process_2018', 'GO_Biological_Process_AutoRIF', 'KEGG_2019', 'WikiPathways_2018'))
  result$cluster <- celltype
  return(result)
}) %>% bind_rows()

filename <- here("output/analysis/D10025_yap1_dataset/enrichR/D10025_yap1_dataset_markers_L3_celltype_top50_enrichR.qs2")
qs_save(object = markers_results_all, file = filename)


# Run for: DEGs (top 50 upregulated only) ----
# no running IL_VEC because its only mutant
degs_up_results_all <- lapply( c("preLEC", "mVEC",  "hmVEC", "cVEC", "LEC"), function(celltype){
  input <- degs %>% filter(., p_val_adj < 0.05 & group == celltype & avg_log2FC > 0) %>%
    top_n(., n=50, wt = -avg_log2FC) %>% pull(Gene)
  result <- run_enrichr(markers = input, databases = c('GO_Biological_Process_2018', 'GO_Biological_Process_AutoRIF', 'KEGG_2019', 'WikiPathways_2018'))
  result$cluster <- celltype
  return(result)
}) %>% bind_rows()

filename <- here("output/analysis/D10025_yap1_dataset/enrichR/D10025_yap1_dataset_DEGs_L3_celltype_mutVSwt_top50_UP_enrichR.qs2")
qs_save(object = degs_up_results_all, file = filename)

# Run for: DEGs (top 50 downregulated only) ----
degs_down_results_all <- lapply( c("preLEC", "mVEC",  "hmVEC", "cVEC", "LEC"), function(celltype){
  input <- degs %>% filter(., p_val_adj < 0.05 & group == celltype & avg_log2FC < 0) %>%
    top_n(., n=50, wt = avg_log2FC) %>% pull(Gene)
  result <- run_enrichr(markers = input, databases = c('GO_Biological_Process_2018', 'GO_Biological_Process_AutoRIF', 'KEGG_2019', 'WikiPathways_2018'))
  result$cluster <- celltype
  return(result)
}) %>% bind_rows()

filename <- here("output/analysis/D10025_yap1_dataset/enrichR/D10025_yap1_dataset_DEGs_L3_celltype_mutVSwt_top50_DOWN_enrichR.qs2")
qs_save(object = degs_down_results_all, file = filename)




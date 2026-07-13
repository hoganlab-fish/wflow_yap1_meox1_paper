#! /usr/bin/Rscript

################################################################################
#     TESTING NEW GENE MODULE SCORE                                            #
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

## read in genes ----
#genes now made here: code/scRNAseq/Grimm_etal/EMBO_L2_LEC_VEC_markers.R, using "EMBO method"
genes <- read_csv(filename <- here("data/genesets/EMBO_L2_VECLEC_markers_EMBOmethod_3_5dpf_genesforscoring_only.csv"))
lec_genes <- genes %>% filter(., direction == "LEC") %>% pull(Gene)
vec_genes <- genes %>% filter(., direction == "VEC") %>% pull(Gene)

## write function ----
run_embo_method_scoring <- function(seurat){
  #LEC score
  seurat <- AddModuleScore(object = seurat,
                           features = list(lec_genes),
                           assay = 'RNA',
                           seed = 100,
                           name = 'LEC_new_score_', #function will add a number after, underscore keeps it cleaner
                           search = TRUE)

  #VEC score
  seurat <- AddModuleScore(object = seurat,
                           features = list(vec_genes),
                           assay = 'RNA',
                           seed = 100,
                           name = 'VEC_new_score_', #function will add a number after, underscore keeps it cleaner
                           search = TRUE)

  # transition
  seurat$transient_score_new <- seurat$LEC_new_score_1 - seurat$VEC_new_score_1

  return(seurat)

}

## run on objects ----

list_objects <- lapply(list_objects, run_embo_method_scoring)

## trial some plots ----
run_plots <- function(seurat){

  plot_lec <- FeaturePlot(seurat,
                      features = "LEC_new_score_1",
                      cols =  c("#d9d9d9", "darkgreen"),
                      # pt.size = point_size,
                      min.cutoff = "q1",
                      max.cutoff =  "q99",
                      order = T)

  plot_vec <- FeaturePlot(seurat,
                          features = "VEC_new_score_1",
                          cols =  c("#d9d9d9", "#406880"),
                          # pt.size = point_size,
                          min.cutoff = "q1",
                          max.cutoff =  "q99",
                          order = T)

  #transition both
  md <- seurat[[]]
  coords <- Embeddings(seurat[["umap"]])
  md <- cbind(md, coords)

  midpoint <- mean(range(md$transient_score_new))

  # Order by absolute difference from the midpoint
  md <- md[order(abs(md$transient_score_new - midpoint)), ]

  # plot
  trans_score_legend <- ggplot(md, aes(x = umap_1, y = umap_2, color = transient_score_new)) +
    geom_point(size = 1) +
    scale_color_gradient2(low = "#406880", high = "darkgreen", mid = "azure2") + #not keeping original green bc its too light
    # viridis::scale_fill_viridis(option = "E", direction = 1) +
    theme_classic()

  plot_combined <- plot_lec + plot_vec + trans_score_legend
  return(plot_combined)

}


plot_list <- lapply(list_objects, run_plots)






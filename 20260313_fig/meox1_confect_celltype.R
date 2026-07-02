# !/usr/bin/Rscript
# merge VECs and LECs into groups for comparison
# 
# context:
#   confects is complementary to pvalues
#     deg asks: what changed?
#     confects asks: what changed substantially?
# important notes:
#   confects is not a replacement for DGE
#   confects runs on top of DGE, and inherits its assumptions/power etc
# steps:
#   here we pseudobulked and ran limma, then topconfects
# output:
# ../Saki_data/Revision_analysis/queue_plots/
#   Level_03_TopConfects_MergedGroups_mutVSwt_FDR0.05.csv
#   Level_03_TopConfects_L3_celltype_merged_VEC_merged_LEC_hippo_targets_MergedGroups_mutVSwt_FDR0.05.csv

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
library(limma)
library(edgeR)
library(topconfects)

# Write function: topconfects for custom groups ----
run_topconfects_custom_groups <- function(seurat,
                                          meta_column,      
                                          custom_groups,    
                                          save_dir,
                                          geno_col = "Genotype",
                                          mut_id = "meox1_mutant",
                                          wt_id = "wildtype",
                                          project_name = "meox1_dataset",
                                          fdr_target = 0.05) {  # topconfects uses FDR instead of p_val
  
  p_name <- if (!is.null(seurat@misc$name)) seurat@misc$name else project_name
  
  all_results <- lapply(names(custom_groups), function(new_group_name) {
    
    target_identities <- custom_groups[[new_group_name]]
    cells_to_keep <- colnames(seurat)[seurat@meta.data[, meta_column] %in% target_identities]
    
    if (length(cells_to_keep) == 0) return(NULL)
    
    seurat_subset <- subset(seurat, cells = cells_to_keep)
    
    # Factor the genotype (WT must be the first level to act as the reference/intercept)
    seurat_subset[[geno_col]] <- factor(seurat_subset@meta.data[, geno_col], levels = c(wt_id, mut_id))
    
    # Ensure sufficient cells for a valid model
    if (all(table(seurat_subset[[geno_col]]) > 5)) {
      cat(sprintf("Running limma + topconfects on group: %s...\n", new_group_name))
      
      # 1. Extract raw counts
      # Note: If using Seurat v5, you may need layer = "counts". For v4, slot = "counts" is standard.
      counts_matrix <- GetAssayData(seurat_subset, assay = "RNA", layer = "counts")
      
      # 2. Build the design matrix for limma
      design <- model.matrix(~ seurat_subset@meta.data[[geno_col]])
      # The coefficient for the mutant will be the 2nd column
      target_coef <- 2 
      
      # 3. Run limma-voom pipeline (Treating cells as replicates)
      dge <- DGEList(counts = counts_matrix)
      dge <- calcNormFactors(dge)
      v <- voom(dge, design, plot = FALSE)
      fit <- lmFit(v, design)
      fit <- eBayes(fit)
      
      # 4. Run topconfects on the limma fit
      confects_res <- limma_confects(fit, coef = target_coef, fdr = fdr_target)
      
      # 5. Extract the results table and tag it with our custom group name
      # topconfects returns a specialized list; the dataframe we want is in $table
      res_df <- confects_res$table %>%
        mutate(group = new_group_name)
      
      return(res_df)
    }
  }) %>% bind_rows()

  # Export the compiled table
  file_name <- sprintf("%s/%s_TopConfects_MergedGroups_mutVSwt_FDR%0.2f.csv", save_dir, p_name, fdr_target)
  write.csv(file = file_name, x = all_results, row.names = FALSE)

  return(all_results)
}

run_topconfects_selected_genes_custom_groups <- function(seurat,
                                                         meta_column,
                                                         custom_groups,
                                                         features_in,
                                                         export_name,
                                                         save_dir,
                                                         geno_col = "Genotype",
                                                         mut_id = "meox1_mutant",
                                                         wt_id = "wildtype",
                                                         project_name = "meox1_dataset",
                                                         fdr_target = 0.05) {
  
  p_name <- if (!is.null(seurat@misc$name)) seurat@misc$name else project_name
  
  all_results <- lapply(names(custom_groups), function(new_group_name) {
    
    target_identities <- custom_groups[[new_group_name]]
    cells_to_keep <- colnames(seurat)[seurat@meta.data[, meta_column] %in% target_identities]
    if (length(cells_to_keep) == 0) return(NULL)
    
    seurat_subset <- subset(seurat, cells = cells_to_keep)
    seurat_subset[[geno_col]] <- factor(seurat_subset@meta.data[, geno_col], levels = c(wt_id, mut_id))
    
    if (all(table(seurat_subset[[geno_col]]) > 5)) {
      cat(sprintf("Running limma + topconfects (Selected Genes) on: %s...\n", new_group_name))
      
      counts_matrix <- GetAssayData(seurat_subset, assay = "RNA", layer = "counts")
      valid_features <- features_in[features_in %in% rownames(counts_matrix)]
      counts_matrix <- counts_matrix[valid_features, ]

      keep_cells <- colSums(counts_matrix) > 0
      counts_matrix <- counts_matrix[, keep_cells]
      seurat_subset <- seurat_subset[, keep_cells]

      design <- model.matrix(~ seurat_subset@meta.data[[geno_col]])
      dge <- DGEList(counts = counts_matrix)
      dge <- calcNormFactors(dge)
      v <- voom(dge, design, plot = FALSE)
      fit <- lmFit(v, design)
      fit <- eBayes(fit)
      
      confects_res <- limma_confects(fit, coef = 2, fdr = fdr_target)
      
      return(confects_res$table %>% mutate(group = new_group_name))
    }
  }) %>% bind_rows()

  file_name <- sprintf("%s/%s_TopConfects_%s_MergedGroups_mutVSwt_FDR%0.2f.csv", 
                       save_dir, p_name, export_name, fdr_target)
  write.csv(file = file_name, x = all_results, row.names = FALSE)

  return(all_results)
}

infile_path <- "../paper/D10051_meox1_dataset_Level_03_annotated.qs2"
data <- qs_read(infile_path)
outfile_dir <- "../Revision_analysis/queue_plots/"

custom_groups <- list(
    "hmVEC__mVEC" = c("hmVEC", "mVEC"),
    "LEC__pre_muLEC" = c("LEC", "pre_muLEC")
)

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
    "bmp4"
)

confect_results <- run_topconfects_custom_groups(
  seurat = data,
  meta_column = "L3_celltype",
  custom_groups = custom_groups,
  save_dir = outfile_dir,
  fdr_target = 0.05
)

confect_hippo <- run_topconfects_selected_genes_custom_groups(
  seurat = data,
  meta_column = "L3_celltype",
  custom_groups = custom_groups,
  features_in = hippo_targets,
  export_name = "L3_celltype_merged_VEC_merged_LEC_hippo_targets",
  save_dir = outfile_dir
)

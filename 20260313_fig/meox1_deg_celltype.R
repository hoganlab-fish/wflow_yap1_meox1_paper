# !/usr/bin/Rscript
# merge VECs and LECs into groups for comparison
# ref: D10051_meox1_dataset_analysis_degs.R
# ../Saki_data/Revision_analysis/queue_plots/
#   Level_03_DEG_mutVSwt_L3_celltype_merged_VEC_merged_LEC_fc0.00_minpct_0.01.csv
#   Level_03_DEG_mutVSwt_L3_celltype_merged_VEC_merged_LEC_hippo_targets_fc0.00_minpct_0.00.csv
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

run_degs_mut_wt_custom_groups <- function(seurat,
                            group,
                            save_dir,
                            custom_groups = NULL,
                            log2_t = 0,
                            min_pct = 0.01) {
  # use custom groups if passed
  if (is.null(custom_groups)) {
    unique_ids <- unique(seurat@meta.data[, group])
    eval_list <- setNames(as.list(unique_ids), unique_ids)
  } else {
    eval_list <- custom_groups
  }

  # parse named blocks
  all_results <- lapply(names(eval_list), function(block_name) {
    target_identities <- eval_list[[block_name]]
    
    # merge multiple idents
    seurat_subset <- subset(seurat, cells = colnames(seurat)[seurat@meta.data[, group] %in% target_identities])
    
    # if no match, skip
    if (ncol(seurat_subset) == 0) return(NULL)
    
    seurat_subset$Genotype <- factor(seurat_subset$Genotype, levels = c("wildtype", "meox1_mutant"))
    
    if (all(table(seurat_subset$Genotype) > 5)) {
      Idents(seurat_subset) <- "Genotype"
      FindMarkers(object = seurat_subset,
                  group.by = "Genotype",
                  ident.1 = "meox1_mutant",
                  ident.2 = "wildtype",
                  min.pct = min_pct,
                  assay = "RNA",
                  logfc.threshold = log2_t) %>%
        rownames_to_column("Gene") %>%
        mutate(group = block_name) # Labels output row with your custom block name
    }
  }) %>% bind_rows()

  # tag filename if custom grouping was deployed
  file_suffix <- ifelse(is.null(custom_groups), group, paste0(group, "_merged_VEC_merged_LEC"))
  file_name <- sprintf("%s/%s_DEG_mutVSwt_%s_fc%0.2f_minpct_%0.2f.csv", save_dir, seurat@misc$name, file_suffix, log2_t, min_pct)
  write.csv(file = file_name, x = all_results, row.names = FALSE)

  return(all_results)
}

run_degs_mut_wt_selected_genes_custom_groups <- function(seurat,
                                                         group,   
                                                         custom_groups, 
                                                         save_dir,
                                                         features_in,   
                                                         export_name,   
                                                         geno_col = "Genotype",
                                                         mut_id = "meox1_mutant",
                                                         wt_id = "wildtype",
                                                         project_name = "meox1_dataset",
                                                         log2_t = 0,
                                                         min_pct = 0) {
                                                         
  p_name <- if (!is.null(seurat@misc$name)) seurat@misc$name else project_name
  
  # Iterate over the NAMES of your custom merged list
  all_results <- lapply(names(custom_groups), function(new_group_name) {
    
    # Extract the vector of original celltypes to merge
    target_identities <- custom_groups[[new_group_name]]
    
    # Grab cells matching ANY of the target identities
    cells_to_keep <- colnames(seurat)[seurat@meta.data[, group] %in% target_identities]
    
    # Skip if none of these cells exist
    if (length(cells_to_keep) == 0) return(NULL)
    
    seurat_subset <- subset(seurat, cells = cells_to_keep)
    
    # Factor the genotype column dynamically
    seurat_subset[[geno_col]] <- factor(
        seurat_subset@meta.data[, geno_col], levels = c(wt_id, mut_id)
        )
    
    # Only run if we have enough cells per genotype in this merged block
    if (all(table(seurat_subset[[geno_col]]) > 5)) {
      Idents(seurat_subset) <- geno_col
      
      FindMarkers(object = seurat_subset,
                  ident.1 = mut_id,
                  ident.2 = wt_id,
                  min.pct = min_pct,
                  assay = "RNA",
                  features = features_in,         # Restrict test to your selected genes
                  logfc.threshold = log2_t) %>%
        rownames_to_column("Gene") %>%
        mutate(group = new_group_name)        # Tag with your custom merged name
    }
  }) %>% bind_rows()

  # Create a distinct filename using the export_name (e.g., "hippo_targets")
  file_name <- sprintf("%s/%s_DEG_mutVSwt_%s_fc%0.2f_minpct_%0.2f.csv", 
                       save_dir, seurat@misc$name, export_name, log2_t, min_pct)
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

# all celltypes
degs <- run_degs_mut_wt_custom_groups(
    seurat = data, 
    group = "L3_celltype", 
    custom_groups = custom_groups, 
    save_dir = outfile_dir
    )

# all celltypes - hippo genes only
degs_hippo <- run_degs_mut_wt_selected_genes_custom_groups(
    seurat = data,
    log2_t = 0,
    min_pct= 0,
    features_in = hippo_targets,
    custom_groups = custom_groups,
    group = "L3_celltype",
    save_dir = outfile_dir,
    export_name = "L3_celltype_merged_VEC_merged_LEC_hippo_targets"
    )

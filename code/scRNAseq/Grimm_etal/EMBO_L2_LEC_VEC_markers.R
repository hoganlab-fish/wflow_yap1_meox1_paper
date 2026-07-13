#! /usr/bin/Rscript

#   First author : Michelle Meier

################################################################################
# LOAD FUNCTIONS
################################################################################
# source seurat functions
library(tidyverse)
library(here)
library(Seurat)

################################################################################
# LOAD DATA
################################################################################
filename <- here("output/data/Datsets/paper/EMBO_JOURNAL_Level_02_seuratObject_LM.RDS")
Level_02_seuratObject <- readRDS(filename)

################################################################################
# GET MARKERS (Level 2)
################################################################################
# Level_02_seuratObject <- SetIdent(Level_02_seuratObject, value = 'Stage')
# Level_02_seuratObject_no40 <- subset(Level_02_seuratObject, idents = "Stage_40hpf", invert = TRUE)
# Level_02_seuratObject_no40 <- SetIdent(Level_02_seuratObject_no40, value = 'L2_Predicted_phenotype')
# markers_lec_embo <- FindMarkers(object = Level_02_seuratObject_no40,
#                                 ident.1 = 'LEC',
#                                 assay = "RNA",
#                                 logfc.threshold = 0)
# markers_vec_embo <- FindMarkers(object = Level_02_seuratObject_no40,
#                                 ident.1 = 'VEC',
#                                 assay = "RNA",
#                                 logfc.threshold = 0)
#
# fileName <- '/team_folders/hogan_lab/Hogan_Lab_Projects/yap1_project/scrnaseq/output/Dataset_Level_03/Markers/EMBO_L2_VECLEC_markers_3_5dpf.RDS'
# saveRDS(object = list(VEC = markers_vec_embo, LEC = markers_lec_embo), file = fileName)


################################################################################
# GET MARKERS (Level 2),for each timepoint
################################################################################
# this is based on the analysis performed for https://link.springer.com/article/10.15252/embj.2022112590
Level_02_seuratObject <- SetIdent(Level_02_seuratObject, value = 'L2_Phenotype_Stage')

# ---- write a function so we can speed things up a bit ----
run_LEC_VEC_markers <- function(seurat = Level_02_seuratObject,
                                timepoint){
  lec_idents <- c(sprintf("LEC_Stage_%idpf", timepoint), sprintf("LEC_prox1a_low_Stage_%idpf", timepoint))
  vec_idents <- c(sprintf("VEC_Stage_%idpf", timepoint))

  tmp_markers <- FindMarkers(object = seurat,
                             ident.1 = lec_idents,
                             ident.2 = vec_idents,
                             assay = "RNA",
                             logfc.threshold = 0.25, #threshold from EMBO paper
                             min.pct = 0.1) #threshold from EMBO paper
  #add expression ratio
  tmp_markers <- tmp_markers %>%
    mutate(., ExpnRatio = ifelse(avg_log2FC < 0, -(1/(pct.1/pct.2)), pct.1/pct.2),
              ExpnRatio = case_when(is.infinite(ExpnRatio) & pct.1 == 0 ~ 10000*pct.2, #replace infinite values with something as they're actually useful
                                    is.infinite(ExpnRatio) & pct.2 == 0 ~ 10000*pct.1,
                                    TRUE ~ ExpnRatio),
              direction = ifelse(avg_log2FC < 0, "VEC", "LEC"),
              timepoint = paste(timepoint, "dpf")) %>%
    rownames_to_column("Gene")

  #filter based on thresholds in EMBO paper
  genes <- tmp_markers %>%
    filter(., abs(ExpnRatio) > 1.5) #both up and down

  return(genes)

}

# ---- run for 3,4,5 dpf ----
list_all_timepoints <- lapply(c(3,4,5), function(TP) run_LEC_VEC_markers(timepoint = TP)) %>%
  bind_rows()


# ---- Get genes for LEC/VEC across all timepoints - to give genes labels ----
# note that this will ignore fold changes etc, this is for marking the genes as LEC or VEC markers

lec_vec <- list_all_timepoints %>%
  group_by(Gene) %>%
  mutate(., always_same_direction = ifelse(all(avg_log2FC > 0)|all(avg_log2FC <0), "yes", "no")) %>% #make sure the directionality is always lec or always vec
  filter(., always_same_direction == "yes") %>%
  select(., Gene, direction) %>%
  distinct() %>%
  filter(., !(grepl("si:|CABZ|zgc:|im:" ,Gene))) #EMBO paper also cleans up more and removes unannotated genes


lec_vec$direction %>% table()
# LEC  VEC
# 2583  321

# ---- Get genes for LEC/VEC across all timepoints - for LEC/VEC scoring ----
lec_vec_scoring_genes <- list_all_timepoints %>%
  filter(., Gene %in% lec_vec$Gene & p_val_adj < 0.05) %>% # use genes from above + add significance threshold
  group_by(Gene) %>%
  reframe(., sum_log2fc = sum(avg_log2FC),
                direction = direction) %>%  #this is just for ranking them
  distinct() %>%
  group_by(direction) %>%
  top_n(., n = 25, wt = abs(sum_log2fc))


# ---- Export for plotting later ----
#only genes
filename <- here("data/genesets/EMBO_L2_VECLEC_markers_EMBOmethod_3_5dpf_genesforlabels_only.csv")
write_csv(file = filename, x = lec_vec, col_names = T)

#for scoring
filename <- here("data/genesets/EMBO_L2_VECLEC_markers_EMBOmethod_3_5dpf_genesforscoring_only.csv")
write_csv(file = filename, x = lec_vec_scoring_genes, col_names = T)

#everything
filename <- here("data/genesets/EMBO_L2_VECLEC_markers_EMBOmethod_3_5dpf.csv")
write_csv(file = filename, x = list_all_timepoints, col_names = T)







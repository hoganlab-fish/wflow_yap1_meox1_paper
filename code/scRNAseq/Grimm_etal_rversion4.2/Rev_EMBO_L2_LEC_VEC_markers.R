#! /usr/bin/Rscript

#   First author : Michelle Meier

#   Last edited October 2024

################################################################################
# LOAD FUNCTIONS
################################################################################
# ---- Load 00_seurat_function ----
# source seurat functions
source("/team_folders/hogan_lab/Hogan_Lab_Scripts/hogan_lab_bitbucket/Hogan_Lab_Scripts/R_scripts/00_ArchR_functions_R4.2.0.Core.R")
library(scales)		# for rescaling color in bar plots
library(Seurat)
library(tidyr)
library(dplyr)

################################################################################
# LOAD DATA
################################################################################
filename <- '/team_folders/hogan_lab/Hogan_Lab_Projects/yap1_project/scrnaseq/output/R_Data_files/EMBO_JOURNAL_Level_02_seuratObject_LM.RDS'
Level_02_seuratObject <- readRDS(filename)

################################################################################
# GET MARKERS (Level 2)
################################################################################
Level_02_seuratObject <- SetIdent(Level_02_seuratObject, value = 'Stage')
Level_02_seuratObject_no40 <- subset(Level_02_seuratObject, idents = "Stage_40hpf", invert = TRUE)
Level_02_seuratObject_no40 <- SetIdent(Level_02_seuratObject_no40, value = 'L2_Predicted_phenotype')
markers_lec_embo <- FindMarkers(object = Level_02_seuratObject_no40, 
                                ident.1 = 'LEC', 
                                assay = "RNA", 
                                logfc.threshold = 0)
markers_vec_embo <- FindMarkers(object = Level_02_seuratObject_no40, 
                                ident.1 = 'VEC', 
                                assay = "RNA",
                                logfc.threshold = 0)

fileName <- '/team_folders/hogan_lab/Hogan_Lab_Projects/yap1_project/scrnaseq/output/Dataset_Level_03/Markers/EMBO_L2_VECLEC_markers_3_5dpf.RDS'
saveRDS(object = list(VEC = markers_vec_embo, LEC = markers_lec_embo), file = fileName)


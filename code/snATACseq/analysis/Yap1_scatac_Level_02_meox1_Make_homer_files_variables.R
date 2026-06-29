#! /usr/bin/Rscript

#   First author : Michelle Meier

#   Last edited November 2023

################################################################################
# START UP RENV + LOAD FUNCTIONS
################################################################################
# ---- Start up renv ----
renv_path <- '/team_folders/hogan_lab/Hogan_Lab_People/Michelle_Meier/scripts/publications/SK_yap1_meox1/scatacseq'
setwd(renv_path)
renv::activate(renv_path)
renv::restore()

# ---- Load ArchR functions ----
# source("/team_folders/hogan_lab/Hogan_Lab_Scripts/hogan_lab_bitbucket/Hogan_Lab_Scripts/R_scripts/00_ArchR_functions_R4.2.0.Core.R")

# ---- Load libraries ----
library(dplyr)
library(tidyr)
# library(ArchR)

library(BSgenome.Drerio.UCSC.danRer11)
library(org.Dr.eg.db)
library(biomaRt)
library(GO.db)
library(httr)

# ---- set threads ----
# addArchRThreads(threads = 8) 

################################################################################
# IN-SCRIPT FUNCTIONS
################################################################################

################################################################################
# GET DATA
################################################################################


################################################################################
# SET VARIABLES
################################################################################
project.name <- 'Yap1_atac_Level_02'
save.dir <- "/hogan_lab/Hogan_Lab_Projects/yap1_project/scatacseq/Dataset_Level_02_Endothelial_cells_resubset_recalcPeaks/PeakFiles"; dir.create(save.dir, showWarnings = F)

################################################################################
# GET DAPS
################################################################################
# calculated in /team_folders/hogan_lab/Hogan_Lab_People/Michelle_Meier/scripts/publications/SK_yap1_meox1/scatacseq/Yap1_scatac_DAP.R
dap.table.all <- read.table('/team_folders/hogan_lab/Hogan_Lab_Projects/yap1_project/scatacseq/Dataset_Level_02_Endothelial_cells_resubset_recalcPeaks/DAP_analysis/231018_VEC_Mutant_vs_VEC_WT_all_DAPs_filteredBy_rawPval_Log2FC.txt',
                            sep = '\t', header = T)


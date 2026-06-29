#! /usr/bin/Rscript

#   First author : Michelle Meier

#   Last edited July 2023

################################################################################
# START UP RENV + LOAD FUNCTIONS
################################################################################
# ---- Start up renv ----
.libPaths('/team_folders/hogan_lab/Hogan_Lab_People/Michelle_Meier/scripts/publications/SK_yap1_meox1/scatacseq/renv/library/R-4.2/x86_64-pc-linux-gnu/')
renv_path <- '/team_folders/hogan_lab/Hogan_Lab_People/Michelle_Meier/scripts/publications/SK_yap1_meox1/scatacseq'
setwd(renv_path)
library(renv)
renv::activate(renv_path)
renv::restore()

# ---- Load ArchR functions ----
source("/team_folders/hogan_lab/Hogan_Lab_Scripts/hogan_lab_bitbucket/Hogan_Lab_Scripts/R_scripts/00_ArchR_functions_R4.2.0.Core.R")
source("/team_folders/hogan_lab/Hogan_Lab_People/Michelle_Meier/scripts/publications/SK_yap1_meox1/scatacseq/DRAFT_AddModuleScore_atac.R")

# ---- Load libraries ----
library(ArchR)
library(plyr)
library(Seurat)
library(parallel)
library(extrafont)
library(GenomicRanges)
library(BSgenome.Drerio.UCSC.danRer11)
library(dplyr)

################################################################################
# IN-SCRIPT FUNCTIONS 
################################################################################


################################################################################
# GET DATA 
################################################################################
# ---- load projects ----
path.project.4dpf.level02 <- '/hogan_lab/Hogan_Lab_Projects/meox1_project/output/scatacseq/Dataset_4dpf_Level_02_proxpaper/'
path.project.yap1.level02 <- '/hogan_lab/Hogan_Lab_Projects/yap1_project/scatacseq/Dataset_Level_02_Endothelial_cells_resubset_recalcPeaks'
project.4dpf.level02  <- loadArchRProject(path.project.4dpf.level02, showLogo=FALSE)
project.yap1.level02 <- loadArchRProject(path.project.yap1.level02, showLogo=FALSE)


# ---- general paths ----
save.dir <- "/hogan_lab/Hogan_Lab_Projects/meox1_project/output/scatacseq/Dataset_4dpf_Level_02_proxpaper/MarkerPeaks"



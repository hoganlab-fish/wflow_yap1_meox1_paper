#!/usr/bin/Rscript

# To load Macs2
# module load  macs/2.2.7.1

# To load ArchR
# cd /team_folders/hogan_lab/Hogan_Lab_People/Lizzie_Mason/r_environments/4.2.0.Core_atac/
# module load 4.2.0_Core_atac
# R

# srun -p prod_med -n 4 --time=0-03:00 --mem=64gb "$@" --pty -u bash -i

#===============================================================================================

# Load libraries 

library(ArchR)
library(plyr)
library(Seurat)
library(parallel)
suppressWarnings(suppressMessages(library(ArchR, quietly=TRUE, verbose=FALSE)))

#===============================================================================================
#   DEFINE ALL FUNCTIONS

#===============================================================================================
#   Hardcoded variables

#   20171 = Mutant
#   20172 = WT

meox_sample_names <- c("Mutant", "WT")

markerGenes <- c(
  "cdh5",
  "dll4",
  "flt1",
  "gata6",
  "kdr",
  "kdrl",
  "lyve1b",
  "prox1a",
  "dab2",
  "gpr182",
  "hand2",
  "ccn2a",
  "hapln3",
  "wwtr1",
  "cdh6",
  "meox1")

marker_genes <- c(
  "ccn2a",
  "prox1a",
  "hapln3",
  "wwtr1",
  "cdh6",
  "meox1",
  "cdh1",
  "cdh5",
  "flt4",
  "cldn11b",
  "dab2",
  "dll4",
  "esama",
  "fn1a",
  "hand2",
  "kdrl",
  "mafa",
  "mafba",
  "mafb",
  "osr2",
  "vegfc",
  "mki67",
  "angpt2a",
  "pcna",
  "mrc1a",
  "plk1",
  "cdk1",
  "flt1",
  "lyve1b",
  "hexa")

#===============================================================================================

# Load ArchR functions
source("../../../functions/00_ArchR_functions_R4.2.0.Core.R")

# Load archr project
ArchR::addArchRThreads(threads=16)

# Ensure separate input and output directories
project_path_out <- "/team_folders/hogan_lab/Hogan_Lab_Projects/yap1_project/scatacseq/Dataset_Level_02_Endothelial_cells_resubset"
project <- loadArchRProject("/team_folders/hogan_lab/Hogan_Lab_Projects/yap1_project/scatacseq/Dataset_Level_01/", showLogo=FALSE)

#===============================================================================================

#---- Subset from Level 01 to Level 02
#---- Recompute clustering and peak calling

# Set variables that are used in these steps
clusters_to_keep      <- c("C24","C25")
data_level            <- paste("02")
clustering_resolution <- 0.8

# Subset the clusters of interest from Level 01

  idxSample <- BiocGenerics::which(project$Clusters %in% clusters_to_keep)
  cellsSample <- project$cellNames[idxSample]
  project <- project[cellsSample, ]

  project <- ArchR::saveArchRProject(
    ArchRProj = project, outputDirectory = project_path_out, load = TRUE
  )

# Perform dimension reduction
project <- addIterativeLSI(
    ArchRProj = project,
    useMatrix = "TileMatrix", 
    name = paste("Level_",data_level,"_IterativeLSI",sep=""), 
    iterations = 2, 
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.2), 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:30
)

# Perform clustering on the new UMAP/tSNE object (reduced dims)
project <- addClusters(
    input = project,
    reducedDims = paste("Level_",data_level,"_IterativeLSI",sep=""),
    method = "Seurat",
    name = "Clusters",
    resolution = clustering_resolution,
    force=TRUE
)

# Perform peak calling again as this is done on clusters
pathToMacs2 <- findMacs2()

# Select the new "clusters" and make sure they are different to Level 01
project <- addReproduciblePeakSet(
    ArchRProj = project, 
    groupBy = "Clusters", 
    pathToMacs2 = pathToMacs2
)

# Add imputation weights

project <- addImputeWeights(project)

project <- ArchR::saveArchRProject(ArchRProj = project, outputDirectory = project_path_out, load = TRUE)

#===============================================================================================

#----- UMAPS

p1 <- plotEmbedding(ArchRProj = project, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = project, colorBy = "cellColData", name = "ClustersPhenotype", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = project, colorBy = "cellColData", name = "Sample", embedding = "UMAP")

ggAlignPlots(p1, p2, p3, type = "h")

plotPDF(p1,p2,p3, name = "Plot-UMAP-Clusters-Phenotype-Genotype.pdf", ArchRProj = project, addDOC = FALSE, width = 5, height = 5)

#   Plot key markers of fate for GeneScore
plot_marker_genes(
  ArchRProject=project,
  MatrixToPlot="GeneScoreMatrix",
  geneList=markerGenes,
  colourScheme="horizonExtra")



#----- Track Plots

#   Plot key markers of fate for Track plots per Cluster

 p <- plotBrowserTrack(
    ArchRProj = project,
    groupBy = "Clusters",
    geneSymbol = marker_genes,
    upstream = 50000,
    downstream = 50000,
    loops = getPeak2GeneLinks(project)
  )
  plotPDF(
    plotList = p, name = "Plot-Tracks-Marker-Genes-with-Peak2GeneLinks.pdf",
    ArchRProj = project, addDOC = FALSE, width = 5, height = 5
  )

sessionInfo()

#!/usr/bin/Rscript

# To load ArchR
# module load R/4.2.0.Core
# cd /team_folders/hogan_lab/Hogan_Lab_People/Lizzie_Mason/r_environments/4.2.0.Core_atac/
# R

# srun -p prod_med -n 4 --time=0-04:00 --mem=64gb "$@" --pty -u bash -i

#===============================================================================================

# Load libraries 

library(ArchR)
library(plyr)
library(Seurat)
library(parallel)
library(ggpubr)
suppressWarnings(suppressMessages(library(ArchR, quietly=TRUE, verbose=FALSE)))

# Load ArchR functions
source("../../../functions/00_ArchR_functions_R4.2.0.Core.R")

#===============================================================================================
#   DEFINE ALL FUNCTIONS

#===============================================================================================
#   Hardcoded variables

GA_colours <- c("#d9d9d9", "#7a0177")

# Endothelial specific genes
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

# Level 01 cell phenotype classification genes
# Use these for a large dotplot in the Supplementary information
marker_genes <- c(
"cdh1",
"epcam",
"krt4",
"tp63",
"nova2",
"mdka",
"pax6a",
"dll4",
"tnc",
"notch3",
"rhbg",
"her15.2",
"ascl1b",
"nrnx1a",
"prox1a",
"crx",
"fli1a",
"gata1a",
"myb",
"runx1",
"hand2",
"runx1",
"scfd2",
"fli1",
"twist1a",
"dmn3a",
"fn1a",
"dmn3a",
"flt1",
"kdrl",
"esama",
"cdh5",
"dll4",
"pdgfra",
"lyve1a",
"lyve1b",
"mrc1a",
"cdh6",
"flt4")

level_01_dotPlot_genes <- c(
  "cdh6",
  "prox1a",
  "stab2",
  "lyve1b",
  "mrc1a",
  "flt1",
  "dll4",
  "esm1",
  "hlx1",
  "hand2",
  "fn1a",
  "gata2b",
  "runx1",
  "alas2",
  "blvrb",
  "spi1b",
  "mpx",
  "neurod4",
  "nova2",
  "cdh1",
  "epcam",
  "myod1",
  "myog",
  "pmela",
  "dct",
  "lum",
  "col1a1a",
  "krt4",
  "cyt1",
  "cryba1b",
  "crybb1",
  "ednrba",
  "pnp4a",
  "gch2",
  "aox5")

#===============================================================================================

# Load objects

# Load archr project
ArchR::addArchRThreads(threads=16)
project_path_out <- "/team_folders/hogan_lab/Hogan_Lab_Projects/yap1_project/scatacseq/Dataset_Level_01"
project <- loadArchRProject("/team_folders/hogan_lab/Hogan_Lab_People/Tyrone_Chen/yap/results/archr_level_1_ensembl_newgenes/", showLogo=FALSE)

#===============================================================================================

# Add imputation weights

project <- addImputeWeights(project)

project <- ArchR::saveArchRProject(ArchRProj = project, outputDirectory = project_path_out, load = TRUE)

#===============================================================================================

#   Plot key markers of fate for GeneScore

project <- addImputeWeights(project)
  
# Add cluster metadata to archr project
mapping_clusterid <- c(
  "C1" = "Neuronal_01",
  "C2" = "Neuronal_02",
  "C3" = "Neuronal_03",
  "C4" = "Neuronal_04",
  "C5" = "Neuronal_05",
  "C6" = "Neuronal_06",
  "C7" = "Unknown_01",
  "C8" = "Neuronal_07",
  "C9" = "Neuronal_08",
  "C10" = "Neuronal_09",
  "C11" = "Fibroblast_01",
  "C12" = "Fibroblast_02",
  "C13" = "Unknown_02",
  "C14" = "Unknown_03",
  "C15" = "Unknown_04",
  "C16" = "Epidermis_01",
  "C17" = "Epidermis_02",
  "C18" = "Epidermis_03",
  "C19" = "Epidermis_04",
  "C20" = "Epidermis_05",
  "C21" = "Epidermis_06",
  "C22" = "Hematopoietic_01",
  "C23" = "Hematopoietic_02",
  "C24" = "Endothelial_01",
  "C25" = "Endocardium_01")

mapping_phenotype <- c(
  "C1" = "Neuronal",
  "C2" = "Neuronal",
  "C3" = "Neuronal",
  "C4" = "Neuronal",
  "C5" = "Neuronal",
  "C6" = "Neuronal",
  "C7" = "Unknown",
  "C8" = "Neuronal",
  "C9" = "Neuronal",
  "C10" = "Neuronal",
  "C11" = "Fibroblast",
  "C12" = "Fibroblast",
  "C13" = "Unknown",
  "C14" = "Unknown",
  "C15" = "Unknown",
  "C16" = "Epidermis",
  "C17" = "Epidermis",
  "C18" = "Epidermis",
  "C19" = "Epidermis",
  "C20" = "Epidermis",
  "C21" = "Epidermis",
  "C22" = "Hematopoietic",
  "C23" = "Hematopoietic",
  "C24" = "Endothelial",
  "C25" = "Endocardium")
  
   project@cellColData$ClustersID <- project@cellColData$Clusters
   project@cellColData$ClustersPhenotype <- project@cellColData$Clusters
  
   project@cellColData$ClustersID <- revalue(
     project@cellColData$ClustersID, mapping_clusterid
   )
   project@cellColData$ClustersPhenotype <- revalue(
     project@cellColData$ClustersPhenotype, mapping_phenotype
   )

  project <- ArchR::saveArchRProject(
    ArchRProj = project, outputDirectory = project_path_out, load = TRUE
  )


yap_sample_names <- c(
  "mut" = "Mutant",
  "wt" = "WT")

  project@cellColData$Genotype <- as.character(project@cellColData$Sample)
  
  project@cellColData$Genotype <- revalue(
     project@cellColData$Genotype, yap_sample_names
   )

  project@cellColData$Phenotype_Genotype <- paste(project@cellColData$ClustersPhenotype, project@cellColData$Genotype, sep="_")

  project <- ArchR::saveArchRProject(
    ArchRProj = project, outputDirectory = project_path_out, load = TRUE
  )


#----- DAP and DAG analysis

DAG_Endothelial_open_in_mutant <- dag_analysis(
  project=project,
  groupBy="Phenotype_Genotype",
  useGroup="Endothelial_Mutant",
  bgdGroup="Endothelial_WT",
  exportFilePath = project_path_out,
  exportFileName = "Level_01_DAG_Endothelial_Mutant_vs_WT",
  chrStatus = "open",
  FDRFilter = FALSE,
  FDRThreshold = 0.05,
  RawPvalFilter = TRUE,
  RawPvalThreshold = 0.05)

DAG_Endothelial_closed_in_mutant <- dag_analysis(
  project=project,
  groupBy="Phenotype_Genotype",
  useGroup="Endothelial_Mutant",
  bgdGroup="Endothelial_WT",
  exportFileName = "Level_01_DAG_Endothelial_Mutant_vs_WT",
  chrStatus = "closed",
  FDRFilter = FALSE,
  FDRThreshold = 0.05,
  RawPvalFilter = TRUE,
  RawPvalThreshold = 0.05)


DAP_Endothelial_open_in_mutant <- dap_analysis(
  project=project,
  groupBy="Phenotype_Genotype",
  useGroup="Endothelial_Mutant",
  bgdGroup="Endothelial_WT",
  exportFileName = "Level_01_DAP_Endothelial_Mutant_vs_WT",
  chrStatus = "open",
  FDRFilter = FALSE,
  FDRThreshold = 0.05,
  RawPvalFilter = TRUE,
  RawPvalThreshold = 0.05)

DAP_Endothelial_closed_in_mutant <- dap_analysis(
  project=project,
  groupBy="Phenotype_Genotype",
  useGroup="Endothelial_Mutant",
  bgdGroup="Endothelial_WT",
  exportFileName = "Level_01_DAP_Endothelial_Mutant_vs_WT",
  chrStatus = "closed",
  FDRFilter = FALSE,
  FDRThreshold = 0.05,
  RawPvalFilter = TRUE,
  RawPvalThreshold = 0.05)

#----- UMAPS

p1 <- plotEmbedding(ArchRProj = project, colorBy = "cellColData", name = "ClustersID", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = project, colorBy = "cellColData", name = "ClustersPhenotype", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = project, colorBy = "cellColData", name = "Sample", embedding = "UMAP")

ggAlignPlots(p1, p2, p3, type = "h")

plotPDF(p1,p2,p3, name = "Plot-UMAP-Clusters-Phenotype-Genotype.pdf", ArchRProj = project, addDOC = FALSE, width = 5, height = 5)

#----- Heatmap of marker genes per phenotype
# Identify markers genes per phenotype
#markersGS <- ArchR::getMarkerFeatures(
#    ArchRProj = project,
#    useMatrix = "GeneScoreMatrix",
#    groupBy = "ClustersPhenotype",
#    bias = c("TSSEnrichment", "log10(nFrags)"),
#    testMethod = "wilcoxon"
#  )
#
#heatmapPeaks <- plotMarkerHeatmap(
#  seMarker = markersGS, 
#  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
#  transpose = TRUE
#)
#
#plotPDF(
#    heatmapPeaks,
#    name = "Plot-Heatmap-markerGenes.pdf",
#    ArchRProj = project,
#    addDOC = FALSE,
#    width = 5,
#    height = 5
#  )
#
##----- Heatmap of marker peaks per phenotype
## Imputed gene score matrix
#project  <- addImputeWeights(project)
#
## Identify marker peaks per phenotype
#markersPeaks <- getMarkerFeatures(
#    ArchRProj = project, 
#    useMatrix = "PeakMatrix", 
#    groupBy = "ClustersPhenotype",
#    bias = c("TSSEnrichment", "log10(nFrags)"),
#    testMethod = "wilcoxon"
#)
#
#markerList    <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
#
#heatmapPeaks <- plotMarkerHeatmap(
#  seMarker = markersPeaks, 
#  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
#  transpose = TRUE
#)
#
#plotPDF(
#    heatmapPeaks,
#    name = "Plot-Heatmap-markerPeaks.pdf",
#    ArchRProj = project,
#    addDOC = FALSE,
#    width = 5,
#    height = 5
#  )


# Make UMAPs for gene accessibility (gene score)
plot_marker_genes(
  ArchRProject=project,
  MatrixToPlot="GeneScoreMatrix",
  geneList=markerGenes,
  colourScheme="horizonExtra")


# Plot key markers of fate for Track plots per Cluster

 p <- plotBrowserTrack(
    ArchRProj = project,
    groupBy = "ClustersPhenotype",
    geneSymbol = markerGenes,
    upstream = 50000,
    downstream = 50000,
    loops = getPeak2GeneLinks(project)
  )
  plotPDF(
    plotList = p, name = "Plot-Tracks-Marker-Genes-with-Peak2GeneLinks.pdf",
    ArchRProj = project, addDOC = FALSE, width = 5, height = 5
  )


  project <- ArchR::saveArchRProject(
    ArchRProj = project, outputDirectory = project_path_out, load = TRUE
  )


#----- FIGURES FOR PAPER - SUPPLEMENTARY INFORMATION AT LEVEL 01

#   Define the cluster markers using gene accessibiity
markerGenes_ClustersID        <- makeMarkerGeneTable(
                                        project = project,
                                        groupBy="ClustersID",
                                        FDRFilter = TRUE,
                                        FDRThreshold = 0.05,
                                        exportFileName=c("Level_01_DAG_cluster_markers_ClustersID"))

markerGenes_ClustersPhenotype <- makeMarkerGeneTable(
                                        project = project,
                                        groupBy="ClustersPhenotype",
                                        FDRFilter = TRUE,
                                        FDRThreshold = 0.05,
                                        exportFileName=c("Level_01_DAG_cluster_markers_ClustersPhenotype"))

#   Define data frame for plotting custom UMAPS
umap                    <- project@embeddings$UMAP$df
umap$Genotype           <- project$Genotype
umap$ClustersID         <- project$ClustersID
umap$ClustersPhenotype  <- project$ClustersPhenotype
umap$Phenotype_Genotype <- project$Phenotype_Genotype
colnames(umap)          <- c("dim1", "dim2", "Genotype", "ClustersID", "ClustersPhenotype", "Phenotype_Genotype")

# Reorder factors for plotting
ClustersPhenotype_ordered <- c(
  "Endothelial",
  "Endocardium",
  "Hematopoietic",
  "Neuronal",
  "Fibroblast",
  "Epidermis",
  "Unknown")

ClustersID_ordered <- c(
  "Endothelial_01",
  "Endocardium_01",
  "Hematopoietic_01",
  "Hematopoietic_02",
  "Neuronal_01",
  "Neuronal_02",
  "Neuronal_03",
  "Neuronal_04",
  "Neuronal_05",
  "Neuronal_06",
  "Neuronal_07",
  "Neuronal_08",
  "Neuronal_09",
  "Fibroblast_01",
  "Fibroblast_02",
  "Epidermis_01",
  "Epidermis_02",
  "Epidermis_03",
  "Epidermis_04",
  "Epidermis_05",
  "Epidermis_06",
  "Unknown_01",
  "Unknown_02",
  "Unknown_03",
  "Unknown_04"
  )

umap$ClustersID         <- factor(umap$ClustersID, levels=ClustersID_ordered)
umap$ClustersPhenotype  <- factor(umap$ClustersPhenotype, levels=ClustersPhenotype_ordered)
umap$Genotype           <- factor(umap$Genotype, levels=c("Mutant", "WT"))

cluster_colours     <- ggplotColours(length(levels(umap$ClustersID)))
phenotype_colours   <- ggplotColours(length(levels(umap$ClustersPhenotype)))


#-----LEVEL 01 UMAPS
umap_clusters   <- ggplot(umap, aes(x=dim1, y=dim2, color=ClustersID)) + 
                    geom_point(size=0.5) + scale_color_manual(values=cluster_colours) + 
                    guides(colour = guide_legend(override.aes = list(size=2))) +
                    theme_classic() +
                    ggtitle(label="ArchR clusters")

umap_phenotype  <- ggplot(umap, aes(x=dim1, y=dim2, color=ClustersPhenotype)) + 
                    geom_point(size=0.5) + scale_color_manual(values=phenotype_colours) +
                    guides(colour = guide_legend(override.aes = list(size=2))) +
                    theme_classic() +
                    ggtitle(label="Cell Phenotype")


umap_genotype   <- ggplot(umap, aes(x=dim1, y=dim2, color=Genotype)) + 
                    geom_point(size=0.5) + scale_color_manual(values=c("#647D97", "#232324")) +
                    guides(colour = guide_legend(override.aes = list(size=2))) +
                    theme_classic() +
                    ggtitle(label="Genotype")

figure_1_umaps  <- ggarrange(umap_clusters,
                            umap_phenotype,
                            umap_genotype,
                            labels = c("a","b","c"),
                            ncol = 3,
                            nrow = 1,
                            align="v")

ggsave(paste(project_path_out,"/Plots/Supplementary_Figure_UMAP_Level_01_cluster_phenotype.EPS",sep=""),
     figure_1_umaps,
     device="eps",
     width=70,
     height=15,
     units="cm",
     scale=1)

#-----LEVEL 01 GENE ACCESSIBILITY DOT PLOT

make_dot_plot_ATAC(
  project=project,
  groupBy=c("ClustersPhenotype"),
  useGroups=c("Endothelial", "Endocardium", "Hematopoietic", "Neuronal", "Fibroblast", "Epidermis", "Unknown"),
  groupName=c("Endothelium", "Endocardium", "Hematopoietic", "Neuronal", "Fibroblast", "Epidermis", "Unknown"),
  markerGenes=level_01_dotPlot_genes,
  exportFileName=c("Supp_figure_01_Level_01_marker_genes"),
  exportFilePath = c("/team_folders/hogan_lab/Hogan_Lab_Projects/yap1_project/scatacseq/Dataset_Level_01/Figures_yap_paper"),
  accessibilityColour = GA_colours,
  plotWidth = 10,
  plotHeight = 20,
  rotateaxis = TRUE)

sessionInfo()


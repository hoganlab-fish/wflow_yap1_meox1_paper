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


# ---- set threads ----
# addArchRThreads(threads = 8)

################################################################################
# IN-SCRIPT FUNCTIONS
################################################################################
createGeneAnnotation_adapted <- function(
    genome = NULL,
    TxDb = NULL,
    OrgDb = NULL,
    genes = NULL,
    exons = NULL,
    TSS = NULL){
  
  
  if(is.null(genes) | is.null(exons) | is.null(TSS)){
    
    inGenes <- genes
    inExons <- exons
    inTSS <- TSS
    
    
    ###########################
    message("Getting Genes..")
    genes <- GenomicFeatures::genes(TxDb)
    #This is where we adapt the ArchR function to not use AnnotationDbi for conversion but biomart!
    
    message("This has been adapted for zebrafish, will not work with any other genomes. Please use default ArchR funtion for other organisms")
    
    ###########################
    # biomart is very slow, so we precalculate the gene id --> zfin symbol translations
    translation.table <- readRDS('/team_folders/hogan_lab/genomes/DRERIO_geneMap_ensembl_gene_id_to_zfin_id_symbol.RDS')
    # make sure trhey are in the same order
    merged.genes <- merge(names(genes) %>% as.data.frame(), translation.table, by.x = '.', by.y = 'ensembl_gene_id', all = T, sort = F)
    mcols(genes)$symbol <- merged.genes$zfin_id_symbol
    names(genes) <- NULL 
    genes <- sort(sortSeqlevels(genes), ignore.strand = TRUE) 
    
    ###########################
    message("Getting Exons..")
    exons <- unlist(GenomicFeatures::exonsBy(TxDb, by = "tx"))
    exons$tx_id <- names(exons)
    mcols(exons)$gene_id <- suppressMessages(AnnotationDbi::select(TxDb, keys = paste0(mcols(exons)$tx_id), column = "GENEID", keytype = "TXID")[, "GENEID"])
    exons <- exons[!is.na(mcols(exons)$gene_id), ]
    # merge and keep order
    merged.exons <- merge(mcols(exons), translation.table, by.x = 'gene_id', by.y = 'ensembl_gene_id', all = T, sort = F)
    mcols(exons)$symbol <- merged.exons$zfin_id_symbol
    names(exons) <- NULL
    mcols(exons)$exon_id <- NULL
    mcols(exons)$exon_name <- NULL
    mcols(exons)$exon_rank <- NULL
    mcols(exons)$tx_id <- NULL
    exons <- sort(sortSeqlevels(exons), ignore.strand = TRUE)
    
    ###########################
    message("Getting TSS..")
    TSS <- unique(GenomicRanges::resize(GenomicFeatures::transcripts(TxDb), width = 1, fix = "start"))
    
    if(!is.null(inGenes)){
      genes <- .validGRanges(inGenes)
    }
    
    if(!is.null(inExons)){
      exons <- .validGRanges(inExons)
    }
    
    if(!is.null(inTSS)){
      TSS <- .validGRanges(inTSS)
    }
    
  }else{
    
    genes <- .validGRanges(genes)
    exons <- .validGRanges(exons)
    TSS <- unique(.validGRanges(TSS))
    
  }
  
  SimpleList(genes = genes, exons = exons, TSS = TSS)
  
}

make_dotplot <- function(project,
                         groupBy,
                         useGroups,
                         groupName,
                         markerGenes,
                         accessibilityColour = c("#d9d9d9","#d7301f"),
                         rotateaxis = FALSE){
  
  #		Arguments:
  #		project = ArchrProject
  #		groupBy = The column name in cellColData used for grouping cells together for marker feature identification (e.g. "Level_02_Phenotype").
  #		useGroups = A character vector that is used to select a subset of groups by name from the designated groupBy column in cellColData (e.g. c("LEC","VEC","AEC")).
  #		groupName = A character vector that contains axis labels in the dot plot (the order of the group name should be corresponding to the useGroups).
  #		markerGenes = A character vector that contains a list of gene names.
  #		exportFilePath = An alternative directory path to save the dot plot to. Default is ./Plots in outputDirectory of the ArchRProject.
  #		exportFileName = A character vector that is used for the output file name (e.g. exportFileName = "Level_02_marker_genes", the output file name would be "DotPlot_Level_02_marker_genes.pdf") ).
  #		accessibilityColour = A character vector that specifies the colour range for accessibility score.
  #		plotWidth = A numeric vector that specifies the width of the dot plot. 
  #		plotHeight = A numeric vector that specifies the height of the dot plot. 
  #		rotateaxis = A boolean vector indicating whether to rotate the axis (default x axis is genes and y axis is groups).
  
  #		Example:	dotPlot 	<-			make_dot_plot_ATAC (project = project,
  #							  									groupBy = "Level_02_Phenotype)",
  #							   									useGroups = c("LEC","VEC","AEC"),
  #							   									groupName = c("Lymphatics", "Veins", "Arterials"),
  #							   									exportFileName = "Level_02_marker_genes"
  #																)
  
  #======================================================================================================
  
  # check if all the gene names are correct
  geneNamesError			<-	NULL
  
  for (genes in markerGenes) {
    if (genes %ni% project@geneAnnotation$genes$symbol) {
      geneNamesError		<-	paste(geneNamesError, genes, sep = ",")
    }
  }
  
  if (!is.null(geneNamesError)) {
    geneNamesError			<-	substring(geneNamesError, 2)
    stop(geneNamesError, " are not found in gene annotations.", "\n")
  }
  
  # check if all the useGroups are correct
  groupNamesError			<-	NULL
  
  for (groups in useGroups) {
    if (groups %ni% as.vector(project@cellColData[[groupBy]])) {
      groupNamesError		<-	paste(groupNamesError, groups, sep = ",")
    }
  }
  
  if (!is.null(groupNamesError)) {
    groupNamesError			<-	substring(groupNamesError, 2)
    stop(groupNamesError, " are not found in ", groupBy, ".", "\n")
  }
  
  
  # === Start function ===
  # get matrix of accessibility score for each cell & gene
  gene_score_matrix 			<- 	getMatrixFromProject(ArchRProj = project, useMatrix = "GeneScoreMatrix")
  matrix <- assays(gene_score_matrix)$GeneScoreMatrix %>% as.matrix() # make sure it is a matrix
  rownames(matrix) <- rowData(gene_score_matrix)$name # use gene names as rownames
  
  # Get relevant metadata as a dataframe
  metadata_df <- colData(gene_score_matrix)[which(colData(gene_score_matrix)[[groupBy]] %in% useGroups), groupBy] %>%
    as.data.frame() %>% # we need ensure this is a dataframe because the line above creates a character vector
    dplyr::rename(., GroupBy = '.') %>% # rename the column to GroupBy
    mutate(., rownames = rownames(colData(gene_score_matrix))[which(colData(gene_score_matrix)[[groupBy]] %in% useGroups)]) #add in the barcodes as a column called rownames
  
  if (length(markerGenes) == 1){ # for some reason, if only one gene the matrix is already transposed
    # Get relevant gene expression
    filtered_gex <- matrix[markerGenes,] %>% # only keep the genes we are interested in
      as.data.frame() %>% #transpose matrix and then convert into dataframe because the functions bellow only work on dfs
      merge(., metadata_df, by.x = 0, by.y = 'rownames') %>% #merge metadata dataframe into this dataframe so we can access GroupBy information
      dplyr::select(., -`Row.names`) %>% # Remove "Row.names" column because we don't need it
      tidyr::gather(., key = 'Gene', 'GeneScore', -GroupBy) %>% #convert dataframe from a wide to a long format
      dplyr::mutate(., Gene = markerGenes) # gene names are removed when there's only one column
    
  } else {
    # Get relevant gene expression
    filtered_gex <- matrix[markerGenes,] %>% # only keep the genes we are interested in
      t() %>% as.data.frame() %>% #transpose matrix and then convert into dataframe because the functions bellow only work on dfs
      merge(., metadata_df, by.x = 0, by.y = 'rownames') %>% #merge metadata dataframe into this dataframe so we can access GroupBy information
      dplyr::select(., -`Row.names`) %>% # Remove "Row.names" column because we don't need it
      tidyr::gather(., key = 'Gene', 'GeneScore', -GroupBy) #convert dataframe from a wide to a long format
  }
  
  
  # Get mean gene expression
  filtered_gex_means <- filtered_gex %>%
    dplyr::group_by(., GroupBy, Gene) %>% # group by GroupBy variable and Gene, allows us to then apply the code bellow on each grouping
    dplyr::summarise(., mean_genescore = mean(GeneScore)) #get mean gene score
  
  # Calculate proportions with > 0
  overall_n <- metadata_df %>% 
    dplyr::count(., GroupBy) %>% #count how many times each factor in GroupBy occurs
    dplyr::rename(., n_total = n)  #rename column to n_total
  
  percentages <- filtered_gex %>% 
    dplyr::group_by(., GroupBy, Gene) %>% # group by GroupBy variable and Gene, allows us to then apply the code bellow on each grouping
    dplyr::mutate(., gene_score_yes = ifelse(GeneScore == 0, 0, 1)) %>% #check if GeneScore is bigger than 0
    dplyr::filter(., gene_score_yes == 1) %>% # only keep the instances where GeneScore was > 0
    dplyr::count(., gene_score_yes) %>% # count how often GeneScore was > 0, will add in a column called n
    merge(., overall_n, by = 'GroupBy') %>% # merge with overall_n dataframe which gives us total number per group in GroupBy
    mutate(., percent = (n/n_total)*100) %>% # get percentages
    dplyr::select(., GroupBy, Gene, percent)
  
  # Make marker_accessibility_table from above 
  marker_accessibility_table <- filtered_gex_means %>%
    merge(., percentages, by = c('GroupBy', 'Gene')) %>%
    dplyr::rename(., Genes = Gene, Groups = GroupBy, Score = mean_genescore, Proportion=percent)
  
  #scale expression score across genes
  marker_accessibility_table <- transform(marker_accessibility_table, norm = ave(Score, Genes, FUN = scale))
  
  
  # change Group and Gene columns into factor
  if (isTRUE(rotateaxis)) {
    marker_accessibility_table$Groups 		<- 	factor(marker_accessibility_table$Groups, levels = groupName)
    marker_accessibility_table$Genes 		<- 	factor(marker_accessibility_table$Genes, levels = rev(markerGenes))
  }	else {
    marker_accessibility_table$Groups 		<- 	factor(marker_accessibility_table$Groups, levels = rev(groupName))
    marker_accessibility_table$Genes 		<- 	factor(marker_accessibility_table$Genes, levels = markerGenes)
  }
  
  
  # make dot plot
  p 		<- ggplot(marker_accessibility_table) +
    geom_point(aes(x= Genes, y = Groups, color = norm, size = Proportion)) + 
    scale_colour_gradient(low = accessibilityColour[1],
                          high = accessibilityColour[2],
                          na.value = "black",
                          guide = "colourbar",
                          aesthetics = "colour") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"))
  # + theme(aspect.ratio=1/1)
  
  # rotate x and y axis if required
  if (isTRUE(rotateaxis)) {
    p <- p + coord_flip()
  }
  
  return(p)
  
}


#	Specifying the ggplot color scheme
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  
  # Make the same spectrum of colours as GGPLOT
  
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

# make trackplots 
make_p3p4_trackplots <- function(project,
                                 groupBy,
                                 useGroups,
                                 name,
                                 up = 50000,
                                 down = 50000){
  region <- getPeakSet(ArchRProj=project)
  region_p4 <- region[region@seqnames == chr & region@ranges@start > p4_limits[1]  &  region@ranges@start< p4_limits[2]]
  #get P3 region
  region_p3 <- region[region@seqnames == chr & region@ranges@start > p3_limits[1]  &  region@ranges@start< p3_limits[2]]
  #combined regions 
  region_filtered <- c(region_p4, region_p3)
  
  #plot
  p_meox1 <- plotBrowserTrack(ArchRProj = project, 
                              geneSymbol = 'meox1', 
                              groupBy = groupBy, 
                              useGroups = useGroups,
                              features = region_filtered,
                              upstream = up,
                              downstream = down)
  
  plotName <- sprintf('%s_PlotTracksP3P4_down%i_up%i_meox1.pdf', name,down, up)
  plotPDF(plotList = p_meox1, 
          name = plotName, 
          ArchRProj = project, 
          addDOC = FALSE, width = 5, height = 5)
  return(p_meox1)
}


make_meox_heatmaps <- function(project,
                               markersPeaks.forplotting = markers.yap2,
                               useGroups,
                               name, 
                               order){
  
  region <- getPeakSet(ArchRProj=project)
  region_filtered <- region[region@elementMetadata@listData$nearestGene == 'meox1']
  
  #filter
  starts <- region_filtered@ranges@start
  markersPeaks.forplotting.filtered <-markersPeaks.forplotting[markersPeaks.forplotting@elementMetadata@listData$start %in% starts & markersPeaks.forplotting@elementMetadata@listData$seqnames == chr] 
  rownames(markersPeaks.forplotting.filtered) <- paste(markersPeaks.forplotting.filtered@elementMetadata@listData$seqnames,paste(markersPeaks.forplotting.filtered@elementMetadata@listData$start, markersPeaks.forplotting.filtered@elementMetadata@listData$end, sep = "-"),sep = ":")
  
  
  #pull out data to make a pretty heatmap 
  plotting.df <- assays(markersPeaks.forplotting.filtered)[['Mean']] 
  plotting.df$PeakID <- rownames(plotting.df)
  plotting.df.long <- plotting.df %>% 
    tidyr::gather(., key = 'Cluster', value = 'Score', -PeakID) %>% 
    dplyr::filter(., Cluster %in% useGroups) 
  
  plotting.df.long <- transform(plotting.df.long, Score = ave(Score, PeakID, FUN = scale))
  small.df <- plotting.df.long[plotting.df.long$Cluster %in% c('VEC_WT', 'VEC_wt'), ]
  levels <- small.df$PeakID[order(small.df$Score, decreasing = T)]
  plotting.df.long$PeakID = factor(plotting.df.long$PeakID, levels = levels)  
  
  
  heatmap.delta <- ggplot(plotting.df.long, aes(x=PeakID, y=Cluster, fill=Score)) +
    geom_tile() + theme_bw() + coord_equal() +
    scale_fill_distiller(palette="PRGn", direction=-1) +
    labs(title = "Heatmap, selected DAPs", fill = 'Scaled Score\nfrom ArchR function') +
    xlab("") + ylab("") +
    theme(axis.title.x=element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.ticks.x=element_blank())
  
  FileName <- sprintf('%s/%s_meox1_ChromatinAccessiblity_heatmap.EPS', figure.dir.out, name)
  print(FileName)
  ggsave(FileName,
         heatmap.delta,
         device="eps",
         width=15,
         height=15,
         units="cm",
         scale=1)
  
  #only show p3 and p4 + promotor and super accessible peak
  region <- getPeakSet(ArchRProj=project)
  region_filtered_p4 <- region[region@seqnames == "chr12" & region@ranges@start > p4_limits[1]  &  region@ranges@start<  p4_limits[2]]
  region_filtered_p3 <- region[region@seqnames == "chr12" & region@ranges@start > p3_limits[1]  &  region@ranges@start< p3_limits[2]]
  region_filtered_promotor <- region[region@seqnames == "chr12" & region@ranges@start > 27461900  &  region@ranges@start < 27461950]
  region_filtered_superaccess <- region[region@seqnames == "chr12" & region@ranges@start > 27469700  &  region@ranges@start<  27469780]
  
  #filter
  starts <- c(region_filtered_p3@ranges@start, region_filtered_p4@ranges@start, region_filtered_promotor@ranges@start, region_filtered_superaccess@ranges@start)
  markersPeaks.forplotting.filtered <-markersPeaks.forplotting[markersPeaks.forplotting@elementMetadata@listData$start %in% starts & markersPeaks.forplotting@elementMetadata@listData$seqnames == chr]
  rownames(markersPeaks.forplotting.filtered) <- paste(markersPeaks.forplotting.filtered@elementMetadata@listData$seqnames,paste(markersPeaks.forplotting.filtered@elementMetadata@listData$start, markersPeaks.forplotting.filtered@elementMetadata@listData$end, sep = "-"),sep = ":")
  
  
  #pull out data to make a pretty heatmap
  plotting.df <- assays(markersPeaks.forplotting.filtered)[['Mean']]
  plotting.df$PeakID <- rownames(plotting.df)
  plotting.df.long <- plotting.df %>%
    tidyr::gather(., key = 'Cluster', value = 'Score', -PeakID) %>%
    dplyr::filter(., Cluster %in% useGroups)
  plotting.df.long <- transform(plotting.df.long, Score = ave(Score, PeakID, FUN = scale))
  plotting.df.long$PeakID = factor(plotting.df.long$PeakID, levels = order)  
  
  heatmap.delta <- ggplot(plotting.df.long, aes(x=PeakID, y=Cluster, fill=Score)) +
    geom_tile() + theme_bw() + coord_equal() +
    scale_fill_distiller(palette="PRGn", direction=-1) +
    labs(title = "Heatmap, selected DAPs", fill = 'Scaled Score\nfrom ArchR function') +
    xlab("") + ylab("") +
    theme(axis.title.x=element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.ticks.x=element_blank())
  
  FileName <- sprintf('%s/%s_P3P4_meox1_ChromatinAccessiblity_heatmap.EPS', figure.dir.out, name)
  print(FileName)
  ggsave(FileName,
         heatmap.delta,
         device="eps",
         width=15,
         height=15,
         units="cm",
         scale=1)
  
}

plot_GEX_umaps <- function(project,
                           title_addition = "",
                           gene, 
                           col_gradient = c("#d9d9d9", "#7a0177"), 
                           exportFileName, 
                           exportFilePath, 
                           size_w = 15,
                           size_h = 15){
  p1 <- plotEmbedding(
    ArchRProj = project, 
    colorBy = "GeneScoreMatrix", 
    name = gene,
    embedding = "UMAP",
    size = 1, 
    imputeWeights = getImputeWeights(project))
  
  x_table <- ggplot_build(p1)$data[[1]]
  umap <- ggplot(x_table, aes(x=x, y=y)) + 
    geom_point(aes(color=value), size=1) +
    scale_colour_gradientn(colours=col_gradient) +
    theme_classic() +
    ggtitle(label=sprintf("%s expression\n%s", gene, title_addition))
  
  umap_nolegend <- umap + NoLegend() + NoAxes() + ggtitle('')
  
  fileName <- sprintf('%s/UMAP_%s_%s.pdf', exportFilePath, exportFileName, gene)
  ggsave(filename =fileName, 
         plot =  umap, 
         device = 'pdf', 
         width = size_w, 
         height = size_h,
         units="cm",
         scale=1)
  fileName <- sprintf('%s/UMAP_%s_%s_nolegend.pdf', exportFilePath, exportFileName, gene)
  ggsave(filename =fileName, 
         plot =  umap_nolegend, 
         device = 'pdf', 
         width = size_w, 
         height = size_h,
         units="cm",
         scale=1)
  return(umap)
  
}


################################################################################
# GET DATA 
################################################################################
# ---- load projects ----
path.project.yap1.level01 <- '/hogan_lab/Hogan_Lab_Projects/yap1_project/scatacseq/Dataset_Level_01'
path.project.yap1.level02 <- '/hogan_lab/Hogan_Lab_Projects/yap1_project/scatacseq/Dataset_Level_02_Endothelial_cells_resubset_recalcPeaks'
path.project.4dpf.level02 <- '/hogan_lab/Hogan_Lab_Projects/meox1_project/output/scatacseq/Dataset_4dpf_Level_02_proxpaper/'
project.yap1.level01  <- loadArchRProject(path.project.yap1.level01, showLogo=FALSE)
project.yap1.level02 <- loadArchRProject(path.project.yap1.level02, showLogo=FALSE)
project.yap1.level02 <- addImputeWeights(project.yap1.level02, reducedDims = 'IterativeLSI_Level_02') 
project.4dpf.level02  <- loadArchRProject(path.project.4dpf.level02, showLogo=FALSE)


# ---- set general variables ----
figure.dir.out <- '/team_folders/hogan_lab/Hogan_Lab_Projects/yap1_project/scatacseq/Figures_yap_paper_revision/'; dir.create(figure.dir.out, showWarnings = F, recursive = T)
p4_limits <- c(27429000, 27430200)
p3_limits <- c(27425800, 27426500)
chr <- 'chr12'

path.markers.yap2 <- "/hogan_lab/Hogan_Lab_Projects/yap1_project/scatacseq/Dataset_Level_02_Endothelial_cells_resubset_recalcPeaks/MarkerGenes/Level_02_markers_Phenotype_Genotype_Level_02.RDS"
markers.yap2  <- readRDS(path.markers.yap2)


path.fimo.results <- '/team_folders/hogan_lab/Hogan_Lab_Projects/yap1_project/scatacseq/Motif_analysis/MEME/231127_meox1_enhancers_TEAD_family_1e2/fimo.tsv'
fimo.results <- read.table(path.fimo.results, 
                           sep = '\t', 
                           header = TRUE)
# ---- get DAPs ----
# dap.table.all <- read.table('/team_folders/hogan_lab/Hogan_Lab_Projects/yap1_project/scatacseq/Dataset_Level_02_Endothelial_cells_resubset_recalcPeaks/DAP_analysis/',
#                             sep = '\t', header = T)

# ---- dotplot ordered ----
ClustersID_ordered_Level_01_dotplot <- c(
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
  'Fibroblast_01',
  'Fibroblast_02', 
  "Epidermis_01",
  "Epidermis_02",
  "Epidermis_03",
  "Epidermis_04",
  "Epidermis_05",
  "Epidermis_06",
  "Unknown_01",
  "Unknown_02",
  "Unknown_03",
  "Unknown_04")

ClustersID_ordered_Level_02_dotplot <- c(
  "Specified_LEC_01",
  "Specified_LEC_02",
  "VEC_01",
  "AEC_01",
  "AEC_02",
  "AEC_03",
  "Endothelial_01",
  "Endothelial_02",
  "Endocardium_01")

# ---- Colours ordered ----
# 
ClustersID_ordered_Level_02 <- c(
  "Specified_LEC_01",
  "Specified_LEC_02",
  "VEC_01",
  "AEC_01",
  "AEC_02",
  "AEC_03",
  "Endothelial_01",
  "Endothelial_02",
  "Endocardium_01")

L2_Phenotype_ordered_Level_02 <- c(
  "Specified_LEC",
  "VEC",
  "AEC",
  "Endothelial",
  "Endocardium")


ClustersID_ordered_cols_Level_02 <- c(
  "#a1d99b",
  '#74c476',
  "#6baed6",
  "#fee0d2",
  "#fcbba1",
  "#fc9272",
  '#ffffd4',
  '#fee391',
  "#d9d9d9"
)


genotype_cols <- c(
  "#dbe2c6",
  "#657c95"
)

genotype_order <- c(
  'WT',
  'Mutant'
)

colours.databases <- c(
  'grey20', #encode
  'grey80' #homer
)
names(colours.databases) <- c('encode', 'homer')

open.close.colours <- c(
  '#c2a5cf', #open
  '#5aae61' #closed
)
names(open.close.colours) <- c('open', 'closed')

# ---- genes ----
marker_genes <- c('meox1', 
                  'ccn2a',
                  'lyve1b', 
                  'prox1a', 
                  'flt1', 
                  'hand2')


hippo_targets_nomeox <- list(hippo = c(
  "yap1",
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
  "bmp4"))

dotplot_marker_genes_L1 <- c(
  'cdh5', #Endothelial
  'kdrl', #Endothelial
  'hand2', #Endocardium
  'fn1a', #Endocardium
  'gata2b', #Hematopoietic
  'runx1', #Hematopoietic
  'nova2', #Neuronal
  'neurod4', #Neuronal
  'lum', #Fibroblast
  'col1a1a', #Fibroblast
  'pdgfrb', #Fibroblast
  'krt4', #Epidermis
  'cyt1', #Epidermis
  'epcam' #Epidermis
)

dotplot_marker_genes_L2 <- c(
  'prox1a', #LEC
  'cdh6', #VEC/LEC
  'stab2', #LEC/VEC
  'mrc1a', #LEC/VEC
  'lyve1b', #VEC/VEC
  'ccn2a',
  'meox1',
  'flt1',#AEC, mAEC
  'dll4', #AEC, mAEC
  'hand2',#Endocardium
  'fn1a' #Endocardium
)





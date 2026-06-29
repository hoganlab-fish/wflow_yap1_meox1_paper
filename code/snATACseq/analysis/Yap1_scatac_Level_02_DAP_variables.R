#! /usr/bin/Rscript

#   First author : Michelle Meier

#   Last edited June 2023

################################################################################
# START UP RENV + LOAD FUNCTIONS
################################################################################
# ---- Start up renv ----
.libPaths('/team_folders/hogan_lab/Hogan_Lab_People/Michelle_Meier/scripts/publications/SK_yap1_meox1/scatacseq/renv/library/R-4.2/x86_64-pc-linux-gnu/')
library(renv)
renv_path <- '/team_folders/hogan_lab/Hogan_Lab_People/Michelle_Meier/scripts/publications/SK_yap1_meox1/scatacseq'
setwd(renv_path)
renv::activate(renv_path)
renv::restore()

# ---- Load ArchR functions ----
source("/team_folders/hogan_lab/Hogan_Lab_Scripts/hogan_lab_bitbucket/Hogan_Lab_Scripts/R_scripts/00_ArchR_functions_R4.2.0.Core.R")

# ---- Load libraries ----
library(ArchR)
library(BSgenome.Drerio.UCSC.danRer11)
library(org.Dr.eg.db)
library(biomaRt)
library(GO.db)
library(httr)
library(ggpubr)

# ---- In script function ----
createGeneAnnotation_adapted <- function(
    genome = NULL,
    TxDb = NULL,
    OrgDb = NULL,
    genes = NULL,
    exons = NULL,
    TSS = NULL
){
  
  
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

make_snakeplot <- function(dap_table_genotype,
                           CT, 
                           gene.list = hippo_targets, 
                           pval = 0.05){
  #closed regions
  dap.table.closed <- dap_table_genotype %>% 
    dplyr::filter(., Status_Group_01 == 'closed' & RawPval < pval) %>%
    mutate(., ranking = rank(-Log2FC, ties.method = 'first'))
  dap.table.open <- dap_table_genotype %>% 
    dplyr::filter(., Status_Group_01 == 'open' & RawPval < pval) %>%
    mutate(ranking = rank(Log2FC, ties.method = 'first'))
  #combine
  dap.table.plot <- rbind(dap.table.closed, dap.table.open)
  
  #label df
  label.peaks <- dap.table.plot %>%
    dplyr::filter(., nearestGene %in% gene.list) 
  
  snakeplot_genotype <- ggplot(dap.table.plot, aes(x = ranking, y = Log2FC)) +
    geom_point(data = dap.table.closed,
               colour =open.close.colours['closed'], 
               size = 2) + 
    geom_point(data = dap.table.open,
               colour =open.close.colours['open'], 
               size = 2) +
    geom_point(data = label.peaks,
               colour ="grey50",
               size = 3) +
    theme_bw() + xlab('Rank') + ylab('Log2(FC)') +
    geom_hline(yintercept = 0) +
    ggrepel::geom_label_repel(data = label.peaks, # Add labels last to appear as the top layer
                              aes(label = nearestGene),
                              force = 2,
                              nudge_y = 1,
                              nudge_x = 0.25) +
    ggtitle(sprintf('%s mutant vs %s WT, highlighted hippo genes (RawPval < 0.05) \nRanked by Log2FC', CT, CT))
  return(snakeplot_genotype)
  
}

# ---- Database ----
txdb <- AnnotationDbi::loadDb(file="/team_folders/hogan_lab/genomes/txdb.sqlite")
seqlevels(txdb) <- paste0("chr", seqlevels(txdb))
seqlevels(txdb) <- paste0("chr", c(seq(1,25)))
macs2_path      <-  "/config/binaries/macs/2.1.1/bin/macs2"


################################################################################
# GET DATA
################################################################################
#general paths: EDIT 18/10/2023 --> COMMENT OUT THIS PATH
# project_path  <- c("/team_folders/hogan_lab/Hogan_Lab_Projects/yap1_project/scatacseq/Dataset_Level_02_Endothelial_cells_resubset")
# EDIT 18/10/2023 --> USE ARCHR PROJECT WITH PEAK SETS CALCULATED USING Phenotype_Genotype_Level_02 
project_path  <- c("/team_folders/hogan_lab/Hogan_Lab_Projects/yap1_project/scatacseq/Dataset_Level_02_Endothelial_cells_resubset_recalcPeaks/")
project_path_genotype <- "/team_folders/hogan_lab/Hogan_Lab_Projects/yap1_project/scatacseq/Dataset_Level_02_Endothelial_cells_resubset_recalcPeaks_genotype" 
project           <- loadArchRProject(project_path, showLogo=FALSE)
figure.dir.out <- '/team_folders/hogan_lab/Hogan_Lab_Projects/yap1_project/scatacseq/Dataset_Level_02_Endothelial_cells_resubset_recalcPeaks/DAP_analysis/'

# ---- Genes ----
hippo_targets <- c(
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
  "bmp4",
  "meox1")

# ---- Colours ----
open.close.colours <- c(
  '#c2a5cf', #open
  '#5aae61' #closed
)
names(open.close.colours) <- c('open', 'closed')

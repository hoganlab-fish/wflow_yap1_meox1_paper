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
source("../functions/00_ArchR_functions_R4.2.0.Core.R")

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


################################################################################
# GET DATA 
################################################################################
# ---- load projects ----
old.path <- '/hogan_lab/Hogan_Lab_People/Oliver_Yu/analysis/prox1a_enhancer/output/scATAC_4dpf_wt_Level_02'
new.path <- '/hogan_lab/Hogan_Lab_Projects/meox1_project/output/scatacseq/Dataset_4dpf_Level_02_proxpaper/'

project <- copy_ArchR_Project(oldProjectPath = old.path, newProjectPath = new.path)
################################################################################
# GET GENOME AND GENE ANNOTATION
################################################################################
macs2_path      <-  "/config/binaries/macs/2.1.1/bin/macs2"

# ---- Database ----
txdb <- AnnotationDbi::loadDb(file="/team_folders/hogan_lab/genomes/txdb.sqlite")
seqlevels(txdb) <- paste0("chr", seqlevels(txdb))
seqlevels(txdb) <- paste0("chr", c(seq(1,25)))

# Load gene and genome annotations
genomeAnnotation <- ArchR::createGenomeAnnotation(
  genome=BSgenome.Drerio.UCSC.danRer11
)

geneAnnotation <- createGeneAnnotation_adapted(
  TxDb=txdb, OrgDb=org.Dr.eg.db
)

################################################################################
# MAKE PSEUDOBULK REPLICATES with Phenotype
################################################################################
project <- ArchR::addGroupCoverages(
  ArchRProj = project, groupBy = "Phenotype", force = TRUE
)

project <- ArchR::saveArchRProject(
  ArchRProj = project, load = TRUE
)

################################################################################
# ADD PEAKSET FOR PROX1   
################################################################################
project <- addReproduciblePeakSet(
  ArchRProj = project,
  groupBy = "Phenotype",
  genomeAnnotation = genomeAnnotation,
  geneAnnotation = geneAnnotation,
  pathToMacs2 = macs2_path,
  genomeSize = 1368780147
)

project <- ArchR::saveArchRProject(
  ArchRProj = project, load = TRUE
)

# ---- Add peak matrix  ----
project <- ArchR::addPeakMatrix(project, force = TRUE)

project <- ArchR::saveArchRProject(
  ArchRProj = project, load = TRUE
)


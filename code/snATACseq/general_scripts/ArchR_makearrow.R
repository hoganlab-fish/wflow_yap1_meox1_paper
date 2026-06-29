#!/usr/bin/Rscript

# run archr pipeline from cellranger-atac 10x fragment data

# ====== Load r environment =====
 warning(
    "This script assumes renv is in working directory"
  )
#renv::activate()
renv::restore()

# ====== Load libraries =====
suppressWarnings(suppressMessages(library(ArchR, quietly=TRUE, verbose=FALSE)))
library(BiocManager)
library(BSgenome.Drerio.UCSC.danRer11)
library(org.Dr.eg.db)
library(chromVARmotifs)
library(JASPAR2020)
library(argparse)

# ====== Helper functions =====
parse_argv <- function(){
  parser <- ArgumentParser(description = "Run ArchR on data." )
  
  # Add command line arguments
  parser$add_argument('-a', '--arrow_path', required = TRUE,  help = 'path to arrow file')
  parser$add_argument('-i', '--infile_paths', required = TRUE,  help = 'Comma seperated paths to fragments (must match order of names in --sample_names)')
  parser$add_argument('-c', '--ncpus', required = TRUE,  help = 'number of cpus',  type="integer")
  parser$add_argument('-s', '--sample_names', required = TRUE,  help = 'Comma seperated names of the data files (must match order in arrow/infile paths)')
  parser$add_argument('-d', '--doublet_scores', required = FALSE,  help = 'filter doublets', default=FALSE)  

  # Parse the command line arguments
  argv <- parser$parse_args()

  # Do work based on the passed arguments
  return(argv)
}

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

build_ref <- function(){
  
  genomeAnnotation <- ArchR::createGenomeAnnotation(
    genome=BSgenome.Drerio.UCSC.danRer11)
  
  # need use custom annotations from ensembl
  txdb <- AnnotationDbi::loadDb(file="/team_folders/hogan_lab/genomes/txdb.sqlite")
  seqlevels(txdb) <- paste0("chr", seqlevels(txdb))
  seqlevels(txdb) <- paste0("chr", c(seq(1,25)))
  geneAnnotation <- createGeneAnnotation_adapted(
    TxDb=txdb, OrgDb=org.Dr.eg.db
  )
  
  return(c(genomeAnnotation=genomeAnnotation, geneAnnotation=geneAnnotation))
}

# ====== Prevent H5AD locking  =====
RHDF5_USE_FILE_LOCKING=FALSE

# ====== ArchR function =====
# from cellranger output --> arrow files
run_archr <- function(){
  sessionInfo()
  args <- parse_argv()
  print(args)

  ## parse arguments
  sample_names <- strsplit(args$sample_names, ",")[[1]]
  arrow_path <- args$arrow_path
  print(arrow_path)
  infile_paths <- strsplit(args$infile_paths, ",")[[1]]
  doublet_scores <- args$doublet_scores

  #adapt arrow_path to get all paths with sample name
  arrow_path_all <- paste0(arrow_path, '/ArrowFiles/'); dir.create(arrow_path_all, showWarnings = F)
  arrow_paths <- paste0(arrow_path_all,sample_names)

  ## Set threads
  ArchR::addArchRThreads(threads=args$ncpus)
  ArchR::getArchRThreads()

  ## build references
  # NOTE: reference assembly and annotations are hardcoded to danRer11
  annotations <- build_ref()


  ## make arrow files
  print("Loading input files:")
  print(infile_paths)
  arrow_paths <- ArchR::createArrowFiles(
    inputFiles = infile_paths,
    sampleNames = sample_names,
    outputNames = arrow_paths,
    minTSS = 4, #Dont set this too high, you can always increase later
    minFrags = 1000,
    addTileMat = TRUE,
    addGeneScoreMat = TRUE,
    geneAnnotation=annotations$geneAnnotation,
    genomeAnnotation=annotations$genomeAnnotation,
    QCDir=paste(arrow_path, "QualityControl", sep="/"), 
    subThreading = F #needed in slurm
    )
 
  ## doublet inference
  if (doublet_scores) {
    doublet_scores <- ArchR::addDoubletScores(
      input = arrow_paths,
      k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
      knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
      LSIMethod = 1, 
      threads = 1 ##needed in slurm
    )
  }

}

run_archr()


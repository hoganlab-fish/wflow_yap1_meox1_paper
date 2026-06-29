#!/usr/bin/Rscript
# run archr pipeline from cellranger-atac 10x fragment data
suppressWarnings(suppressMessages(library(ArchR, quietly=TRUE, verbose=FALSE)))
# library(pheatmap)

build_ref <- function() {
  if (!require("BiocManager")) {
    install.packages("BiocManager")
  }
  if (!require(BSgenome.Drerio.UCSC.danRer11)) {
    BiocManager::install("BSgenome.Drerio.UCSC.danRer11")
    library(BSgenome.Drerio.UCSC.danRer11)
  }
  genomeAnnotation <- ArchR::createGenomeAnnotation(
    genome=BSgenome.Drerio.UCSC.danRer11
  )
  if (!require(TxDb.Drerio.UCSC.danRer11.refGene) & !require(org.Dr.eg.db)) {
    BiocManager::install("org.Dr.eg.db")
    # BiocManager::install("TxDb.Drerio.UCSC.danRer11.refGene")
    # library(TxDb.Drerio.UCSC.danRer11.refGene)
    ## need use custom annotations from ensembl, these are the steps taken:
    # txdb <- GenomicFeatures::makeTxDbFromEnsembl(organism="Danio rerio", release=102)
    # AnnotationDbi::saveDb(txdb, "txdb.sqlite")
    library(org.Dr.eg.db)
  }
  warning(
    "Danio rerio (Zebrafish) sequences as provided by UCSC (danRer11, May 2017",
    "). If your version is different, very bad things will happen."
  )
  ## need use custom annotations from ensembl
  txdb <- AnnotationDbi::loadDb(file="../data/txdb.sqlite")
  # https://github.com/GreenleafLab/ArchR/issues/112
  # geneAnnotation <- ArchR::createGeneAnnotation(
  #   TxDb=TxDb.Drerio.UCSC.danRer11.refGene, OrgDb=org.Dr.eg.db
  # )
  # txdb <- AnnotationDbi::loadDb("../data/txdb.sqlite")
  seqlevels(txdb) <- paste0("chr", seqlevels(txdb))
  seqlevels(txdb) <- paste0("chr", c(seq(1,25)))
  geneAnnotation <- ArchR::createGeneAnnotation(
    TxDb=txdb, OrgDb=org.Dr.eg.db, annoStyle="ENSEMBL"
  )
  if (!require(chromVARmotifs)) {
    devtools::install_github("GreenleafLab/chromVARmotifs")
    library(chromVARmotifs)
  }
  if (!require(JASPAR2020)) {
    BiocManager::install("JASPAR2020")
    library(JASPAR2020)
  }
  return(c(genomeAnnotation=genomeAnnotation, geneAnnotation=geneAnnotation))
}

make_plots <- function(project) {
  # ridge
  p1 <- ArchR::plotGroups(
    ArchRProj = project, groupBy = "Sample", colorBy = "cellColData",
    name = "TSSEnrichment", plotAs = "ridges"
  )
  # violin
  p2 <- ArchR::plotGroups(
    ArchRProj = project, groupBy = "Sample", colorBy = "cellColData",
    name = "TSSEnrichment", plotAs = "violin", alpha = 0.4, addBoxPlot = TRUE
  )
  # ridge log10
  p3 <- ArchR::plotGroups(
    ArchRProj = project, groupBy = "Sample", colorBy = "cellColData",
    name = "log10(nFrags)", plotAs = "ridges"
  )
  # violin log10
  p4 <- ArchR::plotGroups(
    ArchRProj = project, groupBy = "Sample", colorBy = "cellColData",
    name = "log10(nFrags)", plotAs = "violin", alpha = 0.4, addBoxPlot = TRUE
   )
  p5 <- ArchR::plotFragmentSizes(ArchRProj = project)
  p6 <- ArchR::plotTSSEnrichment(ArchRProj = project)
  ArchR::plotPDF(
    p1, p2, p3, p4, p5, p6,
    name="qc_plots.pdf", ArchRProj=project, addDOC=FALSE, width=4, height=4
  )
  return(c(p1=p1, p2=p2, p3=p3, p4=p4, p5=p5, p6=p6))
}

confuse <- function(project) {
  cM <- ArchR::confusionMatrix(paste0(project$Clusters), paste0(project$Sample))
  cM <- cM / Matrix::rowSums(cM)
  p <- pheatmap::pheatmap(
      mat = as.matrix(cM),
      color = paletteContinuous("whiteBlue"),
      border_color = "black"
  )
  return(p)
}

parse_argv <- function() {
  p <- argparser::arg_parser(
    "Run ArchR on data."
  )
  # Add command line arguments
  p <- argparser::add_argument(
    p, "--json", type="character", nargs=1, default=NA,
    help="pass args as json file instead of command line args (overrides args!)"
  )
  p <- argparser::add_argument(
    p, "--arrow_paths", type="character", nargs=Inf, default=NULL,
    help="paths to arrow files (must match order of names in --sample_names)."
  )
  p <- argparser::add_argument(
    p, "--cluster_identities", type="character", nargs=Inf, default=NULL,
    help="cluster identities (user defined, if blank it will be Cluster 1, 2..)."
  )
  p <- argparser::add_argument(
    p, "--doublet_filter", flag=TRUE,
    help="filter doublets."
  )
  p <- argparser::add_argument(
    p, "--footprint_motif", type="character", nargs=Inf, default=NULL,
    help="user defined motifs for annotating plots and making tracks."
  )
  p <- argparser::add_argument(
    p, "--genes_marker", type="character", nargs=Inf, default=NULL,
    help="user defined marker genes for annotating plots and making tracks."
  )
  p <- argparser::add_argument(
    p, "--infile_paths", type="character", nargs=Inf, default=NULL,
    help="paths to fragments (must match order of names in --sample_names)."
  )
  p <- argparser::add_argument(
    p, "--lock_clusters", flag=TRUE,
    help="do not override existing clusters."
  )
  p <- argparser::add_argument(
    p, "--macs2_path", type="character", nargs=Inf, default=NULL,
    help="path to macs2."
  )
  p <- argparser::add_argument(
    p, "--ncpus", help="number of cpus", type="integer", default=1,
  )
  p <- argparser::add_argument(
    p, "--obtain_clusters", type="character", nargs=Inf, default=NULL,
    help="subset clusters based on cluster names (e.g. C1, C2, ...)."
  )
  p <- argparser::add_argument(
    p, "--project_path", type="character", nargs=1, default=FALSE,
    help="path to project"
  )
  p <- argparser::add_argument(
    p, "--sample_names", type="character", nargs=Inf, default=NULL,
    help="names of the data files (must match order in arrow/infile paths)."
  )
  p <- argparser::add_argument(
    p, "--use_existing", type="character", nargs=1, default=NULL,
    help="copies an existing project directory (into project path)."
  )
  # Parse the command line arguments
  argv <- argparser::parse_args(p)

  # Do work based on the passed arguments
  return(argv)
}

run_archr <- function() {
  sessionInfo()
  argv <- parse_argv()

  if (!is.na(argv$json)) {
    print("Json file passed to pipeline, override all other command input!")
    print(argv$json)
    if (!require("rjson")) {
      install.packages("rjson")
    }
    argv <- rjson::fromJSON(file=argv$json)
  }
  print(argv)

  ## parse arguments
  arrow_paths <- argv$arrow_paths
  cluster_id <- argv$cluster_identities
  doublet_filter <- argv$doublet_filter
  infile_paths <- argv$infile_paths
  lock_clusters <- argv$lock_clusters
  macs2_path <- argv$macs2_path
  marker_genes <- toupper(argv$genes_marker)
  motifs <- argv$footprint_motif
  obtain_clusters <- argv$obtain_clusters
  project_path <- argv$project_path
  sample_names <- argv$sample_names
  use_existing <- argv$use_existing

  if (!is.null(arrow_paths)) {names(arrow_paths) <- sample_names}
  if (!is.null(infile_paths)) {names(infile_paths) <- sample_names}

  ArchR::addArchRThreads(threads=argv$ncpus)
  ArchR::getArchRThreads()

  ## build references
  # NOTE: reference assembly and annotations are hardcoded to danRer11
  annotations <- build_ref()

  ## load existing archr project and copy it to project_path
  if (!is.null(use_existing) & !is.null(project_path)) {
    project <- loadArchRProject(use_existing, showLogo=FALSE)
    project <- ArchR::saveArchRProject(
      ArchRProj = project, outputDirectory = project_path, load = TRUE
    )
  }

  ## make arrow files
  if (!is.null(infile_paths) & !is.null(arrow_paths)) {
    print("Loading input files:")
    print(infile_paths)
    arrow_paths <- ArchR::createArrowFiles(
      inputFiles = infile_paths,
      sampleNames = names(infile_paths),
      outputNames = arrow_paths,
      filterTSS = 4, #Dont set this too high, you can always increase later
      filterFrags = 1000,
      addTileMat = TRUE,
      addGeneScoreMat = TRUE,
      geneAnnotation=annotations$geneAnnotation,
      genomeAnnotation=annotations$genomeAnnotation
      # QCDir = paste(project_path, "QualityControl", sep="/")
    )
  } else {
    warning("No input files provided, will look for existing arrow files.")
  }

  # if (!is.null(read_arrows)) {
  #   print("Loading arrow files:")
  #   arrow_paths <- read_arrows
  #   print(arrow_paths)
  # } else {
  #   # TODO: setup logic for loading ArchRProject directly downstream
  #   warning("No arrow files provided, will look for existing ArchRProject.")
  # }

  ## doublet inference
  if (is.null(use_existing) & doublet_filter) {
    doublet_scores <- ArchR::addDoubletScores(
      input = arrow_paths,
      k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
      knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
      LSIMethod = 1
    )
  }

  ## init the project
  if (is.null(use_existing)) {
    print("Creating new project")
    project <- ArchR::ArchRProject(
      ArrowFiles = arrow_paths,
      outputDirectory = project_path,
      copyArrows = TRUE,
      geneAnnotation = annotations$geneAnnotation,
      genomeAnnotation = annotations$genomeAnnotation
    )
  } else {
    print("Use existing project")
    project <- ArchR::loadArchRProject(project_path, showLogo=FALSE)
  }

  # periodically save the project again as checkpoint events
  project <- ArchR::saveArchRProject(
    ArchRProj = project, outputDirectory = project_path, load = TRUE
  )

  ## first batch of qc plots are similar to cellranger output
  plots <- make_plots(project)

  # this uses the previously calculated doublet scores above
  if (is.null(use_existing)) {project <- ArchR::filterDoublets(project)}

  project <- ArchR::saveArchRProject(
    ArchRProj = project, outputDirectory = project_path, load = TRUE
  )

  ## pick clusters (assume you already know what you want)
  if (!is.null(obtain_clusters)) {
    print(paste("Selecting clusters:", obtain_clusters))
    idxSample <- BiocGenerics::which(project$Clusters %in% obtain_clusters)
    cellsSample <- project$cellNames[idxSample]
    project <- project[cellsSample, ]
    project <- ArchR::saveArchRProject(
      ArchRProj=project, outputDirectory=project_path, load=TRUE
    )
  }

  ## reduce dimensions
  # this is non-deterministic, if you get a nice plot make sure to save it!
  if (lock_clusters == FALSE) {
    project <- ArchR::addIterativeLSI(
      ArchRProj = project,
      useMatrix = "TileMatrix",
      name = "IterativeLSI",
      iterations = 2,
      clusterParams = list( #See Seurat::FindClusters
          resolution = c(0.2),
          sampleCells = 10000,
          n.start = 10
      ),
      varFeatures = 25000,
      dimsToUse = 1:30,
      force = TRUE
    )
    project <- ArchR::saveArchRProject(
      ArchRProj=project, outputDirectory=project_path, load=TRUE
    )
  }

  ## do clustering
  # this is deterministic, so all good
  project <- ArchR::addClusters(
    input = project,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    resolution = 0.8,
    force = TRUE
  )
  table(project$Clusters)

  # if confusion matrix matches (cluster identity) the embedding in umap, good
  confusion <- confuse(project)
  ArchR::plotPDF(
    confusion, name="confusion.pdf", ArchRProj=project, addDOC=FALSE,
    width=8, height=8
  )

  ## umaps. tsne also possible but not implemented here
  # self explanatory
  project <- ArchR::addUMAP(
    ArchRProj = project,
    reducedDims = "IterativeLSI",
    name = "UMAP",
    nNeighbors = 30,
    minDist = 0.5,
    metric = "cosine",
    force = TRUE
  )
  umap_1 <- ArchR::plotEmbedding(
    ArchRProj = project, colorBy = "cellColData", name = "Sample",
    embedding = "UMAP",
  )
  umap_2 <- ArchR::plotEmbedding(
    ArchRProj = project, colorBy = "cellColData", name = "Clusters",
    embedding = "UMAP",
  )
  ArchR::plotPDF(
    umap_1, umap_2,
    name="umaps.pdf", ArchRProj=project, addDOC=FALSE, width=8, height=8
  )
  project <- ArchR::saveArchRProject(
    ArchRProj = project, outputDirectory = project_path, load = TRUE
  )

  ## pass in a user defined set of marker genes here
  # get a subset of marker genes per cluster
  markersGS <- ArchR::getMarkerFeatures(
    ArchRProj = project,
    useMatrix = "GeneScoreMatrix",
    groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
  )
  # markerList <- ArchR::getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
  heatmapGS <- ArchR::plotMarkerHeatmap(
    seMarker=markersGS, cutOff="FDR <= 0.01 & Log2FC >= 1.25", transpose=TRUE,
    labelMarkers=marker_genes
  )
  # heatmapGS <- ComplexHeatmap::draw(
  #   heatmapGS, heatmap_legend_side="bot", annotation_legend_side="bot"
  # )
  ArchR::plotPDF(
    heatmapGS, name="heatmap.pdf",
    ArchRProj=project, addDOC=FALSE, width=8, height=8
  )

  # theres an option to make the tracks plot interactive if you run locally
  tracks <- ArchR::plotBrowserTrack(
      ArchRProj = project,
      groupBy = "Clusters",
      geneSymbol = marker_genes,
      upstream = 50000,
      downstream = 50000,
  )
  ArchR::plotPDF(
    tracks, name="tracks.pdf", ArchRProj=project, addDOC=FALSE,
    width=8, height=8
  )

  ## label clusters
  if (!is.null(cluster_id)) {
    project$Clusters_id <- mapLabels(project$Clusters, newLabels=cluster_id)
    umap_cluster_1 <- plotEmbedding(
      project, colorBy="cellColData", name="Clusters_id", embedding = "UMAP"
    )
    umap_cluster_2 <- plotEmbedding(
      project, colorBy="GeneScoreMatrix", name=marker_genes, embedding = "UMAP"
    )
    ArchR::plotPDF(
      umap_cluster_1, umap_cluster_2, name="umap_clusters.pdf",
      ArchRProj = project, addDOC = FALSE, width = 5, height = 5
    )
  } else {
    umap_cluster_1 <- plotEmbedding(
      project, colorBy="cellColData", name="Clusters", embedding = "UMAP"
    )
    umap_cluster_2 <- plotEmbedding(
      project, colorBy="GeneScoreMatrix", name=marker_genes, embedding = "UMAP"
    )
    ArchR::plotPDF(
      umap_cluster_1, umap_cluster_2, name="umap_clusters.pdf",
      ArchRProj = project, addDOC = FALSE,  width = 5, height = 5
    )
  }
  project <- ArchR::saveArchRProject(
    ArchRProj = project, outputDirectory = project_path, load = TRUE
  )

  ## make pseudo bulk replicates
  project <- ArchR::addGroupCoverages(
    ArchRProj = project, groupBy = "Clusters", force = TRUE
  )
  project <- ArchR::saveArchRProject(
    ArchRProj = project, outputDirectory = project_path, load = TRUE
  )

  ## peak calling with macs2
  # est. genomeSize = ./faCount genome.fa : -> 1373471384 - 4691237 = 1368780147
  if (!is.null(macs2_path)) {
    project <- ArchR::addReproduciblePeakSet(
      ArchRProj = project, groupBy = "Clusters", pathToMacs2 = macs2_path,
      genomeSize = 1368780147, force = TRUE
    )
  } else {
    project <- ArchR::addReproduciblePeakSet(
      ArchRProj = project, groupBy = "Clusters", pathToMacs2 = findMacs2(),
      genomeSize = 1368780147, force = TRUE
    )
  }
  project <- ArchR::addPeakMatrix(project, force = TRUE)
  project <- ArchR::saveArchRProject(
    ArchRProj = project, outputDirectory = project_path, load = TRUE
  )

  ## identifying marker peaks
  # TODO: add custom genes
  markersPeaks <- ArchR::getMarkerFeatures(
    ArchRProj = project, useMatrix = "PeakMatrix", groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"), testMethod = "wilcoxon",
    useGroups = NULL, bgdGroups = NULL
  )

  markerList <- ArchR::getMarkers(
    markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1"
  )

  heatmapPeaks <- ArchR::plotMarkerHeatmap(
    seMarker = markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
    transpose = TRUE#, labelMarkers=marker_genes
  )
  # heatmapPeaks <- ArchR::markerHeatmap(
  #   seMarker = markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  #   transpose = TRUE#, labelMarkers=marker_genes
  # )
  # heatmapPeaks <- ComplexHeatmap::draw(
  #   heatmapPeaks, heatmap_legend_side="bot", annotation_legend_side="bot"
  # )
  ArchR::plotPDF(
    heatmapPeaks, name = "peak_marker_heatmap.pdf", width = 8, height = 6,
    ArchRProj = project, addDOC = FALSE
  )

  pma <- ArchR::plotMarkers(
    seMarker = markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 1",
    name = colnames(markersPeaks), plotAs = "MA"
  )
  pv <- ArchR::plotMarkers(
    seMarker = markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 1",
    name = colnames(markersPeaks), plotAs = "Volcano"
  )
  ArchR::plotPDF(
    pma, pv, name = "ma_volcano.pdf", width = 5, height = 5,
    ArchRProj = project, addDOC = FALSE
  )

  # TODO: add custom genes for a specific cluster
  # track <- plotBrowserTrack(
  #   ArchRProj = project, groupBy = "Clusters", geneSymbol = c("GENE"),
  #   features = getMarkers(
  #     markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 1", returnGR = TRUE
  #   )["C1"],
  #   upstream = 50000,
  #   downstream = 50000
  # )
  # plotPDF(
  #   track, name = "tracks_features.pdf", width = 5, height = 5,
  #   ArchRProj = project, addDOC = FALSE
  # )

  ## pairwise comparisons
  # TODO: add custom genes
  markerTest <- ArchR::getMarkerFeatures(
    ArchRProj = project, useMatrix = "PeakMatrix", groupBy = "Clusters",
    testMethod = "wilcoxon", bias = c("TSSEnrichment", "log10(nFrags)"),
    useGroups = NULL, bgdGroups = NULL
  )

  ## motif enrichment
  # https://github.com/GreenleafLab/ArchR/issues/774
  # this block adds the homer motif annotations from chromVARmotifs
  data("homer_pwms")
  homer_motifs <- homer_pwms
  # this block adds the custom prox motif
  prox <- TFBSTools::getMatrixSet(
    x=JASPAR2020,
    opts=list(all_versions=FALSE,species=9606,collection="CORE",matrixtype="PWM")
  )["MA0794.1"]
  # converts motif ids to gene names
  names(prox) <- TFBSTools::name(prox)
  homer_motifs <- c(homer_motifs, prox)
  project <- ArchR::addMotifAnnotations(
    project, motifPWMs=homer_motifs, force=TRUE,
    species="BSgenome.Drerio.UCSC.danRer11"
  )
  # project <- ArchR::addMotifAnnotations(
  #   ArchRProj = project, motifSet = "homer", name = "Motif",
  #   species="BSgenome.Drerio.UCSC.danRer11", force = TRUE
  # )

  project <- ArchR::saveArchRProject(
    ArchRProj = project, outputDirectory = project_path, load = TRUE
  )

  motifsUp <- ArchR::peakAnnoEnrichment(
    seMarker = markerTest,
    ArchRProj = project,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )

  df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
  df <- df[order(df$mlog10Padj, decreasing = TRUE),]
  df$rank <- seq_len(nrow(df))

  ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) +
    geom_point(size = 1) + ggrepel::geom_label_repel(
      data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF),
      size = 1.5, nudge_x = 2, color = "black"
    ) + theme_ArchR() +
    ylab("-log10(P-adj) Motif Enrichment") +
    xlab("Rank Sorted TFs Enriched") +
    scale_color_gradientn(colors = paletteContinuous(set = "comet")
  )
  ggsave(plot=ggUp, filename=paste(project_path, "Plots/tf_up.pdf", sep="/"))

  motifsDo <- peakAnnoEnrichment(
    seMarker = markerTest,
    ArchRProj = project,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC <= -0.5"
  )

  df <- data.frame(TF = rownames(motifsDo), mlog10Padj = assay(motifsDo)[,1])
  df <- df[order(df$mlog10Padj, decreasing = TRUE),]
  df$rank <- seq_len(nrow(df))

  ggDo <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) +
    geom_point(size = 1) +
    ggrepel::geom_label_repel(
      data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF),
      size = 1.5, nudge_x = 2, color = "black"
    ) + theme_ArchR() +
    ylab("-log10(FDR) Motif Enrichment") +
    xlab("Rank Sorted TFs Enriched") +
    scale_color_gradientn(colors = paletteContinuous(set = "comet")
  )
  ggsave(plot=ggDo, filename=paste(project_path, "Plots/tf_do.pdf", sep="/"))

  project <- ArchR::saveArchRProject(
    ArchRProj = project, outputDirectory = project_path, load = TRUE
  )

  enrichMotifs <- ArchR::peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = project,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1"
  )
  heatmapEM <- ArchR::plotEnrichHeatmap(enrichMotifs, n = 7, transpose = TRUE)
  ComplexHeatmap::draw(
    heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot"
  )
  ArchR::plotPDF(
    heatmapEM, name = "motifs_enriched_marker_heatmap.pdf",
    width = 8, height = 6, ArchRProj = project, addDOC = FALSE
  )

  project <- ArchR::saveArchRProject(
    ArchRProj = project, outputDirectory = project_path, load = TRUE
  )

  ## chromVar motif deviations
  project <- addBgdPeaks(project, force = TRUE)
  project <- addDeviationsMatrix(
    ArchRProj = project, peakAnnotation = "Motif", force = TRUE
  )
  plotVarDev <- getVarDeviations(project, name = "MotifMatrix", plot = TRUE)
  plotPDF(
    plotVarDev, name = "variable_motif_deviation_scores.pdf",
    width = 5, height = 5, ArchRProj = project, addDOC = FALSE
  )

  markerMotifs <- getFeatures(
    project, select = paste(marker_genes, collapse="|"),
    useMatrix = "MotifMatrix"
  )
  markerMotifs <- grep("z:", markerMotifs, value = TRUE)

  p <- plotGroups(
    ArchRProj = project,
    groupBy = "Clusters",
    colorBy = "MotifMatrix",
    name = markerMotifs
  )
  plotPDF(
    p, name = "plot_groups_deviations_no_imputation.pdf", width = 5, height = 5,
    ArchRProj = project, addDOC = FALSE
  )

  p <- plotEmbedding(
    ArchRProj = project,
    colorBy = "MotifMatrix",
    name = sort(markerMotifs),
    embedding = "UMAP",
  )
  plotPDF(
    p, name = "umap_groups_deviations_no_imputation.pdf", width = 5, height = 5,
    ArchRProj = project, addDOC = FALSE
  )

  project <- ArchR::saveArchRProject(
    ArchRProj = project, outputDirectory = project_path, load = TRUE
  )

  ## motif footprinting
  if (!is.null(motifs)) {
    motifPositions <- getPositions(project)
    markerMotifs <- unlist(
      lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE))
    )
    markerMotifs <- motifPositions[markerMotifs]

    # we previously grouped cells into pseudobulk
    seFoot <- getFootprints(
      ArchRProj = project, positions = markerMotifs, groupBy = "Clusters"
    )

    plotFootprints(
      seFoot = seFoot, ArchRProj = project, normMethod = "Subtract",
      plotName = "Footprints-Subtract-Bias", addDOC = FALSE, smoothWindow = 5
    )
    plotFootprints(
      seFoot = seFoot, ArchRProj = project, normMethod = "Divide",
      plotName = "Footprints-Divide-Bias", addDOC = FALSE, smoothWindow = 5
    )
  } else {
    warnings("No motif footprinting performed. Provide with --footprint_motif")
  }

  # feature footprinting
  # FIXME: doesnt work on zebrafish genome
  # seTSS <- getFootprints(
  #   ArchRProj = project, positions = GRangesList(TSS = getTSS(project)),
  #   groupBy = "Clusters", flank = 2000
  # )
  #
  # plotFootprints(
  #   seFoot = seTSS, ArchRProj = project, normMethod = "None",
  #   plotName = "TSS-No-Normalization", addDOC = FALSE, flank = 2000,
  #   flankNorm = 100
  # )

  project <- ArchR::saveArchRProject(
    ArchRProj = project, outputDirectory = project_path, load = TRUE
  )

  ## co-accessibility of peaks
  project <- addCoAccessibility(
    ArchRProj = project, reducedDims = "IterativeLSI"
  )

  cA <- getCoAccessibility(
    ArchRProj = project, corCutOff = 0.5, resolution = 1000, returnLoops = TRUE
  )

  p <- plotBrowserTrack(
    ArchRProj = project, groupBy = "Clusters", geneSymbol = marker_genes,
    upstream = 50000, downstream = 50000, loops = getCoAccessibility(project)
  )
  plotPDF(
    plotList = p, name = "Plot-Tracks-Marker-Genes-with-CoAccessibility.pdf",
    ArchRProj = project, addDOC = FALSE, width = 5, height = 5
  )
  rm(p)

  project <- ArchR::saveArchRProject(
    ArchRProj = project, outputDirectory = project_path, load = TRUE
  )

  # peak to gene linkage
  # TODO: this will only work if scRNASeq is available
  # project <- addPeak2GeneLinks(
  #   ArchRProj = project, reducedDims = "IterativeLSI"
  # )
  # p2g <- getPeak2GeneLinks(
  #   ArchRProj = project, corCutOff = 0.5, resolution = 1, returnLoops = FALSE
  # )
  # identify TF regulators
  # seGroupMotif <- getGroupSE(
  #   ArchRProj = project, useMatrix = "MotifMatrix", groupBy = "Clusters"
  # )
  # seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]
  # rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
  #   rowMaxs(assay(seZ) - assay(seZ)[,x])
  # }) %>% Reduce("cbind", .) %>% rowMaxs

  ## trajectory analysis
  # TODO: one day
  # project <- ArchR::saveArchRProject(
  #   ArchRProj = project, outputDirectory = project_path, load = TRUE
  # )

  print(paste0("Memory Size = ", round(object.size(project)/10^6,3), " MB"))
  ArchR::getAvailableMatrices(project)

}

run_archr()

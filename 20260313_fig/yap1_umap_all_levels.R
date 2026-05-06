#!/usr/bin/Rscript
# shared drive: /Revision_analysis/queue_plots/
#   yap1_L01_UMAP_meox1.pdf
#   yap1_L01_UMAP_STRIPPED_meox1.pdf
#   yap1_L02_UMAP_meox1.pdf
#   yap1_L02_UMAP_STRIPPED_meox1.pdf
#   yap1_L03_UMAP_meox1.pdf
#   yap1_L03_UMAP_STRIPPED_meox1.pdf


libraries <- c(
    "chisq.posthoc.test",
    "clustree",
    "ggpubr",
    "ggplot2",   
    "pagoda2",
    "patchwork",
    "qs2",
    "scater",
    "scDblFinder",
    "Seurat",
    "SeuratWrappers",
    "tidyverse"
)

load_packages <- function(packages) {
    suppressPackageStartupMessages({
        cat(paste("Loading", packages, sep=" "), sep="\n")
        invisible(sapply(
            packages, library, character.only = TRUE, quietly = TRUE
            ))
    })
}
invisible(Sys.setenv(OPENBLAS_NUM_THREADS="1"))

make_feature_umap <- function(seurat,
                          feature,
                          width,
                          height,
                          path,
                          name,
                          point_size = 0.5,
                          colours = expression_colours) {
  plot <- FeaturePlot(seurat,
                      features = feature,
                      cols = colours,
                      pt.size = point_size)

  #with everything
  filename <- paste0(name, "_UMAP_", feature, ".pdf")
  ggsave(plot = plot,
         filename = filename,
         path = path,
         device = "pdf",
         height = height+2,
         width = width + 2)

  #stipped
  plot_strip <- plot + NoAxes() + NoLegend() + theme(plot.title = element_blank())
  filename <- paste0(name, "_UMAP_STRIPPED_", feature, ".pdf")
  ggsave(plot = plot_strip,
         filename = filename,
         path = path,
         device = "pdf",
         height = height,
         width = width)

  return(plot)
}

load_packages(libraries)

expression_colours <- c("#d9d9d9", "#40004b")

infile_path <- "../../Saki_data/paper/D10025_yap1_dataset_Level_01_annotated.qs2"
outfile_dir <- "../Revision_analysis/queue_plots/"

data <- qs_read(infile_path)

make_feature_umap(
    seurat=data,
    feature="meox1",
    width=5,
    height=5,
    path=outfile_dir,
    name="yap1_L01",
    point_size = 0.5,
    colours = expression_colours
    )

infile_path <- "../../Saki_data/paper/D10025_yap1_dataset_Level_02_annotated.qs2"
outfile_dir <- "../Revision_analysis/queue_plots/"

data <- qs_read(infile_path)

make_feature_umap(
    seurat=data,
    feature="meox1",
    width=5,
    height=5,
    path=outfile_dir,
    name="yap1_L02",
    point_size = 0.5,
    colours = expression_colours
    )

infile_path <- "../../Saki_data/paper/D10025_yap1_dataset_Level_03_annotated.qs2"
outfile_dir <- "../Revision_analysis/queue_plots/"

data <- qs_read(infile_path)

make_feature_umap(
    seurat=data,
    feature="meox1",
    width=5,
    height=5,
    path=outfile_dir,
    name="yap1_L03",
    point_size = 0.5,
    colours = expression_colours
    )
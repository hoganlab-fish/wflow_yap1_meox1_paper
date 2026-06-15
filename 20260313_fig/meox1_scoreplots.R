#!/usr/bin/Rscript
# meox1_scoreplots_L03_UMAP_ap1_1.pdf 
# meox1_scoreplots_L03_UMAP_hippo_1.pdf
# meox1_scoreplots_L03_UMAP_mapk_1.pdf 
# meox1_scoreplots_L03_UMAP_wnt_1.pdf
# meox1_scoreplots_L03_UMAP_STRIPPED_ap1_1.pdf 
# meox1_scoreplots_L03_UMAP_STRIPPED_hippo_1.pdf 
# meox1_scoreplots_L03_UMAP_STRIPPED_mapk_1.pdf 
# meox1_scoreplots_L03_UMAP_STRIPPED_wnt_1.pdf 

libraries <- c(
    "tidyverse",
    "here",
    "qs2",
    "Seurat",
    "ggplot2",
    "ggpubr",
    "patchwork"
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
load_packages(libraries)

expression_colours <- c("#d9d9d9", "#40004b")

infile_path <- "../../Saki_data/paper/D10051_meox1_dataset_Level_03_annotated.qs2"
outfile_dir <- "../Revision_analysis/queue_plots/tmp/"
data_D10051 <- qs_read(infile_path)
genescores_to_plot <-c("hippo_1", "mapk_1", "wnt_1", "ap1_1")


make_genescore_umap <- function(seurat,
                              feature,
                              width,
                              height,
                              path,
                              name,
                              point_size = 0.5,
                              colours = expression_colours){
  plot <- FeaturePlot(seurat,
                      features = feature,
                      cols = colours,
                      pt.size = point_size,
                      min.cutoff = "q1",
                      max.cutoff =  "q99",
                      order = T)

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


plot_list <- lapply(genescores_to_plot, function(genescore) {
    make_genescore_umap(
        seurat = data_D10051, 
        feature = genescore,
        width = 5, height = 5,
        path = outfile_dir, 
        name = "meox1_scoreplots_L03"
    )
})
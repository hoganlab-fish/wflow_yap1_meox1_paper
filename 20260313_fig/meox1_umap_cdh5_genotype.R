#!/usr/bin/Rscript
# shared drive: /Revision_analysis/queue_plots/
#   meox1_UMAP_cdh5_split_Genotype.pdf
#   meox1_UMAP_STRIPPED_cdh5_split_Genotype.pdf
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
load_packages(libraries)

expression_colours <- c("#d9d9d9", "#40004b")

infile_path <- "../../Saki_data/paper/D10051_meox1_dataset_Level_03_annotated.qs2"
outfile_dir <- "../Revision_analysis/queue_plots/"

# shared drive: /Revision_analysis/queue_plots/
#   meox1_UMAP_cdh5_split_Genotype.pdf
#   meox1_UMAP_STRIPPED_cdh5_split_Genotype.pdf

data <- qs_read(infile_path)
table(data@meta.data$Genotype)
group.by <- "Genotype"

recode_D10051 <- c(
    "wildtype"="wildtype", 
    "meox1_mutant"="meox1 mutant"
)
order_D10051 <- c(
    "wildtype",
    "meox1 mutant"
)
col_map_D10051 <- c(
    "wildtype" = "#dbe2c6", #wildtype
    "meox1 mutant" = "#657c95" #mutant
)

if (!is.null(recode_D10051)) {
    data$Genotype <- recode(data$Genotype, !!!recode_D10051)
}
if (!is.null(order_D10051)) {
    data$Genotype <- factor(data$Genotype, levels=order_D10051)
}
umap_cols <- if (!is.null(col_map_D10051)) col_map_D10051 else NULL
if (is.null(col_map_D10051)) {
    legend_data <- data.frame(
        x=1, y=1,
        Genotype=factor(unique(data$Genotype))
    )
} else {
    legend_data <- data.frame(
        x=1, y=1,
        Genotype=factor(names(col_map_D10051), levels=names(col_map_D10051))
    )
}

make_feature_split_umap <- function(seurat,
                              feature,
                              width,
                              height,
                              path,
                              name,
                              split,
                              point_size = 0.5,
                              colours = expression_colours) {
  plot <- FeaturePlot(seurat,
                      features = feature,
                      cols = colours,
                      split.by = split,
                      pt.size = point_size)

  #with everything
  filename <- paste0(name, "_UMAP_", feature, "_split_", split ,".pdf")
  ggsave(plot = plot,
         filename = filename,
         path = path,
         device = "pdf",
         height = height+2,
         width = width + 2)

  #stipped
  plot_strip <- plot & NoAxes() & NoLegend() & theme(plot.title = element_blank())
  filename <- paste0(name, "_UMAP_STRIPPED_", feature, "_split_", split ,".pdf")
  ggsave(plot = plot_strip,
         filename = filename,
         path = path,
         device = "pdf",
         height = height,
         width = width)

  return(plot)
}

make_feature_split_umap(
    seurat=data,
    feature="cdh5",
    width=10,
    height=5,
    path="../Revision_analysis/queue_plots/", # meox1_umap_cdh5_genotype.pdf
    name="meox1",
    split="Genotype",
    point_size=0.5,
    colours=expression_colours
    )

#!/usr/bin/Rscript
# shared drive: /Revision_analysis/queue_plots/
#   yap1_UMAP_cdh5_split_Genotype.pdf
#   yap1_UMAP_STRIPPED_cdh5_split_Genotype.pdf
#   yap1_UMAP_LEGEND_cdh5_split_Genotype.pdf
#   yap1_UMAP_LEGEND_ONLY_cdh5_split_Genotype.pdf

libraries <- c(
    "chisq.posthoc.test",
    "clustree",
    "cowplot",
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

infile_path <- "../../Saki_data/paper/D10025_yap1_dataset_Level_03_annotated.qs2"
outfile_dir <- "../../wflow_yap1_meox1_paper/Revision_analysis/queue_plots/"

data <- qs_read(infile_path)
table(data@meta.data$Genotype)
group.by <- "Genotype"

recode_D10025 <- c(
    "wildtype"="wildtype", 
    "yap1_mutant"="yap1 mutant"
)
order_D10025 <- c(
    "wildtype",
    "yap1 mutant"
)
col_map_D10025 <- c(
    "wildtype" = "#dbe2c6", #wildtype
    "yap1 mutant" = "#657c95" #mutant
)

if (!is.null(recode_D10025)) {
    data$Genotype <- recode(data$Genotype, !!!recode_D10025)
}
if (!is.null(order_D10025)) {
    data$Genotype <- factor(data$Genotype, levels=order_D10025)
}
umap_cols <- if (!is.null(col_map_D10025)) col_map_D10025 else NULL
if (is.null(col_map_D10025)) {
    legend_data <- data.frame(
        x=1, y=1,
        Genotype=factor(unique(data$Genotype))
    )
} else {
    legend_data <- data.frame(
        x=1, y=1,
        Genotype=factor(names(col_map_D10025), levels=names(col_map_D10025))
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
  # prioritise legend
  plot_legend <- FeaturePlot(seurat,
                      features = feature,
                      cols = colours,
                      split.by = split,
                      pt.size = point_size) & 
                      theme(legend.position="right")
  
  filename <- paste0(name, "_UMAP_LEGEND_", feature, "_split_", split ,".pdf")
  cat(filename, "\n")
  ggsave(plot = plot_legend,
         filename = filename,
         path = path,
         device = "pdf",
         height = height,
         width = width + 2)

  filename <- paste0(name, "_UMAP_LEGEND_ONLY_", feature, "_split_", split ,".pdf")
  cat(filename, "\n")
  legend <- get_legend(plot_legend)
  ggsave(plot = legend, 
        filename = filename, 
        path = path, 
        device = "pdf",
        width = 1, 
        height = 2)

  plot <- FeaturePlot(seurat,
                      features = feature,
                      cols = colours,
                      split.by = split,
                      pt.size = point_size)

  #with everything
  filename <- paste0(name, "_UMAP_", feature, "_split_", split ,".pdf")
  cat(filename, "\n")
  ggsave(plot = plot,
         filename = filename,
         path = path,
         device = "pdf",
         height = height+2,
         width = width + 2)

  #stipped
  plot_strip <- plot & NoAxes() & NoLegend() & theme(plot.title = element_blank())
  filename <- paste0(name, "_UMAP_STRIPPED_", feature, "_split_", split ,".pdf")
  cat(filename, "\n")
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
    path=outfile_dir, # yap1_umap_cdh5_genotype.pdf
    name="yap1",
    split="Genotype",
    point_size=0.5,
    colours=expression_colours
    )

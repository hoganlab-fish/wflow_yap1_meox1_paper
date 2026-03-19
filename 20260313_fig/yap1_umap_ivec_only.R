#!/usr/bin/Rscript
# velocyto redo on L03 yap1 scrnaseq data (with intron retention on cellranger)
# 3 high-resolution plots provided separately for custom composing via Adobe Illustrator
# For each of these plots, legends are provided separately for the same reason
# axis also provided
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

infile_path <- "../../Saki_data/paper/D10025_yap1_dataset_Level_03_annotated.qs2"
outfile_dir <- "../Revision_analysis/queue_plots/"

# shared drive: /Revision_analysis/queue_plots/
#   Yap combined: yap1_umap_ivec_only.pdf
#   Yap combined legend: yap1_umap_ivec_only.legend.pdf

outfile_path <- paste(
    outfile_dir, "yap1_umap_ivec_only.pdf", sep="/"
    )
outfile_path_blank <- paste(
    outfile_dir, "yap1_umap_ivec_only_blank.pdf", sep="/"
    )
outfile_legend <- paste(
    outfile_dir, "yap1_umap_ivec_only.legend.pdf", sep="/"
    )
outfile_obj <- paste(
    outfile_dir, "yap1_umap_ivec_only.D10025.qs2", sep="/"
    )

data <- qs_read(infile_path)
table(data@meta.data$L3_cluster_id)
group.by <- "L3_cluster_id"

recode_D10025 <- c(
    "LEC_01"="LEC", 
    "preLEC_01"="preLEC", 
    "hmVEC_01"="hmsVEC",
    "mVEC_01"="msVEC", 
    "cVEC_01"="cvpVEC", 
    "iVEC_01"="IL_VEC"
)
order_D10025 <- c(
    "LEC", 
    "preLEC", 
    "hmsVEC",
    "msVEC", 
    "cvpVEC", 
    "IL_VEC"
)
col_map_D10025 <- c(
    "LEC"     = "#d9d9d9",
    "preLEC"  = "#d9d9d9",
    "hmsVEC"  = "#d9d9d9",
    "msVEC"   = "#d9d9d9",
    "cvpVEC"  = "#d9d9d9",
    "IL_VEC"  = "#084594"
)                   

if (!is.null(recode_D10025)) {
    data$L3_cluster_id <- recode(data$L3_cluster_id, !!!recode_D10025)
}
if (!is.null(order_D10025)) {
    data$L3_cluster_id <- factor(data$L3_cluster_id, levels=order_D10025)
}
umap_cols <- if (!is.null(col_map_D10025)) col_map_D10025 else NULL
if (is.null(col_map_D10025)) {
    legend_data <- data.frame(
        x=1, y=1,
        L3_cluster_id=factor(unique(data$L3_cluster_id))
    )
} else {
    legend_data <- data.frame(
        x=1, y=1,
        L3_cluster_id=factor(names(col_map_D10025), levels=names(col_map_D10025))
    )
}

gg_legend <- ggplot(legend_data, aes(x=x, y=y, colour=L3_cluster_id)) +
    geom_point(shape=NA) +
    scale_colour_manual(values=col_map_D10025) +
    guides(colour=guide_legend(override.aes=list(shape=16, size=4))) +
    theme_void() +
    theme(
        legend.title=element_blank(),
        legend.background=element_blank(),
        legend.margin=margin(r=200)
    )
ggsave(outfile_legend, gg_legend)

# adapted from wflow_yap1_meox1_paper/code/general_scripts/yap1_meox1_paper_settings.R
plot <- DimPlot(
    data,
    group.by = group.by,
    shuffle = T,
    cols = umap_cols,
    pt.size = 0.5
)

# with everything
ggsave(
    plot = plot,
    filename = outfile_path,
    device = "pdf",
    height = 7,
    width = 7
)

# with less thing
plot_strip <- plot + 
    NoAxes() + 
    NoLegend() + 
    theme(plot.title = element_blank())

ggsave( 
    plot = plot_strip,
    filename = outfile_path_blank,
    device = "pdf",
    height = 5,
    width = 5
)
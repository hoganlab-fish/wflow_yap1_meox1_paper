#!/usr/bin/Rscript
# shared drive: /Revision_analysis/queue_plots/
#   yap1_dotplot_L01.pdf
#   yap1_dotplot_L02.pdf
#   yap1_dotplot_L03.pdf

# please review dotplots and note that:
# - labels are different (but correspond) to names in previous dotplots
# - number of clusters have also changed compared to previous dotplots
# - axes and gaps between dots differ from previous dotplots
# let me know if changes need to be made on my end to any of above
library(ggplot2)
library(Seurat)
library(qs2)

features <- c("meox1")
rotate_x <- theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
expression_colours <- c("#d9d9d9", "#40004b")

#### L01 ####

infile_path <- "../../Saki_data/paper/D10025_yap1_dataset_Level_01_annotated.qs2"
outfile_path <- "../../wflow_yap1_meox1_paper/Revision_analysis/queue_plots/yap1_dotplot_L01.pdf"
level_01 <- qs_read(infile_path)

Idents(level_01) <- level_01@meta.data$L1_celltype

dotplot <- DotPlot(
    object = level_01,
    features = features,
    cluster.idents = F, 
    cols = expression_colours
    ) + coord_flip() + rotate_x + ylab("") + xlab("") +
    guides(
        size = guide_legend(
            title = "Percent Expressed",
            direction = "horizontal"
        ),
        colour = guide_colorbar(
            title = "Average Expression",
            direction = "horizontal",
            barwidth = 10,
            barheight = 1
        )
    ) +
    theme(
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.title = element_text(vjust = 1)
    )

ggsave(
    outfile_path,
    dotplot,
    device="pdf",
    width=12,
    height=4
)

#### L02 ####

infile_path <- "../../Saki_data/paper/D10025_yap1_dataset_Level_02_annotated.qs2"
outfile_path <- "../../wflow_yap1_meox1_paper/Revision_analysis/queue_plots/yap1_dotplot_L02.pdf"
level_02 <- qs_read(infile_path)

Idents(level_02) <- level_02@meta.data$L2_celltype

dotplot <- DotPlot(
    object = level_02,
    features = features,
    cluster.idents = F, 
    cols = expression_colours
    ) + coord_flip() + rotate_x + ylab("") + xlab("") +
    guides(
        size = guide_legend(
            title = "Percent Expressed",
            direction = "horizontal"
        ),
        colour = guide_colorbar(
            title = "Average Expression",
            direction = "horizontal",
            barwidth = 10,
            barheight = 1
        )
    ) +
    theme(
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.title = element_text(vjust = 1)
    )

ggsave(
    outfile_path,
    dotplot,
    device="pdf",
    width=12,
    height=4
)

#### L03 ####

infile_path <- "../../Saki_data/paper/D10025_yap1_dataset_Level_03_annotated.qs2"
outfile_path <- "../../wflow_yap1_meox1_paper/Revision_analysis/queue_plots/yap1_dotplot_L03.pdf"
level_03 <- qs_read(infile_path)

Idents(level_03) <- level_03@meta.data$L3_celltype

dotplot <- DotPlot(
    object = level_03,
    features = features,
    cluster.idents = F, 
    cols = expression_colours
    ) + coord_flip() + rotate_x + ylab("") + xlab("") +
    guides(
        size = guide_legend(
            title = "Percent Expressed",
            direction = "horizontal"
        ),
        colour = guide_colorbar(
            title = "Average Expression",
            direction = "horizontal",
            barwidth = 10,
            barheight = 1
        )
    ) +
    theme(
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.title = element_text(vjust = 1)
    )

ggsave(
    outfile_path,
    dotplot,
    device="pdf",
    width=12,
    height=4
)


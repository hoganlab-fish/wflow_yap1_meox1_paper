#!/usr/bin/Rscript
# Addition of cdh5 genes to ../Revision_analysis/queue_plots/meox1_dataset_Level_03_L3_celltype_genotype_cdh5.pdf

library(ggplot2)
library(Seurat)
library(qs2)

infile_path <- "../../Saki_data/paper/D10051_meox1_dataset_Level_03_annotated.qs2"
outfile_path <- "../Revision_analysis/queue_plots/meox1_dataset_Level_03_L3_celltype_genotype_cdh5.pdf"
# add "cdkn1a", "tp53", "cdkn1bb"
features <- c(
    "cdh5"
    )
rotate_x <- theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# standardised settings across dataset
genotype_colours <- c(
  "#dbe2c6", #wildtype
  "#657c95" #mutant
)
# names(genotype_colours) <- c("WT", "yap1_mutant") -> set name depending on dataset
expression_colours <- c("#d9d9d9", "#40004b")

level_03_celltype_genotype_order <- c(
    "LEC_wildtype",
    "LEC_meox1_mutant",
    "pre_muLEC_wildtype",
    "pre_muLEC_meox1_mutant",
    "hmVEC_wildtype",
    "hmVEC_meox1_mutant",
    "mVEC_wildtype",
    "mVEC_meox1_mutant"
)

level_03 <- qs_read(infile_path)

level_03_relevel <- level_03
level_03_relevel$L3_celltype_genotype <- paste0(
    level_03_relevel$L3_celltype, "_", level_03_relevel$Genotype
    )
level_03_relevel <- SetIdent(
    level_03_relevel, value="L3_celltype_genotype"
    )
Idents(level_03_relevel) <- factor(
    level_03_relevel@active.ident, rev(
        level_03_celltype_genotype_order
        )
    )

dotplot <- DotPlot(
    object = level_03_relevel,
    features = features,
    cluster.idents = F, 
    cols = expression_colours
    ) + rotate_x + ylab("") + xlab("")

ggsave(
    outfile_path,
    dotplot,
    device="pdf",
    width=6,
    height=4
)

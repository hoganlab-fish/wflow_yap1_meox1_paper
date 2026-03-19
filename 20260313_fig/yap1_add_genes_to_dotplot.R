#!/usr/bin/Rscript
# Addition of cdkn1a, tp53 and cdkn1bb genes to yap1_dataset_Level_03_L3_celltype_genotype_cellcyclelec
library(ggplot2)
library(Seurat)
library(qs2)

infile_path <- "../../Saki_data/paper/D10025_yap1_dataset_Level_03_annotated.qs2"
outfile_path <- "../Revision_analysis/queue_plots/yap1_dataset_Level_03_L3_celltype_genotype_cellcyclelec_cdkn1a_tp53_cdkn1bb.pdf"
# add "cdkn1a", "tp53", "cdkn1bb"
features <- c(
    "prox1a", "tbx1", "cdh6", "pcna", "mki67", "cdkn1a", "tp53", "cdkn1bb"
    )
rotate_x <- theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# standardised settings across dataset
genotype_colours <- c(
  "#dbe2c6", #wildtype
  "#657c95" #mutant
)
# names(genotype_colours) <- c("WT", "yap1_mutant") -> set name depending on dataset
expression_colours <- c("#d9d9d9", "#40004b")
cols_cellcycle <- c( '#ffff66', '#cc85ff','#a0d0e0')
names(cols_cellcycle) <- c("G1/G0", "S", "G2M")

level_03_celltype_genotype_order <- c(
    "LEC_wildtype",
    "LEC_yap1_mutant",
    "preLEC_wildtype",
    "preLEC_yap1_mutant",
    "hmVEC_wildtype",
    "hmVEC_yap1_mutant",
    "mVEC_wildtype",
    "mVEC_yap1_mutant",
    "cVEC_wildtype",
    "cVEC_yap1_mutant",
    "iVEC_wildtype",
    "iVEC_yap1_mutant"
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

#!/usr/bin/Rscript
# redo of yap1_dataset_Level_03_ivec_barplot_GOterms so that bar graphs are coloured by enrichR score and y axis is AdjPval -> ii)
library(tidyverse)
library(ggplot2)
library(qs2)
library(Seurat)

outfile_path_original <- "../Revision_analysis/queue_plots/yap1_dataset_Level_03_ivec_barplot_GOterms.OLD.pdf"
outfile_path_requested <- "../Revision_analysis/queue_plots/yap1_dataset_Level_03_ivec_barplot_GOterms.pdf"

outfile_path_sfa <- "../Revision_analysis/queue_plots/yap1_dataset_Level_03_ivec_barplot_GOterms.sfa.pdf"
outfile_path_sba <- "../Revision_analysis/queue_plots/yap1_dataset_Level_03_ivec_barplot_GOterms.sba.pdf"

go_terms <- qs_read("../../Saki_data/analysis/D10025_yap1_dataset/enrichR/D10025_yap1_dataset_markers_L3_celltype_top100_enrichR.qs2")

plot_sizing <- theme_bw() + theme(axis.title = element_text(size = 16),
    axis.text = element_text(size = 14, colour = "black"),
    plot.title = element_text(size = 18),
    legend.text = element_text(size = 14),
    legend.title = element_blank()) 

plot_sizing <- theme_bw() + theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14, colour = "black"),
    plot.title = element_text(size = 18),
    legend.text = element_text(size = 14),
    legend.title = element_blank(),
    # overrides for this plot
    panel.grid = element_blank(),
    legend.position = "top",
    legend.justification = "left"
)

plot_original <- go_terms %>% 
    mutate(., Term = gsub("_WP.*", "", Term)) %>% 
    filter(., cluster == "iVEC" & Adjusted.P.value < 0.05) %>% 
    ggplot(., aes(y = Genes, x = Combined.Score)) +
    geom_bar(stat="identity", fill = "purple4") +
    geom_text(aes(label = Term, x=0), colour = "grey90", size = 6, hjust=-0.01) +
    plot_sizing + xlab("enrichR score") 

plot_requested <- go_terms %>% 
    mutate(Term = gsub("_WP.*", "", Term)) %>% 
    filter(cluster == "iVEC" & Adjusted.P.value < 0.05) %>%
    ggplot(., aes(y = reorder(Genes, -Adjusted.P.value), x = Adjusted.P.value, fill = Combined.Score)) +
    geom_bar(stat="identity") +
    geom_text(aes(label = Term, x=0), colour="green", size=6, hjust=-0.01) +
    geom_text(aes(label = formatC(Adjusted.P.value, digits=3, format="g"),
                  x=Adjusted.P.value,
                  colour="orange"),
                  size=4, hjust=1.1) +
    plot_sizing + theme(legend.title = element_text(size = 14)) +
    scale_fill_viridis_c(name="enrichR score") +
    scale_colour_identity() +
    xlab("AdjPval") + ylab("Genes")

ggsave(
    outfile_path_original,
    plot_original,
    device="pdf",
    width=10,
    height=4
    )

ggsave(
    outfile_path_requested,
    plot_requested,
    device="pdf",
    width=10,
    height=4
    )
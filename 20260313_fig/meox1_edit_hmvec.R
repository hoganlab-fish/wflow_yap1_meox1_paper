#!/usr/bin/Rscript
# v1 20260317: redo of meox1_dataset_Level_03_hmVEC_volcano_pvalue so that y axis is -log (raw Pval) and x axis is log(expression ratio) and DEGs are coloured (using adjusted P val) -> i)
# v2 20260318: do you want me to set the same x and y boundaries as the image you pasted here? -> no its okay, it looks good except could I actually change the colouring so that the "DEGs" are coloured based on rawpval <0.05 rather than Adjpwal <0.05 (old version)
# v3 20260526: hmVEC and mVEC 
# shared drive: /Revision_analysis/queue_plots/
#   meox_dataset_Level_03_hmVEC_volcano_pvalue_RAW.pdf
#   Level_03_DEG_hippo_targets_mutVSwt_L3_celltype_fc0.00_minpct_0.00.csv
#   Level_03_DEG_mutVSwt_L3_celltype_fc0.00_minpct_0.01.csv

library(tidyverse)
library(ggplot2)
library(qs2)
library(Seurat)

infile_path <- "../../Saki_data/paper/D10051_meox1_dataset_Level_03_annotated.qs2"
outfile_path_con <- "../Revision_analysis/queue_plots/meox1_dataset_Level_03_hmVEC_volcano_pvalue.continuous.pdf"
outfile_path_dis <- "../Revision_analysis/queue_plots/meox1_dataset_Level_03_hmVEC_volcano_pvalue_RAW.pdf"

hippo_gene_degs <- read_csv(
    "../../Saki_data/analysis/D10051_meox1_dataset/DEGs/Level_03_DEG_hippo_targets_mutVSwt_L3_celltype_fc0.00_minpct_0.00.csv", 
    show_col_types = F
    )

all_degs <- read_csv(
    "../../Saki_data/analysis/D10051_meox1_dataset/DEGs/Level_03_DEG_mutVSwt_L3_celltype_fc0.00_minpct_0.01.csv", 
    show_col_types = F
    )

hippo_targets <- c(
    "yap1",
    "wwtr1",
    "ccn2a",
    "ccn2b",
    "amotl2a",
    "amotl2b",
    "ccn1",
    "ccn1l2",
    "cavin1b",
    "cavin2a",
    "cavin2b",
    "bmp4",
    "meox1"
    )

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

# add expression ratios
all_degs$ExpnRatio <- all_degs$pct.1 / all_degs$pct.2
all_degs$LogExpnRatio <- log10(all_degs$ExpnRatio)

hmvec_data <- all_degs %>% 
    filter(group == "hmVEC") %>%
    mutate(sig = ifelse(p_val < 0.05, "DEG (RawPVal<0.05)", "non-DEG"))

hippo_mark <- all_degs %>% 
    filter(group == "hmVEC" & Gene %in% hippo_targets) %>%
    mutate(sig = ifelse(p_val < 0.05, "DEG (RawPVal<0.05)", "non-DEG"))

# colour="#525252"
plot_discrete <- hmvec_data %>%
    ggplot(., aes(x = LogExpnRatio, y = -log10(p_val), colour = sig)) +
    geom_point(size=1) +
    geom_point(data = hippo_mark, colour="red", size=2, show.legend=FALSE) +
    ggrepel::geom_text_repel(
        data = hippo_mark, 
        aes(label=Gene, fontface=ifelse(Gene=="meox1", "bold.italic", "italic")), 
        size=4, colour="#ff8a04", show.legend=FALSE
    ) +
    scale_colour_manual(
        values = c(
            "known Hippo target gene" = "red",
            "DEG (RawPVal<0.05)" = "#525252",
            "non-DEG" = "#d9d9d9"
        ),
        breaks = c("known Hippo target gene", "DEG (RawPVal<0.05)", "non-DEG")
    ) +
    geom_point(data=data.frame(LogExpnRatio=NA_real_, p_val=NA_real_, 
                               sig="known Hippo target gene"), aes(colour=sig)) +
    guides(colour = guide_legend(override.aes = list(shape=16, size=5), ncol=1)) +
    geom_text(
        data = data.frame(x=-Inf, y=Inf),
        aes(x=x, y=y, label="DEG analysis"),
        hjust=-0.1, vjust=1.4, size=4, colour="black",
        inherit.aes=FALSE
    ) +
    geom_text(
        data = data.frame(x=-Inf, y=Inf),
        aes(x=x, y=y, label="hmsVEC Mutant vs. WT"),
        hjust=-0.06, vjust=3.0, size=4, colour="black",
        inherit.aes=FALSE
    ) +
    labs(y="-log10(RawPval)", x="log10(ExpnRatio)") +
    theme_bw() +
    theme(
        axis.title = element_text(size=14),
        axis.text = element_text(size=12, colour="black"),
        legend.text = element_text(size=12),
        legend.title = element_blank(),
        legend.position = "top",
        legend.justification = "left",
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.grid = element_blank(),
        legend.spacing.y = unit(0.0, "cm"),
        legend.key.height = unit(0.3, "cm"),
        legend.margin = margin(0, 0, -2, 0),
        aspect.ratio = 1
    )

ggsave(
    outfile_path_dis,
    plot_discrete,
    device="pdf",
    width=5,
    height=5
    )
---
title: "Plot genes of interest for meox1 L3"
subtitle: "Kobayashi et al. 2036"
author: "Tyrone Chen"
date: 'August 05, 2026'
output:
  html_document:
    code_folding: hide
    df_print: paged
    highlight: textmate
    keep_md: TRUE
    number_sections: TRUE
    theme: flatly
    toc: TRUE
    toc_float: TRUE
    toc_title: merged samples
editor_options:
  chunk_output_type: console
params:
  sample_index: 1
  sample_name: "default"
---

# Genes of interest

## Important info

We want to compare ontology across different groups and subgroups of wild type and mutant cells.

- To minimise information overload, a brief description is embedded together with each figure so they can be interpreted in context.

### Quickstart instructions for getting the info you want:

- Select your `GO` analysis of choice in the sidebar (ie 1.5 onwards)
- Suggest starting with the `ORA` plots as the core component.
  - ORA Barplots, dotplots and enrichment map show the same information in different ways
  - ORA enrichment map and category-gene network show higher level pathway information
- FGSEA is plotted as a supplementary complement to ORA and are not always available
  - While useful, they can be safely ignored if absent

## Initial setup

For developers.


```{.r .fold-hide}
# !/usr/bin/Rscript
library(ggplot2)
library(ggpubr)
library(patchwork)
library(qs2)
```

```
## qs2 0.2.1
```

```{.r .fold-hide}
library(Seurat)
```

```
## Loading required package: SeuratObject
```

```
## Loading required package: sp
```

```
## 
## Attaching package: 'SeuratObject'
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, t
```

```{.r .fold-hide}
library(SeuratObject)

outfile_dir <- "../output/figure_extended/dotplots/"

singles_meox_path <- "../../Saki_data/analysis/D10051_meox1_dataset/DEGs/Level_03_DEG_mutVSwt_L3_celltype_fc0.00_minpct_0.01.csv"
grouped_meox_path <- "../output/figure_extended/dge_confects/meox1_Level_03_DEG_mutVSwt_L3_celltype_merged_VEC_merged_LEC_fc0.00_minpct_0.01.csv"

infile_path <- "../../Saki_data/paper/D10051_meox1_dataset_Level_03_annotated.qs2"
data <- qs_read(infile_path)

# make sure wt always on left
data$Genotype <- factor(
    data$Genotype,
    levels = c("wildtype", "meox1_mutant")
)

custom_groups <- list(
    "hmVEC__mVEC" = c("hmVEC", "mVEC"),
    "LEC__pre_muLEC" = c("LEC", "pre_muLEC")
)
```

## Wrapper for downstream functions

For developers. We will reuse this a lot so we make a wrapper function for saving plots.


```{.r .fold-hide}
save_plot_custom <- function(plt, plot_name, width = 8, height = 12) {
    pdf_path <- paste0(outfile_dir, "/", plot_name, ".pdf")
    png_path <- paste0(outfile_dir, "/", plot_name, ".png")
    ggsave(filename = pdf_path, plot = plt, width = width, height = height)
    ggsave(filename = png_path, plot = plt, width = width, height = height)
    cat("*PDF saved to:*", paste0("`", pdf_path, "`"), "\n\n")
}
```

## Plots

For developers. Load in data.


```{.r .fold-hide}
features <- c("cdh5", "cdh6", "flt4", "kdrl", "mafb", "prox1a", "tbx1")
expression_colours <- c("#d9d9d9", "#40004b")

data$celltype_genotype <- paste0(data$L3_celltype, " ", data$Genotype)
data$celltype_genotype <- gsub("pre_muLEC", "pre muLEC", data$celltype_genotype)
data$celltype_genotype <- gsub("meox1_mutant", "Mutant", data$celltype_genotype)
data$celltype_genotype <- gsub("wildtype", "WT", data$celltype_genotype)

plt <- DotPlot(object = data, features = features, cluster.idents = FALSE,
        group.by = "celltype_genotype", cols = expression_colours) +
        Seurat::RotatedAxis() +
        guides(color = guide_colorbar(title = "Expression"))

save_plot_custom(plt, "meox1_L3_dotplot_of_goi")
```

```
## *PDF saved to:* `../output/figure_extended/dotplots//meox1_L3_dotplot_of_goi.pdf`
```

## For developers

Sample run command:


``` bash
Rscript -e "
rmarkdown::render(
  'meox1_dotplot_goi.Rmd',
  output_file = './meox1_dotplot_goi.html'
)
"
```

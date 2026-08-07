---
title: "yap1 UMAP feature plots (cdh5 by genotype, meox1 across levels, iVEC highlight)"
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

# README

Three UMAP panels on the yap1 dataset:

1. `cdh5` expression split by genotype (mirrors `meox1_umap_feature_plots.Rmd`'s cdh5 section).
2. `meox1` expression across all three annotation levels (L01/L02/L03).
3. iVEC-only highlight UMAP (all other clusters greyed out), used to call out the interleukin-responsive VEC
   subcluster specifically.

## Outputs

All under `../output/figure_extended/umap/`:

- `yap1_UMAP_cdh5_split_Genotype.pdf`, `yap1_UMAP_STRIPPED_cdh5_split_Genotype.pdf`,
  `yap1_UMAP_LEGEND_cdh5_split_Genotype.pdf`, `yap1_UMAP_LEGEND_ONLY_cdh5_split_Genotype.pdf`
- `yap1_L01_UMAP_meox1.pdf`, `yap1_L01_UMAP_STRIPPED_meox1.pdf` (and `_L02_`, `_L03_` equivalents)
- `yap1_umap_ivec_only.pdf`, `yap1_umap_ivec_only_blank.pdf`, `yap1_umap_ivec_only.legend.pdf`

## cdh5 UMAP split by genotype

Source: `yap1_umap_cdh5_genotype.R`.


```{.r .fold-hide}
library(cowplot)
library(ggplot2)
library(patchwork)
```

```
## 
## Attaching package: 'patchwork'
```

```
## The following object is masked from 'package:cowplot':
## 
##     align_plots
```

```{.r .fold-hide}
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
library(SeuratWrappers)
```

```
## Warning: package 'SeuratWrappers' was built under R version 4.5.2
```

```{.r .fold-hide}
library(tidyverse)
```

```
## ── Attaching core tidyverse packages ──────────────────────────────────────────────────────────────────────────────── tidyverse 2.0.0 ──
## ✔ dplyr     1.2.1     ✔ readr     2.2.0
## ✔ forcats   1.0.1     ✔ stringr   1.6.0
## ✔ lubridate 1.9.5     ✔ tibble    3.3.1
## ✔ purrr     1.2.2     ✔ tidyr     1.3.2
```

```
## ── Conflicts ────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ dplyr::filter()    masks stats::filter()
## ✖ dplyr::lag()       masks stats::lag()
## ✖ lubridate::stamp() masks cowplot::stamp()
## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors
```

```{.r .fold-hide}
outfile_dir <- "../output/figure_extended/umap/"
expression_colours <- c("#d9d9d9", "#40004b")

data_cdh5 <- qs_read("../../Saki_data/paper/D10025_yap1_dataset_Level_03_annotated.qs2")
table(data_cdh5@meta.data$Genotype)
```

```
## 
##    wildtype yap1_mutant 
##        3895        2662
```

```{.r .fold-hide}
recode_D10025 <- c(
    "wildtype"="wildtype",
    "yap1_mutant"="yap1 mutant"
)
order_D10025 <- c(
    "wildtype",
    "yap1 mutant"
)

data_cdh5$Genotype <- recode(data_cdh5$Genotype, !!!recode_D10025)
data_cdh5$Genotype <- factor(data_cdh5$Genotype, levels=order_D10025)

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

plt <- make_feature_split_umap(
    seurat=data_cdh5,
    feature="cdh5",
    width=10,
    height=5,
    path=outfile_dir,
    name="yap1",
    split="Genotype",
    point_size=0.5,
    colours=expression_colours
    )
```

```
## yap1_UMAP_LEGEND_cdh5_split_Genotype.pdf
```

```
## yap1_UMAP_LEGEND_ONLY_cdh5_split_Genotype.pdf 
## yap1_UMAP_cdh5_split_Genotype.pdf
```

```
## yap1_UMAP_STRIPPED_cdh5_split_Genotype.pdf
```

```{.r .fold-hide}
print(plt)
```

![](yap1_umap_feature_plots_files/figure-html/cdh5_umap-1.png)<!-- -->

## meox1 across all annotation levels

Source: `yap1_umap_all_levels.R`.


```{.r .fold-hide}
make_feature_umap <- function(seurat,
                          feature,
                          width,
                          height,
                          path,
                          name,
                          point_size = 0.5,
                          colours = expression_colours) {
  plot <- FeaturePlot(seurat,
                      features = feature,
                      cols = colours,
                      pt.size = point_size)

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

data_L01 <- qs_read("../../Saki_data/paper/D10025_yap1_dataset_Level_01_annotated.qs2")
plt_L01 <- make_feature_umap(
    seurat=data_L01, feature="meox1", width=5, height=5,
    path=outfile_dir, name="yap1_L01", point_size = 0.5, colours = expression_colours
    )
print(plt_L01)
```

![](yap1_umap_feature_plots_files/figure-html/all_levels_umap-1.png)<!-- -->

```{.r .fold-hide}
data_L02 <- qs_read("../../Saki_data/paper/D10025_yap1_dataset_Level_02_annotated.qs2")
plt_L02 <- make_feature_umap(
    seurat=data_L02, feature="meox1", width=5, height=5,
    path=outfile_dir, name="yap1_L02", point_size = 0.5, colours = expression_colours
    )
print(plt_L02)
```

![](yap1_umap_feature_plots_files/figure-html/all_levels_umap-2.png)<!-- -->

```{.r .fold-hide}
data_L03 <- qs_read("../../Saki_data/paper/D10025_yap1_dataset_Level_03_annotated.qs2")
plt_L03 <- make_feature_umap(
    seurat=data_L03, feature="meox1", width=5, height=5,
    path=outfile_dir, name="yap1_L03", point_size = 0.5, colours = expression_colours
    )
print(plt_L03)
```

![](yap1_umap_feature_plots_files/figure-html/all_levels_umap-3.png)<!-- -->

## iVEC-only highlight

Source: `yap1_umap_ivec_only.R`. All clusters recoloured to grey except the interleukin-responsive VEC
subcluster (`iVEC` → recoded `IL_VEC`), highlighted in dark blue.


```{.r .fold-hide}
data_ivec <- qs_read("../../Saki_data/paper/D10025_yap1_dataset_Level_03_annotated.qs2")
table(data_ivec@meta.data$L3_cluster_id)
```

```
## 
## preLEC_01   mVEC_01  hmVEC_01   cVEC_01    LEC_01   iVEC_01 
##      2461      1385      1209       763       656        83
```

```{.r .fold-hide}
group.by <- "L3_cluster_id"

outfile_legend <- paste(outfile_dir, "yap1_umap_ivec_only.legend.pdf", sep="/")

recode_ivec <- c(
    "LEC_01"="LEC",
    "preLEC_01"="preLEC",
    "hmVEC_01"="hmsVEC",
    "mVEC_01"="msVEC",
    "cVEC_01"="cvpVEC",
    "iVEC_01"="IL_VEC"
)
order_ivec <- c(
    "LEC",
    "preLEC",
    "hmsVEC",
    "msVEC",
    "cvpVEC",
    "IL_VEC"
)
col_map_ivec <- c(
    "LEC"     = "#d9d9d9",
    "preLEC"  = "#d9d9d9",
    "hmsVEC"  = "#d9d9d9",
    "msVEC"   = "#d9d9d9",
    "cvpVEC"  = "#d9d9d9",
    "IL_VEC"  = "#084594"
)

data_ivec$L3_cluster_id <- recode(data_ivec$L3_cluster_id, !!!recode_ivec)
data_ivec$L3_cluster_id <- factor(data_ivec$L3_cluster_id, levels=order_ivec)

legend_data <- data.frame(
    x=1, y=1,
    L3_cluster_id=factor(names(col_map_ivec), levels=names(col_map_ivec))
)

gg_legend <- ggplot(legend_data, aes(x=x, y=y, colour=L3_cluster_id)) +
    geom_point(shape=NA) +
    scale_colour_manual(values=col_map_ivec) +
    guides(colour=guide_legend(override.aes=list(shape=16, size=4))) +
    theme_void() +
    theme(
        legend.title=element_blank(),
        legend.background=element_blank(),
        legend.margin=margin(r=200)
    )
ggsave(outfile_legend, gg_legend)
```

```
## Saving 7 x 7 in image
```

```
## Warning: Removed 6 rows containing missing values or values outside the scale
## range (`geom_point()`).
```

```{.r .fold-hide}
print(gg_legend)
```

```
## Warning: Removed 6 rows containing missing values or values outside the scale
## range (`geom_point()`).
```

![](yap1_umap_feature_plots_files/figure-html/ivec_only-1.png)<!-- -->

```{.r .fold-hide}
# adapted from wflow_yap1_meox1_paper/code/general_scripts/yap1_meox1_paper_settings.R
plot_ivec <- DimPlot(
    data_ivec,
    group.by = group.by,
    shuffle = T,
    cols = col_map_ivec,
    pt.size = 0.5
)
print(plot_ivec)
```

![](yap1_umap_feature_plots_files/figure-html/ivec_only-2.png)<!-- -->

```{.r .fold-hide}
ggsave(
    plot = plot_ivec,
    filename = paste0(outfile_dir, "yap1_umap_ivec_only.pdf"),
    device = "pdf",
    height = 7,
    width = 7
)

plot_ivec_strip <- plot_ivec +
    NoAxes() +
    NoLegend() +
    theme(plot.title = element_blank())

ggsave(
    plot = plot_ivec_strip,
    filename = paste0(outfile_dir, "yap1_umap_ivec_only_blank.pdf"),
    device = "pdf",
    height = 5,
    width = 5
)
```

## For developers

Sample run command:


``` bash
Rscript -e "
rmarkdown::render(
  'yap1_umap_feature_plots.Rmd',
  output_file = './yap1_umap_feature_plots.html'
)
"
```

---
title: "Correlation plots for EMBO LEC VEC correlations"
subtitle: "Kobayashi et al. 2026"
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

# Init


``` r
library(ggplot2)
library(ggrepel)
library(patchwork)
library(qs2)
```

```
## qs2 0.2.1
```

``` r
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

``` r
library(SeuratObject)
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
## ✖ dplyr::filter() masks stats::filter()
## ✖ dplyr::lag()    masks stats::lag()
## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors
```

``` r
save_plot <- function(plt, plot_name, width = 10, height = 8) {
    pdf_path <- paste0(outfile_dir, "/", plot_name, ".pdf")
    png_path <- paste0(outfile_dir, "/", plot_name, ".png")
    ggsave(filename = pdf_path, plot = plt, width = width, height = height)
    ggsave(filename = png_path, plot = plt, width = width, height = height)
    cat("*PDF saved to:*", paste0("`", pdf_path, "`"), "\n\n")
}

infile_path <- "../../Saki_data/paper/D10051_meox1_dataset_Level_03_annotated.qs2"
outfile_dir <- "../output/figure_extended/embo_scores/"
data <- qs_read(infile_path)

data$Genotype <- gsub("meox1_mutant", "meox1_mutant", data$Genotype)
data$Genotype <- gsub("wildtype", "WT", data$Genotype)
data$Genotype <- factor(
    data$Genotype,
    levels = c("WT", "meox1_mutant")
)

all_degs <- "../output/figure_extended/dge_confects/meox1_Level_03_DEG_mutVSwt_L3_celltype_merged_VEC_merged_LEC_fc0.00_minpct_0.01.csv"
all_degs <- read_csv(all_degs, show_col_types = F)

lec_vec_labels <- "../../Saki_data/Revision_analysis/data/genesets/EMBO_L2_VECLEC_markers_EMBOmethod_3_5dpf_genesforlabels_only.csv"
lec_vec_labels <- read_csv(lec_vec_labels)
```

```
## Rows: 2904 Columns: 2
## ── Column specification ────────────────────────────────────────────────────────────────────────────────────────────────────────────────
## Delimiter: ","
## chr (2): Gene, direction
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

``` r
# calculate the scores
embo_data <- "../../Saki_data/paper/EMBO_JOURNAL_Level_02_seuratObject_LM.RDS"
Level_02_seuratObject <- readRDS(embo_data)

genes_labels <- "../../Saki_data/Revision_analysis/data/genesets/EMBO_L2_VECLEC_markers_EMBOmethod_3_5dpf_genesforlabels_only.csv"
genes_labels <- read_csv(genes_labels)
```

```
## Rows: 2904 Columns: 2
## ── Column specification ────────────────────────────────────────────────────────────────────────────────────────────────────────────────
## Delimiter: ","
## chr (2): Gene, direction
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

``` r
genes_scoring <- "../../Saki_data/Revision_analysis/data/genesets/EMBO_L2_VECLEC_markers_EMBOmethod_3_5dpf_genesforscoring_only.csv"
genes_everything <- "../../Saki_data/Revision_analysis/data/genesets/EMBO_L2_VECLEC_markers_EMBOmethod_3_5dpf.csv"

# get markers based on the analysis performed for https://link.springer.com/article/10.15252/embj.2022112590
Level_02_seuratObject <- SetIdent(Level_02_seuratObject, value = 'L2_Phenotype_Stage')

# genes of interest and concordance
genes_of_interest <- genes_labels %>% select(Gene) %>% c()
genes_of_interest <- genes_of_interest[[1]]

data$celltype_genotype <- paste(
    data$L3_celltype, data$Genotype, sep = "_"
    )
data$celltype_broad <- case_when(
    data$L3_celltype %in% c("LEC", "pre_muLEC")   ~ "LEC_combined",
    data$L3_celltype %in% c("hmVEC", "mVEC")      ~ "VEC_combined",
    TRUE ~ as.character(data$L3_celltype)
)
data$celltype_broad_genotype <- paste(
    data$celltype_broad, data$Genotype, sep = "_"
)

# LEC slice only
data$celltype_broad_lec_only <- case_when(
    data$L3_celltype %in% c("hmVEC", "mVEC")      ~ "VEC_combined",
    TRUE ~ as.character(data$L3_celltype)
)
data@meta.data$celltype_broad_lec_only <- paste(
    data$celltype_broad_lec_only, data$Genotype, sep="_"
    )

table(data$celltype_genotype)
```

```
## 
##     hmVEC_meox1_mutant               hmVEC_WT       LEC_meox1_mutant 
##                    168                    337                    231 
##                 LEC_WT      mVEC_meox1_mutant                mVEC_WT 
##                    285                    796                    588 
## pre_muLEC_meox1_mutant           pre_muLEC_WT 
##                     93                     77
```

``` r
table(data$celltype_broad)
```

```
## 
## LEC_combined VEC_combined 
##          686         1889
```

``` r
table(data$celltype_broad_genotype)
```

```
## 
## LEC_combined_meox1_mutant           LEC_combined_WT VEC_combined_meox1_mutant 
##                       324                       362                       964 
##           VEC_combined_WT 
##                       925
```

``` r
table(data$celltype_broad_lec_only)
```

```
## 
##          LEC_meox1_mutant                    LEC_WT    pre_muLEC_meox1_mutant 
##                       231                       285                        93 
##              pre_muLEC_WT VEC_combined_meox1_mutant           VEC_combined_WT 
##                        77                       964                       925
```

``` r
custom_groups <- list(
    "hmVEC__mVEC" = c("hmVEC", "mVEC"),
    "LEC__pre_muLEC" = c("LEC", "pre_muLEC")
)
celltype_broad_genotype <- list(
    "VEC_combined_WT",
    "VEC_combined_meox1_mutant",
    "LEC_combined_WT",
    "LEC_combined_meox1_mutant"
)
```

# All LEC vs All VEC
## Perform differential expression analysis.

All LECs vs All VECs


``` r
# generate additional degs per condition
## WT LEC vs WT VEC // Mut VEC vs Mut LEC
### WT LEC vs WT VEC
Idents(data) <- "celltype_broad_genotype"
data_wt_lec_vec <- FindMarkers(object = data,
                group.by = "celltype_broad_genotype",
                ident.1 = "LEC_combined_WT",
                ident.2 = "VEC_combined_WT",
                min.pct = 0.01,
                assay = "RNA",
                logfc.threshold = 0) %>%
    rownames_to_column("Gene") %>%
    mutate(group = "lec_vec_wt")

### WT VEC vs WT LEC
Idents(data) <- "celltype_broad_genotype"
data_wt_vec_lec <- FindMarkers(object = data,
                group.by = "celltype_broad_genotype",
                ident.1 = "VEC_combined_WT",
                ident.2 = "LEC_combined_WT",
                min.pct = 0.01,
                assay = "RNA",
                logfc.threshold = 0) %>%
    rownames_to_column("Gene") %>%
    mutate(group = "vec_lec_wt")

### Mut VEC vs Mut LEC
Idents(data) <- "celltype_broad_genotype"
data_mut <- FindMarkers(object = data,
                group.by = "celltype_broad_genotype",
                ident.1 = "LEC_combined_meox1_mutant",
                ident.2 = "VEC_combined_meox1_mutant",
                min.pct = 0.01,
                assay = "RNA",
                logfc.threshold = 0) %>%
    rownames_to_column("Gene") %>%
    mutate(group = "lec_vec_mut")

# prep data for vis
data_wt_lec_vec <- data_wt_lec_vec %>%
    select(Gene, avg_log2FC_WT = avg_log2FC)
data_wt_vec_lec <- data_wt_vec_lec %>%
    select(Gene, avg_log2FC_WT = avg_log2FC)
data_mut <- data_mut %>%
    select(Gene, avg_log2FC_meox1_mutant = avg_log2FC)

wt_mut <- inner_join(data_wt_lec_vec, data_mut, by = "Gene")
sum(wt_mut$Gene %in% genes_labels$Gene)
```

```
## [1] 2839
```

## Plot

These are very information-dense plots. How to interpret:

- Individual dots correspond to genes.
- Vertical axis shows _how much more/less is this gene expressed in LEC vs VEC WT_.
- Horizontal axis shows _how much more/less is this gene expressed in LEC vs VEC meox1 mutant_.
- Genes falling in the _top right and bottom left quadrants show concordance_ across both WT and mutant.
- Genes falling in the _top left and bottom right quadrants show discordance_ across both WT and mutant.

Examples for reference:
- top right: genes _up_ in the mutant and wildtype
- bottom left: genes _down_ in the mutant and wildtype
- top left: genes _down_ in mutant but up in wildtype
- bottom right: genes _up_ in the mutant but down in wildtype

_Suggested legend:_ Scatterplot shows the correlation between the average log fold changes of two differential gene expression analyses: [LECs vs VECs wildtype], [LECs vs VECs meox1 mutants]. Individual data points are genes and coloured by a curated list of LEC and VEC markers. The vertical axis answers the question: _how much more/less is this gene expressed in LEC vs VEC WT_, while the horizontal axis shows _how much more/less is this gene expressed in LEC vs VEC meox1 mutants_. Genes falling in the _top right and bottom left quadrants show concordance_ across the wildtype and meox1 mutants, while genes falling in the _top left and bottom right quadrants show discordance_ across both the wildtype and meox1 mutants.


### Plot concordant and discordant markers


``` r
wt_mut <- wt_mut %>%
    mutate(
        highlight = Gene %in% genes_of_interest,
        category = case_when(
            !highlight ~ "Other",
            avg_log2FC_WT > 0 & avg_log2FC_meox1_mutant > 0 ~ "Concordant (up)",
            avg_log2FC_WT < 0 & avg_log2FC_meox1_mutant < 0 ~ "Concordant (down)",
            avg_log2FC_WT * avg_log2FC_meox1_mutant < 0     ~ "Discordant",
            TRUE ~ "Other"  # catches exact-zero edge cases
        )
    )

table(wt_mut$category)
# Concordant (down)   Concordant (up)        Discordant             Other
#              1613               756               470             11267

# colour and plot
category_colours <- c(
    "Concordant (up)"   = "green",   # e.g. red - both conditions upregulated
    "Concordant (down)" = "blue",   # e.g. blue - both conditions downregulated
    "Discordant"        = "purple",   # e.g. orange - opposite direction between conditions
    "Other"             = "grey80"     # background genes, not in your list of interest
)

plt <- ggplot(wt_mut, aes(x = avg_log2FC_meox1_mutant, y = avg_log2FC_WT, color = category)) +
    geom_point(data = filter(wt_mut, !highlight), size = 1, alpha = 0.4) +   # background points first
    geom_point(data = filter(wt_mut, highlight), size = 2.5) +               # highlighted genes on top
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey60") +  # y=x reference line
    scale_color_manual(values = category_colours) +
    labs(
        x = "Av. log2FC (LEC meox1 mutant vs VEC meox1 mutant)",
        y = "Av. log2FC (LEC WT vs VEC WT)",
        color = "Gene category",
        title = "Concordance of DE genes across conditions"
    ) +
    theme_classic()

print(plt)

save_plot(plt, "LECVEC_WT_vs_LECVEC_mut_direction")
```

### Plot LEC VEC scores


``` r
### plot by marker

wt_mut <- wt_mut %>%
    left_join(genes_labels, by = "Gene") %>%
    mutate(
        category = case_when(
            direction == "VEC" ~ "VEC markers",
            direction == "LEC" ~ "LEC markers",
            TRUE ~ "Other"
        )
    )

category_colours <- c(
    "VEC markers" = "#c51b7d",   # red/crimson
    "LEC markers" = "#762a83",   # purple
    "Other"       = "grey80"
)

plt <- ggplot(wt_mut, aes(x = avg_log2FC_meox1_mutant, y = avg_log2FC_WT, color = category)) +
    geom_point(data = filter(wt_mut, category == "Other"), size = 1, alpha = 0.5) +
    geom_point(data = filter(wt_mut, category != "Other"), size = 1) +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
    scale_color_manual(values = category_colours) +
    labs(
        x = "Av. log2FC (LEC meox1 mutant vs VEC meox1 mutant)",
        y = "Av. log2FC (LEC WT vs VEC WT)",
        color = NULL
    ) +
    theme_classic() +
    theme(legend.position = "right")

print(plt)
```

<img src="meox1_L3_correlation_logfc_files/figure-html/vec_lec-1.png" alt="" style="display: block; margin: auto;" />

``` r
save_plot(plt, "LECVEC_WT_vs_LECVEC_mut_markers")
```

```
## *PDF saved to:* `../output/figure_extended/embo_scores//LECVEC_WT_vs_LECVEC_mut_markers.pdf`
```



# Perform differential expression analysis LEC vs all VEC.

Only LEC vs All VECs


``` r
# generate additional degs per condition
## WT VEC vs WT LEC // Mut VEC vs Mut LEC
### WT VEC vs WT LEC
Idents(data) <- "celltype_broad_lec_only"
data_wt_lec_vec <- FindMarkers(object = data,
                group.by = "celltype_broad_lec_only",
                ident.1 = "LEC_WT",
                ident.2 = "VEC_combined_WT",
                min.pct = 0.01,
                assay = "RNA",
                logfc.threshold = 0) %>%
    rownames_to_column("Gene") %>%
    mutate(group = "lec_vec_wt")

### Mut VEC vs Mut LEC
Idents(data) <- "celltype_broad_lec_only"
data_mut <- FindMarkers(object = data,
                group.by = "celltype_broad_lec_only",
                ident.1 = "LEC_meox1_mutant",
                ident.2 = "VEC_combined_meox1_mutant",
                min.pct = 0.01,
                assay = "RNA",
                logfc.threshold = 0) %>%
    rownames_to_column("Gene") %>%
    mutate(group = "lec_vec_mut")

# prep data for vis
data_wt_lec_vec <- data_wt_lec_vec %>%
    select(Gene, avg_log2FC_WT = avg_log2FC)
data_mut <- data_mut %>%
    select(Gene, avg_log2FC_meox1_mutant = avg_log2FC)

wt_mut <- inner_join(data_wt_lec_vec, data_mut, by = "Gene")
sum(wt_mut$Gene %in% genes_labels$Gene)
```

# Only LEC vs All VEC

These are very information-dense plots. How to interpret:

- Individual dots correspond to genes.
- Vertical axis shows _how much more/less is this gene expressed in LEC vs VEC WT_.
- Horizontal axis shows _how much more/less is this gene expressed in LEC vs VEC meox1 mutant_.
- Genes falling in the _top right and bottom left quadrants show concordance_ across both WT and mutant.
- Genes falling in the _top left and bottom right quadrants show discordance_ across both WT and mutant.

Examples for reference:
- top right: genes _up_ in the mutant and wildtype
- bottom left: genes _down_ in the mutant and wildtype
- top left: genes _down_ in mutant but up in wildtype
- bottom right: genes _up_ in the mutant but down in wildtype

_Suggested legend:_ Scatterplot shows the correlation between the average log fold changes of two differential gene expression analyses: [LECs vs VECs wildtype], [LECs vs VECs meox1 mutants]. Individual data points are genes and coloured by a curated list of LEC and VEC markers. The vertical axis answers the question: _how much more/less is this gene expressed in LEC vs VEC WT_, while the horizontal axis shows _how much more/less is this gene expressed in LEC vs VEC meox1 mutants_. Genes falling in the _top right and bottom left quadrants show concordance_ across the wildtype and meox1 mutants, while genes falling in the _top left and bottom right quadrants show discordance_ across both the wildtype and meox1 mutants.

## Plot concordant and discordant markers


``` r
wt_mut <- wt_mut %>%
    mutate(
        highlight = Gene %in% genes_of_interest,
        category = case_when(
            !highlight ~ "Other",
            avg_log2FC_WT > 0 & avg_log2FC_meox1_mutant > 0 ~ "Concordant (up)",
            avg_log2FC_WT < 0 & avg_log2FC_meox1_mutant < 0 ~ "Concordant (down)",
            avg_log2FC_WT * avg_log2FC_meox1_mutant < 0     ~ "Discordant",
            TRUE ~ "Other"  # catches exact-zero edge cases
        )
    )

table(wt_mut$category)
# Concordant (down)   Concordant (up)        Discordant             Other
#              1613               756               470             11267
# Concordant (down)   Concordant (up)        Discordant             Other
#               845              1508               488             11450

# colour and plot
category_colours <- c(
    "Concordant (up)"   = "green",   # e.g. red - both conditions upregulated
    "Concordant (down)" = "blue",   # e.g. blue - both conditions downregulated
    "Discordant"        = "purple",   # e.g. orange - opposite direction between conditions
    "Other"             = "grey80"     # background genes, not in your list of interest
)

plt <- ggplot(wt_mut, aes(x = avg_log2FC_meox1_mutant, y = avg_log2FC_WT, color = category)) +
    geom_point(data = filter(wt_mut, !highlight), size = 1, alpha = 0.4) +   # background points first
    geom_point(data = filter(wt_mut, highlight), size = 2.5) +               # highlighted genes on top
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey60") +  # y=x reference line
    scale_color_manual(values = category_colours) +
    labs(
        x = "Av. log2FC (LEC meox1 mutant vs VEC meox1 mutant)",
        y = "Av. log2FC (LEC WT vs VEC WT)",
        color = "Gene category",
        title = "Concordance of DE genes across conditions"
    ) +
    theme_classic()

print(plt)

save_plot(plt, "LECVEC_WT_vs_LECVEC_mut_direction_leconly")
```

## Plot LEC VEC scores


``` r
### plot by marker

wt_mut <- wt_mut %>%
    left_join(genes_labels, by = "Gene") %>%
    mutate(
        category = case_when(
            direction == "VEC" ~ "VEC markers",
            direction == "LEC" ~ "LEC markers",
            TRUE ~ "Other"
        )
    )

category_colours <- c(
    "VEC markers" = "#c51b7d",   # red/crimson
    "LEC markers" = "#762a83",   # purple
    "Other"       = "grey80"
)

plt <- ggplot(wt_mut, aes(x = avg_log2FC_meox1_mutant, y = avg_log2FC_WT, color = category)) +
    geom_point(data = filter(wt_mut, category == "Other"), size = 1, alpha = 0.5) +
    geom_point(data = filter(wt_mut, category != "Other"), size = 1) +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
    scale_color_manual(values = category_colours) +
    labs(
        x = "Av. log2FC (LEC meox1 mutant vs VEC meox1 mutant)",
        y = "Av. log2FC (LEC WT vs VEC WT)",
        color = NULL
    ) +
    theme_classic() +
    theme(legend.position = "right")

print(plt)

save_plot(plt, "LECVEC_WT_vs_LECVEC_mut_markers_leconly")
```



## Plot LEC VEC scores plus gene tags

Not in use, reference only.


``` r
### plot by marker and show outliers
top10_discordant <- wt_mut %>%
    mutate(divergence = abs(avg_log2FC_meox1_mutant - avg_log2FC_WT)) %>%
    arrange(desc(divergence)) %>%
    slice_head(n = 10)

# concordant: smallest |logFC_x - logFC_y| AMONG genes with matching sign in both axes
# (excludes near-zero/near-origin genes, which would trivially have small divergence
#  without representing a genuine concordant effect in either condition)
top10_concordant <- wt_mut %>%
    filter(sign(avg_log2FC_meox1_mutant) == sign(avg_log2FC_WT), avg_log2FC_meox1_mutant != 0, avg_log2FC_WT != 0) %>%
    mutate(divergence = abs(avg_log2FC_meox1_mutant - avg_log2FC_WT)) %>%
    arrange(divergence) %>%
    slice_head(n = 10)

outlier_genes <- bind_rows(
    top10_discordant %>% mutate(outlier_type = "Top discordant"),
    top10_concordant %>% mutate(outlier_type = "Top concordant")
)

plt_outliers <- plt +
    geom_text_repel(
        data = outlier_genes,
        aes(label = Gene),
        size = 3,
        max.overlaps = Inf,
        show.legend = FALSE,
        color = "black"
    ) +
    geom_point(
        data = outlier_genes,
        shape = 21, size = 3, stroke = 0.8, color = "black", fill = NA
    )

print(plt_outliers)

save_plot(plt, "LECVEC_WT_vs_LECVEC_mut_markers_outliers")
```

## Plot other subtypes

Not in use, reference only.


``` r
#### hold ####
# generate additional degs per condition
## WT LEC vs WT VEC // Mut VEC vs Mut VEC
### Mut VEC vs WT VEC
Idents(data) <- "celltype_broad_genotype"
data_vec <- FindMarkers(object = data,
                group.by = "celltype_broad_genotype",
                ident.1 = "VEC_combined_meox1_mutant",
                ident.2 = "VEC_combined_WT",
                min.pct = 0.01,
                assay = "RNA",
                logfc.threshold = 0) %>%
    rownames_to_column("Gene") %>%
    mutate(group = "vec_mut_wt")

# prep data for vis
data_vec <- data_vec %>%
    select(Gene, avg_log2FC_vec_mut_wt = avg_log2FC)

wt_vs_vec <- inner_join(data_wt_vec_lec, data_vec, by = "Gene")
wt_vs_vec <- wt_vs_vec %>%
    left_join(genes_labels, by = "Gene") %>%
    mutate(
        category = case_when(
            direction == "VEC" ~ "VEC markers",
            direction == "LEC" ~ "LEC markers",
            TRUE ~ "Other"
        )
    )

category_colours <- c(
    "VEC markers" = "#c51b7d",   # red/crimson
    "LEC markers" = "#762a83",   # purple
    "Other"       = "grey80"
)

plt <- ggplot(wt_vs_vec, aes(x = avg_log2FC_vec_mut_wt, y = avg_log2FC_WT, color = category)) +
    geom_point(data = filter(wt_vs_vec, category == "Other"), size = 1, alpha = 0.5) +
    geom_point(data = filter(wt_vs_vec, category != "Other"), size = 1) +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
    scale_color_manual(values = category_colours) +
    labs(
        x = "Av. log2FC (VEC meox1 mutant vs VEC WT)",
        y = "Av. log2FC (VEC WT vs LEC WT)",
        color = NULL
    ) +
    theme_classic() +
    theme(legend.position = "right")

print(plt)
```

<img src="meox1_L3_correlation_logfc_files/figure-html/subtypes-1.png" alt="" style="display: block; margin: auto;" />

``` r
save_plot(plt, "VECLEC_WT_vs_VEC")
```

```
## *PDF saved to:* `../output/figure_extended/embo_scores//VECLEC_WT_vs_VEC.pdf`
```

``` r
## WT VEC vs WT LEC // Mut LEC vs WT LEC
### Mut LEC vs WT LEC
Idents(data) <- "celltype_broad_genotype"
data_lec <- FindMarkers(object = data,
                group.by = "celltype_broad_genotype",
                ident.1 = "LEC_combined_meox1_mutant",
                ident.2 = "LEC_combined_WT",
                min.pct = 0.01,
                assay = "RNA",
                logfc.threshold = 0) %>%
    rownames_to_column("Gene") %>%
    mutate(group = "lec_mut_wt")

# prep data for vis
data_lec <- data_lec %>%
    select(Gene, avg_log2FC_lec_mut_wt = avg_log2FC)

wt_vs_lec <- inner_join(data_wt_lec_vec, data_lec, by = "Gene")
wt_vs_lec <- wt_vs_lec %>%
    left_join(genes_labels, by = "Gene") %>%
    mutate(
        category = case_when(
            direction == "VEC" ~ "VEC markers",
            direction == "LEC" ~ "LEC markers",
            TRUE ~ "Other"
        )
    )

category_colours <- c(
    "VEC markers" = "#c51b7d",   # red/crimson
    "LEC markers" = "#762a83",   # purple
    "Other"       = "grey80"
)

plt <- ggplot(wt_vs_lec, aes(x = avg_log2FC_lec_mut_wt, y = avg_log2FC_WT, color = category)) +
    geom_point(data = filter(wt_vs_lec, category == "Other"), size = 1, alpha = 0.5) +
    geom_point(data = filter(wt_vs_lec, category != "Other"), size = 1) +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
    scale_color_manual(values = category_colours) +
    labs(
        x = "Av. log2FC (LEC meox1 mutant vs LEC WT)",
        y = "Av. log2FC (LEC WT vs VEC WT)",
        color = NULL
    ) +
    theme_classic() +
    theme(legend.position = "right")

print(plt)
```

<img src="meox1_L3_correlation_logfc_files/figure-html/subtypes-2.png" alt="" style="display: block; margin: auto;" />

``` r
save_plot(plt, "LECVEC_WT_vs_LEC")
```

```
## *PDF saved to:* `../output/figure_extended/embo_scores//LECVEC_WT_vs_LEC.pdf`
```

# For developers

Sample run command:


``` bash
Rscript -e "
rmarkdown::render(
  'meox1_L3_correlation_logfc.Rmd',
  output_file = './meox1_L3_correlation_logfc.html'
)
"
```

---
title: "GO analyses for meox1"
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

# GO analyses for meox1

## Important info

We want to compare ontology across different groups and subgroups of wild type and mutant cells.

- This assumes you ran `meox1_deg_confects_prep.Rmd` first, in which we merged VECs and LECs into groups for comparison.
- Two algorithms are used for complementarity (they are conceptually different but functionally similar)
- `FGSEA` can "fail" in some conditions since it makes assumptions which are more consistent with bulk RNA-Seq data. At the same time, expect that `ORA/EGO` will return more results. This is not an error.
- To minimise information overload, a brief description is embedded together with each figure so they can be interpreted in context.

### Quickstart instructions for getting the info you want:

- Select your `GO` analysis of choice in the sidebar (ie 1.5 onwards)
- Suggest starting with the `ORA` plots as the core component.
  - ORA Barplots, dotplots and enrichment map show the same information in different ways
  - ORA enrichment map and category-gene network show higher level pathway information
- FGSEA is plotted as a supplementary complement to ORA and are not always available
  - While useful, they can be safely ignored if absent

### Experimental design information

We perform two different analyses, on two different sample groupings, on two different sample groups. Best described with the tree below

```
DATASET -- CELL STATES (MERGED/SINGLE) -- GO ANALYSIS TYPE

                                 +---- FGSEA ranked analysis
                                 |
              +-- single state --+
              |                  |
              |                  +---- ORA/EGO overrepresentation analysis
meox data ----|
              |                  +---- FGSEA ranked analysis
              |                  |
              +-- merged state --+
                                 |
                                 +---- ORA/EGO overrepresentation analysis

                                 +---- FGSEA ranked analysis
                                 |
              +-- single state --+
              |                  |
              |                  +---- ORA/EGO overrepresentation analysis
yap1 data ----|
              |                  +---- FGSEA ranked analysis
              |                  |
              +-- merged state --+
                                 |
                                 +---- ORA/EGO overrepresentation analysis
```

### Methods information

We compare and contrast two different gene ontology enrichment methods.

| Aspect | GSEA (`gseGO`/fgsea) | ORA (`enrichGO`) |
|---|---|---|
| **Input** | Full ranked gene list (all genes, ranked by `avg_log2FC`) | Binary significant/not-significant gene list (`p_val < 0.05` cutoff) |
| **Statistical test** | Running-sum enrichment score against a ranked list (Kolmogorov-Smirnov-like) | Hypergeometric / Fisher's exact test vs. background |
| **Cutoff dependency** | None — every gene contributes proportionally to its rank | Hard cutoff required — genes just below threshold contribute nothing |
| **Directionality** | Native — `NES > 0` / `NES < 0` splits up- vs down-enriched | Not native — requires re-running ORA separately on up/down subsets (reduces power) |
| **Sensitive to** | Coordinated small shifts across many genes in a pathway | Strong individual gene-level significance |
| **Best suited for** | Well-powered fold-change estimates across full transcriptome/gene set | Clean list of "hits" where direction of each term isn't the focus |
| **Known caveat here** | Single-cell `p_val` is often inflated (large cell counts) — ranking uses `avg_log2FC`, not p-value, to avoid this | Empty results possible for weak/underpowered comparisons (guarded in code) |

Note that `GSEA` methods were originally developed for and inherit the major assumptions of bulk RNA-Seq data. They can still be used for single-cell data, but p-value inflation means their precision can be reduced. To buffer this effect, we rank by fold changes and not by p values. Expect the majority of the output below to work on `ORA` but not `GSEA`.

## Initial setup

For developers.


```{.r .fold-hide}
# !/usr/bin/Rscript
library(clusterProfiler)
```

```
## 
```

```
## clusterProfiler v4.14.6 Learn more at https://yulab-smu.top/contribution-knowledge-mining/
## 
## Please cite:
## 
## T Wu, E Hu, S Xu, M Chen, P Guo, Z Dai, T Feng, L Zhou, W Tang, L Zhan,
## X Fu, S Liu, X Bo, and G Yu. clusterProfiler 4.0: A universal
## enrichment tool for interpreting omics data. The Innovation. 2021,
## 2(3):100141
```

```
## 
## Attaching package: 'clusterProfiler'
```

```
## The following object is masked from 'package:stats':
## 
##     filter
```

```{.r .fold-hide}
library(enrichplot)
```

```
## enrichplot v1.26.6 Learn more at https://yulab-smu.top/contribution-knowledge-mining/
## 
## Please cite:
## 
## G Yu. Thirteen years of clusterProfiler. The Innovation. 2024,
## 5(6):100722
```

```{.r .fold-hide}
library(ggpubr)
```

```
## Loading required package: ggplot2
```

```
## 
## Attaching package: 'ggpubr'
```

```
## The following object is masked from 'package:enrichplot':
## 
##     color_palette
```

```{.r .fold-hide}
library(org.Dr.eg.db)
```

```
## Loading required package: AnnotationDbi
```

```
## Loading required package: stats4
```

```
## Loading required package: BiocGenerics
```

```
## 
## Attaching package: 'BiocGenerics'
```

```
## The following objects are masked from 'package:stats':
## 
##     IQR, mad, sd, var, xtabs
```

```
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
##     get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
##     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
##     Position, rank, rbind, Reduce, rownames, sapply, saveRDS, setdiff,
##     table, tapply, union, unique, unsplit, which.max, which.min
```

```
## Loading required package: Biobase
```

```
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
```

```
## Loading required package: IRanges
```

```
## Loading required package: S4Vectors
```

```
## 
## Attaching package: 'S4Vectors'
```

```
## The following object is masked from 'package:clusterProfiler':
## 
##     rename
```

```
## The following object is masked from 'package:utils':
## 
##     findMatches
```

```
## The following objects are masked from 'package:base':
## 
##     expand.grid, I, unname
```

```
## 
## Attaching package: 'IRanges'
```

```
## The following object is masked from 'package:clusterProfiler':
## 
##     slice
```

```
## 
## Attaching package: 'AnnotationDbi'
```

```
## The following object is masked from 'package:clusterProfiler':
## 
##     select
```

```
## 
```

```{.r .fold-hide}
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
## Attaching package: 'sp'
```

```
## The following object is masked from 'package:IRanges':
## 
##     %over%
```

```
## 
## Attaching package: 'SeuratObject'
```

```
## The following object is masked from 'package:IRanges':
## 
##     intersect
```

```
## The following object is masked from 'package:S4Vectors':
## 
##     intersect
```

```
## The following object is masked from 'package:BiocGenerics':
## 
##     intersect
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, t
```

```{.r .fold-hide}
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
## ✖ lubridate::%within%()    masks IRanges::%within%()
## ✖ dplyr::collapse()        masks IRanges::collapse()
## ✖ dplyr::combine()         masks Biobase::combine(), BiocGenerics::combine()
## ✖ dplyr::desc()            masks IRanges::desc()
## ✖ tidyr::expand()          masks S4Vectors::expand()
## ✖ dplyr::filter()          masks clusterProfiler::filter(), stats::filter()
## ✖ dplyr::first()           masks S4Vectors::first()
## ✖ dplyr::lag()             masks stats::lag()
## ✖ BiocGenerics::Position() masks ggplot2::Position(), base::Position()
## ✖ purrr::reduce()          masks IRanges::reduce()
## ✖ dplyr::rename()          masks S4Vectors::rename(), clusterProfiler::rename()
## ✖ lubridate::second()      masks S4Vectors::second()
## ✖ lubridate::second<-()    masks S4Vectors::second<-()
## ✖ dplyr::select()          masks AnnotationDbi::select(), clusterProfiler::select()
## ✖ purrr::simplify()        masks clusterProfiler::simplify()
## ✖ dplyr::slice()           masks IRanges::slice(), clusterProfiler::slice()
## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors
```

```{.r .fold-hide}
library(UpSetR)

outfile_dir <- "../output/figure_extended/go/"

singles_meox_path <- "../../Saki_data/analysis/D10051_meox1_dataset/DEGs/Level_03_DEG_mutVSwt_L3_celltype_fc0.00_minpct_0.01.csv"
grouped_meox_path <- "../output/figure_extended/dge_confects/meox1_Level_03_DEG_mutVSwt_L3_celltype_merged_VEC_merged_LEC_fc0.00_minpct_0.01.csv"

singles_meox_hippo_path <- "../../Saki_data/analysis/D10051_meox1_dataset/DEGs/analysis/D10051_meox1_dataset/DEGs/Level_03_DEG_hippo_targets_mutVSwt_L3_celltype_fc0.00_minpct_0.00.csv"
grouped_meox_hippo_path <- "../output/figure_extended/dge_confects/meox1_Level_03_DEG_mutVSwt_L3_celltype_merged_VEC_merged_LEC_hippo_targets_fc0.00_minpct_0.00.csv"

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

keytypes(org.Dr.eg.db)
```

```
##  [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS"
##  [6] "ENTREZID"     "ENZYME"       "EVIDENCE"     "EVIDENCEALL"  "GENENAME"    
## [11] "GO"           "GOALL"        "IPI"          "ONTOLOGY"     "ONTOLOGYALL" 
## [16] "PATH"         "PFAM"         "PMID"         "PROSITE"      "REFSEQ"      
## [21] "SYMBOL"       "UNIPROT"      "ZFIN"
```

## Wrapper for downstream functions

For developers. We will reuse this a lot so we make a wrapper function for `FGSEA` and `enrichGO`.


```{.r .fold-hide}
run_go <- function(data, data_name) {
    cat("\n\n## GO analysis:", data_name, "\n\n")
    cat("Input genes:", nrow(data), "\n\n")

    # enrichGO - split into up/down, plus combined
    sig_gene_list      <- data$Gene[data$p_val < 0.05]
    sig_gene_list_up    <- data$Gene[data$p_val < 0.05 & data$avg_log2FC > 0]
    sig_gene_list_down  <- data$Gene[data$p_val < 0.05 & data$avg_log2FC < 0]

    cat("**ORA input genes:** combined =", length(sig_gene_list),
        "&nbsp;&nbsp; up =", length(sig_gene_list_up),
        "&nbsp;&nbsp; down =", length(sig_gene_list_down), "\n\n")

    run_ora <- function(genes) {
        if (length(genes) < 2) return(NULL)
        enrichGO(
            gene = genes,
            keyType = "SYMBOL",
            OrgDb = org.Dr.eg.db,
            ont = "BP",
            pAdjustMethod = "BH",
            pvalueCutoff = 0.05,
            qvalueCutoff = 0.05,
            readable = TRUE
        )
    }

    ego_result      <- run_ora(sig_gene_list)
    ego_result_up   <- run_ora(sig_gene_list_up)
    ego_result_down <- run_ora(sig_gene_list_down)

    n_ora_terms      <- if (!is.null(ego_result))      nrow(ego_result@result)      else 0
    n_ora_up_terms   <- if (!is.null(ego_result_up))   nrow(ego_result_up@result)   else 0
    n_ora_down_terms <- if (!is.null(ego_result_down)) nrow(ego_result_down@result) else 0

    cat("**ORA (combined):**", n_ora_terms, "terms &nbsp;&nbsp; ",
        "**ORA (up):**", n_ora_up_terms, "terms &nbsp;&nbsp; ",
        "**ORA (down):**", n_ora_down_terms, "terms\n\n")

    # export full ORA result tables (all terms, not just significant) - one
    # checkpoint per block, right after each enrichGO() call completes
    export_ora_table <- function(ego_res, n_terms, file_tag) {
        if (n_terms == 0 || is.null(ego_res)) return(invisible(NULL))
        csv_path <- paste0(outfile_dir, "/", data_name, "_ego_", file_tag, "_full_results.csv")
        write.csv(ego_res@result, csv_path, row.names = FALSE)
        cat("*Full ORA (", file_tag, ") table saved to:*", paste0("`", csv_path, "`"), "\n\n")
    }

    export_ora_table(ego_result,      n_ora_terms,      "combined")
    export_ora_table(ego_result_up,   n_ora_up_terms,   "up")
    export_ora_table(ego_result_down, n_ora_down_terms, "down")

    # significant terms - guard against empty results
    ora_sig_terms <- if (n_ora_terms > 0) {
        ego_result@result$ID[ego_result@result$p.adjust < 0.05]
    } else character(0)

    ora_up_sig_terms <- if (n_ora_up_terms > 0) {
        ego_result_up@result$ID[ego_result_up@result$p.adjust < 0.05]
    } else character(0)

    ora_down_sig_terms <- if (n_ora_down_terms > 0) {
        ego_result_down@result$ID[ego_result_down@result$p.adjust < 0.05]
    } else character(0)

    cat("**ORA significant terms:** combined =", length(ora_sig_terms),
        ", up =", length(ora_up_sig_terms), ", down =", length(ora_down_sig_terms), "\n\n")

    # also export just the significant subset separately, for convenience -
    # avoids downstream users having to re-filter the full table by p.adjust
    export_sig_table <- function(ego_res, sig_ids, file_tag) {
        if (length(sig_ids) == 0 || is.null(ego_res)) return(invisible(NULL))
        sig_df <- ego_res@result[ego_res@result$ID %in% sig_ids, ]
        csv_path <- paste0(outfile_dir, "/", data_name, "_ego_", file_tag, "_significant_only.csv")
        write.csv(sig_df, csv_path, row.names = FALSE)
        cat("*Significant-only ORA (", file_tag, ") table saved to:*", paste0("`", csv_path, "`"), "\n\n")
    }

    export_sig_table(ego_result,      ora_sig_terms,      "combined")
    export_sig_table(ego_result_up,   ora_up_sig_terms,   "up")
    export_sig_table(ego_result_down, ora_down_sig_terms, "down")

    # helper to save a ggplot object in both pdf and png, log the pdf path clearly
    save_plot <- function(plt, plot_name, width = 8, height = 12) {
        pdf_path <- paste0(outfile_dir, "/", data_name, "_", plot_name, ".pdf")
        png_path <- paste0(outfile_dir, "/", data_name, "_", plot_name, ".png")
        ggsave(filename = pdf_path, plot = plt, width = width, height = height)
        ggsave(filename = png_path, plot = plt, width = width, height = height)
        cat("*PDF saved to:*", paste0("`", pdf_path, "`"), "\n\n")
    }

    # reusable block for one ORA result (combined / up / down) - avoids
    # triplicating the same five-plot sequence three times over
    run_ora_plots <- function(ego_res, n_terms, block_label, file_tag) {
        if (n_terms == 0 || dim(ego_res)[1] == 0) {
            cat("*No ORA (", block_label, ") terms found for", data_name,
                "- skipping ORA (", block_label, ") plots*\n\n")
            return(invisible(NULL))
        }

        cat("### ORA (", block_label, "): Barplot\n\n")
        cat("Shows the top enriched GO terms ranked by significance, bar length/colour reflecting ",
            "adjusted p-value or gene count. **Expected:** terms biologically consistent with the ",
            "comparison (e.g. vascular/endothelial terms for EC-related contrasts). **Unexpected:** ",
            "dominance by generic housekeeping terms (translation, metabolism) with no thematic ",
            "coherence - may indicate a weak or noisy input gene list rather than a real signal.\n\n")
        plt <- barplot(ego_res, showCategory = 20)
        print(plt)
        save_plot(plt, paste0("ego_", file_tag, "_barplot"))

        cat("### ORA (", block_label, "): Dotplot\n\n")
        cat("Adds gene ratio (fraction of input genes hitting each term) as point size, with colour ",
            "for significance - more informative than the barplot since it shows both effect size and ",
            "confidence simultaneously. **Expected:** larger, more significant dots clustering around a ",
            "coherent biological theme. **Unexpected:** significant terms with very small gene ratios ",
            "(few genes driving a 'significant' result) - treat cautiously as potentially fragile.\n\n")
        plt <- dotplot(ego_res, showCategory = 20)
        print(plt)
        save_plot(plt, paste0("ego_", file_tag, "_dotplot"))

        cat("### ORA (", block_label, "): Enrichment map\n\n")
        cat("Network graph connecting GO terms that share overlapping genes - clusters of connected ",
            "terms suggest a shared underlying biological process rather than independent findings. ",
            "**Expected:** tight clusters of related terms. **Unexpected:** many isolated, unconnected ",
            "terms, suggesting the significant terms are thematically scattered rather than coherent.\n\n")
        ego_res_sim <- pairwise_termsim(ego_res)
        plt <- emapplot(ego_res_sim, showCategory = 20)
        print(plt)
        save_plot(plt, paste0("ego_", file_tag, "_emapplot"), width = 10, height = 10)

        cat("### ORA (", block_label, "): Category-gene network\n\n")
        cat("Shows which specific genes drive which enriched terms, useful for tracing a term back to ",
            "the actual DE genes behind it. **Expected:** recognisable marker genes linking to expected ",
            "terms. **Unexpected:** a term driven entirely by one or two genes rather than a broad gene ",
            "set - weaker evidence than a term supported by many genes.\n\n")
        plt <- cnetplot(ego_res, showCategory = 10)
        print(plt)
        save_plot(plt, paste0("ego_", file_tag, "_cnetplot"), width = 10, height = 10)
    }

    ## EGO plots - combined, then up, then down - each independently guarded
    run_ora_plots(ego_result,      n_ora_terms,      "combined", "combined")
    run_ora_plots(ego_result_up,   n_ora_up_terms,   "up",       "up")
    run_ora_plots(ego_result_down, n_ora_down_terms, "down",     "down")

    cat("**Done:**", data_name, "\n\n---\n\n")

    # return everything for later reuse/inspection
    list(
        ego_result = ego_result,
        ego_result_up = ego_result_up,
        ego_result_down = ego_result_down,
        ora_sig_terms = ora_sig_terms,
        ora_up_sig_terms = ora_up_sig_terms,
        ora_down_sig_terms = ora_down_sig_terms
    )
}

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
grouped_meox <- read.csv(grouped_meox_path)

grouped_lecs <- grouped_meox[grouped_meox$group == "LEC__pre_muLEC", ]
rownames(grouped_lecs) <- grouped_lecs$Gene

grouped_vecs <- grouped_meox[grouped_meox$group == "hmVEC__mVEC", ]
rownames(grouped_vecs) <- grouped_vecs$Gene


singles_meox <- read.csv(singles_meox_path)

singles_LEC <- singles_meox[singles_meox$group == "LEC", ]
rownames(singles_LEC) <- singles_LEC$Gene

singles_pre_muLEC <- singles_meox[singles_meox$group == "pre_muLEC", ]
rownames(singles_pre_muLEC) <- singles_pre_muLEC$Gene

singles_hmVEC <- singles_meox[singles_meox$group == "hmVEC", ]
rownames(singles_hmVEC) <- singles_hmVEC$Gene

singles_mVEC <- singles_meox[singles_meox$group == "mVEC", ]
rownames(singles_mVEC) <- singles_mVEC$Gene
```

We have separate sections for each plot series for easy comparison.


```{.r .fold-hide}
grouped_lecs_out <- run_go(grouped_lecs, data_name="LEC__pre_muLEC")
```



## GO analysis: LEC__pre_muLEC 

Input genes: 15179 

**ORA input genes:** combined = 1624 &nbsp;&nbsp; up = 900 &nbsp;&nbsp; down = 724 

**ORA (combined):** 2548 terms &nbsp;&nbsp;  **ORA (up):** 2077 terms &nbsp;&nbsp;  **ORA (down):** 1803 terms

*Full ORA ( combined ) table saved to:* `../output/figure_extended/go//LEC__pre_muLEC_ego_combined_full_results.csv` 

*Full ORA ( up ) table saved to:* `../output/figure_extended/go//LEC__pre_muLEC_ego_up_full_results.csv` 

*Full ORA ( down ) table saved to:* `../output/figure_extended/go//LEC__pre_muLEC_ego_down_full_results.csv` 

**ORA significant terms:** combined = 40 , up = 75 , down = 12 

*Significant-only ORA ( combined ) table saved to:* `../output/figure_extended/go//LEC__pre_muLEC_ego_combined_significant_only.csv` 

*Significant-only ORA ( up ) table saved to:* `../output/figure_extended/go//LEC__pre_muLEC_ego_up_significant_only.csv` 

*Significant-only ORA ( down ) table saved to:* `../output/figure_extended/go//LEC__pre_muLEC_ego_down_significant_only.csv` 

### ORA ( combined ): Barplot

Shows the top enriched GO terms ranked by significance, bar length/colour reflecting  adjusted p-value or gene count. **Expected:** terms biologically consistent with the  comparison (e.g. vascular/endothelial terms for EC-related contrasts). **Unexpected:**  dominance by generic housekeeping terms (translation, metabolism) with no thematic  coherence - may indicate a weak or noisy input gene list rather than a real signal.

```
## Warning in fortify(object, showCategory = showCategory, by = x, ...): Arguments in `...` must be used.
## ✖ Problematic argument:
## • by = x
## ℹ Did you misspell an argument name?
```

```
## Warning: `aes_string()` was deprecated in ggplot2 3.0.0.
## ℹ Please use tidy evaluation idioms with `aes()`.
## ℹ See also `vignette("ggplot2-in-packages")` for more information.
## ℹ The deprecated feature was likely used in the enrichplot package.
##   Please report the issue at <https://github.com/GuangchuangYu/enrichplot/issues>.
## This warning is displayed once per session.
## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was generated.
```

![](meox1_goanalysis_files/figure-html/grp_lec-1.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//LEC__pre_muLEC_ego_combined_barplot.pdf` 

### ORA ( combined ): Dotplot

Adds gene ratio (fraction of input genes hitting each term) as point size, with colour  for significance - more informative than the barplot since it shows both effect size and  confidence simultaneously. **Expected:** larger, more significant dots clustering around a  coherent biological theme. **Unexpected:** significant terms with very small gene ratios  (few genes driving a 'significant' result) - treat cautiously as potentially fragile.

![](meox1_goanalysis_files/figure-html/grp_lec-2.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//LEC__pre_muLEC_ego_combined_dotplot.pdf` 

### ORA ( combined ): Enrichment map

Network graph connecting GO terms that share overlapping genes - clusters of connected  terms suggest a shared underlying biological process rather than independent findings.  **Expected:** tight clusters of related terms. **Unexpected:** many isolated, unconnected  terms, suggesting the significant terms are thematically scattered rather than coherent.

```
## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
## ℹ Please use `linewidth` instead.
## ℹ The deprecated feature was likely used in the ggtangle package.
##   Please report the issue to the authors.
## This warning is displayed once per session.
## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was generated.
```

![](meox1_goanalysis_files/figure-html/grp_lec-3.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//LEC__pre_muLEC_ego_combined_emapplot.pdf` 

### ORA ( combined ): Category-gene network

Shows which specific genes drive which enriched terms, useful for tracing a term back to  the actual DE genes behind it. **Expected:** recognisable marker genes linking to expected  terms. **Unexpected:** a term driven entirely by one or two genes rather than a broad gene  set - weaker evidence than a term supported by many genes.

![](meox1_goanalysis_files/figure-html/grp_lec-4.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//LEC__pre_muLEC_ego_combined_cnetplot.pdf` 

### ORA ( up ): Barplot

Shows the top enriched GO terms ranked by significance, bar length/colour reflecting  adjusted p-value or gene count. **Expected:** terms biologically consistent with the  comparison (e.g. vascular/endothelial terms for EC-related contrasts). **Unexpected:**  dominance by generic housekeeping terms (translation, metabolism) with no thematic  coherence - may indicate a weak or noisy input gene list rather than a real signal.

```
## Warning in fortify(object, showCategory = showCategory, by = x, ...): Arguments in `...` must be used.
## ✖ Problematic argument:
## • by = x
## ℹ Did you misspell an argument name?
```

![](meox1_goanalysis_files/figure-html/grp_lec-5.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//LEC__pre_muLEC_ego_up_barplot.pdf` 

### ORA ( up ): Dotplot

Adds gene ratio (fraction of input genes hitting each term) as point size, with colour  for significance - more informative than the barplot since it shows both effect size and  confidence simultaneously. **Expected:** larger, more significant dots clustering around a  coherent biological theme. **Unexpected:** significant terms with very small gene ratios  (few genes driving a 'significant' result) - treat cautiously as potentially fragile.

![](meox1_goanalysis_files/figure-html/grp_lec-6.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//LEC__pre_muLEC_ego_up_dotplot.pdf` 

### ORA ( up ): Enrichment map

Network graph connecting GO terms that share overlapping genes - clusters of connected  terms suggest a shared underlying biological process rather than independent findings.  **Expected:** tight clusters of related terms. **Unexpected:** many isolated, unconnected  terms, suggesting the significant terms are thematically scattered rather than coherent.

![](meox1_goanalysis_files/figure-html/grp_lec-7.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//LEC__pre_muLEC_ego_up_emapplot.pdf` 

### ORA ( up ): Category-gene network

Shows which specific genes drive which enriched terms, useful for tracing a term back to  the actual DE genes behind it. **Expected:** recognisable marker genes linking to expected  terms. **Unexpected:** a term driven entirely by one or two genes rather than a broad gene  set - weaker evidence than a term supported by many genes.

![](meox1_goanalysis_files/figure-html/grp_lec-8.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//LEC__pre_muLEC_ego_up_cnetplot.pdf` 

### ORA ( down ): Barplot

Shows the top enriched GO terms ranked by significance, bar length/colour reflecting  adjusted p-value or gene count. **Expected:** terms biologically consistent with the  comparison (e.g. vascular/endothelial terms for EC-related contrasts). **Unexpected:**  dominance by generic housekeeping terms (translation, metabolism) with no thematic  coherence - may indicate a weak or noisy input gene list rather than a real signal.

```
## Warning in fortify(object, showCategory = showCategory, by = x, ...): Arguments in `...` must be used.
## ✖ Problematic argument:
## • by = x
## ℹ Did you misspell an argument name?
```

![](meox1_goanalysis_files/figure-html/grp_lec-9.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//LEC__pre_muLEC_ego_down_barplot.pdf` 

### ORA ( down ): Dotplot

Adds gene ratio (fraction of input genes hitting each term) as point size, with colour  for significance - more informative than the barplot since it shows both effect size and  confidence simultaneously. **Expected:** larger, more significant dots clustering around a  coherent biological theme. **Unexpected:** significant terms with very small gene ratios  (few genes driving a 'significant' result) - treat cautiously as potentially fragile.

![](meox1_goanalysis_files/figure-html/grp_lec-10.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//LEC__pre_muLEC_ego_down_dotplot.pdf` 

### ORA ( down ): Enrichment map

Network graph connecting GO terms that share overlapping genes - clusters of connected  terms suggest a shared underlying biological process rather than independent findings.  **Expected:** tight clusters of related terms. **Unexpected:** many isolated, unconnected  terms, suggesting the significant terms are thematically scattered rather than coherent.

![](meox1_goanalysis_files/figure-html/grp_lec-11.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//LEC__pre_muLEC_ego_down_emapplot.pdf` 

### ORA ( down ): Category-gene network

Shows which specific genes drive which enriched terms, useful for tracing a term back to  the actual DE genes behind it. **Expected:** recognisable marker genes linking to expected  terms. **Unexpected:** a term driven entirely by one or two genes rather than a broad gene  set - weaker evidence than a term supported by many genes.

![](meox1_goanalysis_files/figure-html/grp_lec-12.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//LEC__pre_muLEC_ego_down_cnetplot.pdf` 

**Done:** LEC__pre_muLEC 

---

```{.r .fold-hide}
cat("### Custom plots for downregulated genes in mutant vs wildtype\n\n")
```

### Custom plots for downregulated genes in mutant vs wildtype

```{.r .fold-hide}
# extract all genes associated with
ego_lecs <- grouped_lecs_out$ego_result_down@result

ego_lecs$GeneRatio_numeric <- sapply(
    strsplit(ego_lecs$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])
    )
ego_lecs_top15_terms <- ego_lecs[
    order(-ego_lecs$GeneRatio_numeric), ][1:15,
    ]
ego_lecs_top15_gene_map <- setNames(strsplit(
    ego_lecs_top15_terms$geneID, "/"), ego_lecs_top15_terms$Description
    )
ego_lecs_top15_gene_list <- unique(unlist(unname(ego_lecs_top15_gene_map)))

# strip the si: stuff
ego_lecs_top15_gene_list <- ego_lecs_top15_gene_list[
    grep(":", ego_lecs_top15_gene_list, invert=TRUE)
    ]

ego_lecs_top15_pathways <- data.frame(ego_lecs$ID, ego_lecs$Description, ego_lecs$geneID)
# head(ego_lecs_top15_pathways, n=15)

cat("\n\nDotplots of gene ratio top 15 only.\n\n")
```



Dotplots of gene ratio top 15 only.

```{.r .fold-hide}
cat("Adds gene ratio (fraction of input genes hitting each term) as point size, with colour ",
    "for significance - more informative than the barplot since it shows both effect size and ",
    "confidence simultaneously. **Expected:** larger, more significant dots clustering around a ",
    "coherent biological theme. **Unexpected:** significant terms with very small gene ratios ",
    "(few genes driving a 'significant' result) - treat cautiously as potentially fragile.\n\n")
```

Adds gene ratio (fraction of input genes hitting each term) as point size, with colour  for significance - more informative than the barplot since it shows both effect size and  confidence simultaneously. **Expected:** larger, more significant dots clustering around a  coherent biological theme. **Unexpected:** significant terms with very small gene ratios  (few genes driving a 'significant' result) - treat cautiously as potentially fragile.

```{.r .fold-hide}
cat("\n\n_Suggested figure legend:_ Dotplot shows the effect size and confidence simultaneously for enriched GO terms. Units on the horizontal axis correspond to the fraction of input genes which hit each GO term. Point size corresponds to the quantity of terms, and point colour corresponds to adjusted p values. Top 15 GO terms by gene ratio are shown.\n\n")
```



_Suggested figure legend:_ Dotplot shows the effect size and confidence simultaneously for enriched GO terms. Units on the horizontal axis correspond to the fraction of input genes which hit each GO term. Point size corresponds to the quantity of terms, and point colour corresponds to adjusted p values. Top 15 GO terms by gene ratio are shown.

```{.r .fold-hide}
plt_name <- "top15_LECs_GO_pathways_downregulated_generatio_dotplot"
plt <- dotplot(grouped_lecs_out$ego_result_down, showCategory=15)
print(plt)
```

![](meox1_goanalysis_files/figure-html/grp_lec-13.png)<!-- -->

```{.r .fold-hide}
save_plot_custom(plt, plt_name)
```

*PDF saved to:* `../output/figure_extended/go//top15_LECs_GO_pathways_downregulated_generatio_dotplot.pdf` 

```{.r .fold-hide}
data_lecs <- data[, (data$L3_celltype == "LEC" | data$L3_celltype == "pre_muLEC")]
data_lecs$Genotype <- gsub("meox1_mutant", "sibling", data_lecs$Genotype)
data_lecs$Genotype <- gsub("wildtype", "sibling", data_lecs$Genotype)
data_lecs$Genotype <- factor(
    data_lecs$Genotype,
    levels = c("sibling", "mutant")
)

cat("\n\n**Genotype counts (cells) for LECs grouped:**",
    paste(names(table(data_lecs$Genotype)), table(data_lecs$Genotype), sep = " = ", collapse = ", "),
    "\n\n")
```



**Genotype counts (cells) for LECs grouped:** sibling = 686, mutant = 0 

```{.r .fold-hide}
cat("\n\nDotplots of gene expression (key pathways).\n\n")
```



Dotplots of gene expression (key pathways).

```{.r .fold-hide}
cat("For plot expansion, settled on these:\n\n
1. selected top 15 pathways as a cutoff\n
2. extracted the genes associated with these\n
3. removed genes with the `[a-z]+:` prefix\n
4. selected genes associated with key pathways and added meox1\n
5. plotted dotplots on the remainder\n\n

Attempted heatmaps, not very useful.\n\n")
```

For plot expansion, settled on these:


1. selected top 15 pathways as a cutoff

2. extracted the genes associated with these

3. removed genes with the `[a-z]+:` prefix

4. selected genes associated with key pathways and added meox1

5. plotted dotplots on the remainder



Attempted heatmaps, not very useful.

```{.r .fold-hide}
cat("\n\n_Suggested figure legend:_ Dotplot shows the scaled expression level for genes from enriched GO terms. Units on the horizontal axis correspond to the genotypes. Point size corresponds to the percent of expression in each cell group, and point colour corresponds to average expression z-scores. Genes were selected by membership in key pathways detected by gene enrichment and ontology analyses: lymphoangiogenesis, lymph vessel morphogenesis, lymph vessel development, angiogenesis, and blood vessel morphogenesis. In addition, meox1 is shown for reference.\n\n")
```



_Suggested figure legend:_ Dotplot shows the scaled expression level for genes from enriched GO terms. Units on the horizontal axis correspond to the genotypes. Point size corresponds to the percent of expression in each cell group, and point colour corresponds to average expression z-scores. Genes were selected by membership in key pathways detected by gene enrichment and ontology analyses: lymphoangiogenesis, lymph vessel morphogenesis, lymph vessel development, angiogenesis, and blood vessel morphogenesis. In addition, meox1 is shown for reference.

```{.r .fold-hide}
cat("\n\n_Caution against overinterpretation:_ Not an invalid result, but the scale may exaggerate the magnitude of difference.\n
The dotplots look very clean but theres some important nuance\n
1. the pathway genes are sourced from DEG so these would naturally be different
2. the expression values on the dotplot are z-score values and the relatively low scale indicates that the magnitude of the expression isnt massively different (but still different enough)\n\n")
```



_Caution against overinterpretation:_ Not an invalid result, but the scale may exaggerate the magnitude of difference.

The dotplots look very clean but theres some important nuance

1. the pathway genes are sourced from DEG so these would naturally be different
2. the expression values on the dotplot are z-score values and the relatively low scale indicates that the magnitude of the expression isnt massively different (but still different enough)

```{.r .fold-hide}
# plot gene expr
plt_name <- "top15_LECs_GO_pathways_downregulated_geneexpr_dotplot_curated"
ego_lecs_pathways <- c(
    "lymphangiogenesis",
    "lymph vessel morphogenesis",
    "lymph vessel development",
    "angiogenesis",
    "blood vessel morphogenesis"
    )
ego_lecs_selected_rows <- ego_lecs_top15_pathways[
    ego_lecs_top15_pathways$ego_lecs.Description %in% ego_lecs_pathways,
    ]
ego_lecs_genes_per_pathway <- setNames(
    strsplit(ego_lecs_selected_rows$ego_lecs.geneID, "/"),
    ego_lecs_selected_rows$ego_lecs.Description
    )
ego_lecs_genes_per_pathway <- c(unique(
    unlist(ego_lecs_genes_per_pathway)
    ), "meox1")
cat("\nGene quantity: ", length(ego_lecs_genes_per_pathway))
```


Gene quantity:  34

```{.r .fold-hide}
plt <- DotPlot(data_lecs, features=ego_lecs_genes_per_pathway,
    group.by="Genotype", cols=c("#d9d9d9",  "#40004b")
    ) +
    coord_flip() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
```

```
## Warning: Only one identity present, the expression values will be not scaled
```

```{.r .fold-hide}
print(plt)
```

![](meox1_goanalysis_files/figure-html/grp_lec-14.png)<!-- -->

```{.r .fold-hide}
save_plot_custom(plt, plt_name, height=12, width=4)
```

*PDF saved to:* `../output/figure_extended/go//top15_LECs_GO_pathways_downregulated_geneexpr_dotplot_curated.pdf` 

```{.r .fold-hide}
cat("\n\n_Suggested figure legend:_ Violin plot shows the log2 normalised expression level for differentially expressed genes associated with enriched GO terms. Units on the horizontal axis correspond to the genotypes. Genes were selected by membership in key pathways detected by gene enrichment and ontology analyses: lymphoangiogenesis, lymph vessel morphogenesis, lymph vessel development, angiogenesis, and blood vessel morphogenesis. In addition, meox1 is shown for reference.\n\n")
```



_Suggested figure legend:_ Violin plot shows the log2 normalised expression level for differentially expressed genes associated with enriched GO terms. Units on the horizontal axis correspond to the genotypes. Genes were selected by membership in key pathways detected by gene enrichment and ontology analyses: lymphoangiogenesis, lymph vessel morphogenesis, lymph vessel development, angiogenesis, and blood vessel morphogenesis. In addition, meox1 is shown for reference.

```{.r .fold-hide}
cat("\n_Tyrone: Not every single canonical gene associated with the pathway will be shown. It is a function of the data, only the differentially expressed genes associated with the GO term will be shown._\n\n")
```


_Tyrone: Not every single canonical gene associated with the pathway will be shown. It is a function of the data, only the differentially expressed genes associated with the GO term will be shown._

```{.r .fold-hide}
plt_name <- "top15_LECs_GO_pathways_downregulated_geneexpr_vlnplot_curated"
plt <- VlnPlot(
    data_lecs,
    features = rev(ego_lecs_genes_per_pathway),
    group.by = "Genotype",
    pt.size = 1,
    ncol = 8,
    flip = TRUE,
    cols = c("#dbe2c6", "#657c95")
) + scale_fill_manual(values = c("#dbe2c6", "#657c95")) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
```

```
## Scale for fill is already present.
## Adding another scale for fill, which will replace the existing scale.
```

```{.r .fold-hide}
print(plt)
```

![](meox1_goanalysis_files/figure-html/grp_lec-15.png)<!-- -->

```{.r .fold-hide}
save_plot_custom(plt, plt_name, height=16, width=16)
```

*PDF saved to:* `../output/figure_extended/go//top15_LECs_GO_pathways_downregulated_geneexpr_vlnplot_curated.pdf` 


```{.r .fold-hide}
singles_LEC_out <- run_go(singles_LEC, data_name="LEC")
```



## GO analysis: LEC 

Input genes: 15771 

**ORA input genes:** combined = 1499 &nbsp;&nbsp; up = 944 &nbsp;&nbsp; down = 555 

**ORA (combined):** 2559 terms &nbsp;&nbsp;  **ORA (up):** 2116 terms &nbsp;&nbsp;  **ORA (down):** 1742 terms

*Full ORA ( combined ) table saved to:* `../output/figure_extended/go//LEC_ego_combined_full_results.csv` 

*Full ORA ( up ) table saved to:* `../output/figure_extended/go//LEC_ego_up_full_results.csv` 

*Full ORA ( down ) table saved to:* `../output/figure_extended/go//LEC_ego_down_full_results.csv` 

**ORA significant terms:** combined = 128 , up = 124 , down = 15 

*Significant-only ORA ( combined ) table saved to:* `../output/figure_extended/go//LEC_ego_combined_significant_only.csv` 

*Significant-only ORA ( up ) table saved to:* `../output/figure_extended/go//LEC_ego_up_significant_only.csv` 

*Significant-only ORA ( down ) table saved to:* `../output/figure_extended/go//LEC_ego_down_significant_only.csv` 

### ORA ( combined ): Barplot

Shows the top enriched GO terms ranked by significance, bar length/colour reflecting  adjusted p-value or gene count. **Expected:** terms biologically consistent with the  comparison (e.g. vascular/endothelial terms for EC-related contrasts). **Unexpected:**  dominance by generic housekeeping terms (translation, metabolism) with no thematic  coherence - may indicate a weak or noisy input gene list rather than a real signal.

```
## Warning in fortify(object, showCategory = showCategory, by = x, ...): Arguments in `...` must be used.
## ✖ Problematic argument:
## • by = x
## ℹ Did you misspell an argument name?
```

![](meox1_goanalysis_files/figure-html/sin_LEC-1.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//LEC_ego_combined_barplot.pdf` 

### ORA ( combined ): Dotplot

Adds gene ratio (fraction of input genes hitting each term) as point size, with colour  for significance - more informative than the barplot since it shows both effect size and  confidence simultaneously. **Expected:** larger, more significant dots clustering around a  coherent biological theme. **Unexpected:** significant terms with very small gene ratios  (few genes driving a 'significant' result) - treat cautiously as potentially fragile.

![](meox1_goanalysis_files/figure-html/sin_LEC-2.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//LEC_ego_combined_dotplot.pdf` 

### ORA ( combined ): Enrichment map

Network graph connecting GO terms that share overlapping genes - clusters of connected  terms suggest a shared underlying biological process rather than independent findings.  **Expected:** tight clusters of related terms. **Unexpected:** many isolated, unconnected  terms, suggesting the significant terms are thematically scattered rather than coherent.

![](meox1_goanalysis_files/figure-html/sin_LEC-3.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//LEC_ego_combined_emapplot.pdf` 

### ORA ( combined ): Category-gene network

Shows which specific genes drive which enriched terms, useful for tracing a term back to  the actual DE genes behind it. **Expected:** recognisable marker genes linking to expected  terms. **Unexpected:** a term driven entirely by one or two genes rather than a broad gene  set - weaker evidence than a term supported by many genes.

![](meox1_goanalysis_files/figure-html/sin_LEC-4.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//LEC_ego_combined_cnetplot.pdf` 

### ORA ( up ): Barplot

Shows the top enriched GO terms ranked by significance, bar length/colour reflecting  adjusted p-value or gene count. **Expected:** terms biologically consistent with the  comparison (e.g. vascular/endothelial terms for EC-related contrasts). **Unexpected:**  dominance by generic housekeeping terms (translation, metabolism) with no thematic  coherence - may indicate a weak or noisy input gene list rather than a real signal.

```
## Warning in fortify(object, showCategory = showCategory, by = x, ...): Arguments in `...` must be used.
## ✖ Problematic argument:
## • by = x
## ℹ Did you misspell an argument name?
```

![](meox1_goanalysis_files/figure-html/sin_LEC-5.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//LEC_ego_up_barplot.pdf` 

### ORA ( up ): Dotplot

Adds gene ratio (fraction of input genes hitting each term) as point size, with colour  for significance - more informative than the barplot since it shows both effect size and  confidence simultaneously. **Expected:** larger, more significant dots clustering around a  coherent biological theme. **Unexpected:** significant terms with very small gene ratios  (few genes driving a 'significant' result) - treat cautiously as potentially fragile.

![](meox1_goanalysis_files/figure-html/sin_LEC-6.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//LEC_ego_up_dotplot.pdf` 

### ORA ( up ): Enrichment map

Network graph connecting GO terms that share overlapping genes - clusters of connected  terms suggest a shared underlying biological process rather than independent findings.  **Expected:** tight clusters of related terms. **Unexpected:** many isolated, unconnected  terms, suggesting the significant terms are thematically scattered rather than coherent.

![](meox1_goanalysis_files/figure-html/sin_LEC-7.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//LEC_ego_up_emapplot.pdf` 

### ORA ( up ): Category-gene network

Shows which specific genes drive which enriched terms, useful for tracing a term back to  the actual DE genes behind it. **Expected:** recognisable marker genes linking to expected  terms. **Unexpected:** a term driven entirely by one or two genes rather than a broad gene  set - weaker evidence than a term supported by many genes.

![](meox1_goanalysis_files/figure-html/sin_LEC-8.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//LEC_ego_up_cnetplot.pdf` 

### ORA ( down ): Barplot

Shows the top enriched GO terms ranked by significance, bar length/colour reflecting  adjusted p-value or gene count. **Expected:** terms biologically consistent with the  comparison (e.g. vascular/endothelial terms for EC-related contrasts). **Unexpected:**  dominance by generic housekeeping terms (translation, metabolism) with no thematic  coherence - may indicate a weak or noisy input gene list rather than a real signal.

```
## Warning in fortify(object, showCategory = showCategory, by = x, ...): Arguments in `...` must be used.
## ✖ Problematic argument:
## • by = x
## ℹ Did you misspell an argument name?
```

![](meox1_goanalysis_files/figure-html/sin_LEC-9.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//LEC_ego_down_barplot.pdf` 

### ORA ( down ): Dotplot

Adds gene ratio (fraction of input genes hitting each term) as point size, with colour  for significance - more informative than the barplot since it shows both effect size and  confidence simultaneously. **Expected:** larger, more significant dots clustering around a  coherent biological theme. **Unexpected:** significant terms with very small gene ratios  (few genes driving a 'significant' result) - treat cautiously as potentially fragile.

![](meox1_goanalysis_files/figure-html/sin_LEC-10.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//LEC_ego_down_dotplot.pdf` 

### ORA ( down ): Enrichment map

Network graph connecting GO terms that share overlapping genes - clusters of connected  terms suggest a shared underlying biological process rather than independent findings.  **Expected:** tight clusters of related terms. **Unexpected:** many isolated, unconnected  terms, suggesting the significant terms are thematically scattered rather than coherent.

![](meox1_goanalysis_files/figure-html/sin_LEC-11.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//LEC_ego_down_emapplot.pdf` 

### ORA ( down ): Category-gene network

Shows which specific genes drive which enriched terms, useful for tracing a term back to  the actual DE genes behind it. **Expected:** recognisable marker genes linking to expected  terms. **Unexpected:** a term driven entirely by one or two genes rather than a broad gene  set - weaker evidence than a term supported by many genes.

![](meox1_goanalysis_files/figure-html/sin_LEC-12.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//LEC_ego_down_cnetplot.pdf` 

**Done:** LEC 

---

```{.r .fold-hide}
cat("### Custom plots for downregulated genes in mutant vs wildtype\n\n")
```

### Custom plots for downregulated genes in mutant vs wildtype

```{.r .fold-hide}
ego_lec <- singles_LEC_out$ego_result_down@result
ego_lec$GeneRatio_numeric <- sapply(
    strsplit(ego_lec$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])
    )
ego_lec_top15_terms <- ego_lec[
    order(-ego_lec$GeneRatio_numeric), ][1:15,
    ]
ego_lec_top15_gene_map <- setNames(strsplit(
    ego_lec_top15_terms$geneID, "/"), ego_lec_top15_terms$Description
    )
ego_lec_top15_gene_list <- unique(unlist(unname(ego_lec_top15_gene_map)))

ego_lec_top15_pathways <- data.frame(ego_lec$ID, ego_lec$Description, ego_lec$geneID)

ego_lec_top15_gene_list <- ego_lec_top15_gene_list[
    grep(":", ego_lec_top15_gene_list, invert=TRUE)
    ]

cat("\n\nDotplots of gene ratio top 15 only.\n\n")
```



Dotplots of gene ratio top 15 only.

```{.r .fold-hide}
cat("Adds gene ratio (fraction of input genes hitting each term) as point size, with colour ",
    "for significance - more informative than the barplot since it shows both effect size and ",
    "confidence simultaneously. **Expected:** larger, more significant dots clustering around a ",
    "coherent biological theme. **Unexpected:** significant terms with very small gene ratios ",
    "(few genes driving a 'significant' result) - treat cautiously as potentially fragile.\n\n")
```

Adds gene ratio (fraction of input genes hitting each term) as point size, with colour  for significance - more informative than the barplot since it shows both effect size and  confidence simultaneously. **Expected:** larger, more significant dots clustering around a  coherent biological theme. **Unexpected:** significant terms with very small gene ratios  (few genes driving a 'significant' result) - treat cautiously as potentially fragile.

```{.r .fold-hide}
cat("\n\n_Suggested figure legend:_ Dotplot shows the effect size and confidence simultaneously for enriched GO terms. Units on the horizontal axis correspond to the fraction of input genes which hit each GO term. Point size corresponds to the quantity of terms, and point colour corresponds to adjusted p values. Top 15 GO terms by gene ratio are shown.\n\n")
```



_Suggested figure legend:_ Dotplot shows the effect size and confidence simultaneously for enriched GO terms. Units on the horizontal axis correspond to the fraction of input genes which hit each GO term. Point size corresponds to the quantity of terms, and point colour corresponds to adjusted p values. Top 15 GO terms by gene ratio are shown.

```{.r .fold-hide}
plt_name <- "top15_LEC_GO_pathways_downregulated_generatio_dotplot"
plt <- dotplot(singles_LEC_out$ego_result_down, showCategory=15)
print(plt)
```

![](meox1_goanalysis_files/figure-html/sin_LEC-13.png)<!-- -->

```{.r .fold-hide}
save_plot_custom(plt, plt_name)
```

*PDF saved to:* `../output/figure_extended/go//top15_LEC_GO_pathways_downregulated_generatio_dotplot.pdf` 

```{.r .fold-hide}
data_lec <- data[, data$L3_celltype == "LEC"]
data_lec$Genotype <- gsub("meox1_mutant", "mutant", data_lec$Genotype)
data_lec$Genotype <- gsub("wildtype", "sibling", data_lec$Genotype)
data_lec$Genotype <- factor(
    data_lec$Genotype,
    levels = c("sibling", "mutant")
)

cat("\n\n**Genotype counts (cells) for LEC singles:**",
    paste(names(table(data_lec$Genotype)), table(data_lec$Genotype), sep = " = ", collapse = ", "),
    "\n\n")
```



**Genotype counts (cells) for LEC singles:** sibling = 285, mutant = 231 

```{.r .fold-hide}
cat("\n\nDotplots of gene expression (key pathways).\n\n")
```



Dotplots of gene expression (key pathways).

```{.r .fold-hide}
cat("For plot expansion, settled on these:\n\n
1. selected top 15 pathways as a cutoff\n
2. extracted the genes associated with these\n
3. removed genes with the `[a-z]+:` prefix\n
4. selected genes associated with key pathways and added meox1\n
5. plotted dotplots on the remainder\n\n

Attempted heatmaps, not very useful.\n\n")
```

For plot expansion, settled on these:


1. selected top 15 pathways as a cutoff

2. extracted the genes associated with these

3. removed genes with the `[a-z]+:` prefix

4. selected genes associated with key pathways and added meox1

5. plotted dotplots on the remainder



Attempted heatmaps, not very useful.

```{.r .fold-hide}
cat("\n\n_Suggested figure legend:_ Dotplot shows the scaled expression level for available genes from enriched GO terms. Units on the horizontal axis correspond to the genotypes. Point size corresponds to the percent of expression in each cell group, and point colour corresponds to average expression z-scores. Genes were selected by membership in key pathways detected by gene enrichment and ontology analyses: lymphoangiogenesis, lymph vessel morphogenesis, lymph vessel development, angiogenesis, and blood vessel morphogenesis. In addition, meox1 is shown for reference.\n\n")
```



_Suggested figure legend:_ Dotplot shows the scaled expression level for available genes from enriched GO terms. Units on the horizontal axis correspond to the genotypes. Point size corresponds to the percent of expression in each cell group, and point colour corresponds to average expression z-scores. Genes were selected by membership in key pathways detected by gene enrichment and ontology analyses: lymphoangiogenesis, lymph vessel morphogenesis, lymph vessel development, angiogenesis, and blood vessel morphogenesis. In addition, meox1 is shown for reference.

```{.r .fold-hide}
cat("\n\n_Caution against overinterpretation:_ Not an invalid result, but the scale may exaggerate the magnitude of difference.\n
The dotplots look very clean but theres some important nuance\n
1. the pathway genes are sourced from DEG so these would naturally be different
2. the expression values on the dotplot are z-score values and the relatively low scale indicates that the magnitude of the expression isnt massively different (but still different enough)\n\n")
```



_Caution against overinterpretation:_ Not an invalid result, but the scale may exaggerate the magnitude of difference.

The dotplots look very clean but theres some important nuance

1. the pathway genes are sourced from DEG so these would naturally be different
2. the expression values on the dotplot are z-score values and the relatively low scale indicates that the magnitude of the expression isnt massively different (but still different enough)

```{.r .fold-hide}
# plot gene expr
plt_name <- "top15_LEC_GO_pathways_downregulated_geneexpr_dotplot_curated"
ego_lec_pathways <- c(
    "lymphangiogenesis",
    "lymph vessel morphogenesis",
    "lymph vessel development",
    "angiogenesis",
    "blood vessel morphogenesis"
    )
ego_lec_selected_rows <- ego_lec_top15_pathways[
    ego_lec_top15_pathways$ego_lec.Description %in% ego_lec_pathways,
    ]
ego_lec_genes_per_pathway <- setNames(
    strsplit(ego_lec_selected_rows$ego_lec.geneID, "/"),
    ego_lec_selected_rows$ego_lec.Description
    )
ego_lec_genes_per_pathway <- c(unique(
    unlist(ego_lec_genes_per_pathway)
    ), "meox1")
cat("\nGene quantity: ", length(ego_lec_genes_per_pathway))
```


Gene quantity:  33

```{.r .fold-hide}
plt <- DotPlot(data_lec, features=ego_lec_genes_per_pathway,
    group.by="Genotype", cols=c("#d9d9d9",  "#40004b")
    ) +
    coord_flip() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
```

```
## Warning: Scaling data with a low number of groups may produce misleading
## results
```

```{.r .fold-hide}
print(plt)
```

![](meox1_goanalysis_files/figure-html/sin_LEC-14.png)<!-- -->

```{.r .fold-hide}
save_plot_custom(plt, plt_name, height=12, width=4)
```

*PDF saved to:* `../output/figure_extended/go//top15_LEC_GO_pathways_downregulated_geneexpr_dotplot_curated.pdf` 

```{.r .fold-hide}
cat("\n\n_Suggested figure legend:_ Violin plot shows the log2 normalised expression level for differentially expressed genes associated with enriched GO terms. Units on the horizontal axis correspond to the genotypes. Genes were selected by membership in key pathways detected by gene enrichment and ontology analyses: lymphoangiogenesis, lymph vessel morphogenesis, lymph vessel development, angiogenesis, and blood vessel morphogenesis. In addition, meox1 is shown for reference.\n\n")
```



_Suggested figure legend:_ Violin plot shows the log2 normalised expression level for differentially expressed genes associated with enriched GO terms. Units on the horizontal axis correspond to the genotypes. Genes were selected by membership in key pathways detected by gene enrichment and ontology analyses: lymphoangiogenesis, lymph vessel morphogenesis, lymph vessel development, angiogenesis, and blood vessel morphogenesis. In addition, meox1 is shown for reference.

```{.r .fold-hide}
cat("\n_Tyrone: Not every single canonical gene associated with the pathway will be shown. It is a function of the data, only the differentially expressed genes associated with the GO term will be shown._\n\n")
```


_Tyrone: Not every single canonical gene associated with the pathway will be shown. It is a function of the data, only the differentially expressed genes associated with the GO term will be shown._

```{.r .fold-hide}
plt_name <- "top15_LEC_GO_pathways_downregulated_geneexpr_vlnplot_curated"
plt <- VlnPlot(
    data_lec,
    features = rev(ego_lec_genes_per_pathway),
    group.by = "Genotype",
    pt.size = 1,
    ncol = 8,
    flip = TRUE,
    cols = c("#dbe2c6", "#657c95")
) + scale_fill_manual(values = c("#dbe2c6", "#657c95")) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
```

```
## Scale for fill is already present.
## Adding another scale for fill, which will replace the existing scale.
```

```{.r .fold-hide}
print(plt)
```

![](meox1_goanalysis_files/figure-html/sin_LEC-15.png)<!-- -->

```{.r .fold-hide}
save_plot_custom(plt, plt_name, height=12, width=16)
```

*PDF saved to:* `../output/figure_extended/go//top15_LEC_GO_pathways_downregulated_geneexpr_vlnplot_curated.pdf` 


```{.r .fold-hide}
singles_pre_muLEC_out <- run_go(singles_pre_muLEC, data_name="pre_muLEC")
```



## GO analysis: pre_muLEC 

Input genes: 16252 

**ORA input genes:** combined = 753 &nbsp;&nbsp; up = 426 &nbsp;&nbsp; down = 327 

**ORA (combined):** 2006 terms &nbsp;&nbsp;  **ORA (up):** 1582 terms &nbsp;&nbsp;  **ORA (down):** 1222 terms

*Full ORA ( combined ) table saved to:* `../output/figure_extended/go//pre_muLEC_ego_combined_full_results.csv` 

*Full ORA ( up ) table saved to:* `../output/figure_extended/go//pre_muLEC_ego_up_full_results.csv` 

*Full ORA ( down ) table saved to:* `../output/figure_extended/go//pre_muLEC_ego_down_full_results.csv` 

**ORA significant terms:** combined = 48 , up = 64 , down = 4 

*Significant-only ORA ( combined ) table saved to:* `../output/figure_extended/go//pre_muLEC_ego_combined_significant_only.csv` 

*Significant-only ORA ( up ) table saved to:* `../output/figure_extended/go//pre_muLEC_ego_up_significant_only.csv` 

*Significant-only ORA ( down ) table saved to:* `../output/figure_extended/go//pre_muLEC_ego_down_significant_only.csv` 

### ORA ( combined ): Barplot

Shows the top enriched GO terms ranked by significance, bar length/colour reflecting  adjusted p-value or gene count. **Expected:** terms biologically consistent with the  comparison (e.g. vascular/endothelial terms for EC-related contrasts). **Unexpected:**  dominance by generic housekeeping terms (translation, metabolism) with no thematic  coherence - may indicate a weak or noisy input gene list rather than a real signal.

```
## Warning in fortify(object, showCategory = showCategory, by = x, ...): Arguments in `...` must be used.
## ✖ Problematic argument:
## • by = x
## ℹ Did you misspell an argument name?
```

![](meox1_goanalysis_files/figure-html/sin_pre_muLEC-1.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//pre_muLEC_ego_combined_barplot.pdf` 

### ORA ( combined ): Dotplot

Adds gene ratio (fraction of input genes hitting each term) as point size, with colour  for significance - more informative than the barplot since it shows both effect size and  confidence simultaneously. **Expected:** larger, more significant dots clustering around a  coherent biological theme. **Unexpected:** significant terms with very small gene ratios  (few genes driving a 'significant' result) - treat cautiously as potentially fragile.

![](meox1_goanalysis_files/figure-html/sin_pre_muLEC-2.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//pre_muLEC_ego_combined_dotplot.pdf` 

### ORA ( combined ): Enrichment map

Network graph connecting GO terms that share overlapping genes - clusters of connected  terms suggest a shared underlying biological process rather than independent findings.  **Expected:** tight clusters of related terms. **Unexpected:** many isolated, unconnected  terms, suggesting the significant terms are thematically scattered rather than coherent.

![](meox1_goanalysis_files/figure-html/sin_pre_muLEC-3.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//pre_muLEC_ego_combined_emapplot.pdf` 

### ORA ( combined ): Category-gene network

Shows which specific genes drive which enriched terms, useful for tracing a term back to  the actual DE genes behind it. **Expected:** recognisable marker genes linking to expected  terms. **Unexpected:** a term driven entirely by one or two genes rather than a broad gene  set - weaker evidence than a term supported by many genes.

![](meox1_goanalysis_files/figure-html/sin_pre_muLEC-4.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//pre_muLEC_ego_combined_cnetplot.pdf` 

### ORA ( up ): Barplot

Shows the top enriched GO terms ranked by significance, bar length/colour reflecting  adjusted p-value or gene count. **Expected:** terms biologically consistent with the  comparison (e.g. vascular/endothelial terms for EC-related contrasts). **Unexpected:**  dominance by generic housekeeping terms (translation, metabolism) with no thematic  coherence - may indicate a weak or noisy input gene list rather than a real signal.

```
## Warning in fortify(object, showCategory = showCategory, by = x, ...): Arguments in `...` must be used.
## ✖ Problematic argument:
## • by = x
## ℹ Did you misspell an argument name?
```

![](meox1_goanalysis_files/figure-html/sin_pre_muLEC-5.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//pre_muLEC_ego_up_barplot.pdf` 

### ORA ( up ): Dotplot

Adds gene ratio (fraction of input genes hitting each term) as point size, with colour  for significance - more informative than the barplot since it shows both effect size and  confidence simultaneously. **Expected:** larger, more significant dots clustering around a  coherent biological theme. **Unexpected:** significant terms with very small gene ratios  (few genes driving a 'significant' result) - treat cautiously as potentially fragile.

![](meox1_goanalysis_files/figure-html/sin_pre_muLEC-6.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//pre_muLEC_ego_up_dotplot.pdf` 

### ORA ( up ): Enrichment map

Network graph connecting GO terms that share overlapping genes - clusters of connected  terms suggest a shared underlying biological process rather than independent findings.  **Expected:** tight clusters of related terms. **Unexpected:** many isolated, unconnected  terms, suggesting the significant terms are thematically scattered rather than coherent.

![](meox1_goanalysis_files/figure-html/sin_pre_muLEC-7.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//pre_muLEC_ego_up_emapplot.pdf` 

### ORA ( up ): Category-gene network

Shows which specific genes drive which enriched terms, useful for tracing a term back to  the actual DE genes behind it. **Expected:** recognisable marker genes linking to expected  terms. **Unexpected:** a term driven entirely by one or two genes rather than a broad gene  set - weaker evidence than a term supported by many genes.

![](meox1_goanalysis_files/figure-html/sin_pre_muLEC-8.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//pre_muLEC_ego_up_cnetplot.pdf` 

### ORA ( down ): Barplot

Shows the top enriched GO terms ranked by significance, bar length/colour reflecting  adjusted p-value or gene count. **Expected:** terms biologically consistent with the  comparison (e.g. vascular/endothelial terms for EC-related contrasts). **Unexpected:**  dominance by generic housekeeping terms (translation, metabolism) with no thematic  coherence - may indicate a weak or noisy input gene list rather than a real signal.

```
## Warning in fortify(object, showCategory = showCategory, by = x, ...): Arguments in `...` must be used.
## ✖ Problematic argument:
## • by = x
## ℹ Did you misspell an argument name?
```

![](meox1_goanalysis_files/figure-html/sin_pre_muLEC-9.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//pre_muLEC_ego_down_barplot.pdf` 

### ORA ( down ): Dotplot

Adds gene ratio (fraction of input genes hitting each term) as point size, with colour  for significance - more informative than the barplot since it shows both effect size and  confidence simultaneously. **Expected:** larger, more significant dots clustering around a  coherent biological theme. **Unexpected:** significant terms with very small gene ratios  (few genes driving a 'significant' result) - treat cautiously as potentially fragile.

![](meox1_goanalysis_files/figure-html/sin_pre_muLEC-10.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//pre_muLEC_ego_down_dotplot.pdf` 

### ORA ( down ): Enrichment map

Network graph connecting GO terms that share overlapping genes - clusters of connected  terms suggest a shared underlying biological process rather than independent findings.  **Expected:** tight clusters of related terms. **Unexpected:** many isolated, unconnected  terms, suggesting the significant terms are thematically scattered rather than coherent.

![](meox1_goanalysis_files/figure-html/sin_pre_muLEC-11.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//pre_muLEC_ego_down_emapplot.pdf` 

### ORA ( down ): Category-gene network

Shows which specific genes drive which enriched terms, useful for tracing a term back to  the actual DE genes behind it. **Expected:** recognisable marker genes linking to expected  terms. **Unexpected:** a term driven entirely by one or two genes rather than a broad gene  set - weaker evidence than a term supported by many genes.

![](meox1_goanalysis_files/figure-html/sin_pre_muLEC-12.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//pre_muLEC_ego_down_cnetplot.pdf` 

**Done:** pre_muLEC 

---


```{.r .fold-hide}
grouped_vecs_out <- run_go(grouped_vecs, data_name="hmVEC__mVEC")
```



## GO analysis: hmVEC__mVEC 

Input genes: 14057 

**ORA input genes:** combined = 2818 &nbsp;&nbsp; up = 1575 &nbsp;&nbsp; down = 1243 

**ORA (combined):** 2894 terms &nbsp;&nbsp;  **ORA (up):** 2565 terms &nbsp;&nbsp;  **ORA (down):** 2315 terms

*Full ORA ( combined ) table saved to:* `../output/figure_extended/go//hmVEC__mVEC_ego_combined_full_results.csv` 

*Full ORA ( up ) table saved to:* `../output/figure_extended/go//hmVEC__mVEC_ego_up_full_results.csv` 

*Full ORA ( down ) table saved to:* `../output/figure_extended/go//hmVEC__mVEC_ego_down_full_results.csv` 

**ORA significant terms:** combined = 200 , up = 160 , down = 174 

*Significant-only ORA ( combined ) table saved to:* `../output/figure_extended/go//hmVEC__mVEC_ego_combined_significant_only.csv` 

*Significant-only ORA ( up ) table saved to:* `../output/figure_extended/go//hmVEC__mVEC_ego_up_significant_only.csv` 

*Significant-only ORA ( down ) table saved to:* `../output/figure_extended/go//hmVEC__mVEC_ego_down_significant_only.csv` 

### ORA ( combined ): Barplot

Shows the top enriched GO terms ranked by significance, bar length/colour reflecting  adjusted p-value or gene count. **Expected:** terms biologically consistent with the  comparison (e.g. vascular/endothelial terms for EC-related contrasts). **Unexpected:**  dominance by generic housekeeping terms (translation, metabolism) with no thematic  coherence - may indicate a weak or noisy input gene list rather than a real signal.

```
## Warning in fortify(object, showCategory = showCategory, by = x, ...): Arguments in `...` must be used.
## ✖ Problematic argument:
## • by = x
## ℹ Did you misspell an argument name?
```

![](meox1_goanalysis_files/figure-html/grp_vec-1.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//hmVEC__mVEC_ego_combined_barplot.pdf` 

### ORA ( combined ): Dotplot

Adds gene ratio (fraction of input genes hitting each term) as point size, with colour  for significance - more informative than the barplot since it shows both effect size and  confidence simultaneously. **Expected:** larger, more significant dots clustering around a  coherent biological theme. **Unexpected:** significant terms with very small gene ratios  (few genes driving a 'significant' result) - treat cautiously as potentially fragile.

![](meox1_goanalysis_files/figure-html/grp_vec-2.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//hmVEC__mVEC_ego_combined_dotplot.pdf` 

### ORA ( combined ): Enrichment map

Network graph connecting GO terms that share overlapping genes - clusters of connected  terms suggest a shared underlying biological process rather than independent findings.  **Expected:** tight clusters of related terms. **Unexpected:** many isolated, unconnected  terms, suggesting the significant terms are thematically scattered rather than coherent.

![](meox1_goanalysis_files/figure-html/grp_vec-3.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//hmVEC__mVEC_ego_combined_emapplot.pdf` 

### ORA ( combined ): Category-gene network

Shows which specific genes drive which enriched terms, useful for tracing a term back to  the actual DE genes behind it. **Expected:** recognisable marker genes linking to expected  terms. **Unexpected:** a term driven entirely by one or two genes rather than a broad gene  set - weaker evidence than a term supported by many genes.

![](meox1_goanalysis_files/figure-html/grp_vec-4.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//hmVEC__mVEC_ego_combined_cnetplot.pdf` 

### ORA ( up ): Barplot

Shows the top enriched GO terms ranked by significance, bar length/colour reflecting  adjusted p-value or gene count. **Expected:** terms biologically consistent with the  comparison (e.g. vascular/endothelial terms for EC-related contrasts). **Unexpected:**  dominance by generic housekeeping terms (translation, metabolism) with no thematic  coherence - may indicate a weak or noisy input gene list rather than a real signal.

```
## Warning in fortify(object, showCategory = showCategory, by = x, ...): Arguments in `...` must be used.
## ✖ Problematic argument:
## • by = x
## ℹ Did you misspell an argument name?
```

![](meox1_goanalysis_files/figure-html/grp_vec-5.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//hmVEC__mVEC_ego_up_barplot.pdf` 

### ORA ( up ): Dotplot

Adds gene ratio (fraction of input genes hitting each term) as point size, with colour  for significance - more informative than the barplot since it shows both effect size and  confidence simultaneously. **Expected:** larger, more significant dots clustering around a  coherent biological theme. **Unexpected:** significant terms with very small gene ratios  (few genes driving a 'significant' result) - treat cautiously as potentially fragile.

![](meox1_goanalysis_files/figure-html/grp_vec-6.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//hmVEC__mVEC_ego_up_dotplot.pdf` 

### ORA ( up ): Enrichment map

Network graph connecting GO terms that share overlapping genes - clusters of connected  terms suggest a shared underlying biological process rather than independent findings.  **Expected:** tight clusters of related terms. **Unexpected:** many isolated, unconnected  terms, suggesting the significant terms are thematically scattered rather than coherent.

![](meox1_goanalysis_files/figure-html/grp_vec-7.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//hmVEC__mVEC_ego_up_emapplot.pdf` 

### ORA ( up ): Category-gene network

Shows which specific genes drive which enriched terms, useful for tracing a term back to  the actual DE genes behind it. **Expected:** recognisable marker genes linking to expected  terms. **Unexpected:** a term driven entirely by one or two genes rather than a broad gene  set - weaker evidence than a term supported by many genes.

![](meox1_goanalysis_files/figure-html/grp_vec-8.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//hmVEC__mVEC_ego_up_cnetplot.pdf` 

### ORA ( down ): Barplot

Shows the top enriched GO terms ranked by significance, bar length/colour reflecting  adjusted p-value or gene count. **Expected:** terms biologically consistent with the  comparison (e.g. vascular/endothelial terms for EC-related contrasts). **Unexpected:**  dominance by generic housekeeping terms (translation, metabolism) with no thematic  coherence - may indicate a weak or noisy input gene list rather than a real signal.

```
## Warning in fortify(object, showCategory = showCategory, by = x, ...): Arguments in `...` must be used.
## ✖ Problematic argument:
## • by = x
## ℹ Did you misspell an argument name?
```

![](meox1_goanalysis_files/figure-html/grp_vec-9.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//hmVEC__mVEC_ego_down_barplot.pdf` 

### ORA ( down ): Dotplot

Adds gene ratio (fraction of input genes hitting each term) as point size, with colour  for significance - more informative than the barplot since it shows both effect size and  confidence simultaneously. **Expected:** larger, more significant dots clustering around a  coherent biological theme. **Unexpected:** significant terms with very small gene ratios  (few genes driving a 'significant' result) - treat cautiously as potentially fragile.

![](meox1_goanalysis_files/figure-html/grp_vec-10.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//hmVEC__mVEC_ego_down_dotplot.pdf` 

### ORA ( down ): Enrichment map

Network graph connecting GO terms that share overlapping genes - clusters of connected  terms suggest a shared underlying biological process rather than independent findings.  **Expected:** tight clusters of related terms. **Unexpected:** many isolated, unconnected  terms, suggesting the significant terms are thematically scattered rather than coherent.

![](meox1_goanalysis_files/figure-html/grp_vec-11.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//hmVEC__mVEC_ego_down_emapplot.pdf` 

### ORA ( down ): Category-gene network

Shows which specific genes drive which enriched terms, useful for tracing a term back to  the actual DE genes behind it. **Expected:** recognisable marker genes linking to expected  terms. **Unexpected:** a term driven entirely by one or two genes rather than a broad gene  set - weaker evidence than a term supported by many genes.

![](meox1_goanalysis_files/figure-html/grp_vec-12.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//hmVEC__mVEC_ego_down_cnetplot.pdf` 

**Done:** hmVEC__mVEC 

---

```{.r .fold-hide}
cat("### Custom plots for downregulated genes in mutant vs wildtype\n\n")
```

### Custom plots for downregulated genes in mutant vs wildtype

```{.r .fold-hide}
ego_vecs <- grouped_vecs_out$ego_result_down@result
ego_vecs$GeneRatio_numeric <- sapply(
    strsplit(ego_vecs$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])
    )
ego_vecs_top15_terms <- ego_vecs[
    order(-ego_vecs$GeneRatio_numeric), ][1:15,
    ]
ego_vecs_top15_gene_map <- setNames(strsplit(
    ego_vecs_top15_terms$geneID, "/"), ego_vecs_top15_terms$Description
    )
ego_vecs_top15_gene_list <- unique(unlist(unname(ego_vecs_top15_gene_map)))

ego_vecs_top15_pathways <- data.frame(ego_vecs$ID, ego_vecs$Description, ego_vecs$geneID)

ego_vecs_top15_gene_list <- ego_vecs_top15_gene_list[
    grep(":", ego_vecs_top15_gene_list, invert=TRUE)
    ]

cat("\n\nDotplots of gene ratio top 15 only.\n\n")
```



Dotplots of gene ratio top 15 only.

```{.r .fold-hide}
cat("Adds gene ratio (fraction of input genes hitting each term) as point size, with colour ",
    "for significance - more informative than the barplot since it shows both effect size and ",
    "confidence simultaneously. **Expected:** larger, more significant dots clustering around a ",
    "coherent biological theme. **Unexpected:** significant terms with very small gene ratios ",
    "(few genes driving a 'significant' result) - treat cautiously as potentially fragile.\n\n")
```

Adds gene ratio (fraction of input genes hitting each term) as point size, with colour  for significance - more informative than the barplot since it shows both effect size and  confidence simultaneously. **Expected:** larger, more significant dots clustering around a  coherent biological theme. **Unexpected:** significant terms with very small gene ratios  (few genes driving a 'significant' result) - treat cautiously as potentially fragile.

```{.r .fold-hide}
cat("\n\n_Suggested figure legend:_ Dotplot shows the effect size and confidence simultaneously for differentially expressed genes associated with enriched GO terms. Units on the horizontal axis correspond to the fraction of input genes which hit each GO term. Point size corresponds to the quantity of terms, and point colour corresponds to adjusted p values. Top 15 GO terms by gene ratio are shown.\n\n")
```



_Suggested figure legend:_ Dotplot shows the effect size and confidence simultaneously for differentially expressed genes associated with enriched GO terms. Units on the horizontal axis correspond to the fraction of input genes which hit each GO term. Point size corresponds to the quantity of terms, and point colour corresponds to adjusted p values. Top 15 GO terms by gene ratio are shown.

```{.r .fold-hide}
plt_name <- "top15_VECs_GO_pathways_downregulated_generatio_dotplot"
plt <- dotplot(grouped_vecs_out$ego_result_down, showCategory=15)
print(plt)
```

![](meox1_goanalysis_files/figure-html/grp_vec-13.png)<!-- -->

```{.r .fold-hide}
save_plot_custom(plt, plt_name)
```

*PDF saved to:* `../output/figure_extended/go//top15_VECs_GO_pathways_downregulated_generatio_dotplot.pdf` 

```{.r .fold-hide}
data_vecs <- data[, (data$L3_celltype == "hmVEC" | data$L3_celltype == "mVEC")]
data_vecs$Genotype <- gsub("meox1_mutant", "mutant", data_vecs$Genotype)
data_vecs$Genotype <- gsub("wildtype", "sibling", data_vecs$Genotype)
data_vecs$Genotype <- factor(
    data_vecs$Genotype,
    levels = c("sibling", "mutant")
)

cat("\n\n**Genotype counts (cells) for VECs grouped:**",
    paste(names(table(data_vecs$Genotype)), table(data_vecs$Genotype), sep = " = ", collapse = ", "),
    "\n\n")
```



**Genotype counts (cells) for VECs grouped:** sibling = 925, mutant = 964 

```{.r .fold-hide}
cat("\n\nDotplots of gene expression (key pathways).\n\n")
```



Dotplots of gene expression (key pathways).

```{.r .fold-hide}
cat("For plot expansion, settled on these:\n\n
1. selected top 15 pathways as a cutoff\n
2. extracted the genes associated with these\n
3. removed genes with the `[a-z]+:` prefix\n
4. selected genes associated with key pathways and added meox1\n
5. plotted dotplots on the remainder\n\n

Attempted heatmaps, not very useful.\n\n")
```

For plot expansion, settled on these:


1. selected top 15 pathways as a cutoff

2. extracted the genes associated with these

3. removed genes with the `[a-z]+:` prefix

4. selected genes associated with key pathways and added meox1

5. plotted dotplots on the remainder



Attempted heatmaps, not very useful.

```{.r .fold-hide}
cat("\n\n_Suggested figure legend:_ Dotplot shows the scaled expression level for genes from enriched GO terms. Units on the horizontal axis correspond to the genotypes. Point size corresponds to the percent of expression in each cell group, and point colour corresponds to average expression z-scores. Genes were selected by membership in key pathways detected by gene enrichment and ontology analyses: sprouting angiogenesis. In addition, meox1 is shown for reference.\n\n")
```



_Suggested figure legend:_ Dotplot shows the scaled expression level for genes from enriched GO terms. Units on the horizontal axis correspond to the genotypes. Point size corresponds to the percent of expression in each cell group, and point colour corresponds to average expression z-scores. Genes were selected by membership in key pathways detected by gene enrichment and ontology analyses: sprouting angiogenesis. In addition, meox1 is shown for reference.

```{.r .fold-hide}
cat("\n\n_Caution against overinterpretation:_ Not an invalid result, but the scale may exaggerate the magnitude of difference.\n
The dotplots look very clean but theres some important nuance\n
1. the pathway genes are sourced from DEG so these would naturally be different
2. the expression values on the dotplot are z-score values and the relatively low scale indicates that the magnitude of the expression isnt massively different (but still different enough)\n\n")
```



_Caution against overinterpretation:_ Not an invalid result, but the scale may exaggerate the magnitude of difference.

The dotplots look very clean but theres some important nuance

1. the pathway genes are sourced from DEG so these would naturally be different
2. the expression values on the dotplot are z-score values and the relatively low scale indicates that the magnitude of the expression isnt massively different (but still different enough)

```{.r .fold-hide}
# plot gene expr
plt_name <- "top15_VECs_GO_pathways_downregulated_geneexpr_dotplot_curated"
ego_vecs_pathways <- c(
    "sprouting angiogenesis"
    )
ego_vecs_selected_rows <- ego_vecs_top15_pathways[
    ego_vecs_top15_pathways$ego_vecs.Description %in% ego_vecs_pathways,
    ]
ego_vecs_genes_per_pathway <- setNames(
    strsplit(ego_vecs_selected_rows$ego_vecs.geneID, "/"),
    ego_vecs_selected_rows$ego_vecs.Description
    )
ego_vecs_genes_per_pathway <- c(unique(
    unlist(ego_vecs_genes_per_pathway)
    ), "meox1")
ego_vecs_genes_per_pathway <- ego_vecs_genes_per_pathway[
    grep(":", ego_vecs_genes_per_pathway, invert=TRUE)
    ]
cat("\nGene quantity: ", length(ego_vecs_genes_per_pathway))
```


Gene quantity:  29

```{.r .fold-hide}
plt <- DotPlot(data_vecs, features=ego_vecs_genes_per_pathway,
    group.by="Genotype", cols=c("#d9d9d9",  "#40004b")
    ) +
    coord_flip() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
```

```
## Warning: Scaling data with a low number of groups may produce misleading
## results
```

```{.r .fold-hide}
print(plt)
```

![](meox1_goanalysis_files/figure-html/grp_vec-14.png)<!-- -->

```{.r .fold-hide}
save_plot_custom(plt, plt_name, height=12, width=4)
```

*PDF saved to:* `../output/figure_extended/go//top15_VECs_GO_pathways_downregulated_geneexpr_dotplot_curated.pdf` 

```{.r .fold-hide}
cat("\n\n_Suggested figure legend:_ Violin plot shows the log2 normalised expression level for differentially expressed genes associated with enriched GO terms. Units on the horizontal axis correspond to the genotypes. Genes were selected by membership in key pathways detected by gene enrichment and ontology analyses: sprouting angiogenesis. In addition, meox1 is shown for reference.\n\n")
```



_Suggested figure legend:_ Violin plot shows the log2 normalised expression level for differentially expressed genes associated with enriched GO terms. Units on the horizontal axis correspond to the genotypes. Genes were selected by membership in key pathways detected by gene enrichment and ontology analyses: sprouting angiogenesis. In addition, meox1 is shown for reference.

```{.r .fold-hide}
cat("\n_Tyrone: Not every single canonical gene associated with the pathway will be shown. It is a function of the data, only the differentially expressed genes associated with the GO term will be shown. Also, I trimmed this down to sprouting angiogenesis only to get a manageable gene list, but its worth noting that their associated genes in the context of this dataset are a full subset of angiogenesis and blood vessel morphogenesis pathways._\n\n")
```


_Tyrone: Not every single canonical gene associated with the pathway will be shown. It is a function of the data, only the differentially expressed genes associated with the GO term will be shown. Also, I trimmed this down to sprouting angiogenesis only to get a manageable gene list, but its worth noting that their associated genes in the context of this dataset are a full subset of angiogenesis and blood vessel morphogenesis pathways._

```{.r .fold-hide}
plt_name <- "top15_VECs_GO_pathways_downregulated_geneexpr_vlnplot_curated"
plt <- VlnPlot(
    data_vecs,
    features = rev(ego_vecs_genes_per_pathway),
    group.by = "Genotype",
    pt.size = 1,
    ncol = 8,
    flip = TRUE,
    cols = c("#dbe2c6", "#657c95")
) + scale_fill_manual(values = c("#dbe2c6", "#657c95")) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
```

```
## Scale for fill is already present.
## Adding another scale for fill, which will replace the existing scale.
```

```{.r .fold-hide}
print(plt)
```

![](meox1_goanalysis_files/figure-html/grp_vec-15.png)<!-- -->

```{.r .fold-hide}
save_plot_custom(plt, plt_name, height=12, width=16)
```

*PDF saved to:* `../output/figure_extended/go//top15_VECs_GO_pathways_downregulated_geneexpr_vlnplot_curated.pdf` 


```{.r .fold-hide}
singles_hmVEC_out <- run_go(singles_hmVEC, data_name="hmVEC")
```



## GO analysis: hmVEC 

Input genes: 14308 

**ORA input genes:** combined = 1798 &nbsp;&nbsp; up = 1247 &nbsp;&nbsp; down = 551 

**ORA (combined):** 2711 terms &nbsp;&nbsp;  **ORA (up):** 2416 terms &nbsp;&nbsp;  **ORA (down):** 1789 terms

*Full ORA ( combined ) table saved to:* `../output/figure_extended/go//hmVEC_ego_combined_full_results.csv` 

*Full ORA ( up ) table saved to:* `../output/figure_extended/go//hmVEC_ego_up_full_results.csv` 

*Full ORA ( down ) table saved to:* `../output/figure_extended/go//hmVEC_ego_down_full_results.csv` 

**ORA significant terms:** combined = 92 , up = 78 , down = 41 

*Significant-only ORA ( combined ) table saved to:* `../output/figure_extended/go//hmVEC_ego_combined_significant_only.csv` 

*Significant-only ORA ( up ) table saved to:* `../output/figure_extended/go//hmVEC_ego_up_significant_only.csv` 

*Significant-only ORA ( down ) table saved to:* `../output/figure_extended/go//hmVEC_ego_down_significant_only.csv` 

### ORA ( combined ): Barplot

Shows the top enriched GO terms ranked by significance, bar length/colour reflecting  adjusted p-value or gene count. **Expected:** terms biologically consistent with the  comparison (e.g. vascular/endothelial terms for EC-related contrasts). **Unexpected:**  dominance by generic housekeeping terms (translation, metabolism) with no thematic  coherence - may indicate a weak or noisy input gene list rather than a real signal.

```
## Warning in fortify(object, showCategory = showCategory, by = x, ...): Arguments in `...` must be used.
## ✖ Problematic argument:
## • by = x
## ℹ Did you misspell an argument name?
```

![](meox1_goanalysis_files/figure-html/sin_hmVEC-1.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//hmVEC_ego_combined_barplot.pdf` 

### ORA ( combined ): Dotplot

Adds gene ratio (fraction of input genes hitting each term) as point size, with colour  for significance - more informative than the barplot since it shows both effect size and  confidence simultaneously. **Expected:** larger, more significant dots clustering around a  coherent biological theme. **Unexpected:** significant terms with very small gene ratios  (few genes driving a 'significant' result) - treat cautiously as potentially fragile.

![](meox1_goanalysis_files/figure-html/sin_hmVEC-2.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//hmVEC_ego_combined_dotplot.pdf` 

### ORA ( combined ): Enrichment map

Network graph connecting GO terms that share overlapping genes - clusters of connected  terms suggest a shared underlying biological process rather than independent findings.  **Expected:** tight clusters of related terms. **Unexpected:** many isolated, unconnected  terms, suggesting the significant terms are thematically scattered rather than coherent.

![](meox1_goanalysis_files/figure-html/sin_hmVEC-3.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//hmVEC_ego_combined_emapplot.pdf` 

### ORA ( combined ): Category-gene network

Shows which specific genes drive which enriched terms, useful for tracing a term back to  the actual DE genes behind it. **Expected:** recognisable marker genes linking to expected  terms. **Unexpected:** a term driven entirely by one or two genes rather than a broad gene  set - weaker evidence than a term supported by many genes.

![](meox1_goanalysis_files/figure-html/sin_hmVEC-4.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//hmVEC_ego_combined_cnetplot.pdf` 

### ORA ( up ): Barplot

Shows the top enriched GO terms ranked by significance, bar length/colour reflecting  adjusted p-value or gene count. **Expected:** terms biologically consistent with the  comparison (e.g. vascular/endothelial terms for EC-related contrasts). **Unexpected:**  dominance by generic housekeeping terms (translation, metabolism) with no thematic  coherence - may indicate a weak or noisy input gene list rather than a real signal.

```
## Warning in fortify(object, showCategory = showCategory, by = x, ...): Arguments in `...` must be used.
## ✖ Problematic argument:
## • by = x
## ℹ Did you misspell an argument name?
```

![](meox1_goanalysis_files/figure-html/sin_hmVEC-5.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//hmVEC_ego_up_barplot.pdf` 

### ORA ( up ): Dotplot

Adds gene ratio (fraction of input genes hitting each term) as point size, with colour  for significance - more informative than the barplot since it shows both effect size and  confidence simultaneously. **Expected:** larger, more significant dots clustering around a  coherent biological theme. **Unexpected:** significant terms with very small gene ratios  (few genes driving a 'significant' result) - treat cautiously as potentially fragile.

![](meox1_goanalysis_files/figure-html/sin_hmVEC-6.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//hmVEC_ego_up_dotplot.pdf` 

### ORA ( up ): Enrichment map

Network graph connecting GO terms that share overlapping genes - clusters of connected  terms suggest a shared underlying biological process rather than independent findings.  **Expected:** tight clusters of related terms. **Unexpected:** many isolated, unconnected  terms, suggesting the significant terms are thematically scattered rather than coherent.

![](meox1_goanalysis_files/figure-html/sin_hmVEC-7.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//hmVEC_ego_up_emapplot.pdf` 

### ORA ( up ): Category-gene network

Shows which specific genes drive which enriched terms, useful for tracing a term back to  the actual DE genes behind it. **Expected:** recognisable marker genes linking to expected  terms. **Unexpected:** a term driven entirely by one or two genes rather than a broad gene  set - weaker evidence than a term supported by many genes.

![](meox1_goanalysis_files/figure-html/sin_hmVEC-8.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//hmVEC_ego_up_cnetplot.pdf` 

### ORA ( down ): Barplot

Shows the top enriched GO terms ranked by significance, bar length/colour reflecting  adjusted p-value or gene count. **Expected:** terms biologically consistent with the  comparison (e.g. vascular/endothelial terms for EC-related contrasts). **Unexpected:**  dominance by generic housekeeping terms (translation, metabolism) with no thematic  coherence - may indicate a weak or noisy input gene list rather than a real signal.

```
## Warning in fortify(object, showCategory = showCategory, by = x, ...): Arguments in `...` must be used.
## ✖ Problematic argument:
## • by = x
## ℹ Did you misspell an argument name?
```

![](meox1_goanalysis_files/figure-html/sin_hmVEC-9.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//hmVEC_ego_down_barplot.pdf` 

### ORA ( down ): Dotplot

Adds gene ratio (fraction of input genes hitting each term) as point size, with colour  for significance - more informative than the barplot since it shows both effect size and  confidence simultaneously. **Expected:** larger, more significant dots clustering around a  coherent biological theme. **Unexpected:** significant terms with very small gene ratios  (few genes driving a 'significant' result) - treat cautiously as potentially fragile.

![](meox1_goanalysis_files/figure-html/sin_hmVEC-10.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//hmVEC_ego_down_dotplot.pdf` 

### ORA ( down ): Enrichment map

Network graph connecting GO terms that share overlapping genes - clusters of connected  terms suggest a shared underlying biological process rather than independent findings.  **Expected:** tight clusters of related terms. **Unexpected:** many isolated, unconnected  terms, suggesting the significant terms are thematically scattered rather than coherent.

![](meox1_goanalysis_files/figure-html/sin_hmVEC-11.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//hmVEC_ego_down_emapplot.pdf` 

### ORA ( down ): Category-gene network

Shows which specific genes drive which enriched terms, useful for tracing a term back to  the actual DE genes behind it. **Expected:** recognisable marker genes linking to expected  terms. **Unexpected:** a term driven entirely by one or two genes rather than a broad gene  set - weaker evidence than a term supported by many genes.

![](meox1_goanalysis_files/figure-html/sin_hmVEC-12.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//hmVEC_ego_down_cnetplot.pdf` 

**Done:** hmVEC 

---


```{.r .fold-hide}
singles_mVEC_out <- run_go(singles_mVEC, data_name="mVEC")
```



## GO analysis: mVEC 

Input genes: 14187 

**ORA input genes:** combined = 1979 &nbsp;&nbsp; up = 884 &nbsp;&nbsp; down = 1095 

**ORA (combined):** 2738 terms &nbsp;&nbsp;  **ORA (up):** 2179 terms &nbsp;&nbsp;  **ORA (down):** 2190 terms

*Full ORA ( combined ) table saved to:* `../output/figure_extended/go//mVEC_ego_combined_full_results.csv` 

*Full ORA ( up ) table saved to:* `../output/figure_extended/go//mVEC_ego_up_full_results.csv` 

*Full ORA ( down ) table saved to:* `../output/figure_extended/go//mVEC_ego_down_full_results.csv` 

**ORA significant terms:** combined = 112 , up = 60 , down = 58 

*Significant-only ORA ( combined ) table saved to:* `../output/figure_extended/go//mVEC_ego_combined_significant_only.csv` 

*Significant-only ORA ( up ) table saved to:* `../output/figure_extended/go//mVEC_ego_up_significant_only.csv` 

*Significant-only ORA ( down ) table saved to:* `../output/figure_extended/go//mVEC_ego_down_significant_only.csv` 

### ORA ( combined ): Barplot

Shows the top enriched GO terms ranked by significance, bar length/colour reflecting  adjusted p-value or gene count. **Expected:** terms biologically consistent with the  comparison (e.g. vascular/endothelial terms for EC-related contrasts). **Unexpected:**  dominance by generic housekeeping terms (translation, metabolism) with no thematic  coherence - may indicate a weak or noisy input gene list rather than a real signal.

```
## Warning in fortify(object, showCategory = showCategory, by = x, ...): Arguments in `...` must be used.
## ✖ Problematic argument:
## • by = x
## ℹ Did you misspell an argument name?
```

![](meox1_goanalysis_files/figure-html/sin_mVEC-1.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//mVEC_ego_combined_barplot.pdf` 

### ORA ( combined ): Dotplot

Adds gene ratio (fraction of input genes hitting each term) as point size, with colour  for significance - more informative than the barplot since it shows both effect size and  confidence simultaneously. **Expected:** larger, more significant dots clustering around a  coherent biological theme. **Unexpected:** significant terms with very small gene ratios  (few genes driving a 'significant' result) - treat cautiously as potentially fragile.

![](meox1_goanalysis_files/figure-html/sin_mVEC-2.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//mVEC_ego_combined_dotplot.pdf` 

### ORA ( combined ): Enrichment map

Network graph connecting GO terms that share overlapping genes - clusters of connected  terms suggest a shared underlying biological process rather than independent findings.  **Expected:** tight clusters of related terms. **Unexpected:** many isolated, unconnected  terms, suggesting the significant terms are thematically scattered rather than coherent.

![](meox1_goanalysis_files/figure-html/sin_mVEC-3.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//mVEC_ego_combined_emapplot.pdf` 

### ORA ( combined ): Category-gene network

Shows which specific genes drive which enriched terms, useful for tracing a term back to  the actual DE genes behind it. **Expected:** recognisable marker genes linking to expected  terms. **Unexpected:** a term driven entirely by one or two genes rather than a broad gene  set - weaker evidence than a term supported by many genes.

![](meox1_goanalysis_files/figure-html/sin_mVEC-4.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//mVEC_ego_combined_cnetplot.pdf` 

### ORA ( up ): Barplot

Shows the top enriched GO terms ranked by significance, bar length/colour reflecting  adjusted p-value or gene count. **Expected:** terms biologically consistent with the  comparison (e.g. vascular/endothelial terms for EC-related contrasts). **Unexpected:**  dominance by generic housekeeping terms (translation, metabolism) with no thematic  coherence - may indicate a weak or noisy input gene list rather than a real signal.

```
## Warning in fortify(object, showCategory = showCategory, by = x, ...): Arguments in `...` must be used.
## ✖ Problematic argument:
## • by = x
## ℹ Did you misspell an argument name?
```

![](meox1_goanalysis_files/figure-html/sin_mVEC-5.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//mVEC_ego_up_barplot.pdf` 

### ORA ( up ): Dotplot

Adds gene ratio (fraction of input genes hitting each term) as point size, with colour  for significance - more informative than the barplot since it shows both effect size and  confidence simultaneously. **Expected:** larger, more significant dots clustering around a  coherent biological theme. **Unexpected:** significant terms with very small gene ratios  (few genes driving a 'significant' result) - treat cautiously as potentially fragile.

![](meox1_goanalysis_files/figure-html/sin_mVEC-6.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//mVEC_ego_up_dotplot.pdf` 

### ORA ( up ): Enrichment map

Network graph connecting GO terms that share overlapping genes - clusters of connected  terms suggest a shared underlying biological process rather than independent findings.  **Expected:** tight clusters of related terms. **Unexpected:** many isolated, unconnected  terms, suggesting the significant terms are thematically scattered rather than coherent.

![](meox1_goanalysis_files/figure-html/sin_mVEC-7.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//mVEC_ego_up_emapplot.pdf` 

### ORA ( up ): Category-gene network

Shows which specific genes drive which enriched terms, useful for tracing a term back to  the actual DE genes behind it. **Expected:** recognisable marker genes linking to expected  terms. **Unexpected:** a term driven entirely by one or two genes rather than a broad gene  set - weaker evidence than a term supported by many genes.

![](meox1_goanalysis_files/figure-html/sin_mVEC-8.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//mVEC_ego_up_cnetplot.pdf` 

### ORA ( down ): Barplot

Shows the top enriched GO terms ranked by significance, bar length/colour reflecting  adjusted p-value or gene count. **Expected:** terms biologically consistent with the  comparison (e.g. vascular/endothelial terms for EC-related contrasts). **Unexpected:**  dominance by generic housekeeping terms (translation, metabolism) with no thematic  coherence - may indicate a weak or noisy input gene list rather than a real signal.

```
## Warning in fortify(object, showCategory = showCategory, by = x, ...): Arguments in `...` must be used.
## ✖ Problematic argument:
## • by = x
## ℹ Did you misspell an argument name?
```

![](meox1_goanalysis_files/figure-html/sin_mVEC-9.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//mVEC_ego_down_barplot.pdf` 

### ORA ( down ): Dotplot

Adds gene ratio (fraction of input genes hitting each term) as point size, with colour  for significance - more informative than the barplot since it shows both effect size and  confidence simultaneously. **Expected:** larger, more significant dots clustering around a  coherent biological theme. **Unexpected:** significant terms with very small gene ratios  (few genes driving a 'significant' result) - treat cautiously as potentially fragile.

![](meox1_goanalysis_files/figure-html/sin_mVEC-10.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//mVEC_ego_down_dotplot.pdf` 

### ORA ( down ): Enrichment map

Network graph connecting GO terms that share overlapping genes - clusters of connected  terms suggest a shared underlying biological process rather than independent findings.  **Expected:** tight clusters of related terms. **Unexpected:** many isolated, unconnected  terms, suggesting the significant terms are thematically scattered rather than coherent.

![](meox1_goanalysis_files/figure-html/sin_mVEC-11.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//mVEC_ego_down_emapplot.pdf` 

### ORA ( down ): Category-gene network

Shows which specific genes drive which enriched terms, useful for tracing a term back to  the actual DE genes behind it. **Expected:** recognisable marker genes linking to expected  terms. **Unexpected:** a term driven entirely by one or two genes rather than a broad gene  set - weaker evidence than a term supported by many genes.

![](meox1_goanalysis_files/figure-html/sin_mVEC-12.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//mVEC_ego_down_cnetplot.pdf` 

**Done:** mVEC 

---


## For developers

Sample run command:


``` bash
Rscript -e "
rmarkdown::render(
  'meox1_goanalysis.Rmd',
  output_file = './meox1_goanalysis.html'
)
"
```

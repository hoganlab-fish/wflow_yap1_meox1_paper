---
title: "GO analyses for yap1"
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

# GO analyses for yap1

## Important info

We want to compare ontology across different groups and subgroups of wild type and mutant cells.

- This assumes you ran `yap1_deg_confects_prep.Rmd` first, in which we merged VECs and LECs into groups for comparison.
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
## S Xu, E Hu, Y Cai, Z Xie, X Luo, L Zhan, W Tang, Q Wang, B Liu, R Wang,
## W Xie, T Wu, L Xie, G Yu. Using clusterProfiler to characterize
## multiomics data. Nature Protocols. 2024, 19(11):3292-3320
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
## S Xu, E Hu, Y Cai, Z Xie, X Luo, L Zhan, W Tang, Q Wang, B Liu, R Wang,
## W Xie, T Wu, L Xie, G Yu. Using clusterProfiler to characterize
## multiomics data. Nature Protocols. 2024, 19(11):3292-3320
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
library(UpSetR)

outfile_dir <- "../output/figure_extended/go/"

singles_yap_path <- "../../Saki_data/analysis/D10025_yap1_dataset/DEGs/Level_03_DEG_mutVSwt_L3_celltype_fc0.00_minpct_0.01.csv"
grouped_yap_path <- "../output/figure_extended/dge_confects/yap1_Level_03_DEG_mutVSwt_L3_celltype_merged_VEC_merged_LEC_fc0.00_minpct_0.01.csv"

singles_yap_hippo_path <- "../../Saki_data/analysis/D10025_yap1_dataset/DEGs/analysis/D10025_yap1_dataset/DEGs/Level_03_DEG_hippo_targets_mutVSwt_L3_celltype_fc0.00_minpct_0.00.csv"
grouped_yap_hippo_path <- "../output/figure_extended/dge_confects/yap1_Level_03_DEG_mutVSwt_L3_celltype_merged_VEC_merged_LEC_hippo_targets_fc0.00_minpct_0.00.csv"

custom_groups <- list(
    "LEC__preLEC" = c("LEC", "preLEC"),
    "cVEC__hmVEC__mVEC" = c("cVEC", "hmVEC", "mVEC")
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
grouped_yap <- read.csv(grouped_yap_path)


custom_groups <- list(
    "LEC__preLEC" = c("LEC", "preLEC"),
    "cVEC__hmVEC__mVEC" = c("cVEC", "hmVEC", "mVEC")
)


grouped_lecs <- grouped_yap[grouped_yap$group == "LEC__preLEC", ]
rownames(grouped_lecs) <- grouped_lecs$Gene

grouped_vecs <- grouped_yap[grouped_yap$group == "cVEC__hmVEC__mVEC", ]
rownames(grouped_vecs) <- grouped_vecs$Gene


singles_yap <- read.csv(singles_yap_path)

singles_LEC <- singles_yap[singles_yap$group == "LEC", ]
rownames(singles_LEC) <- singles_LEC$Gene

singles_pre_muLEC <- singles_yap[singles_yap$group == "preLEC", ]
rownames(singles_pre_muLEC) <- singles_pre_muLEC$Gene

singles_cVEC <- singles_yap[singles_yap$group == "cVEC", ]
rownames(singles_cVEC) <- singles_cVEC$Gene

singles_hmVEC <- singles_yap[singles_yap$group == "hmVEC", ]
rownames(singles_hmVEC) <- singles_hmVEC$Gene

singles_mVEC <- singles_yap[singles_yap$group == "mVEC", ]
rownames(singles_mVEC) <- singles_mVEC$Gene
```

We have separate sections for each plot series for easy comparison.


```{.r .fold-hide}
grouped_lecs_out <- run_go(grouped_lecs, data_name="LEC__preLEC")
```



## GO analysis: LEC__preLEC 

Input genes: 10790 

**ORA input genes:** combined = 7030 &nbsp;&nbsp; up = 5531 &nbsp;&nbsp; down = 1499 

**ORA (combined):** 3103 terms &nbsp;&nbsp;  **ORA (up):** 3068 terms &nbsp;&nbsp;  **ORA (down):** 2513 terms

*Full ORA ( combined ) table saved to:* `../output/figure_extended/go//LEC__preLEC_ego_combined_full_results.csv` 

*Full ORA ( up ) table saved to:* `../output/figure_extended/go//LEC__preLEC_ego_up_full_results.csv` 

*Full ORA ( down ) table saved to:* `../output/figure_extended/go//LEC__preLEC_ego_down_full_results.csv` 

**ORA significant terms:** combined = 503 , up = 413 , down = 40 

*Significant-only ORA ( combined ) table saved to:* `../output/figure_extended/go//LEC__preLEC_ego_combined_significant_only.csv` 

*Significant-only ORA ( up ) table saved to:* `../output/figure_extended/go//LEC__preLEC_ego_up_significant_only.csv` 

*Significant-only ORA ( down ) table saved to:* `../output/figure_extended/go//LEC__preLEC_ego_down_significant_only.csv` 

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

![](yap1_goanalysis_files/figure-html/grp_lec-1.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//LEC__preLEC_ego_combined_barplot.pdf` 

### ORA ( combined ): Dotplot

Adds gene ratio (fraction of input genes hitting each term) as point size, with colour  for significance - more informative than the barplot since it shows both effect size and  confidence simultaneously. **Expected:** larger, more significant dots clustering around a  coherent biological theme. **Unexpected:** significant terms with very small gene ratios  (few genes driving a 'significant' result) - treat cautiously as potentially fragile.

![](yap1_goanalysis_files/figure-html/grp_lec-2.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//LEC__preLEC_ego_combined_dotplot.pdf` 

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

![](yap1_goanalysis_files/figure-html/grp_lec-3.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//LEC__preLEC_ego_combined_emapplot.pdf` 

### ORA ( combined ): Category-gene network

Shows which specific genes drive which enriched terms, useful for tracing a term back to  the actual DE genes behind it. **Expected:** recognisable marker genes linking to expected  terms. **Unexpected:** a term driven entirely by one or two genes rather than a broad gene  set - weaker evidence than a term supported by many genes.

![](yap1_goanalysis_files/figure-html/grp_lec-4.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//LEC__preLEC_ego_combined_cnetplot.pdf` 

### ORA ( up ): Barplot

Shows the top enriched GO terms ranked by significance, bar length/colour reflecting  adjusted p-value or gene count. **Expected:** terms biologically consistent with the  comparison (e.g. vascular/endothelial terms for EC-related contrasts). **Unexpected:**  dominance by generic housekeeping terms (translation, metabolism) with no thematic  coherence - may indicate a weak or noisy input gene list rather than a real signal.

```
## Warning in fortify(object, showCategory = showCategory, by = x, ...): Arguments in `...` must be used.
## ✖ Problematic argument:
## • by = x
## ℹ Did you misspell an argument name?
```

![](yap1_goanalysis_files/figure-html/grp_lec-5.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//LEC__preLEC_ego_up_barplot.pdf` 

### ORA ( up ): Dotplot

Adds gene ratio (fraction of input genes hitting each term) as point size, with colour  for significance - more informative than the barplot since it shows both effect size and  confidence simultaneously. **Expected:** larger, more significant dots clustering around a  coherent biological theme. **Unexpected:** significant terms with very small gene ratios  (few genes driving a 'significant' result) - treat cautiously as potentially fragile.

![](yap1_goanalysis_files/figure-html/grp_lec-6.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//LEC__preLEC_ego_up_dotplot.pdf` 

### ORA ( up ): Enrichment map

Network graph connecting GO terms that share overlapping genes - clusters of connected  terms suggest a shared underlying biological process rather than independent findings.  **Expected:** tight clusters of related terms. **Unexpected:** many isolated, unconnected  terms, suggesting the significant terms are thematically scattered rather than coherent.

![](yap1_goanalysis_files/figure-html/grp_lec-7.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//LEC__preLEC_ego_up_emapplot.pdf` 

### ORA ( up ): Category-gene network

Shows which specific genes drive which enriched terms, useful for tracing a term back to  the actual DE genes behind it. **Expected:** recognisable marker genes linking to expected  terms. **Unexpected:** a term driven entirely by one or two genes rather than a broad gene  set - weaker evidence than a term supported by many genes.

![](yap1_goanalysis_files/figure-html/grp_lec-8.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//LEC__preLEC_ego_up_cnetplot.pdf` 

### ORA ( down ): Barplot

Shows the top enriched GO terms ranked by significance, bar length/colour reflecting  adjusted p-value or gene count. **Expected:** terms biologically consistent with the  comparison (e.g. vascular/endothelial terms for EC-related contrasts). **Unexpected:**  dominance by generic housekeeping terms (translation, metabolism) with no thematic  coherence - may indicate a weak or noisy input gene list rather than a real signal.

```
## Warning in fortify(object, showCategory = showCategory, by = x, ...): Arguments in `...` must be used.
## ✖ Problematic argument:
## • by = x
## ℹ Did you misspell an argument name?
```

![](yap1_goanalysis_files/figure-html/grp_lec-9.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//LEC__preLEC_ego_down_barplot.pdf` 

### ORA ( down ): Dotplot

Adds gene ratio (fraction of input genes hitting each term) as point size, with colour  for significance - more informative than the barplot since it shows both effect size and  confidence simultaneously. **Expected:** larger, more significant dots clustering around a  coherent biological theme. **Unexpected:** significant terms with very small gene ratios  (few genes driving a 'significant' result) - treat cautiously as potentially fragile.

![](yap1_goanalysis_files/figure-html/grp_lec-10.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//LEC__preLEC_ego_down_dotplot.pdf` 

### ORA ( down ): Enrichment map

Network graph connecting GO terms that share overlapping genes - clusters of connected  terms suggest a shared underlying biological process rather than independent findings.  **Expected:** tight clusters of related terms. **Unexpected:** many isolated, unconnected  terms, suggesting the significant terms are thematically scattered rather than coherent.

![](yap1_goanalysis_files/figure-html/grp_lec-11.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//LEC__preLEC_ego_down_emapplot.pdf` 

### ORA ( down ): Category-gene network

Shows which specific genes drive which enriched terms, useful for tracing a term back to  the actual DE genes behind it. **Expected:** recognisable marker genes linking to expected  terms. **Unexpected:** a term driven entirely by one or two genes rather than a broad gene  set - weaker evidence than a term supported by many genes.

![](yap1_goanalysis_files/figure-html/grp_lec-12.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//LEC__preLEC_ego_down_cnetplot.pdf` 

**Done:** LEC__preLEC 

---


```{.r .fold-hide}
singles_LEC_out <- run_go(singles_LEC, data_name="LEC")
```



## GO analysis: LEC 

Input genes: 11463 

**ORA input genes:** combined = 2624 &nbsp;&nbsp; up = 2537 &nbsp;&nbsp; down = 87 

**ORA (combined):** 2889 terms &nbsp;&nbsp;  **ORA (up):** 2877 terms &nbsp;&nbsp;  **ORA (down):** 568 terms

*Full ORA ( combined ) table saved to:* `../output/figure_extended/go//LEC_ego_combined_full_results.csv` 

*Full ORA ( up ) table saved to:* `../output/figure_extended/go//LEC_ego_up_full_results.csv` 

*Full ORA ( down ) table saved to:* `../output/figure_extended/go//LEC_ego_down_full_results.csv` 

**ORA significant terms:** combined = 111 , up = 107 , down = 0 

*Significant-only ORA ( combined ) table saved to:* `../output/figure_extended/go//LEC_ego_combined_significant_only.csv` 

*Significant-only ORA ( up ) table saved to:* `../output/figure_extended/go//LEC_ego_up_significant_only.csv` 

### ORA ( combined ): Barplot

Shows the top enriched GO terms ranked by significance, bar length/colour reflecting  adjusted p-value or gene count. **Expected:** terms biologically consistent with the  comparison (e.g. vascular/endothelial terms for EC-related contrasts). **Unexpected:**  dominance by generic housekeeping terms (translation, metabolism) with no thematic  coherence - may indicate a weak or noisy input gene list rather than a real signal.

```
## Warning in fortify(object, showCategory = showCategory, by = x, ...): Arguments in `...` must be used.
## ✖ Problematic argument:
## • by = x
## ℹ Did you misspell an argument name?
```

![](yap1_goanalysis_files/figure-html/sin_LEC-1.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//LEC_ego_combined_barplot.pdf` 

### ORA ( combined ): Dotplot

Adds gene ratio (fraction of input genes hitting each term) as point size, with colour  for significance - more informative than the barplot since it shows both effect size and  confidence simultaneously. **Expected:** larger, more significant dots clustering around a  coherent biological theme. **Unexpected:** significant terms with very small gene ratios  (few genes driving a 'significant' result) - treat cautiously as potentially fragile.

![](yap1_goanalysis_files/figure-html/sin_LEC-2.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//LEC_ego_combined_dotplot.pdf` 

### ORA ( combined ): Enrichment map

Network graph connecting GO terms that share overlapping genes - clusters of connected  terms suggest a shared underlying biological process rather than independent findings.  **Expected:** tight clusters of related terms. **Unexpected:** many isolated, unconnected  terms, suggesting the significant terms are thematically scattered rather than coherent.

![](yap1_goanalysis_files/figure-html/sin_LEC-3.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//LEC_ego_combined_emapplot.pdf` 

### ORA ( combined ): Category-gene network

Shows which specific genes drive which enriched terms, useful for tracing a term back to  the actual DE genes behind it. **Expected:** recognisable marker genes linking to expected  terms. **Unexpected:** a term driven entirely by one or two genes rather than a broad gene  set - weaker evidence than a term supported by many genes.

![](yap1_goanalysis_files/figure-html/sin_LEC-4.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//LEC_ego_combined_cnetplot.pdf` 

### ORA ( up ): Barplot

Shows the top enriched GO terms ranked by significance, bar length/colour reflecting  adjusted p-value or gene count. **Expected:** terms biologically consistent with the  comparison (e.g. vascular/endothelial terms for EC-related contrasts). **Unexpected:**  dominance by generic housekeeping terms (translation, metabolism) with no thematic  coherence - may indicate a weak or noisy input gene list rather than a real signal.

```
## Warning in fortify(object, showCategory = showCategory, by = x, ...): Arguments in `...` must be used.
## ✖ Problematic argument:
## • by = x
## ℹ Did you misspell an argument name?
```

![](yap1_goanalysis_files/figure-html/sin_LEC-5.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//LEC_ego_up_barplot.pdf` 

### ORA ( up ): Dotplot

Adds gene ratio (fraction of input genes hitting each term) as point size, with colour  for significance - more informative than the barplot since it shows both effect size and  confidence simultaneously. **Expected:** larger, more significant dots clustering around a  coherent biological theme. **Unexpected:** significant terms with very small gene ratios  (few genes driving a 'significant' result) - treat cautiously as potentially fragile.

![](yap1_goanalysis_files/figure-html/sin_LEC-6.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//LEC_ego_up_dotplot.pdf` 

### ORA ( up ): Enrichment map

Network graph connecting GO terms that share overlapping genes - clusters of connected  terms suggest a shared underlying biological process rather than independent findings.  **Expected:** tight clusters of related terms. **Unexpected:** many isolated, unconnected  terms, suggesting the significant terms are thematically scattered rather than coherent.

![](yap1_goanalysis_files/figure-html/sin_LEC-7.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//LEC_ego_up_emapplot.pdf` 

### ORA ( up ): Category-gene network

Shows which specific genes drive which enriched terms, useful for tracing a term back to  the actual DE genes behind it. **Expected:** recognisable marker genes linking to expected  terms. **Unexpected:** a term driven entirely by one or two genes rather than a broad gene  set - weaker evidence than a term supported by many genes.

![](yap1_goanalysis_files/figure-html/sin_LEC-8.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//LEC_ego_up_cnetplot.pdf` 

*No ORA ( down ) terms found for LEC - skipping ORA ( down ) plots*

**Done:** LEC 

---


```{.r .fold-hide}
singles_pre_muLEC_out <- run_go(singles_pre_muLEC, data_name="preLEC")
```



## GO analysis: preLEC 

Input genes: 10717 

**ORA input genes:** combined = 6212 &nbsp;&nbsp; up = 5215 &nbsp;&nbsp; down = 997 

**ORA (combined):** 3084 terms &nbsp;&nbsp;  **ORA (up):** 3050 terms &nbsp;&nbsp;  **ORA (down):** 2264 terms

*Full ORA ( combined ) table saved to:* `../output/figure_extended/go//preLEC_ego_combined_full_results.csv` 

*Full ORA ( up ) table saved to:* `../output/figure_extended/go//preLEC_ego_up_full_results.csv` 

*Full ORA ( down ) table saved to:* `../output/figure_extended/go//preLEC_ego_down_full_results.csv` 

**ORA significant terms:** combined = 439 , up = 375 , down = 16 

*Significant-only ORA ( combined ) table saved to:* `../output/figure_extended/go//preLEC_ego_combined_significant_only.csv` 

*Significant-only ORA ( up ) table saved to:* `../output/figure_extended/go//preLEC_ego_up_significant_only.csv` 

*Significant-only ORA ( down ) table saved to:* `../output/figure_extended/go//preLEC_ego_down_significant_only.csv` 

### ORA ( combined ): Barplot

Shows the top enriched GO terms ranked by significance, bar length/colour reflecting  adjusted p-value or gene count. **Expected:** terms biologically consistent with the  comparison (e.g. vascular/endothelial terms for EC-related contrasts). **Unexpected:**  dominance by generic housekeeping terms (translation, metabolism) with no thematic  coherence - may indicate a weak or noisy input gene list rather than a real signal.

```
## Warning in fortify(object, showCategory = showCategory, by = x, ...): Arguments in `...` must be used.
## ✖ Problematic argument:
## • by = x
## ℹ Did you misspell an argument name?
```

![](yap1_goanalysis_files/figure-html/sin_pre_muLEC-1.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//preLEC_ego_combined_barplot.pdf` 

### ORA ( combined ): Dotplot

Adds gene ratio (fraction of input genes hitting each term) as point size, with colour  for significance - more informative than the barplot since it shows both effect size and  confidence simultaneously. **Expected:** larger, more significant dots clustering around a  coherent biological theme. **Unexpected:** significant terms with very small gene ratios  (few genes driving a 'significant' result) - treat cautiously as potentially fragile.

![](yap1_goanalysis_files/figure-html/sin_pre_muLEC-2.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//preLEC_ego_combined_dotplot.pdf` 

### ORA ( combined ): Enrichment map

Network graph connecting GO terms that share overlapping genes - clusters of connected  terms suggest a shared underlying biological process rather than independent findings.  **Expected:** tight clusters of related terms. **Unexpected:** many isolated, unconnected  terms, suggesting the significant terms are thematically scattered rather than coherent.

![](yap1_goanalysis_files/figure-html/sin_pre_muLEC-3.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//preLEC_ego_combined_emapplot.pdf` 

### ORA ( combined ): Category-gene network

Shows which specific genes drive which enriched terms, useful for tracing a term back to  the actual DE genes behind it. **Expected:** recognisable marker genes linking to expected  terms. **Unexpected:** a term driven entirely by one or two genes rather than a broad gene  set - weaker evidence than a term supported by many genes.

![](yap1_goanalysis_files/figure-html/sin_pre_muLEC-4.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//preLEC_ego_combined_cnetplot.pdf` 

### ORA ( up ): Barplot

Shows the top enriched GO terms ranked by significance, bar length/colour reflecting  adjusted p-value or gene count. **Expected:** terms biologically consistent with the  comparison (e.g. vascular/endothelial terms for EC-related contrasts). **Unexpected:**  dominance by generic housekeeping terms (translation, metabolism) with no thematic  coherence - may indicate a weak or noisy input gene list rather than a real signal.

```
## Warning in fortify(object, showCategory = showCategory, by = x, ...): Arguments in `...` must be used.
## ✖ Problematic argument:
## • by = x
## ℹ Did you misspell an argument name?
```

![](yap1_goanalysis_files/figure-html/sin_pre_muLEC-5.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//preLEC_ego_up_barplot.pdf` 

### ORA ( up ): Dotplot

Adds gene ratio (fraction of input genes hitting each term) as point size, with colour  for significance - more informative than the barplot since it shows both effect size and  confidence simultaneously. **Expected:** larger, more significant dots clustering around a  coherent biological theme. **Unexpected:** significant terms with very small gene ratios  (few genes driving a 'significant' result) - treat cautiously as potentially fragile.

![](yap1_goanalysis_files/figure-html/sin_pre_muLEC-6.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//preLEC_ego_up_dotplot.pdf` 

### ORA ( up ): Enrichment map

Network graph connecting GO terms that share overlapping genes - clusters of connected  terms suggest a shared underlying biological process rather than independent findings.  **Expected:** tight clusters of related terms. **Unexpected:** many isolated, unconnected  terms, suggesting the significant terms are thematically scattered rather than coherent.

![](yap1_goanalysis_files/figure-html/sin_pre_muLEC-7.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//preLEC_ego_up_emapplot.pdf` 

### ORA ( up ): Category-gene network

Shows which specific genes drive which enriched terms, useful for tracing a term back to  the actual DE genes behind it. **Expected:** recognisable marker genes linking to expected  terms. **Unexpected:** a term driven entirely by one or two genes rather than a broad gene  set - weaker evidence than a term supported by many genes.

![](yap1_goanalysis_files/figure-html/sin_pre_muLEC-8.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//preLEC_ego_up_cnetplot.pdf` 

### ORA ( down ): Barplot

Shows the top enriched GO terms ranked by significance, bar length/colour reflecting  adjusted p-value or gene count. **Expected:** terms biologically consistent with the  comparison (e.g. vascular/endothelial terms for EC-related contrasts). **Unexpected:**  dominance by generic housekeeping terms (translation, metabolism) with no thematic  coherence - may indicate a weak or noisy input gene list rather than a real signal.

```
## Warning in fortify(object, showCategory = showCategory, by = x, ...): Arguments in `...` must be used.
## ✖ Problematic argument:
## • by = x
## ℹ Did you misspell an argument name?
```

![](yap1_goanalysis_files/figure-html/sin_pre_muLEC-9.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//preLEC_ego_down_barplot.pdf` 

### ORA ( down ): Dotplot

Adds gene ratio (fraction of input genes hitting each term) as point size, with colour  for significance - more informative than the barplot since it shows both effect size and  confidence simultaneously. **Expected:** larger, more significant dots clustering around a  coherent biological theme. **Unexpected:** significant terms with very small gene ratios  (few genes driving a 'significant' result) - treat cautiously as potentially fragile.

![](yap1_goanalysis_files/figure-html/sin_pre_muLEC-10.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//preLEC_ego_down_dotplot.pdf` 

### ORA ( down ): Enrichment map

Network graph connecting GO terms that share overlapping genes - clusters of connected  terms suggest a shared underlying biological process rather than independent findings.  **Expected:** tight clusters of related terms. **Unexpected:** many isolated, unconnected  terms, suggesting the significant terms are thematically scattered rather than coherent.

![](yap1_goanalysis_files/figure-html/sin_pre_muLEC-11.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//preLEC_ego_down_emapplot.pdf` 

### ORA ( down ): Category-gene network

Shows which specific genes drive which enriched terms, useful for tracing a term back to  the actual DE genes behind it. **Expected:** recognisable marker genes linking to expected  terms. **Unexpected:** a term driven entirely by one or two genes rather than a broad gene  set - weaker evidence than a term supported by many genes.

![](yap1_goanalysis_files/figure-html/sin_pre_muLEC-12.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//preLEC_ego_down_cnetplot.pdf` 

**Done:** preLEC 

---


```{.r .fold-hide}
grouped_vecs_out <- run_go(grouped_vecs, data_name="cVEC__hmVEC__iVEC__mVEC")
```



## GO analysis: cVEC__hmVEC__iVEC__mVEC 

Input genes: 10826 

**ORA input genes:** combined = 6713 &nbsp;&nbsp; up = 5580 &nbsp;&nbsp; down = 1133 

**ORA (combined):** 3114 terms &nbsp;&nbsp;  **ORA (up):** 3093 terms &nbsp;&nbsp;  **ORA (down):** 2453 terms

*Full ORA ( combined ) table saved to:* `../output/figure_extended/go//cVEC__hmVEC__iVEC__mVEC_ego_combined_full_results.csv` 

*Full ORA ( up ) table saved to:* `../output/figure_extended/go//cVEC__hmVEC__iVEC__mVEC_ego_up_full_results.csv` 

*Full ORA ( down ) table saved to:* `../output/figure_extended/go//cVEC__hmVEC__iVEC__mVEC_ego_down_full_results.csv` 

**ORA significant terms:** combined = 476 , up = 367 , down = 48 

*Significant-only ORA ( combined ) table saved to:* `../output/figure_extended/go//cVEC__hmVEC__iVEC__mVEC_ego_combined_significant_only.csv` 

*Significant-only ORA ( up ) table saved to:* `../output/figure_extended/go//cVEC__hmVEC__iVEC__mVEC_ego_up_significant_only.csv` 

*Significant-only ORA ( down ) table saved to:* `../output/figure_extended/go//cVEC__hmVEC__iVEC__mVEC_ego_down_significant_only.csv` 

### ORA ( combined ): Barplot

Shows the top enriched GO terms ranked by significance, bar length/colour reflecting  adjusted p-value or gene count. **Expected:** terms biologically consistent with the  comparison (e.g. vascular/endothelial terms for EC-related contrasts). **Unexpected:**  dominance by generic housekeeping terms (translation, metabolism) with no thematic  coherence - may indicate a weak or noisy input gene list rather than a real signal.

```
## Warning in fortify(object, showCategory = showCategory, by = x, ...): Arguments in `...` must be used.
## ✖ Problematic argument:
## • by = x
## ℹ Did you misspell an argument name?
```

![](yap1_goanalysis_files/figure-html/grp_vec-1.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//cVEC__hmVEC__iVEC__mVEC_ego_combined_barplot.pdf` 

### ORA ( combined ): Dotplot

Adds gene ratio (fraction of input genes hitting each term) as point size, with colour  for significance - more informative than the barplot since it shows both effect size and  confidence simultaneously. **Expected:** larger, more significant dots clustering around a  coherent biological theme. **Unexpected:** significant terms with very small gene ratios  (few genes driving a 'significant' result) - treat cautiously as potentially fragile.

![](yap1_goanalysis_files/figure-html/grp_vec-2.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//cVEC__hmVEC__iVEC__mVEC_ego_combined_dotplot.pdf` 

### ORA ( combined ): Enrichment map

Network graph connecting GO terms that share overlapping genes - clusters of connected  terms suggest a shared underlying biological process rather than independent findings.  **Expected:** tight clusters of related terms. **Unexpected:** many isolated, unconnected  terms, suggesting the significant terms are thematically scattered rather than coherent.

![](yap1_goanalysis_files/figure-html/grp_vec-3.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//cVEC__hmVEC__iVEC__mVEC_ego_combined_emapplot.pdf` 

### ORA ( combined ): Category-gene network

Shows which specific genes drive which enriched terms, useful for tracing a term back to  the actual DE genes behind it. **Expected:** recognisable marker genes linking to expected  terms. **Unexpected:** a term driven entirely by one or two genes rather than a broad gene  set - weaker evidence than a term supported by many genes.

![](yap1_goanalysis_files/figure-html/grp_vec-4.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//cVEC__hmVEC__iVEC__mVEC_ego_combined_cnetplot.pdf` 

### ORA ( up ): Barplot

Shows the top enriched GO terms ranked by significance, bar length/colour reflecting  adjusted p-value or gene count. **Expected:** terms biologically consistent with the  comparison (e.g. vascular/endothelial terms for EC-related contrasts). **Unexpected:**  dominance by generic housekeeping terms (translation, metabolism) with no thematic  coherence - may indicate a weak or noisy input gene list rather than a real signal.

```
## Warning in fortify(object, showCategory = showCategory, by = x, ...): Arguments in `...` must be used.
## ✖ Problematic argument:
## • by = x
## ℹ Did you misspell an argument name?
```

![](yap1_goanalysis_files/figure-html/grp_vec-5.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//cVEC__hmVEC__iVEC__mVEC_ego_up_barplot.pdf` 

### ORA ( up ): Dotplot

Adds gene ratio (fraction of input genes hitting each term) as point size, with colour  for significance - more informative than the barplot since it shows both effect size and  confidence simultaneously. **Expected:** larger, more significant dots clustering around a  coherent biological theme. **Unexpected:** significant terms with very small gene ratios  (few genes driving a 'significant' result) - treat cautiously as potentially fragile.

![](yap1_goanalysis_files/figure-html/grp_vec-6.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//cVEC__hmVEC__iVEC__mVEC_ego_up_dotplot.pdf` 

### ORA ( up ): Enrichment map

Network graph connecting GO terms that share overlapping genes - clusters of connected  terms suggest a shared underlying biological process rather than independent findings.  **Expected:** tight clusters of related terms. **Unexpected:** many isolated, unconnected  terms, suggesting the significant terms are thematically scattered rather than coherent.

![](yap1_goanalysis_files/figure-html/grp_vec-7.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//cVEC__hmVEC__iVEC__mVEC_ego_up_emapplot.pdf` 

### ORA ( up ): Category-gene network

Shows which specific genes drive which enriched terms, useful for tracing a term back to  the actual DE genes behind it. **Expected:** recognisable marker genes linking to expected  terms. **Unexpected:** a term driven entirely by one or two genes rather than a broad gene  set - weaker evidence than a term supported by many genes.

![](yap1_goanalysis_files/figure-html/grp_vec-8.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//cVEC__hmVEC__iVEC__mVEC_ego_up_cnetplot.pdf` 

### ORA ( down ): Barplot

Shows the top enriched GO terms ranked by significance, bar length/colour reflecting  adjusted p-value or gene count. **Expected:** terms biologically consistent with the  comparison (e.g. vascular/endothelial terms for EC-related contrasts). **Unexpected:**  dominance by generic housekeeping terms (translation, metabolism) with no thematic  coherence - may indicate a weak or noisy input gene list rather than a real signal.

```
## Warning in fortify(object, showCategory = showCategory, by = x, ...): Arguments in `...` must be used.
## ✖ Problematic argument:
## • by = x
## ℹ Did you misspell an argument name?
```

![](yap1_goanalysis_files/figure-html/grp_vec-9.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//cVEC__hmVEC__iVEC__mVEC_ego_down_barplot.pdf` 

### ORA ( down ): Dotplot

Adds gene ratio (fraction of input genes hitting each term) as point size, with colour  for significance - more informative than the barplot since it shows both effect size and  confidence simultaneously. **Expected:** larger, more significant dots clustering around a  coherent biological theme. **Unexpected:** significant terms with very small gene ratios  (few genes driving a 'significant' result) - treat cautiously as potentially fragile.

![](yap1_goanalysis_files/figure-html/grp_vec-10.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//cVEC__hmVEC__iVEC__mVEC_ego_down_dotplot.pdf` 

### ORA ( down ): Enrichment map

Network graph connecting GO terms that share overlapping genes - clusters of connected  terms suggest a shared underlying biological process rather than independent findings.  **Expected:** tight clusters of related terms. **Unexpected:** many isolated, unconnected  terms, suggesting the significant terms are thematically scattered rather than coherent.

![](yap1_goanalysis_files/figure-html/grp_vec-11.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//cVEC__hmVEC__iVEC__mVEC_ego_down_emapplot.pdf` 

### ORA ( down ): Category-gene network

Shows which specific genes drive which enriched terms, useful for tracing a term back to  the actual DE genes behind it. **Expected:** recognisable marker genes linking to expected  terms. **Unexpected:** a term driven entirely by one or two genes rather than a broad gene  set - weaker evidence than a term supported by many genes.

![](yap1_goanalysis_files/figure-html/grp_vec-12.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//cVEC__hmVEC__iVEC__mVEC_ego_down_cnetplot.pdf` 

**Done:** cVEC__hmVEC__iVEC__mVEC 

---


```{.r .fold-hide}
singles_cVEC_out <- run_go(singles_cVEC, data_name="cVEC")
```



## GO analysis: cVEC 

Input genes: 11049 

**ORA input genes:** combined = 2923 &nbsp;&nbsp; up = 2741 &nbsp;&nbsp; down = 182 

**ORA (combined):** 2967 terms &nbsp;&nbsp;  **ORA (up):** 2936 terms &nbsp;&nbsp;  **ORA (down):** 1156 terms

*Full ORA ( combined ) table saved to:* `../output/figure_extended/go//cVEC_ego_combined_full_results.csv` 

*Full ORA ( up ) table saved to:* `../output/figure_extended/go//cVEC_ego_up_full_results.csv` 

*Full ORA ( down ) table saved to:* `../output/figure_extended/go//cVEC_ego_down_full_results.csv` 

**ORA significant terms:** combined = 154 , up = 119 , down = 21 

*Significant-only ORA ( combined ) table saved to:* `../output/figure_extended/go//cVEC_ego_combined_significant_only.csv` 

*Significant-only ORA ( up ) table saved to:* `../output/figure_extended/go//cVEC_ego_up_significant_only.csv` 

*Significant-only ORA ( down ) table saved to:* `../output/figure_extended/go//cVEC_ego_down_significant_only.csv` 

### ORA ( combined ): Barplot

Shows the top enriched GO terms ranked by significance, bar length/colour reflecting  adjusted p-value or gene count. **Expected:** terms biologically consistent with the  comparison (e.g. vascular/endothelial terms for EC-related contrasts). **Unexpected:**  dominance by generic housekeeping terms (translation, metabolism) with no thematic  coherence - may indicate a weak or noisy input gene list rather than a real signal.

```
## Warning in fortify(object, showCategory = showCategory, by = x, ...): Arguments in `...` must be used.
## ✖ Problematic argument:
## • by = x
## ℹ Did you misspell an argument name?
```

![](yap1_goanalysis_files/figure-html/sin_cVEC-1.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//cVEC_ego_combined_barplot.pdf` 

### ORA ( combined ): Dotplot

Adds gene ratio (fraction of input genes hitting each term) as point size, with colour  for significance - more informative than the barplot since it shows both effect size and  confidence simultaneously. **Expected:** larger, more significant dots clustering around a  coherent biological theme. **Unexpected:** significant terms with very small gene ratios  (few genes driving a 'significant' result) - treat cautiously as potentially fragile.

![](yap1_goanalysis_files/figure-html/sin_cVEC-2.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//cVEC_ego_combined_dotplot.pdf` 

### ORA ( combined ): Enrichment map

Network graph connecting GO terms that share overlapping genes - clusters of connected  terms suggest a shared underlying biological process rather than independent findings.  **Expected:** tight clusters of related terms. **Unexpected:** many isolated, unconnected  terms, suggesting the significant terms are thematically scattered rather than coherent.

![](yap1_goanalysis_files/figure-html/sin_cVEC-3.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//cVEC_ego_combined_emapplot.pdf` 

### ORA ( combined ): Category-gene network

Shows which specific genes drive which enriched terms, useful for tracing a term back to  the actual DE genes behind it. **Expected:** recognisable marker genes linking to expected  terms. **Unexpected:** a term driven entirely by one or two genes rather than a broad gene  set - weaker evidence than a term supported by many genes.

![](yap1_goanalysis_files/figure-html/sin_cVEC-4.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//cVEC_ego_combined_cnetplot.pdf` 

### ORA ( up ): Barplot

Shows the top enriched GO terms ranked by significance, bar length/colour reflecting  adjusted p-value or gene count. **Expected:** terms biologically consistent with the  comparison (e.g. vascular/endothelial terms for EC-related contrasts). **Unexpected:**  dominance by generic housekeeping terms (translation, metabolism) with no thematic  coherence - may indicate a weak or noisy input gene list rather than a real signal.

```
## Warning in fortify(object, showCategory = showCategory, by = x, ...): Arguments in `...` must be used.
## ✖ Problematic argument:
## • by = x
## ℹ Did you misspell an argument name?
```

![](yap1_goanalysis_files/figure-html/sin_cVEC-5.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//cVEC_ego_up_barplot.pdf` 

### ORA ( up ): Dotplot

Adds gene ratio (fraction of input genes hitting each term) as point size, with colour  for significance - more informative than the barplot since it shows both effect size and  confidence simultaneously. **Expected:** larger, more significant dots clustering around a  coherent biological theme. **Unexpected:** significant terms with very small gene ratios  (few genes driving a 'significant' result) - treat cautiously as potentially fragile.

![](yap1_goanalysis_files/figure-html/sin_cVEC-6.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//cVEC_ego_up_dotplot.pdf` 

### ORA ( up ): Enrichment map

Network graph connecting GO terms that share overlapping genes - clusters of connected  terms suggest a shared underlying biological process rather than independent findings.  **Expected:** tight clusters of related terms. **Unexpected:** many isolated, unconnected  terms, suggesting the significant terms are thematically scattered rather than coherent.

![](yap1_goanalysis_files/figure-html/sin_cVEC-7.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//cVEC_ego_up_emapplot.pdf` 

### ORA ( up ): Category-gene network

Shows which specific genes drive which enriched terms, useful for tracing a term back to  the actual DE genes behind it. **Expected:** recognisable marker genes linking to expected  terms. **Unexpected:** a term driven entirely by one or two genes rather than a broad gene  set - weaker evidence than a term supported by many genes.

![](yap1_goanalysis_files/figure-html/sin_cVEC-8.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//cVEC_ego_up_cnetplot.pdf` 

### ORA ( down ): Barplot

Shows the top enriched GO terms ranked by significance, bar length/colour reflecting  adjusted p-value or gene count. **Expected:** terms biologically consistent with the  comparison (e.g. vascular/endothelial terms for EC-related contrasts). **Unexpected:**  dominance by generic housekeeping terms (translation, metabolism) with no thematic  coherence - may indicate a weak or noisy input gene list rather than a real signal.

```
## Warning in fortify(object, showCategory = showCategory, by = x, ...): Arguments in `...` must be used.
## ✖ Problematic argument:
## • by = x
## ℹ Did you misspell an argument name?
```

![](yap1_goanalysis_files/figure-html/sin_cVEC-9.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//cVEC_ego_down_barplot.pdf` 

### ORA ( down ): Dotplot

Adds gene ratio (fraction of input genes hitting each term) as point size, with colour  for significance - more informative than the barplot since it shows both effect size and  confidence simultaneously. **Expected:** larger, more significant dots clustering around a  coherent biological theme. **Unexpected:** significant terms with very small gene ratios  (few genes driving a 'significant' result) - treat cautiously as potentially fragile.

![](yap1_goanalysis_files/figure-html/sin_cVEC-10.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//cVEC_ego_down_dotplot.pdf` 

### ORA ( down ): Enrichment map

Network graph connecting GO terms that share overlapping genes - clusters of connected  terms suggest a shared underlying biological process rather than independent findings.  **Expected:** tight clusters of related terms. **Unexpected:** many isolated, unconnected  terms, suggesting the significant terms are thematically scattered rather than coherent.

![](yap1_goanalysis_files/figure-html/sin_cVEC-11.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//cVEC_ego_down_emapplot.pdf` 

### ORA ( down ): Category-gene network

Shows which specific genes drive which enriched terms, useful for tracing a term back to  the actual DE genes behind it. **Expected:** recognisable marker genes linking to expected  terms. **Unexpected:** a term driven entirely by one or two genes rather than a broad gene  set - weaker evidence than a term supported by many genes.

![](yap1_goanalysis_files/figure-html/sin_cVEC-12.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//cVEC_ego_down_cnetplot.pdf` 

**Done:** cVEC 

---


```{.r .fold-hide}
singles_hmVEC_out <- run_go(singles_hmVEC, data_name="hmVEC")
```



## GO analysis: hmVEC 

Input genes: 11245 

**ORA input genes:** combined = 4758 &nbsp;&nbsp; up = 4360 &nbsp;&nbsp; down = 398 

**ORA (combined):** 3058 terms &nbsp;&nbsp;  **ORA (up):** 3027 terms &nbsp;&nbsp;  **ORA (down):** 1571 terms

*Full ORA ( combined ) table saved to:* `../output/figure_extended/go//hmVEC_ego_combined_full_results.csv` 

*Full ORA ( up ) table saved to:* `../output/figure_extended/go//hmVEC_ego_up_full_results.csv` 

*Full ORA ( down ) table saved to:* `../output/figure_extended/go//hmVEC_ego_down_full_results.csv` 

**ORA significant terms:** combined = 323 , up = 319 , down = 10 

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

![](yap1_goanalysis_files/figure-html/sin_hmVEC-1.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//hmVEC_ego_combined_barplot.pdf` 

### ORA ( combined ): Dotplot

Adds gene ratio (fraction of input genes hitting each term) as point size, with colour  for significance - more informative than the barplot since it shows both effect size and  confidence simultaneously. **Expected:** larger, more significant dots clustering around a  coherent biological theme. **Unexpected:** significant terms with very small gene ratios  (few genes driving a 'significant' result) - treat cautiously as potentially fragile.

![](yap1_goanalysis_files/figure-html/sin_hmVEC-2.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//hmVEC_ego_combined_dotplot.pdf` 

### ORA ( combined ): Enrichment map

Network graph connecting GO terms that share overlapping genes - clusters of connected  terms suggest a shared underlying biological process rather than independent findings.  **Expected:** tight clusters of related terms. **Unexpected:** many isolated, unconnected  terms, suggesting the significant terms are thematically scattered rather than coherent.

![](yap1_goanalysis_files/figure-html/sin_hmVEC-3.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//hmVEC_ego_combined_emapplot.pdf` 

### ORA ( combined ): Category-gene network

Shows which specific genes drive which enriched terms, useful for tracing a term back to  the actual DE genes behind it. **Expected:** recognisable marker genes linking to expected  terms. **Unexpected:** a term driven entirely by one or two genes rather than a broad gene  set - weaker evidence than a term supported by many genes.

![](yap1_goanalysis_files/figure-html/sin_hmVEC-4.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//hmVEC_ego_combined_cnetplot.pdf` 

### ORA ( up ): Barplot

Shows the top enriched GO terms ranked by significance, bar length/colour reflecting  adjusted p-value or gene count. **Expected:** terms biologically consistent with the  comparison (e.g. vascular/endothelial terms for EC-related contrasts). **Unexpected:**  dominance by generic housekeeping terms (translation, metabolism) with no thematic  coherence - may indicate a weak or noisy input gene list rather than a real signal.

```
## Warning in fortify(object, showCategory = showCategory, by = x, ...): Arguments in `...` must be used.
## ✖ Problematic argument:
## • by = x
## ℹ Did you misspell an argument name?
```

![](yap1_goanalysis_files/figure-html/sin_hmVEC-5.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//hmVEC_ego_up_barplot.pdf` 

### ORA ( up ): Dotplot

Adds gene ratio (fraction of input genes hitting each term) as point size, with colour  for significance - more informative than the barplot since it shows both effect size and  confidence simultaneously. **Expected:** larger, more significant dots clustering around a  coherent biological theme. **Unexpected:** significant terms with very small gene ratios  (few genes driving a 'significant' result) - treat cautiously as potentially fragile.

![](yap1_goanalysis_files/figure-html/sin_hmVEC-6.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//hmVEC_ego_up_dotplot.pdf` 

### ORA ( up ): Enrichment map

Network graph connecting GO terms that share overlapping genes - clusters of connected  terms suggest a shared underlying biological process rather than independent findings.  **Expected:** tight clusters of related terms. **Unexpected:** many isolated, unconnected  terms, suggesting the significant terms are thematically scattered rather than coherent.

![](yap1_goanalysis_files/figure-html/sin_hmVEC-7.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//hmVEC_ego_up_emapplot.pdf` 

### ORA ( up ): Category-gene network

Shows which specific genes drive which enriched terms, useful for tracing a term back to  the actual DE genes behind it. **Expected:** recognisable marker genes linking to expected  terms. **Unexpected:** a term driven entirely by one or two genes rather than a broad gene  set - weaker evidence than a term supported by many genes.

![](yap1_goanalysis_files/figure-html/sin_hmVEC-8.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//hmVEC_ego_up_cnetplot.pdf` 

### ORA ( down ): Barplot

Shows the top enriched GO terms ranked by significance, bar length/colour reflecting  adjusted p-value or gene count. **Expected:** terms biologically consistent with the  comparison (e.g. vascular/endothelial terms for EC-related contrasts). **Unexpected:**  dominance by generic housekeeping terms (translation, metabolism) with no thematic  coherence - may indicate a weak or noisy input gene list rather than a real signal.

```
## Warning in fortify(object, showCategory = showCategory, by = x, ...): Arguments in `...` must be used.
## ✖ Problematic argument:
## • by = x
## ℹ Did you misspell an argument name?
```

![](yap1_goanalysis_files/figure-html/sin_hmVEC-9.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//hmVEC_ego_down_barplot.pdf` 

### ORA ( down ): Dotplot

Adds gene ratio (fraction of input genes hitting each term) as point size, with colour  for significance - more informative than the barplot since it shows both effect size and  confidence simultaneously. **Expected:** larger, more significant dots clustering around a  coherent biological theme. **Unexpected:** significant terms with very small gene ratios  (few genes driving a 'significant' result) - treat cautiously as potentially fragile.

![](yap1_goanalysis_files/figure-html/sin_hmVEC-10.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//hmVEC_ego_down_dotplot.pdf` 

### ORA ( down ): Enrichment map

Network graph connecting GO terms that share overlapping genes - clusters of connected  terms suggest a shared underlying biological process rather than independent findings.  **Expected:** tight clusters of related terms. **Unexpected:** many isolated, unconnected  terms, suggesting the significant terms are thematically scattered rather than coherent.

![](yap1_goanalysis_files/figure-html/sin_hmVEC-11.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//hmVEC_ego_down_emapplot.pdf` 

### ORA ( down ): Category-gene network

Shows which specific genes drive which enriched terms, useful for tracing a term back to  the actual DE genes behind it. **Expected:** recognisable marker genes linking to expected  terms. **Unexpected:** a term driven entirely by one or two genes rather than a broad gene  set - weaker evidence than a term supported by many genes.

![](yap1_goanalysis_files/figure-html/sin_hmVEC-12.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//hmVEC_ego_down_cnetplot.pdf` 

**Done:** hmVEC 

---


```{.r .fold-hide}
singles_mVEC_out <- run_go(singles_mVEC, data_name="mVEC")
```



## GO analysis: mVEC 

Input genes: 10738 

**ORA input genes:** combined = 3347 &nbsp;&nbsp; up = 3197 &nbsp;&nbsp; down = 150 

**ORA (combined):** 2911 terms &nbsp;&nbsp;  **ORA (up):** 2879 terms &nbsp;&nbsp;  **ORA (down):** 898 terms

*Full ORA ( combined ) table saved to:* `../output/figure_extended/go//mVEC_ego_combined_full_results.csv` 

*Full ORA ( up ) table saved to:* `../output/figure_extended/go//mVEC_ego_up_full_results.csv` 

*Full ORA ( down ) table saved to:* `../output/figure_extended/go//mVEC_ego_down_full_results.csv` 

**ORA significant terms:** combined = 258 , up = 241 , down = 6 

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

![](yap1_goanalysis_files/figure-html/sin_mVEC-1.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//mVEC_ego_combined_barplot.pdf` 

### ORA ( combined ): Dotplot

Adds gene ratio (fraction of input genes hitting each term) as point size, with colour  for significance - more informative than the barplot since it shows both effect size and  confidence simultaneously. **Expected:** larger, more significant dots clustering around a  coherent biological theme. **Unexpected:** significant terms with very small gene ratios  (few genes driving a 'significant' result) - treat cautiously as potentially fragile.

![](yap1_goanalysis_files/figure-html/sin_mVEC-2.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//mVEC_ego_combined_dotplot.pdf` 

### ORA ( combined ): Enrichment map

Network graph connecting GO terms that share overlapping genes - clusters of connected  terms suggest a shared underlying biological process rather than independent findings.  **Expected:** tight clusters of related terms. **Unexpected:** many isolated, unconnected  terms, suggesting the significant terms are thematically scattered rather than coherent.

![](yap1_goanalysis_files/figure-html/sin_mVEC-3.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//mVEC_ego_combined_emapplot.pdf` 

### ORA ( combined ): Category-gene network

Shows which specific genes drive which enriched terms, useful for tracing a term back to  the actual DE genes behind it. **Expected:** recognisable marker genes linking to expected  terms. **Unexpected:** a term driven entirely by one or two genes rather than a broad gene  set - weaker evidence than a term supported by many genes.

![](yap1_goanalysis_files/figure-html/sin_mVEC-4.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//mVEC_ego_combined_cnetplot.pdf` 

### ORA ( up ): Barplot

Shows the top enriched GO terms ranked by significance, bar length/colour reflecting  adjusted p-value or gene count. **Expected:** terms biologically consistent with the  comparison (e.g. vascular/endothelial terms for EC-related contrasts). **Unexpected:**  dominance by generic housekeeping terms (translation, metabolism) with no thematic  coherence - may indicate a weak or noisy input gene list rather than a real signal.

```
## Warning in fortify(object, showCategory = showCategory, by = x, ...): Arguments in `...` must be used.
## ✖ Problematic argument:
## • by = x
## ℹ Did you misspell an argument name?
```

![](yap1_goanalysis_files/figure-html/sin_mVEC-5.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//mVEC_ego_up_barplot.pdf` 

### ORA ( up ): Dotplot

Adds gene ratio (fraction of input genes hitting each term) as point size, with colour  for significance - more informative than the barplot since it shows both effect size and  confidence simultaneously. **Expected:** larger, more significant dots clustering around a  coherent biological theme. **Unexpected:** significant terms with very small gene ratios  (few genes driving a 'significant' result) - treat cautiously as potentially fragile.

![](yap1_goanalysis_files/figure-html/sin_mVEC-6.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//mVEC_ego_up_dotplot.pdf` 

### ORA ( up ): Enrichment map

Network graph connecting GO terms that share overlapping genes - clusters of connected  terms suggest a shared underlying biological process rather than independent findings.  **Expected:** tight clusters of related terms. **Unexpected:** many isolated, unconnected  terms, suggesting the significant terms are thematically scattered rather than coherent.

![](yap1_goanalysis_files/figure-html/sin_mVEC-7.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//mVEC_ego_up_emapplot.pdf` 

### ORA ( up ): Category-gene network

Shows which specific genes drive which enriched terms, useful for tracing a term back to  the actual DE genes behind it. **Expected:** recognisable marker genes linking to expected  terms. **Unexpected:** a term driven entirely by one or two genes rather than a broad gene  set - weaker evidence than a term supported by many genes.

![](yap1_goanalysis_files/figure-html/sin_mVEC-8.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//mVEC_ego_up_cnetplot.pdf` 

### ORA ( down ): Barplot

Shows the top enriched GO terms ranked by significance, bar length/colour reflecting  adjusted p-value or gene count. **Expected:** terms biologically consistent with the  comparison (e.g. vascular/endothelial terms for EC-related contrasts). **Unexpected:**  dominance by generic housekeeping terms (translation, metabolism) with no thematic  coherence - may indicate a weak or noisy input gene list rather than a real signal.

```
## Warning in fortify(object, showCategory = showCategory, by = x, ...): Arguments in `...` must be used.
## ✖ Problematic argument:
## • by = x
## ℹ Did you misspell an argument name?
```

![](yap1_goanalysis_files/figure-html/sin_mVEC-9.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//mVEC_ego_down_barplot.pdf` 

### ORA ( down ): Dotplot

Adds gene ratio (fraction of input genes hitting each term) as point size, with colour  for significance - more informative than the barplot since it shows both effect size and  confidence simultaneously. **Expected:** larger, more significant dots clustering around a  coherent biological theme. **Unexpected:** significant terms with very small gene ratios  (few genes driving a 'significant' result) - treat cautiously as potentially fragile.

![](yap1_goanalysis_files/figure-html/sin_mVEC-10.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//mVEC_ego_down_dotplot.pdf` 

### ORA ( down ): Enrichment map

Network graph connecting GO terms that share overlapping genes - clusters of connected  terms suggest a shared underlying biological process rather than independent findings.  **Expected:** tight clusters of related terms. **Unexpected:** many isolated, unconnected  terms, suggesting the significant terms are thematically scattered rather than coherent.

![](yap1_goanalysis_files/figure-html/sin_mVEC-11.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//mVEC_ego_down_emapplot.pdf` 

### ORA ( down ): Category-gene network

Shows which specific genes drive which enriched terms, useful for tracing a term back to  the actual DE genes behind it. **Expected:** recognisable marker genes linking to expected  terms. **Unexpected:** a term driven entirely by one or two genes rather than a broad gene  set - weaker evidence than a term supported by many genes.

![](yap1_goanalysis_files/figure-html/sin_mVEC-12.png)<!-- -->*PDF saved to:* `../output/figure_extended/go//mVEC_ego_down_cnetplot.pdf` 

**Done:** mVEC 

---


## For developers

Sample run command:


``` bash
Rscript -e "
rmarkdown::render(
  'yap1_goanalysis.Rmd',
  output_file = './yap1_goanalysis.html'
)
"
```

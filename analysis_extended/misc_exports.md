---
title: "Misc data exports"
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

Catch-all for small, non-figure data exports that don't belong in any of the other extended documents. Currently
just one: cell composition counts/percentages per `L3_cluster_id` x `Genotype`, formatted for import into
GraphPad Prism.

## Outputs

Under `../output/figure_extended/misc/`:

- `Meox1_Level_03_L3_Seurat_cluster_predicted_phenotype_Genotype_composition_for_prism.csv`

## Export

Source: `export_prism.R`.


```{.r .fold-hide}
library(dplyr)
```

```
## 
## Attaching package: 'dplyr'
```

```
## The following objects are masked from 'package:stats':
## 
##     filter, lag
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

```{.r .fold-hide}
library(qs2)
```

```
## qs2 0.2.1
```

```{.r .fold-hide}
infile_path <- "../../Saki_data/paper/D10051_meox1_dataset_Level_03_annotated.qs2"
outfile_path <- "../output/figure_extended/misc/Meox1_Level_03_L3_Seurat_cluster_predicted_phenotype_Genotype_composition_for_prism.csv"

data <- qs2::qs_read(infile_path)
table(data@meta.data[["L3_cluster_id"]])
```

```
## 
##      mVEC_01       LEC_01     hmVEC_01 pre_muLEC_01 
##         1384          516          505          170
```

```{.r .fold-hide}
for_prism <- data@meta.data %>%
    mutate(Genotype = recode(Genotype,
        "wildtype" = "WT",
        "meox1_mutant" = "Mutant"
    )) %>%
    group_by(Genotype, L3_cluster_id) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(Genotype) %>%
    mutate(percent = count / sum(count) * 100) %>%
    ungroup()

for_prism
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["Genotype"],"name":[1],"type":["chr"],"align":["left"]},{"label":["L3_cluster_id"],"name":[2],"type":["fct"],"align":["left"]},{"label":["count"],"name":[3],"type":["int"],"align":["right"]},{"label":["percent"],"name":[4],"type":["dbl"],"align":["right"]}],"data":[{"1":"Mutant","2":"mVEC_01","3":"796","4":"61.801242"},{"1":"Mutant","2":"LEC_01","3":"231","4":"17.934783"},{"1":"Mutant","2":"hmVEC_01","3":"168","4":"13.043478"},{"1":"Mutant","2":"pre_muLEC_01","3":"93","4":"7.220497"},{"1":"WT","2":"mVEC_01","3":"588","4":"45.687646"},{"1":"WT","2":"LEC_01","3":"285","4":"22.144522"},{"1":"WT","2":"hmVEC_01","3":"337","4":"26.184926"},{"1":"WT","2":"pre_muLEC_01","3":"77","4":"5.982906"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

```{.r .fold-hide}
write.csv(
    for_prism, outfile_path, row.names = FALSE
    )
```

## For developers

Sample run command:


``` bash
Rscript -e "
rmarkdown::render(
  'misc_exports.Rmd',
  output_file = './misc_exports.html'
)
"
```

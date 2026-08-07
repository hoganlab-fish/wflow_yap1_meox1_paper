---
title: 'Saki et al list of extended plots'
subtitle: ""
author: "Tyrone Chen, Michelle Meier"
date: 'August 06, 2026'
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
    toc_title: 'Saki et al list of extended plots'
editor_options:
  chunk_output_type: console
params:
  sample_index: 1
  sample_name: 'Saki et al list of extended plots'
---

<a href="https://opensource.org/licenses/MIT">
    <img src="https://img.shields.io/badge/License-MIT-yellow.svg?style=for-the-badge" alt="MIT License">
</a>

<a href="https://creativecommons.org/licenses/by/4.0/">
    <img src="https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg?style=for-the-badge" alt="CC BY 4.0">
</a>

<a href="https://biomedicalsciences.unimelb.edu.au/sbs-research-groups/anatomy-and-physiology-research/stem-cell-and-developmental-biology/hogan-laboratory-vascular-developmental-genetics-and-cell-biology">
    <img src="https://img.shields.io/badge/Website-hoganlab-4285F4?style=for-the-badge&logo=google-chrome&logoColor=white" alt="Website">
</a>

Copyright © 2026 <a href="https://orcid.org/0000-0002-9207-0385">Tyrone Chen <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>, <a href="https://orcid.org/0009-0005-5595-3882">Michelle Meier <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>

# README

This is a record of a collection of plots made separately from the original workflow. In other words, a figure in the manuscript will be covered either by one of the documents linked below or by `wflow_yap1_meox1_paper/docs/figures_meox1_dataset.html` / `figures_yap1_dataset.html`. 

> **_NOTE:_** You may see extra plots which form part of the exploratory analysis, and are not part of the manuscript. Any final plots in the manuscript may have been realigned, resized or reannotated through Adobe Illustrator, but the content will remain unchanged.

# Overview

To minimise cognitive load we used subdirectory naming conventions which have functional (but not conceptual)
similarity to `workflowr`'s layout:

- `wflow_yap1_meox1_paper/analysis_extended/` contains the code associated with the dataset (this directory).
- `wflow_yap1_meox1_paper/output/figure_extended/` contains the output figures (and, where the source script mixed
  the two, their backing data tables) associated with the dataset.

This is a **side channel, not part of the official `workflowr` build** (`_workflowr.yml` still only points at
`analysis/`) - render each document manually with the command in its own "For developers" section, or via
`rmarkdown::render_site()` scoped to this directory if you set one up later.

## Recommended run order

Several documents depend on CSVs produced by others. Suggested order:

1. `meox1_deg_confects_prep.Rmd`, `yap1_deg_confects_prep.Rmd` - generates the merged-group DEG/topconfects
   tables everything else reads.
2. `meox1_embo_lecvec_scores.Rmd` - regenerates the EMBO marker geneset CSVs consumed by step 3.
3. Everything else (GO analyses, dotplots, volcano plots, UMAPs, correlation plots) - order between these
   doesn't matter, they're independent of each other.
4. `velocyto_analysis.Rmd`, `meox1_slingshot_pseudotime.Rmd` - independent, but by far the most compute-heavy;
   run last / separately, expect long runtimes.

# Figures

## Non-pseudotime

| Document | Topic |
|---|---|
| [meox1_deg_confects_prep.html](meox1_deg_confects_prep.html) | meox1 DEG + topconfects prep (Level 3 merged groups, Level 1 fibroblasts) |
| [yap1_deg_confects_prep.html](yap1_deg_confects_prep.html) | yap1 DEG + topconfects prep (Level 3 merged groups) |
| [meox1_goanalysis.html](meox1_goanalysis.html) | meox1 GO analysis (ORA, raw p-value) |
| [meox1_goanalysis_adjpval.html](meox1_goanalysis_adjpval.html) | meox1 GO analysis (ORA, adjusted p-value variant) |
| [meox1_L1_fibroblast_goanalysis.html](meox1_L1_fibroblast_goanalysis.html) | meox1 GO analysis, L1 fibroblasts (FGSEA + ORA + UpSet) |
| [yap1_goanalysis.html](yap1_goanalysis.html) | yap1 GO analysis (ORA, raw p-value) |
| [meox1_dotplot_goi.html](meox1_dotplot_goi.html) | meox1 genes-of-interest dotplot |
| [meox1_dotplot_panels.html](meox1_dotplot_panels.html) | meox1 cdh5 + cell-cycle/p53 dotplot panels |
| [yap1_dotplot_panels.html](yap1_dotplot_panels.html) | yap1 meox1-across-levels, L3 DEG ratio, cell-cycle/p53 dotplot panels |
| [meox1_volcano_plots.html](meox1_volcano_plots.html) | meox1 hmVEC/mVEC volcano plots |
| [yap1_volcano_plots.html](yap1_volcano_plots.html) | yap1 hmVEC volcano plot + iVEC GO term barplot |
| [meox1_umap_feature_plots.html](meox1_umap_feature_plots.html) | meox1 cdh5-by-genotype UMAP + pathway score UMAPs |
| [yap1_umap_feature_plots.html](yap1_umap_feature_plots.html) | yap1 cdh5-by-genotype, meox1-across-levels, iVEC-highlight UMAPs |
| [meox1_embo_lecvec_scores.html](meox1_embo_lecvec_scores.html) | meox1 EMBO-method LEC/VEC marker scores, split by genotype |
| [meox1_L3_correlation_logfc.html](meox1_L3_correlation_logfc.html) | meox1 LEC-vs-VEC log2FC correlation (WT vs mutant concordance) |
| [misc_exports.html](misc_exports.html) | Misc data exports (Prism composition table) |

<div class="figure" style="text-align: center">
<embed src="../output/figure_extended/embo_scores/embo_transition_scores_split_genotype_dotplot.pdf" title="Example: meox1 EMBO transition score dotplot, split by genotype" type="application/pdf" />
<p class="caption">Example: meox1 EMBO transition score dotplot, split by genotype</p>
</div>

## Pseudotime

| Document | Topic |
|---|---|
| [velocyto_analysis.html](velocyto_analysis.html) | RNA velocity (velocyto.R) for meox1 and yap1, plus unused exon/intron QC comparison |
| [meox1_slingshot_pseudotime.html](meox1_slingshot_pseudotime.html) | meox1 slingshot pseudotime + tradeSeq divergence testing (heaviest document here) |

<div class="figure" style="text-align: center">
<embed src="../output/figure_extended/slingshot/meox1__FeaLabel_Pseudotime_Contours_start_hmVEC_lin1.D10051.pdf" title="Example: meox1 slingshot pseudotime contours, hmVEC lineage 1" type="application/pdf" />
<p class="caption">Example: meox1 slingshot pseudotime contours, hmVEC lineage 1</p>
</div>

# For developers

## Non-pseudotime

Each document above is self-contained and can be rendered independently once its own prerequisites (see
"Recommended run order") have been run.

## Pseudotime

`meox1_slingshot_pseudotime.Rmd` and `velocyto_analysis.Rmd` have no cross-document dependencies but are the
slowest by a wide margin (the former fits a `tradeSeq` GAM per lineage per gene subset, across every candidate
start celltype, for the combined dataset and both samples individually). Budget accordingly and consider running
them on a machine with more cores/memory than a laptop.

## Command line arguments

Rerun everything, in the order from "Recommended run order" above (non-pseudotime documents first, pseudotime
last), with:


```{.bash .fold-hide}
Rscript -e "rmarkdown::render('meox1_deg_confects_prep.Rmd', output_file = 'meox1_deg_confects_prep.html')"
Rscript -e "rmarkdown::render('yap1_deg_confects_prep.Rmd', output_file = 'yap1_deg_confects_prep.html')"
Rscript -e "rmarkdown::render('meox1_embo_lecvec_scores.Rmd', output_file = 'meox1_embo_lecvec_scores.html')"
Rscript -e "rmarkdown::render('meox1_goanalysis.Rmd', output_file = 'meox1_goanalysis.html')"
Rscript -e "rmarkdown::render('meox1_goanalysis_adjpval.Rmd', output_file = 'meox1_goanalysis_adjpval.html')"
Rscript -e "rmarkdown::render('meox1_L1_fibroblast_goanalysis.Rmd', output_file = 'meox1_L1_fibroblast_goanalysis.html')"
Rscript -e "rmarkdown::render('yap1_goanalysis.Rmd', output_file = 'yap1_goanalysis.html')"
Rscript -e "rmarkdown::render('meox1_dotplot_goi.Rmd', output_file = 'meox1_dotplot_goi.html')"
Rscript -e "rmarkdown::render('meox1_dotplot_panels.Rmd', output_file = 'meox1_dotplot_panels.html')"
Rscript -e "rmarkdown::render('yap1_dotplot_panels.Rmd', output_file = 'yap1_dotplot_panels.html')"
Rscript -e "rmarkdown::render('meox1_volcano_plots.Rmd', output_file = 'meox1_volcano_plots.html')"
Rscript -e "rmarkdown::render('yap1_volcano_plots.Rmd', output_file = 'yap1_volcano_plots.html')"
Rscript -e "rmarkdown::render('meox1_umap_feature_plots.Rmd', output_file = 'meox1_umap_feature_plots.html')"
Rscript -e "rmarkdown::render('yap1_umap_feature_plots.Rmd', output_file = 'yap1_umap_feature_plots.html')"
Rscript -e "rmarkdown::render('meox1_L3_correlation_logfc.Rmd', output_file = 'meox1_L3_correlation_logfc.html')"
Rscript -e "rmarkdown::render('misc_exports.Rmd', output_file = 'misc_exports.html')"

# warning - run this separately
Rscript -e "rmarkdown::render('velocyto_analysis.Rmd', output_file = 'velocyto_analysis.html')"
Rscript -e "rmarkdown::render('meox1_slingshot_pseudotime.Rmd', output_file = 'meox1_slingshot_pseudotime.html')"

Rscript -e "rmarkdown::render('README_extended.Rmd', output_file = 'README_extended.html')"
```

## Dependencies

> **_NOTE:_** Software versions are not in perfect alignment across the plots described here to the other plots
> shown in the manuscript due to a series of technical constraints. This applies in particular to `velocyto`
> pseudotime plots, since it has extremely specific install quirks. While little functional impact is expected,
> the software versions used are still recorded for full transparency and reproducibility.

Two conda environment specs, exported from the machines these documents were actually rendered on:

<details>
<summary>non_pseudotime_and_slingshot.yml (527 lines) - click to expand</summary>

Environment `saki` - covers every non-pseudotime document plus `meox1_slingshot_pseudotime.Rmd` (slingshot/tradeSeq's own install quirks are captured in this one environment, not a separate one).

```yaml
name: saki
channels:
  - bioconda
  - conda-forge
dependencies:
  - _openmp_mutex=4.5=20_gnu
  - _r-mutex=1.0.1=anacondar_1
  - argcomplete=3.7.0=pyhcf101f3_0
  - binutils_impl_linux-64=2.45.1=default_hfdba357_101
  - bioconductor-annotationdbi=1.72.0=r45hdfd78af_0
  - bioconductor-assorthead=1.4.0=r45hdfd78af_0
  - bioconductor-beachmat=2.26.0=r45ha27e39d_0
  - bioconductor-biobase=2.70.0=r45h01b2380_0
  - bioconductor-biocgenerics=0.56.0=r45hdfd78af_2
  - bioconductor-biocio=1.20.0=r45hdfd78af_0
  - bioconductor-biocneighbors=2.4.0=r45ha27e39d_0
  - bioconductor-biocparallel=1.44.0=r45ha27e39d_1
  - bioconductor-biocsingular=1.26.1=r45ha27e39d_0
  - bioconductor-biostrings=2.78.0=r45h01b2380_0
  - bioconductor-bluster=1.20.0=r45ha27e39d_0
  - bioconductor-cigarillo=1.0.0=r45h01b2380_0
  - bioconductor-clusterprofiler=4.18.4=r45hdfd78af_0
  - bioconductor-data-packages=20260207=hdfd78af_0
  - bioconductor-delayedarray=0.36.0=r45h01b2380_0
  - bioconductor-delayedmatrixstats=1.32.0=r45hdfd78af_0
  - bioconductor-dose=4.4.0=r45hdfd78af_0
  - bioconductor-edger=4.8.2=r45h01b2380_0
  - bioconductor-enrichplot=1.30.4=r45hdfd78af_0
  - bioconductor-fgsea=1.36.2=r45ha27e39d_0
  - bioconductor-genomeinfodb=1.46.2=r45hdfd78af_0
  - bioconductor-genomicalignments=1.46.0=r45h01b2380_0
  - bioconductor-genomicranges=1.62.1=r45h01b2380_0
  - bioconductor-ggtree=4.0.4=r45hdfd78af_1
  - bioconductor-go.db=3.22.0=r45hdfd78af_0
  - bioconductor-gosemsim=2.36.0=r45ha27e39d_0
  - bioconductor-iranges=2.44.0=r45h01b2380_1
  - bioconductor-keggrest=1.50.0=r45hdfd78af_0
  - bioconductor-limma=3.66.0=r45h01b2380_0
  - bioconductor-matrixgenerics=1.22.0=r45hdfd78af_1
  - bioconductor-metapod=1.18.0=r45ha27e39d_0
  - bioconductor-qvalue=2.42.0=r45hdfd78af_0
  - bioconductor-rhtslib=3.6.0=r45h01b2380_0
  - bioconductor-rsamtools=2.26.0=r45ha27e39d_0
  - bioconductor-rtracklayer=1.70.1=r45h01b2380_0
  - bioconductor-s4arrays=1.10.1=r45h01b2380_0
  - bioconductor-s4vectors=0.48.0=r45h01b2380_1
  - bioconductor-scaledmatrix=1.18.0=r45hdfd78af_0
  - bioconductor-scater=1.38.0=r45hdfd78af_0
  - bioconductor-scdblfinder=1.24.0=r45hdfd78af_0
  - bioconductor-scran=1.38.0=r45ha27e39d_0
  - bioconductor-scuttle=1.20.0=r45ha27e39d_0
  - bioconductor-seqinfo=1.0.0=r45hdfd78af_0
  - bioconductor-singlecellexperiment=1.32.0=r45hdfd78af_0
  - bioconductor-slingshot=2.18.0=r45hdfd78af_0
  - bioconductor-sparsearray=1.10.8=r45h01b2380_0
  - bioconductor-sparsematrixstats=1.22.0=r45ha27e39d_0
  - bioconductor-summarizedexperiment=1.40.0=r45hdfd78af_0
  - bioconductor-tradeseq=1.24.0=r45hdfd78af_0
  - bioconductor-trajectoryutils=1.18.0=r45hdfd78af_0
  - bioconductor-treeio=1.34.0=r45hdfd78af_0
  - bioconductor-ucsc.utils=1.6.1=r45hdfd78af_0
  - bioconductor-xvector=0.50.0=r45h01b2380_0
  - bwidget=1.10.1=ha770c72_1
  - bzip2=1.0.8=hda65f42_9
  - c-ares=1.34.6=hb03c661_0
  - ca-certificates=2026.2.25=hbd8a1cb_0
  - cairo=1.18.4=he90730b_1
  - cmake=4.2.3=hc85cc9f_1
  - curl=8.18.0=hcf29cc6_1
  - font-ttf-dejavu-sans-mono=2.37=hab24e00_0
  - font-ttf-inconsolata=3.000=h77eed37_0
  - font-ttf-source-code-pro=2.038=h77eed37_0
  - font-ttf-ubuntu=0.83=h77eed37_3
  - fontconfig=2.17.1=h27c8c51_0
  - fonts-conda-ecosystem=1=0
  - fonts-conda-forge=1=hc364b38_1
  - fribidi=1.0.16=hb03c661_0
  - gcc_impl_linux-64=15.2.0=he420e7e_18
  - gfortran_impl_linux-64=15.2.0=h281d09f_18
  - glpk=5.0=h445213a_0
  - gmp=6.3.0=hac33072_2
  - graphite2=1.3.14=hecca717_2
  - gsl=2.7=he838d99_0
  - gxx_impl_linux-64=15.2.0=hda75c37_18
  - harfbuzz=12.3.2=h6083320_0
  - icu=78.2=h33c6efd_0
  - kernel-headers_linux-64=5.14.0=he073ed8_3
  - keyutils=1.6.3=hb9d3cd8_0
  - krb5=1.22.2=ha1258a1_0
  - ld_impl_linux-64=2.45.1=default_hbd61a6d_101
  - lerc=4.0.0=h0aef613_1
  - libblas=3.11.0=5_h4a7cf45_openblas
  - libcblas=3.11.0=5_h0358290_openblas
  - libcurl=8.18.0=hcf29cc6_1
  - libdeflate=1.25=h17f619e_0
  - libedit=3.1.20250104=pl5321h7949ede_0
  - libev=4.33=hd590300_2
  - libexpat=2.7.4=hecca717_0
  - libffi=3.5.2=h3435931_0
  - libfreetype=2.14.1=ha770c72_0
  - libfreetype6=2.14.1=h73754d4_0
  - libgcc=15.2.0=he0feb66_18
  - libgcc-devel_linux-64=15.2.0=hcc6f6b0_118
  - libgcc-ng=15.2.0=h69a702a_18
  - libgfortran=15.2.0=h69a702a_18
  - libgfortran-ng=15.2.0=h69a702a_18
  - libgfortran5=15.2.0=h68bc16d_18
  - libgit2=1.9.2=hc20babb_0
  - libglib=2.86.4=h6548e54_1
  - libgomp=15.2.0=he0feb66_18
  - libhwloc=2.12.2=default_hafda6a7_1000
  - libiconv=1.18=h3b78370_2
  - libjpeg-turbo=3.1.2=hb03c661_0
  - liblapack=3.11.0=5_h47877c9_openblas
  - liblzma=5.8.2=hb03c661_0
  - libmpdec=4.0.0=hb03c661_1
  - libnghttp2=1.67.0=had1ee68_0
  - libopenblas=0.3.30=pthreads_h94d23a6_4
  - libpng=1.6.55=h421ea60_0
  - libsanitizer=15.2.0=h90f66d4_18
  - libsqlite=3.52.0=hf4e2dac_0
  - libssh2=1.11.1=hcf80075_0
  - libstdcxx=15.2.0=h934c35e_18
  - libstdcxx-devel_linux-64=15.2.0=hd446a21_118
  - libstdcxx-ng=15.2.0=hdf11a46_18
  - libtiff=4.7.1=h9d88235_1
  - libuuid=2.41.3=h5347b49_0
  - libuv=1.51.0=hb03c661_1
  - libwebp-base=1.6.0=hd42ef1d_0
  - libxcb=1.17.0=h8a09558_0
  - libxgboost=3.2.0=cpu_h2ebb00f_1
  - libxml2=2.15.2=he237659_0
  - libxml2-16=2.15.2=hca6bf5a_0
  - libxml2-devel=2.15.2=he237659_0
  - libzlib=1.3.1=hb9d3cd8_2
  - make=4.4.1=hb9d3cd8_2
  - ncurses=6.5=h2d0b736_3
  - nlopt=2.10.1=np2py314h6477eea_2
  - numpy=2.4.2=py314h2b28147_1
  - openssl=3.6.1=h35e630c_1
  - pandoc=3.9=ha770c72_0
  - pango=1.56.4=hadf4263_0
  - pcre2=10.47=haa7fec5_0
  - pip=26.0.1=pyh145f28c_0
  - pixman=0.46.4=h54a6638_1
  - pthread-stubs=0.4=hb9d3cd8_1002
  - python=3.14.3=h32b2ec7_101_cp314
  - python_abi=3.14=8_cp314
  - pyyaml=6.0.3=py314h67df5f8_1
  - r-abind=1.4_8=r45hc72bb7e_1
  - r-ape=5.8_1=r45h3704496_2
  - r-aplot=0.3.1=r45hc72bb7e_0
  - r-askpass=1.2.1=r45h54b55ab_1
  - r-assertthat=0.2.1=r45hc72bb7e_6
  - r-backports=1.5.0=r45h54b55ab_2
  - r-base=4.5.2=h1fbe982_4
  - r-base64enc=0.1_6=r45h54b55ab_0
  - r-beeswarm=0.4.0=r45h54b55ab_5
  - r-bh=1.90.0_1=r45hc72bb7e_0
  - r-bit=4.6.0=r45h54b55ab_1
  - r-bit64=4.6.0_1=r45h54b55ab_1
  - r-bitops=1.0_9=r45h54b55ab_1
  - r-blob=1.3.0=r45hc72bb7e_0
  - r-boot=1.3_32=r45hc72bb7e_1
  - r-brew=1.0_10=r45hc72bb7e_2
  - r-brio=1.1.5=r45h54b55ab_2
  - r-broom=1.0.12=r45hc72bb7e_0
  - r-bslib=0.10.0=r45hc72bb7e_0
  - r-cachem=1.1.0=r45h54b55ab_2
  - r-cairo=1.7_0=r45h0057c2c_1
  - r-callr=3.7.6=r45hc72bb7e_2
  - r-car=3.1_5=r45hc72bb7e_0
  - r-cardata=3.0_6=r45hc72bb7e_0
  - r-caret=7.0_1=r45h54b55ab_0
  - r-catools=1.18.3=r45h3697838_1
  - r-cellranger=1.1.0=r45hc72bb7e_1008
  - r-checkmate=2.3.3=r45h54b55ab_1
  - r-class=7.3_23=r45h54b55ab_1
  - r-cli=3.6.5=r45h3697838_1
  - r-clipr=0.8.0=r45hc72bb7e_4
  - r-clock=0.7.4=r45h3697838_0
  - r-cluster=2.1.8.2=r45heaba542_0
  - r-clustree=0.5.1=r45hc72bb7e_2
  - r-codetools=0.2_20=r45hc72bb7e_2
  - r-colorspace=2.1_2=r45h54b55ab_0
  - r-commonmark=2.0.0=r45h54b55ab_1
  - r-conflicted=1.2.0=r45h785f33e_3
  - r-conquer=1.3.3=r45h3704496_5
  - r-corrplot=0.95=r45hc72bb7e_1
  - r-cowplot=1.2.0=r45hc72bb7e_2
  - r-cpp11=0.5.3=r45h785f33e_0
  - r-crayon=1.5.3=r45hc72bb7e_2
  - r-credentials=2.0.3=r45hc72bb7e_1
  - r-crosstalk=1.2.2=r45hc72bb7e_1
  - r-crul=1.6.0=r45hc72bb7e_1
  - r-curl=7.0.0=r45h10955f1_1
  - r-data.table=1.17.8=r45h1c8cec4_1
  - r-dbi=1.3.0=r45hc72bb7e_0
  - r-dbplyr=2.5.2=r45hc72bb7e_0
  - r-deldir=2.0_4=r45heaba542_2
  - r-dendsort=0.3.4=r45ha770c72_4
  - r-deriv=4.2.0=r45hc72bb7e_1
  - r-desc=1.4.3=r45hc72bb7e_2
  - r-devtools=2.4.6=r45hc72bb7e_0
  - r-diagram=1.6.5=r45ha770c72_4
  - r-diffobj=0.3.6=r45h54b55ab_1
  - r-digest=0.6.39=r45h3697838_0
  - r-doby=4.7.1=r45hc72bb7e_0
  - r-dotcall64=1.2=r45heaba542_1
  - r-downlit=0.4.5=r45hc72bb7e_0
  - r-dplyr=1.2.0=r45h3697838_0
  - r-dqrng=0.3.2=r45h3697838_2
  - r-drat=0.2.5=r45hc72bb7e_1
  - r-dtplyr=1.3.3=r45hc72bb7e_0
  - r-e1071=1.7_17=r45h3697838_0
  - r-ellipsis=0.3.2=r45h54b55ab_4
  - r-evaluate=1.0.5=r45hc72bb7e_1
  - r-fansi=1.0.7=r45h54b55ab_0
  - r-farver=2.1.2=r45h3697838_2
  - r-fastcluster=1.3.0=r45h3697838_1
  - r-fastdummies=1.7.5=r45hc72bb7e_1
  - r-fastmap=1.2.0=r45h3697838_2
  - r-fastmatch=1.1_8=r45h54b55ab_0
  - r-fitdistrplus=1.2_6=r45hc72bb7e_0
  - r-fnn=1.1.4.1=r45h3697838_2
  - r-fontawesome=0.5.3=r45hc72bb7e_1
  - r-fontbitstreamvera=0.1.1=r45hc72bb7e_1007
  - r-fontliberation=0.1.0=r45hc72bb7e_1007
  - r-fontquiver=0.2.1=r45hc72bb7e_1007
  - r-forcats=1.0.1=r45hc72bb7e_0
  - r-foreach=1.5.2=r45hc72bb7e_4
  - r-forecast=9.0.1=r45hf1899b2_0
  - r-formatr=1.14=r45hc72bb7e_3
  - r-formula=1.2_5=r45hc72bb7e_3
  - r-fracdiff=1.5_3=r45h7b2432b_2
  - r-fs=1.6.6=r45h3697838_1
  - r-futile.logger=1.4.9=r45hc72bb7e_0
  - r-futile.options=1.0.1=r45hc72bb7e_1006
  - r-future=1.70.0=r45hc72bb7e_0
  - r-future.apply=1.20.2=r45hc72bb7e_0
  - r-gargle=1.6.1=r45h785f33e_0
  - r-gbrd=0.4.12=r45hc72bb7e_2
  - r-gdtools=0.5.0=r45ha0c56d8_0
  - r-generics=0.1.4=r45hc72bb7e_1
  - r-gert=2.3.1=r45h5e22a44_0
  - r-getpass=0.2_4=r45h54b55ab_2
  - r-gfonts=0.2.0=r45hc72bb7e_3
  - r-ggbeeswarm=0.7.3=r45hc72bb7e_0
  - r-ggforce=0.5.0=r45h3697838_1
  - r-ggfun=0.2.1=r45hc72bb7e_0
  - r-ggiraph=0.9.5=r45h8ff94db_0
  - r-ggnewscale=0.5.2=r45hc72bb7e_1
  - r-ggplot2=4.0.2=r45h785f33e_0
  - r-ggplotify=0.1.3=r45hc72bb7e_1
  - r-ggpubr=0.6.3=r45hc72bb7e_0
  - r-ggraph=2.2.2=r45h3697838_0
  - r-ggrastr=1.0.2=r45hc72bb7e_3
  - r-ggrepel=0.9.6=r45h3697838_2
  - r-ggridges=0.5.7=r45hc72bb7e_1
  - r-ggsci=4.2.0=r45hc72bb7e_0
  - r-ggsignif=0.6.4=r45hc72bb7e_3
  - r-ggtangle=0.1.2=r45hc72bb7e_0
  - r-gh=1.5.0=r45hc72bb7e_1
  - r-git2r=0.35.0=r45h08ddeb0_1
  - r-gitcreds=0.1.2=r45hc72bb7e_4
  - r-globals=0.19.1=r45hc72bb7e_0
  - r-glue=1.8.0=r45h54b55ab_1
  - r-goftest=1.2_3=r45h54b55ab_4
  - r-googledrive=2.1.2=r45hc72bb7e_1
  - r-googlesheets4=1.1.2=r45h785f33e_1
  - r-gower=1.0.2=r45h54b55ab_0
  - r-gplots=3.3.0=r45hc72bb7e_0
  - r-graphlayouts=1.2.3=r45h3697838_0
  - r-gridextra=2.3=r45hc72bb7e_1007
  - r-gridgraphics=0.5_1=r45hc72bb7e_4
  - r-gson=0.2.0=r45hc72bb7e_0
  - r-gtable=0.3.6=r45hc72bb7e_1
  - r-gtools=3.9.5=r45h54b55ab_2
  - r-hardhat=1.4.2=r45hc72bb7e_1
  - r-haven=2.5.5=r45h6d565e7_1
  - r-here=1.0.2=r45hc72bb7e_0
  - r-hexbin=1.28.5=r45heaba542_1
  - r-highr=0.11=r45hc72bb7e_2
  - r-hms=1.1.4=r45hc72bb7e_0
  - r-htmltools=0.5.9=r45h3697838_0
  - r-htmlwidgets=1.6.4=r45h785f33e_4
  - r-httpcode=0.3.0=r45ha770c72_5
  - r-httpuv=1.6.16=r45h6d565e7_1
  - r-httr=1.4.8=r45hc72bb7e_0
  - r-httr2=1.2.2=r45hc72bb7e_0
  - r-ica=1.0_3=r45hc72bb7e_4
  - r-ids=1.0.1=r45hc72bb7e_5
  - r-igraph=2.1.4=r45hf411e2a_2
  - r-ini=0.3.1=r45hc72bb7e_1007
  - r-ipred=0.9_15=r45h54b55ab_2
  - r-irlba=2.3.7=r45h0e4624f_0
  - r-isoband=0.3.0=r45h3697838_0
  - r-iterators=1.0.14=r45hc72bb7e_4
  - r-jquerylib=0.1.4=r45hc72bb7e_4
  - r-jsonlite=2.0.0=r45h54b55ab_1
  - r-kernsmooth=2.23_26=r45ha0a88a1_1
  - r-knitr=1.51=r45hc72bb7e_0
  - r-labeling=0.4.3=r45hc72bb7e_2
  - r-lambda.r=1.2.4=r45hc72bb7e_5
  - r-later=1.4.7=r45h3697838_0
  - r-lattice=0.22_9=r45h54b55ab_0
  - r-lava=1.8.2=r45hc72bb7e_0
  - r-lazyeval=0.2.2=r45h54b55ab_6
  - r-leidenbase=0.1.36=r45ha11a66c_0
  - r-lifecycle=1.0.5=r45hc72bb7e_0
  - r-listenv=0.10.1=r45hc72bb7e_0
  - r-lme4=1.1_38=r45h3697838_0
  - r-lmtest=0.9_40=r45heaba542_4
  - r-locfit=1.5_9.12=r45h54b55ab_1
  - r-lubridate=1.9.5=r45h54b55ab_0
  - r-magrittr=2.0.4=r45h54b55ab_0
  - r-mass=7.3_65=r45h54b55ab_0
  - r-matrix=1.7_4=r45h0e4624f_1
  - r-matrixmodels=0.5_4=r45hc72bb7e_1
  - r-matrixstats=1.5.0=r45h54b55ab_1
  - r-memoise=2.0.1=r45hc72bb7e_4
  - r-mgcv=1.9_4=r45h0e4624f_0
  - r-microbenchmark=1.5.0=r45h54b55ab_1
  - r-mime=0.13=r45h54b55ab_1
  - r-miniui=0.1.2=r45hc72bb7e_1
  - r-minqa=1.2.8=r45ha36cffa_2
  - r-modelmetrics=1.2.2.2=r45h3697838_5
  - r-modelr=0.1.11=r45hc72bb7e_3
  - r-munsell=0.5.1=r45hc72bb7e_2
  - r-n2r=1.0.4=r45h3697838_0
  - r-nlme=3.1_168=r45heaba542_1
  - r-nloptr=2.2.1=r45h8ae9fae_1
  - r-nnet=7.3_20=r45h54b55ab_1
  - r-numderiv=2016.8_1.1=r45hc72bb7e_7
  - r-openssl=2.3.5=r45h68c19f5_0
  - r-otel=0.2.0=r45hc72bb7e_1
  - r-pagoda2=1.0.14=r45h3704496_0
  - r-parallelly=1.46.1=r45h54b55ab_0
  - r-patchwork=1.3.2=r45hc72bb7e_1
  - r-pbapply=1.7_4=r45hc72bb7e_1
  - r-pbkrtest=0.5.5=r45hc72bb7e_1
  - r-pbmcapply=1.5.1=r45h54b55ab_5
  - r-pheatmap=1.0.13=r45hc72bb7e_1
  - r-pillar=1.11.1=r45hc72bb7e_0
  - r-pkgbuild=1.4.8=r45hc72bb7e_1
  - r-pkgconfig=2.0.3=r45hc72bb7e_5
  - r-pkgdown=2.2.0=r45hc72bb7e_0
  - r-pkgload=1.5.0=r45hc72bb7e_0
  - r-plogr=0.2.0=r45hc72bb7e_1007
  - r-plotly=4.12.0=r45hc72bb7e_0
  - r-plyr=1.8.9=r45h3697838_3
  - r-png=0.1_8=r45h6b2d295_3
  - r-polyclip=1.10_7=r45h3697838_1
  - r-polynom=1.4_1=r45hc72bb7e_4
  - r-praise=1.0.0=r45hc72bb7e_1009
  - r-prettyunits=1.2.0=r45hc72bb7e_2
  - r-princurve=2.1.6=r45h3697838_4
  - r-proc=1.19.0.1=r45h3697838_1
  - r-processx=3.8.6=r45h54b55ab_1
  - r-prodlim=2026.03.11=r45h3697838_0
  - r-profvis=0.4.0=r45h54b55ab_1
  - r-progress=1.2.3=r45hc72bb7e_2
  - r-progressr=0.18.0=r45hc72bb7e_0
  - r-promises=1.5.0=r45hc72bb7e_1
  - r-proxy=0.4_29=r45h54b55ab_0
  - r-ps=1.9.1=r45h54b55ab_1
  - r-purrr=1.2.1=r45h54b55ab_0
  - r-qs2=0.1.7=r45h1a6caf7_0
  - r-quadprog=1.5_8=r45ha0a88a1_7
  - r-quantmod=0.4.28=r45hc72bb7e_1
  - r-quantreg=6.1=r45h11cdb10_1
  - r-r.methodss3=1.8.2=r45hc72bb7e_4
  - r-r.oo=1.27.1=r45hc72bb7e_1
  - r-r.utils=2.13.0=r45hc72bb7e_1
  - r-r6=2.6.1=r45hc72bb7e_1
  - r-ragg=1.5.0=r45h9f1dc4d_1
  - r-rann=2.6.2=r45h3697838_2
  - r-rappdirs=0.3.4=r45h54b55ab_0
  - r-rbibutils=2.4.1=r45h54b55ab_0
  - r-rcmdcheck=1.4.0=r45h785f33e_4
  - r-rcolorbrewer=1.1_3=r45h785f33e_4
  - r-rcpp=1.1.1=r45h3697838_0
  - r-rcppannoy=0.0.23=r45h3697838_0
  - r-rcpparmadillo=15.2.3_1=r45h3704496_0
  - r-rcppeigen=0.3.4.0.2=r45h3704496_1
  - r-rcpphnsw=0.6.0=r45h3697838_2
  - r-rcppml=0.3.7=r45h3697838_4
  - r-rcppparallel=5.1.11_2=r45h0d96847_0
  - r-rcppprogress=0.4.2=r45hc72bb7e_5
  - r-rcppspdlog=0.0.27=r45h3697838_0
  - r-rcpptoml=0.2.3=r45h3697838_1
  - r-rcurl=1.98_1.17=r45h46721d4_2
  - r-rdpack=2.6.6=r45hc72bb7e_0
  - r-readr=2.2.0=r45h3697838_0
  - r-readxl=1.4.5=r45h10e25cc_1
  - r-recipes=1.3.1=r45hc72bb7e_1
  - r-reformulas=0.4.4=r45hc72bb7e_0
  - r-rematch=2.0.0=r45hc72bb7e_2
  - r-rematch2=2.1.2=r45hc72bb7e_5
  - r-remotes=2.5.0=r45hc72bb7e_2
  - r-reprex=2.1.1=r45hc72bb7e_2
  - r-reshape2=1.4.5=r45h3697838_0
  - r-restfulr=0.0.16=r45hf7ecca6_1
  - r-reticulate=1.45.0=r45h3697838_0
  - r-rglpk=0.6_5.1=r45hcbc0784_2
  - r-rjson=0.2.23=r45h3697838_1
  - r-rlang=1.1.7=r45h3697838_0
  - r-rmarkdown=2.30=r45hc72bb7e_0
  - r-rmtstat=0.3.1=r45ha770c72_4
  - r-rocr=1.0_12=r45hc72bb7e_0
  - r-rook=1.2=r45h54b55ab_3
  - r-roxygen2=7.3.3=r45h3697838_1
  - r-rpart=4.1.24=r45h54b55ab_1
  - r-rprojroot=2.1.1=r45hc72bb7e_1
  - r-rspectra=0.16_2=r45h3704496_2
  - r-rsqlite=3.53.1=r45h3697838_0
  - r-rstatix=0.7.3=r45hc72bb7e_0
  - r-rstudioapi=0.18.0=r45hc72bb7e_0
  - r-rsvd=1.0.5=r45hc72bb7e_4
  - r-rtsne=0.17=r45hf1899b2_3
  - r-rversions=3.0.0=r45hc72bb7e_0
  - r-rvest=1.0.5=r45hc72bb7e_1
  - r-s7=0.2.1=r45h54b55ab_0
  - r-sass=0.4.10=r45h3697838_1
  - r-scales=1.4.0=r45hc72bb7e_1
  - r-scattermore=1.2=r45h3697838_5
  - r-scatterpie=0.2.6=r45hc72bb7e_1
  - r-sccore=1.0.6=r45h3697838_2
  - r-sctransform=0.4.3=r45h3704496_0
  - r-selectr=0.5_1=r45hc72bb7e_0
  - r-sessioninfo=1.2.3=r45hc72bb7e_1
  - r-seurat=5.4.0=r45h3697838_0
  - r-seuratobject=5.3.0=r45h3697838_0
  - r-shape=1.4.6.1=r45ha770c72_2
  - r-shiny=1.13.0=r45h785f33e_0
  - r-sitmo=2.0.2=r45h3697838_4
  - r-slam=0.1_55=r45h7b2432b_1
  - r-snow=0.4_4=r45hc72bb7e_4
  - r-sourcetools=0.1.7_1=r45h3697838_3
  - r-sp=2.2_1=r45h54b55ab_0
  - r-spam=2.11_3=r45h2ddecb4_0
  - r-sparsem=1.84_2=r45heaba542_1
  - r-sparsevctrs=0.3.6=r45h54b55ab_0
  - r-spatstat.data=3.1_9=r45hc72bb7e_0
  - r-spatstat.explore=3.7_0=r45h54b55ab_0
  - r-spatstat.geom=3.7_0=r45h54b55ab_0
  - r-spatstat.random=3.4_4=r45h3697838_0
  - r-spatstat.sparse=3.1_0=r45h54b55ab_2
  - r-spatstat.univar=3.1_6=r45h54b55ab_0
  - r-spatstat.utils=3.2_1=r45h54b55ab_0
  - r-squarem=2021.1=r45hc72bb7e_4
  - r-statmod=1.5.1=r45hb1d0f04_0
  - r-stringfish=0.18.0=r45h2a3d9df_0
  - r-stringi=1.8.7=r45h3d52c89_2
  - r-stringr=1.6.0=r45h785f33e_0
  - r-survival=3.8_6=r45h54b55ab_0
  - r-sys=3.4.3=r45h54b55ab_1
  - r-systemfonts=1.3.1=r45h74f4acd_0
  - r-tensor=1.5.1=r45hc72bb7e_1
  - r-testthat=3.3.2=r45h3697838_0
  - r-textshaping=1.0.4=r45h74f4acd_0
  - r-tibble=3.3.1=r45h54b55ab_0
  - r-tidydr=0.0.6=r45hc72bb7e_1
  - r-tidygraph=1.3.0=r45h3697838_2
  - r-tidyr=1.3.2=r45h3697838_0
  - r-tidyselect=1.2.1=r45hc72bb7e_2
  - r-tidytree=0.4.8=r45hc72bb7e_0
  - r-tidyverse=2.0.0=r45h785f33e_3
  - r-timechange=0.4.0=r45h3697838_0
  - r-timedate=4052.112=r45hc72bb7e_0
  - r-tinytex=0.58=r45hc72bb7e_0
  - r-triebeard=0.4.1=r45h3697838_4
  - r-tseries=0.10_60=r45ha0a88a1_0
  - r-ttr=0.24.4=r45h54b55ab_2
  - r-tweenr=2.0.3=r45h3697838_2
  - r-tzdb=0.5.0=r45h3697838_2
  - r-upsetr=1.4.1=r45hc72bb7e_0
  - r-urca=1.3_4=r45heaba542_2
  - r-urlchecker=1.0.1=r45hc72bb7e_4
  - r-urltools=1.7.3.1=r45h3697838_1
  - r-usethis=3.2.1=r45hc72bb7e_1
  - r-utf8=1.2.6=r45h54b55ab_1
  - r-uuid=1.2_2=r45h54b55ab_0
  - r-uwot=0.2.4=r45h3697838_0
  - r-vctrs=0.7.1=r45h3697838_0
  - r-vipor=0.4.7=r45hc72bb7e_2
  - r-viridis=0.6.5=r45hc72bb7e_2
  - r-viridislite=0.4.3=r45hc72bb7e_0
  - r-vroom=1.7.0=r45h3697838_0
  - r-waldo=0.6.2=r45hc72bb7e_1
  - r-whisker=0.4.1=r45hc72bb7e_3
  - r-withr=3.0.2=r45hc72bb7e_1
  - r-workflowr=1.7.2=r45hc72bb7e_1
  - r-xfun=0.56=r45h3697838_0
  - r-xgboost=3.2.0=cpu_r45h2ebb00f_1
  - r-xml=3.99_0.22=r45hf705907_0
  - r-xml2=1.5.2=r45he78afff_0
  - r-xopen=1.0.1=r45hc72bb7e_2
  - r-xtable=1.8_8=r45hc72bb7e_0
  - r-xts=0.14.2=r45h54b55ab_0
  - r-yaml=2.3.12=r45h54b55ab_0
  - r-yulab.utils=0.2.4=r45hc72bb7e_0
  - r-zip=2.3.3=r45h54b55ab_1
  - r-zoo=1.8_15=r45h54b55ab_0
  - readline=8.3=h853b02a_0
  - rhash=1.4.6=hb9d3cd8_1
  - sed=4.9=h6688a6e_0
  - sysroot_linux-64=2.34=h087de78_3
  - tbb=2022.3.0=hb700be7_2
  - tbb-devel=2022.3.0=h51de99f_2
  - tk=8.6.13=noxft_h366c992_103
  - tktable=2.10=h8d826fa_7
  - tomlkit=0.15.0=pyha770c72_0
  - tzdata=2025c=hc9c84f9_1
  - xmltodict=1.0.4=pyhcf101f3_0
  - xorg-libice=1.1.2=hb9d3cd8_0
  - xorg-libsm=1.2.6=he73a12e_0
  - xorg-libx11=1.8.13=he1eb515_0
  - xorg-libxau=1.0.12=hb03c661_1
  - xorg-libxdmcp=1.1.5=hb03c661_1
  - xorg-libxext=1.3.7=hb03c661_0
  - xorg-libxrender=0.9.12=hb9d3cd8_0
  - xorg-libxt=1.3.1=hb9d3cd8_0
  - xorg-xextproto=7.3.0=hb9d3cd8_1004
  - yaml=0.2.5=h280c20c_3
  - yq=4.1.1=pyhcf101f3_0
  - zstd=1.5.7=hb78ec9c_6
```

</details>

<details>
<summary>velocyto.yml (366 lines) - click to expand</summary>

Environment `velocyto` - covers `velocyto_analysis.Rmd` only. Kept separate deliberately: `velocyto.R` has its own old, narrow set of compatible dependency versions that don't mix cleanly with the environment above.

```yaml
name: velocyto
channels:
  - bioconda
  - conda-forge
dependencies:
  - _openmp_mutex=4.5=20_gnu
  - _r-mutex=1.0.1=anacondar_1
  - binutils_impl_linux-64=2.45.1=default_hfdba357_102
  - bioconductor-biobase=2.66.0=r44h3df3fcb_0
  - bioconductor-biocgenerics=0.52.0=r44hdfd78af_3
  - bioconductor-pcamethods=1.98.0=r44he5774e6_1
  - boost-cpp=1.85.0=h3c6214e_4
  - bwidget=1.10.1=ha770c72_1
  - bzip2=1.0.8=hda65f42_9
  - c-ares=1.34.6=hb03c661_0
  - ca-certificates=2026.5.20=hbd8a1cb_0
  - cairo=1.18.4=h3394656_0
  - curl=8.20.0=hcf29cc6_0
  - font-ttf-dejavu-sans-mono=2.37=hab24e00_0
  - font-ttf-inconsolata=3.000=h77eed37_0
  - font-ttf-source-code-pro=2.038=h77eed37_0
  - font-ttf-ubuntu=0.83=h77eed37_3
  - fontconfig=2.18.1=h27c8c51_0
  - fonts-conda-ecosystem=1=0
  - fonts-conda-forge=1=hc364b38_1
  - freetype=2.14.3=ha770c72_0
  - fribidi=1.0.16=hb03c661_0
  - gcc_impl_linux-64=15.2.0=he0086c7_19
  - gfortran_impl_linux-64=15.2.0=h281d09f_19
  - glpk=5.0=h445213a_0
  - gmp=6.3.0=hac33072_2
  - graphite2=1.3.15=hecca717_0
  - gsl=2.7=he838d99_0
  - gxx_impl_linux-64=15.2.0=hda75c37_19
  - harfbuzz=12.2.0=h15599e2_0
  - hdf5=1.14.6=nompi_h19486de_109
  - icu=75.1=he02047a_0
  - kernel-headers_linux-64=5.14.0=he073ed8_3
  - keyutils=1.6.3=hb9d3cd8_0
  - krb5=1.22.2=ha1258a1_0
  - ld_impl_linux-64=2.45.1=default_hbd61a6d_102
  - lerc=4.1.0=hdb68285_0
  - libaec=1.1.5=h088129d_0
  - libblas=3.11.0=8_h4a7cf45_openblas
  - libboost=1.85.0=h0ccab89_4
  - libboost-devel=1.85.0=h00ab1b0_4
  - libboost-headers=1.85.0=ha770c72_4
  - libcblas=3.11.0=8_h0358290_openblas
  - libcurl=8.20.0=hcf29cc6_0
  - libdeflate=1.25=h17f619e_0
  - libedit=3.1.20250104=pl5321h7949ede_0
  - libev=4.33=hd590300_2
  - libexpat=2.8.1=hecca717_0
  - libffi=3.5.2=h3435931_0
  - libfreetype=2.14.3=ha770c72_0
  - libfreetype6=2.14.3=h73754d4_0
  - libgcc=15.2.0=he0feb66_19
  - libgcc-devel_linux-64=15.2.0=hcc6f6b0_119
  - libgcc-ng=15.2.0=h69a702a_19
  - libgfortran=15.2.0=h69a702a_19
  - libgfortran-ng=15.2.0=h69a702a_19
  - libgfortran5=15.2.0=h68bc16d_19
  - libgit2=1.9.4=hc20babb_0
  - libglib=2.88.1=h0d30a3d_2
  - libgomp=15.2.0=he0feb66_19
  - libhwloc=2.12.2=default_hafda6a7_1000
  - libiconv=1.18=h3b78370_2
  - libjpeg-turbo=3.1.4.1=hb03c661_0
  - liblapack=3.11.0=8_h47877c9_openblas
  - liblzma=5.8.3=hb03c661_0
  - liblzma-devel=5.8.3=hb03c661_0
  - libnghttp2=1.68.1=h877daf1_0
  - libopenblas=0.3.33=pthreads_h94d23a6_0
  - libpng=1.6.58=h421ea60_0
  - libsanitizer=15.2.0=h90f66d4_19
  - libssh2=1.11.1=hcf80075_0
  - libstdcxx=15.2.0=h934c35e_19
  - libstdcxx-devel_linux-64=15.2.0=hd446a21_119
  - libstdcxx-ng=15.2.0=hdf11a46_19
  - libtiff=4.7.1=h9d88235_1
  - libuuid=2.42.1=h5347b49_0
  - libuv=1.52.1=h280c20c_0
  - libwebp-base=1.6.0=hd42ef1d_0
  - libxcb=1.17.0=h8a09558_0
  - libxml2=2.15.1=h26afc86_0
  - libxml2-16=2.15.1=ha9997c6_0
  - libzlib=1.3.2=h25fd6f3_2
  - make=4.4.1=hb9d3cd8_2
  - ncurses=6.6=hdb14827_0
  - openssl=3.6.2=h35e630c_0
  - pandoc=3.9.0.2=ha770c72_0
  - pango=1.56.4=hadf4263_0
  - pcre2=10.47=haa7fec5_0
  - pixman=0.46.4=h54a6638_1
  - pthread-stubs=0.4=hb9d3cd8_1002
  - r-abind=1.4_8=r44hc72bb7e_1
  - r-askpass=1.2.1=r44h54b55ab_1
  - r-assertthat=0.2.1=r44hc72bb7e_6
  - r-backports=1.5.1=r44h54b55ab_0
  - r-base=4.4.3=h835929b_7
  - r-base64enc=0.1_6=r44h54b55ab_0
  - r-bh=1.87.0_1=r44hc72bb7e_1
  - r-bit=4.6.0=r44h54b55ab_1
  - r-bit64=4.8.2=r44h54b55ab_0
  - r-bitops=1.0_9=r44h54b55ab_1
  - r-blob=1.3.0=r44hc72bb7e_0
  - r-brew=1.0_10=r44hc72bb7e_2
  - r-brio=1.1.5=r44h54b55ab_2
  - r-broom=1.0.13=r44hc72bb7e_0
  - r-bslib=0.11.0=r44hc72bb7e_0
  - r-cachem=1.1.0=r44h54b55ab_2
  - r-callr=3.7.6=r44hc72bb7e_2
  - r-catools=1.18.3=r44h3697838_1
  - r-cellranger=1.1.0=r44hc72bb7e_1008
  - r-cli=3.6.6=r44h3697838_0
  - r-cliapp=0.1.2=r44hc72bb7e_2
  - r-clipr=0.8.1=r44hc72bb7e_0
  - r-cluster=2.1.8.2=r44heaba542_0
  - r-codetools=0.2_20=r44hc72bb7e_2
  - r-colorspace=2.1_2=r44h54b55ab_0
  - r-commonmark=2.0.0=r44h54b55ab_1
  - r-conflicted=1.2.0=r44h785f33e_3
  - r-cowplot=1.2.0=r44hc72bb7e_2
  - r-cpp11=0.5.5=r44h785f33e_0
  - r-crayon=1.5.3=r44hc72bb7e_2
  - r-credentials=2.0.3=r44hc72bb7e_1
  - r-crosstalk=1.2.2=r44hc72bb7e_1
  - r-curl=7.1.0=r44h10955f1_0
  - r-data.table=1.17.8=r44h1c8cec4_1
  - r-dbi=1.3.0=r44hc72bb7e_0
  - r-dbplyr=2.5.2=r44hc72bb7e_0
  - r-deldir=2.0_4=r44heaba542_2
  - r-dendsort=0.3.4=r44ha770c72_4
  - r-desc=1.4.3=r44hc72bb7e_2
  - r-devtools=2.5.2=r44hc72bb7e_0
  - r-diffobj=0.3.6=r44h54b55ab_1
  - r-digest=0.6.39=r44h3697838_0
  - r-dotcall64=1.2=r44heaba542_1
  - r-downlit=0.4.5=r44hc72bb7e_0
  - r-dplyr=1.2.1=r44h3697838_0
  - r-dqrng=0.3.2=r44h3697838_2
  - r-drat=0.2.5=r44hc72bb7e_1
  - r-dtplyr=1.3.3=r44hc72bb7e_0
  - r-ellipsis=0.3.3=r44h54b55ab_0
  - r-evaluate=1.0.5=r44hc72bb7e_1
  - r-fansi=1.0.7=r44h54b55ab_0
  - r-farver=2.1.2=r44h3697838_2
  - r-fastcluster=1.3.0=r44h3697838_1
  - r-fastdummies=1.7.6=r44hc72bb7e_0
  - r-fastmap=1.2.0=r44h3697838_2
  - r-filelock=1.0.3=r44h54b55ab_2
  - r-fitdistrplus=1.2_6=r44hc72bb7e_0
  - r-fnn=1.1.4.1=r44h3697838_2
  - r-fontawesome=0.5.3=r44hc72bb7e_1
  - r-forcats=1.0.1=r44hc72bb7e_0
  - r-fs=2.1.0=r44h3697838_0
  - r-future=1.70.0=r44hc72bb7e_0
  - r-future.apply=1.20.2=r44hc72bb7e_0
  - r-gargle=1.6.1=r44h785f33e_0
  - r-generics=0.1.4=r44hc72bb7e_1
  - r-gert=2.3.1=r44h5e22a44_0
  - r-ggplot2=4.0.3=r44h785f33e_0
  - r-ggrepel=0.9.8=r44h3697838_0
  - r-ggridges=0.5.7=r44hc72bb7e_1
  - r-gh=1.6.0=r44hc72bb7e_0
  - r-gitcreds=0.1.2=r44hc72bb7e_4
  - r-globals=0.19.1=r44hc72bb7e_0
  - r-glue=1.8.1=r44h54b55ab_0
  - r-goftest=1.2_3=r44h54b55ab_4
  - r-googledrive=2.1.2=r44hc72bb7e_1
  - r-googlesheets4=1.1.2=r44h785f33e_1
  - r-gplots=3.3.0=r44hc72bb7e_0
  - r-gridextra=2.3=r44hc72bb7e_1007
  - r-gtable=0.3.6=r44hc72bb7e_1
  - r-gtools=3.9.5=r44h54b55ab_2
  - r-haven=2.5.5=r44h6d565e7_1
  - r-hdf5r=1.3.12=r44h39a46f8_2
  - r-here=1.0.2=r44hc72bb7e_0
  - r-hexbin=1.28.5=r44heaba542_1
  - r-highr=0.12=r44hc72bb7e_0
  - r-hms=1.1.4=r44hc72bb7e_0
  - r-htmltools=0.5.9=r44h3697838_0
  - r-htmlwidgets=1.6.4=r44h785f33e_4
  - r-httpuv=1.6.17=r44h6d565e7_0
  - r-httr=1.4.8=r44hc72bb7e_0
  - r-httr2=1.2.2=r44hc72bb7e_0
  - r-ica=1.0_3=r44hc72bb7e_4
  - r-ids=1.0.1=r44hc72bb7e_5
  - r-igraph=2.3.2=r44hf411e2a_0
  - r-ini=0.3.1=r44hc72bb7e_1007
  - r-irlba=2.3.7=r44h0e4624f_0
  - r-isoband=0.3.0=r44h3697838_0
  - r-jquerylib=0.1.4=r44hc72bb7e_4
  - r-jsonlite=2.0.0=r44h54b55ab_1
  - r-kernsmooth=2.23_26=r44ha0a88a1_1
  - r-knitr=1.51=r44hc72bb7e_0
  - r-labeling=0.4.3=r44hc72bb7e_2
  - r-later=1.4.8=r44h3697838_0
  - r-lattice=0.22_9=r44h54b55ab_0
  - r-lazyeval=0.2.3=r44h54b55ab_0
  - r-leidenbase=0.1.36=r44ha11a66c_0
  - r-lifecycle=1.0.5=r44hc72bb7e_0
  - r-listenv=0.10.1=r44hc72bb7e_0
  - r-lmtest=0.9_40=r44heaba542_4
  - r-lpsolve=5.6.23=r44h54b55ab_1
  - r-lubridate=1.9.5=r44h54b55ab_0
  - r-magrittr=2.0.5=r44h54b55ab_0
  - r-mass=7.3_65=r44h54b55ab_0
  - r-matrix=1.7_5=r44h0e4624f_0
  - r-matrixstats=1.5.0=r44h54b55ab_1
  - r-memoise=2.0.1=r44hc72bb7e_4
  - r-mgcv=1.9_4=r44h0e4624f_0
  - r-mime=0.13=r44h54b55ab_1
  - r-miniui=0.1.2=r44hc72bb7e_1
  - r-modelr=0.1.11=r44hc72bb7e_3
  - r-munsell=0.5.1=r44hc72bb7e_2
  - r-n2r=1.0.5=r44h3697838_0
  - r-nlme=3.1_169=r44heaba542_0
  - r-openssl=2.4.1=r44h68c19f5_0
  - r-otel=0.2.0=r44hc72bb7e_1
  - r-pagoda2=1.0.15=r44h3704496_0
  - r-pak=0.9.5=r44hc72bb7e_0
  - r-parallelly=1.47.0=r44h54b55ab_0
  - r-patchwork=1.3.2=r44hc72bb7e_1
  - r-pbapply=1.7_4=r44hc72bb7e_1
  - r-pbmcapply=1.5.1=r44h54b55ab_5
  - r-pillar=1.11.1=r44hc72bb7e_0
  - r-pkgbuild=1.4.8=r44hc72bb7e_1
  - r-pkgcache=2.2.4=r44hc72bb7e_1
  - r-pkgconfig=2.0.3=r44hc72bb7e_5
  - r-pkgdown=2.2.0=r44hc72bb7e_0
  - r-pkgload=1.5.2=r44hc72bb7e_0
  - r-plotly=4.12.0=r44hc72bb7e_0
  - r-plyr=1.8.9=r44h3697838_3
  - r-png=0.1_9=r44haf2892b_0
  - r-polyclip=1.10_7=r44h3697838_1
  - r-praise=1.0.0=r44hc72bb7e_1009
  - r-prettycode=1.1.0=r44hc72bb7e_5
  - r-prettyunits=1.2.0=r44hc72bb7e_2
  - r-proc=1.19.0.1=r44h3697838_1
  - r-processx=3.9.0=r44h54b55ab_0
  - r-profvis=0.4.0=r44h54b55ab_1
  - r-progress=1.2.3=r44hc72bb7e_2
  - r-progressr=0.19.0=r44hc72bb7e_0
  - r-promises=1.5.0=r44hc72bb7e_1
  - r-ps=1.9.3=r44h54b55ab_0
  - r-purrr=1.2.2=r44h54b55ab_0
  - r-qs2=0.2.1=r44h1a6caf7_0
  - r-r.methodss3=1.8.2=r44hc72bb7e_4
  - r-r.oo=1.27.1=r44hc72bb7e_1
  - r-r.utils=2.13.0=r44hc72bb7e_1
  - r-r6=2.6.1=r44hc72bb7e_1
  - r-ragg=1.5.2=r44h9f1dc4d_0
  - r-rann=2.6.2=r44h3697838_2
  - r-rappdirs=0.3.4=r44h54b55ab_0
  - r-rcmdcheck=1.4.0=r44h785f33e_4
  - r-rcolorbrewer=1.1_3=r44h785f33e_4
  - r-rcpp=1.1.1_1.1=r44h3697838_0
  - r-rcppannoy=0.0.23=r44h3697838_0
  - r-rcpparmadillo=15.2.7_1=r44h3704496_0
  - r-rcppeigen=0.3.4.0.2=r44h3704496_1
  - r-rcpphnsw=0.7.0=r44h3697838_0
  - r-rcppparallel=5.1.11_2=r44h0d96847_0
  - r-rcppprogress=0.4.2=r44hc72bb7e_5
  - r-rcppspdlog=0.0.28=r44h3697838_0
  - r-rcpptoml=0.2.3=r44h3697838_1
  - r-readr=2.2.0=r44h3697838_0
  - r-readxl=1.5.0=r44h10e25cc_0
  - r-rematch=2.0.0=r44hc72bb7e_2
  - r-rematch2=2.1.2=r44hc72bb7e_5
  - r-remotes=2.5.0=r44hc72bb7e_2
  - r-reprex=2.1.1=r44hc72bb7e_2
  - r-reshape2=1.4.5=r44h3697838_0
  - r-reticulate=1.46.0=r44h3697838_0
  - r-rjson=0.2.23=r44h3697838_1
  - r-rlang=1.2.0=r44h3697838_0
  - r-rmarkdown=2.31=r44hc72bb7e_0
  - r-rmtstat=0.3.1=r44ha770c72_4
  - r-rocr=1.0_12=r44hc72bb7e_0
  - r-rook=1.2=r44h54b55ab_3
  - r-roxygen2=8.0.0=r44h3697838_0
  - r-rprojroot=2.1.1=r44hc72bb7e_1
  - r-rspectra=0.16_2=r44h3704496_2
  - r-rstudioapi=0.18.0=r44hc72bb7e_0
  - r-rtsne=0.17=r44hf1899b2_3
  - r-rversions=3.0.0=r44hc72bb7e_0
  - r-rvest=1.0.5=r44hc72bb7e_1
  - r-s7=0.2.2=r44h54b55ab_0
  - r-sass=0.4.10=r44h3697838_1
  - r-scales=1.4.0=r44hc72bb7e_1
  - r-scattermore=1.2=r44h3697838_5
  - r-sccore=1.0.7=r44h3697838_0
  - r-sctransform=0.4.3=r44h3704496_0
  - r-selectr=0.5_1=r44hc72bb7e_0
  - r-sessioninfo=1.2.3=r44hc72bb7e_1
  - r-seurat=5.5.0=r44h3697838_0
  - r-seuratobject=5.4.0=r44h3697838_0
  - r-shiny=1.13.0=r44h785f33e_0
  - r-sitmo=2.0.2=r44h3697838_4
  - r-sourcetools=0.1.7_2=r44h3697838_0
  - r-sp=2.2_1=r44h54b55ab_0
  - r-spam=2.11_4=r44h2ddecb4_0
  - r-spatstat.data=3.1_9=r44hc72bb7e_0
  - r-spatstat.explore=3.8_0=r44h54b55ab_0
  - r-spatstat.geom=3.8_1=r44h54b55ab_0
  - r-spatstat.random=3.4_5=r44h3697838_0
  - r-spatstat.sparse=3.2_0=r44h54b55ab_0
  - r-spatstat.univar=3.2_0=r44h54b55ab_0
  - r-spatstat.utils=3.2_3=r44h54b55ab_0
  - r-stringfish=0.19.0=r44h2a3d9df_0
  - r-stringi=1.8.7=r44h2dae267_1
  - r-stringr=1.6.0=r44h785f33e_0
  - r-survival=3.8_6=r44h54b55ab_0
  - r-sys=3.4.3=r44h54b55ab_1
  - r-systemfonts=1.3.2=r44h74f4acd_0
  - r-tensor=1.5.1=r44hc72bb7e_1
  - r-testthat=3.3.2=r44h3697838_0
  - r-textshaping=1.0.4=r44h74f4acd_0
  - r-tibble=3.3.1=r44h54b55ab_0
  - r-tidyr=1.3.2=r44h3697838_0
  - r-tidyselect=1.2.1=r44hc72bb7e_2
  - r-tidyverse=2.0.0=r44h785f33e_3
  - r-timechange=0.4.0=r44h3697838_0
  - r-tinytex=0.59=r44hc72bb7e_0
  - r-triebeard=0.4.1=r44h3697838_4
  - r-tzdb=0.5.0=r44h3697838_2
  - r-urlchecker=1.0.1=r44hc72bb7e_4
  - r-urltools=1.7.3.1=r44h3697838_1
  - r-usethis=3.2.1=r44hc72bb7e_1
  - r-utf8=1.2.6=r44h54b55ab_1
  - r-uuid=1.2_2=r44h54b55ab_0
  - r-uwot=0.2.4=r44h3697838_0
  - r-vctrs=0.7.3=r44h3697838_0
  - r-velocyto.r=0.6=r44h503566f_9
  - r-viridislite=0.4.3=r44hc72bb7e_0
  - r-vroom=1.7.1=r44h3697838_0
  - r-waldo=0.6.2=r44hc72bb7e_1
  - r-whisker=0.4.1=r44hc72bb7e_3
  - r-withr=3.0.2=r44hc72bb7e_1
  - r-xfun=0.57=r44h3697838_0
  - r-xml2=1.5.2=r44he78afff_0
  - r-xopen=1.0.1=r44hc72bb7e_2
  - r-xtable=1.8_8=r44hc72bb7e_0
  - r-yaml=2.3.12=r44h54b55ab_0
  - r-zip=2.3.3=r44h54b55ab_1
  - r-zoo=1.8_15=r44h54b55ab_0
  - readline=8.3=h853b02a_0
  - sed=4.10=h19d0853_0
  - sysroot_linux-64=2.34=h087de78_3
  - tbb=2022.3.0=hb700be7_2
  - tbb-devel=2022.3.0=h51de99f_2
  - tk=8.6.13=noxft_h366c992_103
  - tktable=2.10=h5a7a40f_8
  - tzdata=2025c=hc9c84f9_1
  - xorg-libice=1.1.2=hb9d3cd8_0
  - xorg-libsm=1.2.6=he73a12e_0
  - xorg-libx11=1.8.13=he1eb515_0
  - xorg-libxau=1.0.12=hb03c661_1
  - xorg-libxdmcp=1.1.5=hb03c661_1
  - xorg-libxext=1.3.7=hb03c661_0
  - xorg-libxrender=0.9.12=hb9d3cd8_0
  - xorg-libxt=1.3.1=hb9d3cd8_0
  - xz=5.8.3=ha02ee65_0
  - xz-gpl-tools=5.8.3=ha02ee65_0
  - xz-tools=5.8.3=hb03c661_0
  - zstd=1.5.7=hb78ec9c_6
```

</details>

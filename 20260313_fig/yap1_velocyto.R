#!/usr/bin/Rscript
# velocyto redo on L03 yap1 scrnaseq data (with intron retention on cellranger)
# 3 high-resolution plots provided separately for custom composing via Adobe Illustrator
# For each of these plots, legends are provided separately for the same reason
# axis also provided
libraries <- c(
    "tidyverse",
    "velocyto.R",
    "pagoda2",
    "qs2",
    "Seurat",
    "SeuratWrappers",
    "ggplot2"
)

load_packages <- function(packages) {
    suppressPackageStartupMessages({
        cat(paste("Loading", packages, sep=" "), sep="\n")
        invisible(sapply(
            packages, library, character.only = TRUE, quietly = TRUE
            ))
    })
}
invisible(Sys.setenv(OPENBLAS_NUM_THREADS="1"))
load_packages(libraries)

detect_sample_prefixes <- function(loom_cellnames, n_samples) {
    prefixes <- unique(sub(":.*", "", loom_cellnames))
    if (length(prefixes) != n_samples) {
        warning(sprintf("Expected %d samples but found %d prefixes: %s", 
                        n_samples, length(prefixes), paste(prefixes, collapse=", ")))
    }
    return(sort(prefixes))  # sort so sample 1 < sample 2
}

run_velocyto <- function(
    data_sub, loom_cellnames, ldat, outfile_path, outfile_legend, outfile_obj, label,
    recode_ids=NULL,      # named vector: c("old"="new")
    col_map=NULL,         # named colour vector
    level_order=NULL      # character vector of level order
    ) {                   # -- subset loom assays to matched cells --
    # -- optional recoding --
    if (!is.null(recode_ids)) {
        data_sub$L3_cluster_id <- recode(data_sub$L3_cluster_id, !!!recode_ids)
    }

    # -- optional reordering --
    if (!is.null(level_order)) {
        data_sub$L3_cluster_id <- factor(data_sub$L3_cluster_id, levels=level_order)
    }
    
    umap_cols <- if (!is.null(col_map)) col_map else NULL

    list_all_assays_filtered <- lapply(names(x=ldat), function(name){
        assay <- ldat[[name]]
        colnames(assay) <- loom_cellnames
        assay <- assay[, intersect(colnames(data_sub), colnames(assay))]
        # assay <- assay[, colnames(data_sub)]
    })
    names(list_all_assays_filtered) <- names(x=ldat)

    emat <- list_all_assays_filtered$spliced
    nmat <- list_all_assays_filtered$unspliced

    emb <- Embeddings(data_sub, reduction="umap")
    cell.dist <- as.dist(1 - armaCor(t(emb)))

    gg <- UMAPPlot(
        data_sub, 
        group.by="L3_cluster_id", 
        cols=if (!is.null(umap_cols)) umap_cols else waiver()) +
        guides(colour=guide_legend(override.aes=list(shape=16, size=4))) +
        theme_void() +
        theme(
            legend.position=c(1, 1),
            legend.justification=c(1, 1),
            legend.title=element_blank(),
            legend.background=element_blank(),
            legend.key=element_blank()
        )

    colors <- as.list(ggplot_build(gg)$data[[1]]$colour)
    names(colors) <- rownames(emb)

    if (is.null(col_map)) {
        legend_data <- data.frame(
            x=1, y=1,
            L3_cluster_id=factor(unique(data_sub$L3_cluster_id))
        )
    } else {
        legend_data <- data.frame(
            x=1, y=1,
            L3_cluster_id=factor(names(col_map), levels=names(col_map))
        )
    }

    gg_legend <- ggplot(legend_data, aes(x=x, y=y, colour=L3_cluster_id)) +
        geom_point(shape=NA) +
        scale_colour_manual(values=col_map) +
        guides(colour=guide_legend(override.aes=list(shape=16, size=4))) +
        theme_void() +
        theme(
            legend.title=element_blank(),
            legend.background=element_blank(),
            legend.margin=margin(r=200)
        )
    ggsave(outfile_legend, gg_legend)

    cell_identity <- data_sub$L3_cluster_id

    emat_filtered <- filter.genes.by.cluster.expression(
        emat, cell_identity, min.max.cluster.average=0.1
    )
    nmat_filtered <- filter.genes.by.cluster.expression(
        nmat, cell_identity, min.max.cluster.average=0.05
    )
    cat(sprintf("emat_filtered: %d genes x %d cells\n", nrow(emat_filtered), ncol(emat_filtered)))
    cat(sprintf("nmat_filtered: %d genes x %d cells\n", nrow(nmat_filtered), ncol(nmat_filtered)))

    rvel <- gene.relative.velocity.estimates(
        emat_filtered,
        nmat_filtered,
        deltaT=deltaT,
        kCells=kCells,
        cell.dist=cell.dist,
        fit.quantile=fit.quantile,
        n.cores=n.cores
    )

    n <- 400
    kCells <- 10
    deltaT <- 3
    fit.quantile <- 0.02

    pdf(outfile_path)
    par(bty="n")
    show.velocity.on.embedding.cor(
        emb, rvel, n=n, scale='sqrt',
        cell.colors=ac(colors, alpha=alpha),
        cex=cex, arrow.scale=arrow.scale,
        show.grid.flow=show.grid.flow,
        min.grid.cell.mass=min.grid.cell.mass,
        grid.n=grid.n, arrow.lwd=arrow.lwd,
        do.par=F, cell.border.alpha=cell.border.alpha,
        n.cores=n.cores,
        main=sprintf(
            "Cell Velocity: %s\nfitquantile=%.2f, deltaT=%.1f, kCells=%i, n=%i",
            label, fit.quantile, deltaT, kCells, n
        )
    )
    dev.off()

    qs_save(rvel, outfile_obj)
    cat("Saved: ", outfile_obj, "\n")
}

# plot parameters
alpha <- 0.5
cex <<- 0.7
arrow.scale <<- 2
show.grid.flow <<- TRUE
min.grid.cell.mass <<- 1
grid.n <<- 50
arrow.lwd <<- 1
cell.border.alpha <<- 0
n.cores <<- 1

n <- 400
kCells <<- 10
deltaT <<- 3
fit.quantile <<- 0.02


infile_path <- "../../Saki_data/paper/D10025_yap1_dataset_Level_03_annotated.qs2"
infile_loom <- "Dataset_10025_Yap_scRNAseq.20260311.intron.loom"
outfile_dir <- "../Revision_analysis/queue_plots/"

# shared drive: /Revision_analysis/queue_plots/
#   Yap combined: yap_scrna_Level_03_L3_velocyto_introns.D10025.pdf
#   Yap combined legend: yap_scrna_Level_03_L3_velocyto_introns.D10025.legend.pdf
#   Yap (Mut sample): yap_scrna_Level_03_L3_velocyto_introns.S20127.pdf
#   Yap (Mut sample) legend: yap_scrna_Level_03_L3_velocyto_introns.S20127.legend.pdf
#   Yap (WT sample): yap_scrna_Level_03_L3_velocyto_introns.S20128.pdf
#   Yap (WT sample) legend: yap_scrna_Level_03_L3_velocyto_introns.S20128.legend.pdf
#   mini axis for umaps if needed: mini_axis.pdf

outfile_path_D10025 <- paste(
    outfile_dir, "yap_scrna_Level_03_L3_velocyto_introns.D10025.pdf", sep="/"
    )
outfile_legend_D10025 <- paste(
    outfile_dir, "yap_scrna_Level_03_L3_velocyto_introns.D10025.legend.pdf", sep="/"
    )
outfile_obj_D10025 <- paste(
    outfile_dir, "velocyto.D10025.qs2", sep="/"
    )

outfile_path_S20127 <- paste(
    outfile_dir, "yap_scrna_Level_03_L3_velocyto_introns.S20127.pdf", sep="/"
    )
outfile_legend_S20127 <- paste(
    outfile_dir, "yap_scrna_Level_03_L3_velocyto_introns.S20127.legend.pdf", sep="/"
    )
outfile_obj_S20127 <- paste(
    outfile_dir, "velocyto.S20127.qs2", sep="/"
    )

outfile_path_S20128 <- paste(
    outfile_dir, "yap_scrna_Level_03_L3_velocyto_introns.S20128.pdf", sep="/"
    )
outfile_legend_S20128 <- paste(
    outfile_dir, "yap_scrna_Level_03_L3_velocyto_introns.S20128.legend.pdf", sep="/"
    )
outfile_obj_S20128 <- paste(
    outfile_dir, "velocyto.S20128.qs2", sep="/"
    )

data <- qs_read(infile_path)
table(data@meta.data$L3_cluster_id)
group.by <- "L3_cluster_id"
ldat <- ReadVelocity(infile_loom)

# -- cell name matching --
loom_cellnames <- colnames(ldat$spliced)
sample_names <- detect_sample_prefixes(loom_cellnames, n_samples = 2)
cat("Detected prefixes:", paste(sample_names, collapse=", "), "\n")

for (i in 1:length(sample_names)) {
    current_name <- sprintf("%s:(.*)x", sample_names[i])
    loom_cellnames <- gsub(current_name, sprintf("\\1-1_%i", i), loom_cellnames)
}

# -- subset Seurat to cells present in loom --
n_missing <- sum(!colnames(data) %in% loom_cellnames)
if (n_missing > 0) {
    cat(sprintf("Note: %d / %d Seurat cells missing from loom, subsetting\n", 
                n_missing, ncol(data)))
}
data_sub <- data[, colnames(data) %in% loom_cellnames]
cat(sprintf("Proceeding with %d cells\n", ncol(data_sub)))

data_S20127 <- data_sub[, grepl("_1$", colnames(data_sub))]
data_S20128 <- data_sub[, grepl("_2$", colnames(data_sub))]
loom_S20127 <- loom_cellnames[grepl("_1$", loom_cellnames)]
loom_S20128 <- loom_cellnames[grepl("_2$", loom_cellnames)]

recode_D10025 <- c(
    "LEC_01"="LEC", 
    "preLEC_01"="preLEC", 
    "hmVEC_01"="hmsVEC",
    "mVEC_01"="msVEC", 
    "cVEC_01"="cvpVEC", 
    "iVEC_01"="IL_VEC"
)
order_D10025 <- c(
    "LEC", 
    "preLEC", 
    "hmsVEC",
    "msVEC", 
    "cvpVEC", 
    "IL_VEC"
)
col_map_D10025 <- c(
    "LEC"     = "#41ae76",
    "preLEC"  = "#a1d99b",
    "hmsVEC"  = "#9ecae1",
    "msVEC"   = "#6baed6",
    "cvpVEC"  = "#4292c6",
    "IL_VEC"  = "#084594"
)                   

# combined - all samples, with full colour map and ordering
run_velocyto(
    data_sub=data_sub, 
    loom_cellnames=loom_cellnames, 
    ldat=ldat,
    outfile_path=outfile_path_D10025, 
    outfile_legend=outfile_legend_D10025, 
    outfile_obj=outfile_obj_D10025,
    label="D10025",
    recode_ids=recode_D10025,
    col_map=col_map_D10025,
    level_order=order_D10025
)

recode_S20127 <- recode_D10025
order_S20127 <- order_D10025
col_map_S20127 <- col_map_D10025
loom_S20127 <- loom_cellnames

run_velocyto(
    data_sub=data_S20127, 
    loom_cellnames=loom_S20127, 
    ldat=ldat,
    outfile_path=outfile_path_S20127, 
    outfile_legend=outfile_legend_S20127, 
    outfile_obj=outfile_obj_S20127,
    label="S20027_MUT",
    recode_ids=recode_S20127,
    col_map=col_map_S20127,
    level_order=order_S20127
)

recode_S20128 <- recode_D10025
order_S20128 <- order_D10025
col_map_S20128 <- col_map_D10025
loom_S20128 <- loom_cellnames

run_velocyto(
    data_sub=data_S20128, 
    loom_cellnames=loom_S20128, 
    ldat=ldat,
    outfile_path=outfile_path_S20128, 
    outfile_legend=outfile_legend_S20128, 
    outfile_obj=outfile_obj_S20128,
    label="S20028_WT",
    recode_ids=recode_S20128,
    col_map=col_map_S20128,
    level_order=order_S20128
)

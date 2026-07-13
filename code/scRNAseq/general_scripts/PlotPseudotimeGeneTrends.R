#! /usr/bin/Rscript

#' @title Plot Gene Expression Trends Along a Pseudotime Trajectory
#' @description
#' Plots gene expression along a pseudotime trajectory for any Seurat object
#' that has a numeric pseudotime column in its metadata.
#'
#' Combined plot: all requested genes plotted together along pseudotime,
#'   points and a GAM-smoothed trend line coloured by gene.
#' Per-gene plots (one per gene): a single gene's expression along pseudotime,
#'   points coloured by cell type/phenotype, with a single black GAM-smoothed
#'   trend line, and the cell-type legend placed directly under the gene-name
#'   title.
#'
#' The GAM-smoothing approach is heavily inspired by the
#' \href{https://github.com/huayc09/SeuratExtend}{SeuratExtend} package, but is
#' re-implemented here in a self-contained way (no external dependency on that
#' package) so this script can be dropped into any project.
#'
#' @note Written with the assistance of Claude (Anthropic).
#'
#' @param seurat.object A Seurat object.
#' @param features Character vector of >= 1 gene names to plot. The combined
#'   plot shows all of them together; one per-gene plot is produced per gene.
#' @param pseudotime.col Name of a numeric column in `seurat.object@meta.data`
#'   holding pseudotime values for the trajectory of interest (e.g. a single
#'   Slingshot lineage such as `"slingPseudotime_1_nostart"`, a Palantir
#'   pseudotime column, or any other numeric pseudotime).
#' @param group.by Name of the metadata column giving the cell type / cell
#'   state grouping to colour the per-gene plots by (e.g. `"L2_Predicted_phenotype"`).
#' @param group.colours Optional named character vector of colours for the
#'   levels of `group.by`. Defaults to a `scales::hue_pal()` palette.
#' @param gene.colours Optional named vector (named by `features`) of light
#'   "fill" colours for the points in the combined plot. Defaults to lightened
#'   versions of `gene.line.colours`.
#' @param gene.line.colours Optional named vector (named by `features`) of
#'   darker colours for the combined plot's smoothed trend lines (and point
#'   outlines). Defaults to a `scales::hue_pal()` palette.
#' @param assay Optional assay to switch to before fetching expression data
#'   (default: keep whatever is currently the default assay).
#' @param slot/layer Which slot/layer to pull expression values from
#'   (default `"data"`, i.e. normalised expression).
#' @param pt.alpha,pt.size Point transparency/size for all panels.
#' @param line.size Base line width for the smoothed trend line(s).
#' @param legend.text.size Font size for legend text (and, for the combined
#'   plot, the legend title) across all panels. Default 12.
#' @param n.pseudotime.points Number of points used for the GAM prediction grid.
#' @param save.dir Optional directory to save PDFs/PNGs of the panels to. If
#'   `NULL` (default), nothing is written to disk.
#' @param project.name Prefix used for saved file names.
#' @param width,height,units,device Passed to `ggplot2::ggsave()`.
#'
#' @return A list with:
#'   \item{combined_plot}{ggplot object, all genes combined, coloured by gene.}
#'   \item{gene_panels}{named list of ggplot objects, one per gene, coloured by `group.by`.}
#'   \item{combined_figure}{a single arranged figure (via `ggpubr::ggarrange`) with the combined plot on top and the per-gene plots below.}
#'   \item{data}{the tidy long-format data frame used for plotting.}
#'   \item{predictions}{the GAM-smoothed prediction data frame used for the trend lines.}
#'
#' @examples
#' \dontrun{
#' res <- PlotPseudotimeGeneTrends(
#'   seurat.object   = seu,
#'   features        = c("gene1", "gene2"),
#'   pseudotime.col  = "slingPseudotime_1",
#'   group.by        = "cell_type",
#'   group.colours   = my_cell_type_colours,
#'   save.dir        = "figures/",
#'   project.name    = "my_project"
#' )
#' res$combined_figure
#' }
#'
#' @export
PlotPseudotimeGeneTrends <- function(seurat.object,
                                      features,
                                      pseudotime.col,
                                      group.by,
                                      group.colours = NULL,
                                      gene.colours = NULL,
                                      gene.line.colours = NULL,
                                      assay = NULL,
                                      slot = "data",
                                      pt.alpha = 0.7,
                                      pt.size = 1,
                                      line.size = 0.5,
                                      legend.text.size = 12,
                                      n.pseudotime.points = 200,
                                      save.dir = NULL,
                                      project.name = "pseudotime_trends",
                                      width = 18,
                                      height = 7,
                                      units = "cm",
                                      device = "pdf") {

  # ---- dependencies ----
  requireNamespace("Seurat")
  requireNamespace("ggplot2")
  requireNamespace("mgcv")
  requireNamespace("reshape2")
  requireNamespace("scales")
  requireNamespace("ggpubr")

  # ---- validate inputs ----
  if (!inherits(seurat.object, "Seurat")) {
    stop("`seurat.object` must be a Seurat object.")
  }
  if (length(features) < 1) {
    stop("Please provide at least one gene in `features`.")
  }
  if (!pseudotime.col %in% colnames(seurat.object@meta.data)) {
    stop(sprintf(
      "'%s' was not found in seurat.object@meta.data. Provide the name of a numeric metadata column holding pseudotime values for the trajectory of interest.",
      pseudotime.col
    ))
  }
  if (!group.by %in% colnames(seurat.object@meta.data)) {
    stop(sprintf("'%s' was not found in seurat.object@meta.data.", group.by))
  }
  if (!is.numeric(seurat.object@meta.data[[pseudotime.col]])) {
    stop(sprintf("'%s' must be a numeric metadata column.", pseudotime.col))
  }

  if (!is.null(assay)) {
    Seurat::DefaultAssay(seurat.object) <- assay
  }

  # ---- fetch expression data (compatible with both older `slot` and Seurat v5 `layer`) ----
  expr.mat <- tryCatch(
    Seurat::FetchData(seurat.object, vars = features, layer = slot),
    error = function(e) Seurat::FetchData(seurat.object, vars = features, slot = slot)
  )

  # ---- assemble tidy long-format data frame ----
  meta.df <- seurat.object@meta.data[, c(pseudotime.col, group.by), drop = FALSE]
  colnames(meta.df) <- c("Pseudotime", "CellType")

  # preserve existing factor level order for the grouping variable if present
  if (is.factor(seurat.object@meta.data[[group.by]])) {
    group.levels <- levels(droplevels(seurat.object@meta.data[[group.by]]))
  } else {
    group.levels <- unique(as.character(meta.df$CellType))
  }
  meta.df$CellType <- factor(as.character(meta.df$CellType), levels = group.levels)

  # make sure expression rows line up with metadata rows regardless of any
  # internal reordering FetchData() may perform
  expr.mat <- expr.mat[rownames(meta.df), , drop = FALSE]

  df <- cbind(meta.df, expr.mat)

  df.long <- reshape2::melt(
    df,
    id.vars = c("Pseudotime", "CellType"),
    measure.vars = features,
    variable.name = "Gene",
    value.name = "Expression"
  )
  df.long$Gene <- factor(df.long$Gene, levels = features)
  df.long <- df.long[!is.na(df.long$Pseudotime) & !is.na(df.long$Expression), ]

  # ---- colours ----
  .lighten_colour <- function(colour, amount = 0.55) {
    rgb.col <- grDevices::col2rgb(colour) / 255
    mixed <- rgb.col * (1 - amount) + amount
    grDevices::rgb(mixed[1, ], mixed[2, ], mixed[3, ])
  }

  if (is.null(group.colours)) {
    group.colours <- stats::setNames(scales::hue_pal()(length(group.levels)), group.levels)
  } else if (is.null(names(group.colours))) {
    names(group.colours) <- group.levels[seq_along(group.colours)]
  }

  if (is.null(gene.line.colours)) {
    gene.line.colours <- stats::setNames(scales::hue_pal()(length(features)), features)
  } else if (is.null(names(gene.line.colours))) {
    names(gene.line.colours) <- features[seq_along(gene.line.colours)]
  }

  if (is.null(gene.colours)) {
    gene.colours <- stats::setNames(
      vapply(gene.line.colours[features], .lighten_colour, character(1)),
      features
    )
  } else if (is.null(names(gene.colours))) {
    names(gene.colours) <- features[seq_along(gene.colours)]
  }

  # ---- GAM smoothing (self-contained equivalent of perform_gam_analysis()) ----
  .fit_gam_trend <- function(sub.df, n = n.pseudotime.points) {
    model <- mgcv::gam(Expression ~ s(Pseudotime), data = sub.df)
    grid.df <- data.frame(
      Pseudotime = seq(min(sub.df$Pseudotime), max(sub.df$Pseudotime), length.out = n)
    )
    grid.df$Predicted <- as.numeric(stats::predict(model, newdata = grid.df))
    grid.df
  }

  predictions <- do.call(rbind, lapply(features, function(gene) {
    sub.df <- df.long[df.long$Gene == gene, ]
    pred <- .fit_gam_trend(sub.df)
    pred$Gene <- gene
    pred
  }))
  predictions$Gene <- factor(predictions$Gene, levels = features)

  # ---- combined plot: all genes combined, coloured by gene ----
  combined_plot <- ggplot2::ggplot() +
    ggplot2::labs(x = "Pseudotime", y = "Gene expression", fill = "Gene", colour = "Gene") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = legend.text.size),
      legend.title = ggplot2::element_text(size = legend.text.size)
    ) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::geom_point(
      data = df.long,
      ggplot2::aes(x = Pseudotime, y = Expression, fill = Gene, colour = Gene),
      alpha = pt.alpha, size = pt.size, stroke = 0, shape = 21
    ) +
    ggplot2::scale_fill_manual(values = gene.colours) +
    ggplot2::geom_line(
      data = predictions,
      ggplot2::aes(x = Pseudotime, y = Predicted, colour = Gene),
      linewidth = line.size * 1.5
    ) +
    ggplot2::scale_colour_manual(values = gene.line.colours)

  # ---- per-gene plots: one per gene, coloured by group.by ----
  gene.panels <- lapply(features, function(gene) {
    ggplot2::ggplot() +
      ggplot2::labs(x = "Pseudotime", y = "Gene expression", title = gene,
                    fill = NULL, colour = NULL) +
      ggplot2::theme_classic() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "italic", hjust = 0.5),
        legend.position = "top",
        legend.title = ggplot2::element_blank(),
        legend.key.size = ggplot2::unit(0.4, "lines"),
        legend.spacing.x = ggplot2::unit(0.15, "cm"),
        legend.spacing.y = ggplot2::unit(0.05, "cm"),
        legend.margin = ggplot2::margin(t = 0, b = 0),
        legend.box.spacing = ggplot2::unit(0.1, "cm"),
        legend.text = ggplot2::element_text(size = legend.text.size)
      ) +
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::geom_point(
        data = df.long[df.long$Gene == gene, ],
        ggplot2::aes(x = Pseudotime, y = Expression, fill = CellType, colour = CellType),
        alpha = pt.alpha, size = pt.size, stroke = 0, shape = 21
      ) +
      ggplot2::scale_fill_manual(values = group.colours) +
      ggplot2::scale_colour_manual(values = group.colours) +
      ggplot2::guides(
        fill = ggplot2::guide_legend(override.aes = list(size = 3, alpha = 1, stroke = 0)),
        colour = ggplot2::guide_legend(override.aes = list(size = 3, alpha = 1, stroke = 0))
      ) +
      ggplot2::geom_line(
        data = predictions[predictions$Gene == gene, ],
        ggplot2::aes(x = Pseudotime, y = Predicted),
        colour = "black", linewidth = line.size
      )
  })
  names(gene.panels) <- features

  # ---- combined figure: combined plot on top, one gene panel per gene below ----
  bottom_row <- ggpubr::ggarrange(
    plotlist = gene.panels,
    ncol = length(gene.panels), nrow = 1
  )
  combined_figure <- ggpubr::ggarrange(
    combined_plot, bottom_row,
    nrow = 2
  )

  # ---- optional saving ----
  if (!is.null(save.dir)) {
    if (!dir.exists(save.dir)) dir.create(save.dir, recursive = TRUE)

    ggplot2::ggsave(
      filename = file.path(save.dir, paste0(project.name, "_pseudotime_trends_combined.", device)),
      plot = combined_figure, device = device, width = width, height = height * 2, units = units
    )
    ggplot2::ggsave(
      filename = file.path(save.dir, paste0(project.name, "_pseudotime_trends_combined_plot.", device)),
      plot = combined_plot, device = device, width = width, height = height, units = units
    )
    for (gene in features) {
      ggplot2::ggsave(
        filename = file.path(save.dir, paste0(project.name, "_pseudotime_trends_", gene, ".", device)),
        plot = gene.panels[[gene]], device = device,
        width = width / max(1, length(features)), height = height, units = units
      )
    }
  }

  return(list(
    combined_plot = combined_plot,
    gene_panels = gene.panels,
    combined_figure = combined_figure,
    data = df.long,
    predictions = predictions
  ))
}

#' Plot a Bonsai tree, optionally colored by a grouping variable
#'
#' A thin, opinionated wrapper around \code{ape::plot.phylo(type = "fan")}
#' tuned for single-cell scale (tens of thousands of tips):
#' \itemize{
#'   \item \code{no.margin = TRUE} so \code{ape} doesn't reserve space for
#'     tip labels it isn't drawing, which otherwise leaves the tree
#'     off-center with a lot of wasted canvas.
#'   \item Branches de-emphasized (thin, light grey) so tip colors read
#'     clearly against them.
#'   \item An automatically generated color palette that maximizes hue
#'     separation across groups, rather than something like
#'     \code{RColorBrewer::brewer.pal(n, "Paired")}, which deliberately
#'     pairs up similar hues -- the wrong choice when every group needs to
#'     look different from every other one.
#'   \item A legend that doesn't overlap the tree by default.
#'   \item Direct \code{file} output, since interactively rendering tens of
#'     thousands of points in R's graphics device can be very slow.
#' }
#'
#' @param bonsai_tree A \code{bonsai_tree} object, as returned by
#'   \code{\link{bonsai_read_tree}}, or a plain \code{ape::phylo} object.
#' @param tip_groups Optional named character/factor vector giving a group
#'   label for each cell, named by cell ID (matching
#'   \code{bonsai_tree$phylo$tip.label}, e.g. Seurat cell barcodes -- build
#'   one with e.g. \code{setNames(seurat_obj$Zone, colnames(seurat_obj))}).
#'   If \code{NULL} (default), tips are plotted uncolored. Tips missing
#'   from \code{tip_groups} (or with \code{NA}) are colored grey.
#' @param colors Optional named character vector mapping group labels (from
#'   \code{tip_groups}) to hex colors. If \code{NULL} (default) and
#'   \code{tip_groups} is given, colors are generated automatically via
#'   evenly-spaced HCL hues.
#' @param file Optional file path to save the plot to (\code{.png} or
#'   \code{.pdf}, inferred from the extension). If \code{NULL} (default),
#'   plots to the current graphics device. Recommended at real single-cell
#'   scale.
#' @param width,height,res Passed to \code{png()}/\code{pdf()} when
#'   \code{file} is given (\code{width}/\code{height} in pixels for
#'   \code{.png}, converted to inches via \code{res} for \code{.pdf}).
#'   Defaults \code{3000}, \code{3000}, \code{300}.
#' @param legend Logical. Whether to draw a legend for \code{tip_groups}.
#'   Default \code{TRUE} (ignored if \code{tip_groups} is \code{NULL}).
#' @param legend_position Passed as \code{legend()}'s first argument (e.g.
#'   \code{"topleft"}, \code{"bottomright"}). Default \code{"topleft"}.
#' @param legend_ncol Number of columns in the legend. Default \code{2},
#'   since a single column can run tall with many groups.
#' @param legend_cex Text size in the legend. Default \code{1} (base R's
#'   own default) -- reduce this (e.g. \code{0.6}) to shrink the legend box
#'   so it covers less of the tree, especially with many groups.
#' @param legend_pt_cex Point/symbol size in the legend, independent of
#'   \code{legend_cex}. Default \code{1}.
#' @param type Passed to \code{ape::plot.phylo()}. Default \code{"fan"};
#'   \code{"unrooted"} is another option worth trying.
#' @param tip_cex Point size for tips. Default \code{0.25} -- reduce this
#'   for a very large number of tips.
#' @param edge_width,edge_color De-emphasize branches relative to tip
#'   colors. Defaults \code{0.1} and \code{"grey80"}.
#' @param ... Additional arguments passed to \code{ape::plot.phylo()}.
#'
#' @return Invisibly, the group-to-color mapping used (a named character
#'   vector, or \code{NULL} if \code{tip_groups} was not given) -- reuse it
#'   to keep a matching Seurat \code{DimPlot()} (or any other plot) colored
#'   consistently with the tree.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' tree <- bonsai_read_tree(result, benv)
#' zone_by_cell <- setNames(seurat_obj$Zone, colnames(seurat_obj))
#' colors <- bonsai_plot_tree(tree, tip_groups = zone_by_cell, file = "tree.png")
#'
#' # Reuse the same colors for a matching UMAP:
#' Seurat::DimPlot(seurat_obj, group.by = "Zone", cols = colors)
#' }
bonsai_plot_tree <- function(bonsai_tree,
                              tip_groups = NULL,
                              colors = NULL,
                              file = NULL,
                              width = 3000, height = 3000, res = 300,
                              legend = TRUE,
                              legend_position = "topleft",
                              legend_ncol = 2,
                              legend_cex = 1,
                              legend_pt_cex = 1,
                              type = "fan",
                              tip_cex = 0.25,
                              edge_width = 0.1,
                              edge_color = "grey80",
                              ...) {

  phylo <- if (inherits(bonsai_tree, "bonsai_tree")) {
    bonsai_tree$phylo
  } else if (inherits(bonsai_tree, "phylo")) {
    bonsai_tree
  } else {
    cli::cli_abort("{.arg bonsai_tree} must be a {.cls bonsai_tree} object (from {.fn bonsai_read_tree}) or an {.cls phylo} object.")
  }

  tip_colors <- NULL
  color_map <- NULL

  if (!is.null(tip_groups)) {
    if (is.null(names(tip_groups))) {
      cli::cli_abort("{.arg tip_groups} must be a named vector (names = cell IDs matching tree tip labels).")
    }

    tip_group_vals <- as.character(tip_groups[phylo$tip.label])
    group_names <- sort(unique(stats::na.omit(tip_group_vals)))

    if (is.null(colors)) {
      n <- length(group_names)
      color_map <- stats::setNames(
        grDevices::hcl(h = seq(15, 375, length.out = n + 1)[seq_len(n)], c = 100, l = 60),
        group_names
      )
    } else {
      missing_groups <- setdiff(group_names, names(colors))
      if (length(missing_groups) > 0) {
        cli::cli_abort("{.arg colors} is missing entries for group(s): {.val {missing_groups}}.")
      }
      color_map <- colors
    }

    tip_colors <- color_map[tip_group_vals]
    tip_colors[is.na(tip_colors)] <- "grey70"

    n_matched <- sum(!is.na(tip_group_vals))
    if (n_matched == 0) {
      cli::cli_warn("None of {.arg tip_groups}'s names matched {.code bonsai_tree$phylo$tip.label} -- check that its names are cell IDs, not e.g. row indices.")
    } else if (n_matched < length(phylo$tip.label)) {
      cli::cli_alert_info("{n_matched}/{length(phylo$tip.label)} tips matched a group in {.arg tip_groups}; the rest are plotted grey.")
    }
  }

  plot_it <- function() {
    withr::local_par(mar = c(0, 0, 0, 0))
    ape::plot.phylo(
      phylo, type = type, show.tip.label = FALSE, no.margin = TRUE,
      edge.width = edge_width, edge.color = edge_color, ...
    )
    if (!is.null(tip_colors)) {
      ape::tiplabels(pch = 20, col = tip_colors, cex = tip_cex)
    }
    if (!is.null(color_map) && legend) {
      graphics::legend(legend_position, legend = names(color_map), pch = 16,
                       col = color_map, ncol = legend_ncol, xpd = NA,
                       cex = legend_cex, pt.cex = legend_pt_cex)
    }
  }

  if (!is.null(file)) {
    ext <- tolower(fs::path_ext(file))
    if (ext == "png") {
      grDevices::png(file, width = width, height = height, res = res)
    } else if (ext == "pdf") {
      grDevices::pdf(file, width = width / res, height = height / res)
    } else {
      cli::cli_abort("{.arg file} must end in {.val .png} or {.val .pdf}, got {.val {file}}.")
    }
    on.exit(grDevices::dev.off(), add = TRUE)
    plot_it()
    cli::cli_alert_success("Tree plot written to {file}")
  } else {
    plot_it()
  }

  invisible(color_map)
}

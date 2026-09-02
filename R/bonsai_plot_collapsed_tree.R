#' Plot a collapsed (cluster-level) Bonsai tree
#'
#' A companion to \code{\link{bonsai_plot_tree}} for the small, cluster-level
#' trees \code{\link{bonsai_collapse_tree}} produces. At that scale a
#' color-coded legend is unnecessary and just adds another element to
#' overlap the tree -- with few enough tips, each one gets a direct text
#' label instead. Point size is scaled to how many original tips collapsed
#' into each cluster by default, so you can see at a glance which clusters
#' are large and which are small even though that information isn't in the
#' topology itself.
#'
#' @param collapsed_tree A \code{bonsai_collapsed_tree} object, as returned
#'   by \code{\link{bonsai_collapse_tree}}.
#' @param colors Optional named character vector mapping cluster labels to
#'   hex colors. If \code{NULL} (default), colors are generated
#'   automatically via evenly-spaced HCL hues (same approach as
#'   \code{\link{bonsai_plot_tree}}).
#' @param size_by_n Logical. If \code{TRUE} (default), point size reflects
#'   \code{collapsed_tree$sizes} (the number of original tips each cluster
#'   represents), scaled by square root so a 100x larger cluster isn't
#'   drawn 100x bigger. If \code{FALSE}, all points are the same size
#'   (\code{mean(cex_range)}).
#' @param cex_range Numeric length-2 vector, the point size range used when
#'   \code{size_by_n = TRUE}. Default \code{NULL}, which shrinks the upper
#'   bound as the number of clusters grows beyond ~10 -- each additional
#'   cluster shares less of the circle's circumference in \code{type =
#'   "fan"}, so a size tuned by eye for a handful of clusters starts
#'   overlapping its neighbors once there are many (verified). Pass an
#'   explicit \code{c(min, max)} to override this scaling.
#' @param label_cex Text size for cluster labels. Default \code{NULL}, which
#'   shrinks it the same way and for the same reason as \code{cex_range}.
#' @param label_offset Gap between each point and its label, in the same
#'   units as the tree's branch lengths. Default \code{NULL}, which picks a
#'   value proportional to the tree's own depth so it scales sensibly
#'   regardless of your data's branch-length units.
#' @param file Optional file path to save the plot to (\code{.png} or
#'   \code{.pdf}, inferred from the extension). If \code{NULL} (default),
#'   plots to the current graphics device.
#' @param width,height,res Passed to \code{png()}/\code{pdf()} when
#'   \code{file} is given. \code{width}/\code{height} default to \code{NULL},
#'   which scales the canvas up with the number of clusters -- point size in
#'   \code{cex} units is roughly a fixed physical size regardless of canvas
#'   size, so a bigger canvas is what actually gives more clusters more room
#'   before they start overlapping (verified: a fixed \code{2000} stayed
#'   crowded at 10 clusters even after shrinking point/label size, while
#'   doubling the canvas at the same point size fixed it directly). A fixed
#'   \code{2000} was tuned by eye for a handful of clusters and is kept as
#'   the floor. \code{res} defaults to \code{300}.
#' @param type Passed to \code{ape::plot.phylo()}. Default \code{"fan"}.
#' @param edge_width,edge_color De-emphasize branches relative to points.
#'   Defaults \code{1} and \code{"grey50"} -- less extreme than
#'   \code{bonsai_plot_tree()}'s defaults, since a collapsed tree has few
#'   enough branches that thin/light ones would be hard to see at all.
#' @param ... Additional arguments passed to \code{ape::plot.phylo()}.
#'
#' @return Invisibly, the cluster-to-color mapping used (a named character
#'   vector) -- reuse it to keep other plots of the same clusters
#'   consistently colored.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' reduced <- bonsai_collapse_tree(tree, tip_groups)
#' bonsai_plot_collapsed_tree(reduced, file = "reduced_tree.png")
#' }
bonsai_plot_collapsed_tree <- function(collapsed_tree,
                                        colors = NULL,
                                        size_by_n = TRUE,
                                        cex_range = NULL,
                                        label_cex = NULL,
                                        label_offset = NULL,
                                        file = NULL,
                                        width = NULL, height = NULL, res = 300,
                                        type = "fan",
                                        edge_width = 1,
                                        edge_color = "grey50",
                                        ...) {

  if (!inherits(collapsed_tree, "bonsai_collapsed_tree")) {
    cli::cli_abort("{.arg collapsed_tree} must be a {.cls bonsai_collapsed_tree} object, as returned by {.fn bonsai_collapse_tree}.")
  }
  phylo <- collapsed_tree$phylo
  sizes <- collapsed_tree$sizes[phylo$tip.label]
  n_tips <- length(phylo$tip.label)

  # A fixed cex_range/label_cex tuned by eye for a handful of clusters
  # overlaps badly once there are many more of them sharing the same 360
  # degrees (verified: cex_range = c(0.8, 4) at 10 clusters produced
  # touching/overlapping points and collided labels). Shrink both
  # proportionally once there are more than a comfortable ~10 clusters,
  # unless the caller passed an explicit value.
  crowding_factor <- min(1, 10 / n_tips)
  if (is.null(cex_range)) {
    cex_range <- c(0.8, max(1.5, 4 * crowding_factor))
  }
  if (is.null(label_cex)) {
    label_cex <- max(0.5, 0.8 * crowding_factor)
  }
  if (is.null(width)) width <- max(2000, 450 * n_tips)
  if (is.null(height)) height <- max(2000, 450 * n_tips)

  if (is.null(colors)) {
    n <- length(phylo$tip.label)
    color_map <- stats::setNames(
      grDevices::hcl(h = seq(15, 375, length.out = n + 1)[seq_len(n)], c = 100, l = 60),
      phylo$tip.label
    )
  } else {
    missing_colors <- setdiff(phylo$tip.label, names(colors))
    if (length(missing_colors) > 0) {
      cli::cli_abort("{.arg colors} is missing entries for cluster(s): {.val {missing_colors}}.")
    }
    color_map <- colors
  }
  tip_colors <- color_map[phylo$tip.label]

  if (size_by_n) {
    w <- sqrt(sizes)
    rng <- range(w)
    tip_cex <- if (diff(rng) == 0) {
      rep(mean(cex_range), length(w))
    } else {
      (w - rng[1]) / diff(rng) * diff(cex_range) + cex_range[1]
    }
  } else {
    tip_cex <- rep(mean(cex_range), length(phylo$tip.label))
  }

  max_depth <- max(ape::node.depth.edgelength(phylo))
  if (is.null(label_offset)) {
    # Scaled by the largest point actually drawn, not just tree depth --
    # a fixed small offset is fine for small/uniform points, but a big
    # size_by_n = TRUE point (a highly populous cluster) can be large
    # enough in its own right to overlap the start of its label if the
    # offset doesn't grow with it too (verified -- this happened with a
    # depth-only offset).
    label_offset <- max_depth * 0.02 * max(tip_cex)
  }

  # Unlike bonsai_plot_tree(), deliberately not forcing no.margin = TRUE
  # here: that setting removes all reserved space around the plot, which is
  # correct when there are no tip labels to draw (bonsai_plot_tree()'s
  # case) but wrong here, where labels are the whole point.
  #
  # ape::plot.phylo()'s own margin auto-sizing (no.margin = FALSE, the
  # default) turns out not to be enough on its own for type = "fan": labels
  # near the circle's edges are drawn at rotated angles, and the auto-sized
  # margin under-accounts for how much horizontal/vertical room a rotated
  # label actually needs, so long labels (e.g. after renaming) still get
  # clipped at the device edge (verified -- this happened even with
  # no.margin = FALSE alone). Explicitly padding x.lim/y.lim well beyond
  # the tree's natural radius is the standard fix for this in ape.
  #
  # A fixed multiplier on max_depth alone isn't enough once labels get
  # longer than a short generic "cl_N" -- e.g. a bonsai_dominant_labels()
  # result like "Hepatocytes (61%)" got clipped at the plot edge with a
  # fixed 1.4x pad (verified). Scale the extra padding by the longest
  # label's character count too, so it grows with the actual text.
  longest_label <- max(nchar(phylo$tip.label))
  pad <- max_depth * 1.4 + label_offset + max_depth * 0.09 * label_cex * longest_label
  xylim <- c(-pad, pad)

  plot_it <- function() {
    ape::plot.phylo(
      phylo, type = type, show.tip.label = TRUE,
      edge.width = edge_width, edge.color = edge_color,
      cex = label_cex, label.offset = label_offset, tip.color = tip_colors,
      x.lim = xylim, y.lim = xylim,
      ...
    )
    ape::tiplabels(pch = 16, col = tip_colors, cex = tip_cex, adj = 0.5)
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
    cli::cli_alert_success("Collapsed tree plot written to {file}")
  } else {
    plot_it()
  }

  invisible(color_map)
}

#' Label collapsed-tree clusters by their dominant annotation
#'
#' A geometry-based clustering (e.g. \code{\link{bonsai_cluster_tree}}'s
#' \code{cl_0}, \code{cl_1}, ...) tells you nothing about what each cluster
#' actually *is* biologically. If you separately have a real annotation for
#' each cell (a cell-type call, a Zone, anything), this finds the majority
#' value of that annotation within each cluster's original members and
#' returns labels built from it -- feed the result straight into
#' \code{\link{bonsai_rename_clusters}}.
#'
#' Each label includes the dominant value's share of its cluster by default
#' (e.g. \code{"Hepatocytes (94\%)"}), since that purity is informative on
#' its own -- a cluster that's 94\% one cell type is a very different result
#' from one that's 52\%, and both would otherwise show up as an identical
#' label. If two different clusters happen to share the same dominant
#' annotation, each is disambiguated with its original cluster label so
#' \code{\link{bonsai_rename_clusters}} doesn't reject the result for
#' having duplicate names.
#'
#' Matches cells to clusters via \code{collapsed_tree$tip_groups}, which
#' records the grouping \code{\link{bonsai_collapse_tree}} was built from --
#' keyed by \code{collapsed_tree}'s *current* tip labels. Call this before
#' renaming the tree by any other means, while those still are the tree's
#' original cluster identifiers.
#'
#' @param collapsed_tree A \code{bonsai_collapsed_tree} object, as returned
#'   by \code{\link{bonsai_collapse_tree}}.
#' @param annotation Named character/factor vector, cell ID -> annotation
#'   label (e.g. a cell-type call), covering (at least) the same cells
#'   \code{collapsed_tree} was built from.
#' @param show_pct Logical. If \code{TRUE} (default), append the dominant
#'   value's percentage share of the cluster, e.g. \code{"LSEC (87\%)"}.
#'
#' @return A named character vector suitable for
#'   \code{\link{bonsai_rename_clusters}}'s \code{rename} argument: names
#'   are \code{collapsed_tree}'s current cluster labels, values are the
#'   dominant-annotation-based labels.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' reduced <- bonsai_collapse_tree(tree, tip_groups)
#' dominant <- bonsai_dominant_labels(reduced, cluster_by_cell)
#' reduced <- bonsai_rename_clusters(reduced, dominant)
#' bonsai_plot_collapsed_tree(reduced, file = "reduced_by_cell_type.png")
#' }
bonsai_dominant_labels <- function(collapsed_tree, annotation, show_pct = TRUE) {
  if (!inherits(collapsed_tree, "bonsai_collapsed_tree")) {
    cli::cli_abort("{.arg collapsed_tree} must be a {.cls bonsai_collapsed_tree} object, as returned by {.fn bonsai_collapse_tree}.")
  }
  if (is.null(names(annotation))) {
    cli::cli_abort("{.arg annotation} must be a named vector (names = cell IDs).")
  }

  tip_groups <- collapsed_tree$tip_groups
  cell_ids <- names(tip_groups)
  group_vals <- as.character(tip_groups)
  ann_vals <- as.character(annotation[cell_ids])

  current_labels <- collapsed_tree$phylo$tip.label
  matched <- current_labels %in% unique(stats::na.omit(group_vals))
  if (!any(matched)) {
    cli::cli_abort(c(
      "None of {.arg collapsed_tree}'s current tip labels match its own recorded cluster assignments.",
      "i" = "Call {.fn bonsai_dominant_labels} on a freshly-collapsed tree, before renaming it any other way."
    ))
  }

  n_ann_matched <- sum(!is.na(ann_vals))
  if (n_ann_matched == 0) {
    cli::cli_warn("None of {.arg annotation}'s names matched {.arg collapsed_tree}'s cell IDs -- check that its names are cell IDs, not e.g. row indices.")
  } else if (n_ann_matched < length(cell_ids)) {
    cli::cli_alert_info("{n_ann_matched}/{length(cell_ids)} cells matched a value in {.arg annotation}; the rest are excluded from the majority vote.")
  }

  dominant <- stats::setNames(character(length(current_labels)), current_labels)
  for (g in current_labels) {
    member_ann <- stats::na.omit(ann_vals[group_vals == g])
    if (length(member_ann) == 0) {
      dominant[g] <- g
      next
    }
    tbl <- sort(table(member_ann), decreasing = TRUE)
    top_type <- names(tbl)[1]
    pct <- round(100 * tbl[1] / sum(tbl))
    dominant[g] <- if (show_pct) sprintf("%s (%d%%)", top_type, pct) else top_type
  }

  dup <- duplicated(dominant) | duplicated(dominant, fromLast = TRUE)
  if (any(dup)) {
    dominant[dup] <- paste0(dominant[dup], " [", names(dominant)[dup], "]")
  }

  dominant
}

#' Collapse a Bonsai tree to one representative tip per group
#'
#' Reduces a (potentially huge) cell-level tree down to one tip per group --
#' e.g. per cluster from \code{\link{bonsai_cluster_tree}} or
#' \code{\link{bonsai_cluster_by_annotation}}, or any other grouping -- by
#' keeping one representative tip per group and dropping the rest.
#' \code{ape::drop.tip()} recomputes branch lengths for what remains, so the
#' reduced tree's topology and distances between group representatives stay
#' phylogenetically meaningful; this isn't just a random subset.
#'
#' This assumes each group is monophyletic (its members form one connected
#' clade in the tree) -- true by construction for
#' \code{\link{bonsai_cluster_tree}}'s output (it cuts subtrees directly),
#' but not guaranteed for \code{\link{bonsai_cluster_by_annotation}}'s
#' annotation-optimized clustering, which is free to group cells that
#' aren't monophyletic. If you use this on a non-monophyletic grouping, the
#' chosen representative is arbitrary (the first member found) and the
#' reduced tree's shape can be misleading for that group.
#'
#' @param bonsai_tree A \code{bonsai_tree} object, as returned by
#'   \code{\link{bonsai_read_tree}}, or a plain \code{ape::phylo} object.
#' @param tip_groups Named character/factor vector giving a group label per
#'   tip, named by tip label (e.g. cell ID). Tips missing from
#'   \code{tip_groups} (or with \code{NA}) are dropped from the reduced
#'   tree entirely, rather than kept as their own singleton group.
#' @param labels Optional named vector mapping group values (as they appear
#'   in \code{tip_groups}) to display names used as the reduced tree's tip
#'   labels. Default \code{NULL} (uses the group values themselves as-is).
#'   Can also be changed later via \code{\link{bonsai_rename_clusters}}
#'   without recomputing the collapse.
#'
#' @return An S3 object of class \code{bonsai_collapsed_tree}, a list with
#'   elements:
#'   \describe{
#'     \item{phylo}{The reduced \code{ape::phylo} tree, one tip per group.}
#'     \item{sizes}{Named integer vector, group label -> number of original
#'       tips that collapsed into it.}
#'     \item{tip_groups}{The original \code{tip_groups} argument, kept for
#'       reference.}
#'   }
#'   Pass this to \code{\link{bonsai_plot_collapsed_tree}} to plot it, or
#'   \code{\link{bonsai_rename_clusters}} to rename groups afterward.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' tree_clusters <- bonsai_cluster_tree(result, n_clusters = 10, bonsai_env = benv)
#' reduced <- bonsai_collapse_tree(
#'   tree,
#'   setNames(tree_clusters$annot_bnsi_cluster_n10, rownames(tree_clusters))
#' )
#' bonsai_plot_collapsed_tree(reduced)
#' }
bonsai_collapse_tree <- function(bonsai_tree, tip_groups, labels = NULL) {

  phylo <- if (inherits(bonsai_tree, "bonsai_tree")) {
    bonsai_tree$phylo
  } else if (inherits(bonsai_tree, "phylo")) {
    bonsai_tree
  } else {
    cli::cli_abort("{.arg bonsai_tree} must be a {.cls bonsai_tree} object (from {.fn bonsai_read_tree}) or an {.cls phylo} object.")
  }
  if (is.null(names(tip_groups))) {
    cli::cli_abort("{.arg tip_groups} must be a named vector (names = tip labels).")
  }

  group_vals <- as.character(tip_groups[phylo$tip.label])
  n_na <- sum(is.na(group_vals))
  if (n_na > 0) {
    cli::cli_alert_info("{n_na}/{length(group_vals)} tips have no group in {.arg tip_groups} and will be dropped from the reduced tree.")
  }

  group_names <- sort(unique(stats::na.omit(group_vals)))
  if (length(group_names) == 0) {
    cli::cli_abort("No tips have a non-NA group in {.arg tip_groups}.")
  }

  if (!is.null(labels)) {
    missing_labels <- setdiff(group_names, names(labels))
    if (length(missing_labels) > 0) {
      cli::cli_abort("{.arg labels} is missing entries for group(s): {.val {missing_labels}}.")
    }
    display_names <- unname(labels[group_names])
  } else {
    display_names <- group_names
  }
  if (anyDuplicated(display_names)) {
    cli::cli_abort("Resulting tip labels have duplicates ({.val {display_names[duplicated(display_names)]}}) -- {.arg labels} must map to distinct names.")
  }

  sizes <- vapply(group_names, function(g) sum(group_vals == g, na.rm = TRUE), integer(1))
  keep_tips <- vapply(group_names, function(g) phylo$tip.label[which(group_vals == g)][1], character(1))

  collapsed_phylo <- ape::drop.tip(phylo, setdiff(phylo$tip.label, keep_tips))
  collapsed_phylo$tip.label <- display_names[match(collapsed_phylo$tip.label, keep_tips)]

  names(sizes) <- display_names

  cli::cli_alert_success("Collapsed {length(phylo$tip.label)} tips to {length(display_names)} groups.")

  structure(
    list(phylo = collapsed_phylo, sizes = sizes, tip_groups = tip_groups),
    class = "bonsai_collapsed_tree"
  )
}

#' @export
print.bonsai_collapsed_tree <- function(x, ...) {
  cli::cli_h3("bonsai_collapsed_tree")
  cli::cli_bullets(c(
    "*" = "{ape::Ntip(x$phylo)} clusters, from {sum(x$sizes)} original tips",
    "*" = "cluster sizes range {min(x$sizes)}-{max(x$sizes)}"
  ))
  invisible(x)
}

#' Rename clusters in a collapsed tree
#'
#' Renames one or more clusters in a \code{bonsai_collapsed_tree}, updating
#' both the tree's tip labels and the cluster-size tracking consistently --
#' safer than editing \code{collapsed_tree$phylo$tip.label} directly, which
#' would silently desynchronize it from \code{collapsed_tree$sizes}.
#'
#' @param collapsed_tree A \code{bonsai_collapsed_tree} object, as returned
#'   by \code{\link{bonsai_collapse_tree}}.
#' @param rename Named character vector: names are the cluster's *current*
#'   labels, values are the new labels. Only clusters you name are
#'   affected -- this doesn't need to cover every cluster.
#'
#' @return The updated \code{bonsai_collapsed_tree}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' reduced <- bonsai_rename_clusters(reduced, c(
#'   "cl_0" = "T cells",
#'   "cl_1" = "Monocytes"
#' ))
#' }
bonsai_rename_clusters <- function(collapsed_tree, rename) {
  if (!inherits(collapsed_tree, "bonsai_collapsed_tree")) {
    cli::cli_abort("{.arg collapsed_tree} must be a {.cls bonsai_collapsed_tree} object, as returned by {.fn bonsai_collapse_tree}.")
  }
  if (is.null(names(rename))) {
    cli::cli_abort("{.arg rename} must be a named vector (names = current cluster labels, values = new labels).")
  }

  current <- collapsed_tree$phylo$tip.label
  missing <- setdiff(names(rename), current)
  if (length(missing) > 0) {
    cli::cli_abort("{.arg rename} names not found among current cluster labels: {.val {missing}}.")
  }

  new_labels <- current
  idx <- match(names(rename), current)
  new_labels[idx] <- unname(rename)
  if (anyDuplicated(new_labels)) {
    cli::cli_abort("Renaming would produce duplicate cluster labels: {.val {new_labels[duplicated(new_labels)]}}.")
  }

  collapsed_tree$phylo$tip.label <- new_labels
  names(collapsed_tree$sizes) <- new_labels

  collapsed_tree
}

#' Attach a Bonsai tree reconstruction back onto a Seurat object
#'
#' Stores the reconstructed tree, per-cell posterior LTQ estimates, and
#' (optionally) cellstates cluster assignments on a Seurat object, keyed
#' and reordered to match the object's own cell order, so downstream
#' Seurat-based analysis continues to work normally.
#'
#' Cell identity matching between Bonsai's tree (leaf names) and the
#' Seurat object's cell barcodes is done \strong{by name, not position} --
#' this function explicitly checks for and reports any mismatches rather
#' than assuming the two are aligned, since silent misalignment here would
#' corrupt every downstream analysis without any error being raised.
#'
#' @param seurat_obj A \code{Seurat} object -- normally the same one
#'   originally passed to \code{\link{bonsai_write_sanity_input}}, though
#'   this isn't checked; only cell-name matching is checked.
#' @param bonsai_tree A \code{bonsai_tree} object, as returned by
#'   \code{\link{bonsai_read_tree}}.
#' @param cellstates_output A \code{cellstates_output} object, as returned
#'   by \code{\link{run_cellstates}}, or \code{NULL} (default) to skip
#'   attaching cellstate cluster labels.
#' @param on_mismatch One of \code{"warn"} (default) or \code{"error"}.
#'   Controls behavior when the tree's leaf names and the Seurat object's
#'   cell names don't match exactly (e.g. if cells were filtered from one
#'   but not the other after the pipeline ran). With \code{"warn"}, only
#'   the intersection of cells is annotated and everything else gets
#'   \code{NA}; with \code{"error"}, the function stops.
#'
#' @return The input \code{seurat_obj}, with:
#'   \describe{
#'     \item{\code{seurat_obj@misc$bonsai}}{A list containing the full
#'       \code{phylo} tree object, \code{edge_info}, \code{vert_info}, and
#'       \code{metadata} (including internal, non-cell nodes).}
#'     \item{\code{seurat_obj@misc$bonsai$leaf_posterior_ltqs}}{A genes x
#'       cells matrix of posterior LTQ point estimates, subset to leaf
#'       (cell) vertices only, reordered to match
#'       \code{colnames(seurat_obj)}. \code{NULL} if
#'       \code{bonsai_tree$posterior_ltqs} was not loaded.}
#'     \item{A new metadata column \code{bonsai_cellstate}}{ (only if
#'       \code{cellstates_output} is supplied), giving each cell's
#'       cellstate cluster label.}
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' seu <- bonsai_to_seurat(seu, tree, cellstates_output = cellstates_out)
#' }
bonsai_to_seurat <- function(seurat_obj,
                              bonsai_tree,
                              cellstates_output = NULL,
                              on_mismatch = c("warn", "error")) {

  on_mismatch <- match.arg(on_mismatch)

  if (!inherits(seurat_obj, "Seurat")) {
    cli::cli_abort("{.arg seurat_obj} must be a Seurat object.")
  }
  if (!inherits(bonsai_tree, "bonsai_tree")) {
    cli::cli_abort("{.arg bonsai_tree} must be a {.cls bonsai_tree} object, as returned by {.fn bonsai_read_tree}.")
  }

  seurat_cells <- colnames(seurat_obj)
  tree_leaves <- bonsai_tree$phylo$tip.label

  missing_from_tree <- setdiff(seurat_cells, tree_leaves)
  missing_from_seurat <- setdiff(tree_leaves, seurat_cells)

  if (length(missing_from_tree) > 0 || length(missing_from_seurat) > 0) {
    msg <- c(
      "Cell name mismatch between the Seurat object and the Bonsai tree's leaves.",
      "x" = "{length(missing_from_tree)} Seurat cells are not in the tree.",
      "x" = "{length(missing_from_seurat)} tree leaves are not in the Seurat object.",
      "i" = "This can happen if cells were filtered from one but not the other after the pipeline ran."
    )
    if (on_mismatch == "error") {
      cli::cli_abort(msg)
    } else {
      cli::cli_warn(c(msg, "i" = "Proceeding with only the intersection annotated; everything else gets NA."))
    }
  }

  common_cells <- intersect(seurat_cells, tree_leaves)
  if (length(common_cells) == 0) {
    cli::cli_abort("No cell names in common between the Seurat object and the Bonsai tree -- cannot proceed.")
  }

  # ---- Attach the tree object + metadata to misc, in full (including internal nodes) ----
  seurat_obj@misc$bonsai <- list(
    phylo = bonsai_tree$phylo,
    edge_info = bonsai_tree$edge_info,
    vert_info = bonsai_tree$vert_info,
    metadata = bonsai_tree$metadata,
    leaf_posterior_ltqs = NULL
  )

  # ---- Subset + reorder posterior LTQs to leaf (cell) vertices only, matching Seurat's cell order ----
  if (!is.null(bonsai_tree$posterior_ltqs) && !is.null(bonsai_tree$vert_info)) {
    vert_names <- bonsai_tree$vert_info$vert_name
    if (is.null(vert_names)) {
      cli::cli_warn("bonsai_tree$vert_info has no 'vert_name' column -- cannot subset posterior LTQs to cells. Skipping.")
    } else {
      leaf_row_idx <- match(seurat_cells, vert_names)
      n_found <- sum(!is.na(leaf_row_idx))
      if (n_found == 0) {
        cli::cli_warn("None of the Seurat object's cell names were found in vert_info$vert_name -- skipping posterior LTQ attachment.")
      } else {
        # posterior_ltqs is vertices x genes (per Bonsai's "vertByGene" naming);
        # transpose to genes x cells to match Seurat's storage convention.
        ltq_subset <- bonsai_tree$posterior_ltqs[leaf_row_idx, , drop = FALSE]
        ltq_subset <- t(ltq_subset)
        colnames(ltq_subset) <- seurat_cells
        seurat_obj@misc$bonsai$leaf_posterior_ltqs <- ltq_subset
        cli::cli_alert_info("Attached posterior LTQs for {n_found}/{length(seurat_cells)} cells (genes x cells matrix in seurat_obj@misc$bonsai$leaf_posterior_ltqs).")
      }
    }
  }

  # ---- Optionally attach cellstates cluster labels, matched by cell ID ----
  if (!is.null(cellstates_output)) {
    if (!inherits(cellstates_output, "cellstates_output")) {
      cli::cli_abort("{.arg cellstates_output} must be a {.cls cellstates_output} object, as returned by {.fn run_cellstates}.")
    }
    if (!fs::file_exists(cellstates_output$cellids_path) || !fs::file_exists(cellstates_output$clusters_path)) {
      cli::cli_warn("cellstates output files not found ({cellstates_output$cellids_path}, {cellstates_output$clusters_path}) -- skipping cellstate attachment.")
    } else {
      cs_cell_ids <- readLines(cellstates_output$cellids_path)
      cs_clusters <- readLines(cellstates_output$clusters_path)
      if (length(cs_cell_ids) != length(cs_clusters)) {
        cli::cli_warn("cellstates CellID.txt ({length(cs_cell_ids)} lines) and optimized_clusters.txt ({length(cs_clusters)} lines) have different lengths -- skipping cellstate attachment, this indicates something is wrong with the cellstates output.")
      } else {
        cluster_by_cell <- stats::setNames(cs_clusters, cs_cell_ids)
        matched <- cluster_by_cell[seurat_cells]
        n_matched <- sum(!is.na(matched))
        seurat_obj$bonsai_cellstate <- matched
        cli::cli_alert_info("Attached cellstate labels for {n_matched}/{length(seurat_cells)} cells as seurat_obj$bonsai_cellstate.")
      }
    }
  }

  cli::cli_alert_success("Bonsai results attached to Seurat object (seurat_obj@misc$bonsai).")

  seurat_obj
}

#' Cut a Bonsai tree into clusters by minimizing within-cluster distance
#'
#' Calls \code{get_min_pdists_clustering_from_nwk_str()} from
#' \code{downstream_analyses/get_clusters_max_diameter.py} in-process via
#' \code{reticulate}, to cut a reconstructed tree into groups of cells that
#' minimize summed pairwise distance within each group -- a purely
#' tree-geometry-based clustering, using no external annotation (unlike
#' \code{\link{bonsai_cluster_by_annotation}}). Roughly analogous to
#' \code{stats::cutree()} on a hierarchical clustering, but respecting the
#' tree's actual branch lengths rather than just topology.
#'
#' \strong{Why this calls a Python function directly instead of a CLI
#' script:} \code{get_clusters_max_diameter.py}'s own \code{if __name__ ==
#' "__main__":} block is leftover development/test code, not a general-
#' purpose CLI -- it hardcodes an absolute path from the original
#' developer's machine, expects a specific test filename
#' (\code{test_binary_tree_8_leafs.nwk}) that won't exist for a real Bonsai
#' output tree, and even calls a function
#' (\code{get_max_diam_clustering_from_nwk_file}) that isn't defined
#' anywhere in the module. Wrapping that CLI would produce a function that
#' always fails. The underlying \code{get_min_pdists_clustering_from_nwk_str()}
#' function it was presumably meant to demonstrate is itself clean, general,
#' and already used this way internally (by
#' \code{bonsai_scout/bonsai_scout_preprocess.py}), so this function calls
#' it the same way, matching this package's stated architecture for
#' lightweight, single-process downstream analyses.
#'
#' \code{get_min_pdists_clustering_from_nwk_str()} internally computes
#' clusterings at \emph{every} cut level from few clusters up to
#' \code{n_clusters} as a byproduct of its greedy-splitting algorithm; this
#' function returns all of them (one column per level) rather than trying
#' to guess which single level you want.
#'
#' @param bonsai_result A \code{bonsai_result} object, as returned by
#'   \code{\link{run_bonsai}}, or a path to a \code{final_bonsai_*}
#'   directory.
#' @param n_clusters Integer. Maximum number of clusters to cut down to;
#'   clusterings at every level up to this are returned (see Details).
#' @param bonsai_env A \code{bonsai_env} object, as returned by
#'   \code{\link{bonsai_install}}.
#' @param cell_ids Optional character vector restricting clustering to
#'   these cell/node IDs. Default \code{NULL} (all leaves).
#' @param verbose Logical, passed through to the underlying Python
#'   function. Default \code{TRUE}.
#'
#' @return A data frame, one row per cell/node ID (row names = IDs), with
#'   one column per clustering level, each giving that cell's cluster
#'   label (e.g. \code{"cl_0"}, \code{"cl_1"}, ...) at that level.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' clusters <- bonsai_cluster_tree(result, n_clusters = 20, bonsai_env = benv)
#' table(clusters[[ncol(clusters)]])  # finest-level clustering's cluster sizes
#' }
bonsai_cluster_tree <- function(bonsai_result,
                                 n_clusters,
                                 bonsai_env,
                                 cell_ids = NULL,
                                 verbose = TRUE) {

  if (!inherits(bonsai_env, "bonsai_env")) {
    cli::cli_abort("{.arg bonsai_env} must be a {.cls bonsai_env} object, as returned by {.fn bonsai_install}.")
  }
  final_dir <- if (inherits(bonsai_result, "bonsai_result")) {
    bonsai_result$final_dir
  } else if (is.character(bonsai_result) && fs::dir_exists(bonsai_result)) {
    bonsai_result
  } else {
    cli::cli_abort("{.arg bonsai_result} must be a {.cls bonsai_result} object or a path to a final_bonsai_* directory.")
  }
  if (is.null(final_dir)) {
    cli::cli_abort("{.arg bonsai_result} has no final tree directory available.")
  }

  nwk_files <- fs::dir_ls(final_dir, glob = "*.nwk")
  if (length(nwk_files) == 0) {
    cli::cli_abort("No .nwk file found in {final_dir}.")
  }
  if (length(nwk_files) > 1) {
    cli::cli_warn("Multiple .nwk files found in {final_dir}; using {nwk_files[[1]]}.")
  }
  nwk_str <- readLines(as.character(nwk_files[[1]]), n = 1, warn = FALSE)

  bonsai_use_conda(bonsai_env)
  reticulate::use_condaenv(bonsai_env$env_name, required = TRUE)
  # import_from_path(), not import(): this module lives in the cloned
  # Bonsai repo, not an installed package, and its own module-level code
  # does absolute imports (e.g. "from bonsai.bonsai_helpers import
  # mp_print") assuming the repo root is on sys.path -- so the repo root,
  # not the downstream_analyses/ subdirectory, is the path to add.
  mod <- reticulate::import_from_path(
    "downstream_analyses.get_clusters_max_diameter",
    path = as.character(bonsai_env$bonsai_repo),
    convert = TRUE
  )

  result <- mod$get_min_pdists_clustering_from_nwk_str(
    tree_nwk_str = nwk_str,
    n_clusters = as.integer(n_clusters),
    cell_ids = if (!is.null(cell_ids)) as.list(cell_ids) else NULL,
    verbose = verbose
  )
  all_clusterings <- result[[1]]

  cl_df <- mod$get_cluster_assignments(all_clusterings = all_clusterings)
  cl_df
}

#' Cut a Bonsai tree into clusters that best match a known annotation
#'
#' Shells out to \code{downstream_analyses/clustering_to_maximize_NMI.py}
#' to cut a reconstructed tree into discrete clusters, choosing cut points
#' that maximize Normalized Mutual Information (NMI) against a given
#' cell-level annotation (e.g. cell type, treatment group, tissue zone).
#' Also computes, for comparison, a purely tree-geometry-based clustering
#' (minimizing within-cluster pairwise distance, ignoring the annotation
#' entirely) and reports both clusterings' NMI against your annotation --
#' a high NMI for the annotation-based clustering just confirms the method
#' worked; a high NMI for the \emph{distance-based} clustering (which never
#' saw your annotation) is the more informative number, since it's a
#' direct, quantitative measure of how well tree structure alone recovers
#' your labels.
#'
#' @param bonsai_result A \code{bonsai_result} object, as returned by
#'   \code{\link{run_bonsai}}, or a path to a \code{final_bonsai_*}
#'   directory.
#' @param annotation A named character or factor vector giving one
#'   annotation label per cell, named by cell ID (matching tree leaf
#'   names/tip labels, e.g. \code{setNames(seurat_obj$Zone,
#'   colnames(seurat_obj))}). Cells with \code{NA} are dropped.
#' @param bonsai_env A \code{bonsai_env} object, as returned by
#'   \code{\link{bonsai_install}}.
#' @param output_dir Character. Directory to write inputs to and read
#'   results from.
#' @param prohibit_small_clusters Logical. If \code{TRUE}, disallows
#'   splitting off many small clusters. Default \code{FALSE}.
#' @param cutting_tol Numeric. How large the NMI improvement must be before
#'   a cluster is split further. Default \code{1e-4}.
#' @param greedy Logical. If \code{FALSE} (default), uses an MCMC-like
#'   optimization scheme (slower, better results); if \code{TRUE}, greedily
#'   cuts subtrees until no move improves the score (faster, somewhat
#'   worse).
#' @param overwrite Logical. If \code{FALSE} (default) and output already
#'   exists in \code{output_dir}, skip re-running.
#'
#' @return An S3 object of class \code{bonsai_nmi_clustering}, a list with
#'   elements:
#'   \describe{
#'     \item{clustering_results}{Data frame, one row per cell (row names =
#'       cell IDs), with one column per clustering method (annotation-based
#'       and distance-based cluster assignments).}
#'     \item{nmi_scores}{Data frame comparing the two methods' NMI against
#'       your annotation, and their runtime.}
#'     \item{output_dir}{Path to the output directory.}
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' zone_by_cell <- setNames(seurat_obj$Zone, colnames(seurat_obj))
#' clustering <- bonsai_cluster_by_annotation(result, zone_by_cell, benv, "nmi_clustering")
#' clustering$nmi_scores
#' }
bonsai_cluster_by_annotation <- function(bonsai_result,
                                          annotation,
                                          bonsai_env,
                                          output_dir,
                                          prohibit_small_clusters = FALSE,
                                          cutting_tol = 1e-4,
                                          greedy = FALSE,
                                          overwrite = FALSE) {

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
  if (is.null(names(annotation))) {
    cli::cli_abort("{.arg annotation} must be a named vector (names = cell IDs matching tree leaf names).")
  }

  nwk_files <- fs::dir_ls(final_dir, glob = "*.nwk")
  if (length(nwk_files) == 0) {
    cli::cli_abort("No .nwk file found in {final_dir}.")
  }
  if (length(nwk_files) > 1) {
    cli::cli_warn("Multiple .nwk files found in {final_dir}; using {nwk_files[[1]]}.")
  }
  nwk_file <- nwk_files[[1]]

  script <- fs::path(bonsai_env$bonsai_repo, "downstream_analyses", "clustering_to_maximize_NMI.py")
  if (!fs::file_exists(script)) {
    cli::cli_abort("Could not find {script} -- check the Bonsai repo layout at {bonsai_env$bonsai_repo}.")
  }

  cell_ids <- names(annotation)
  ann_vals <- as.character(annotation)
  keep <- !is.na(ann_vals)
  if (any(!keep)) {
    cli::cli_alert_info("{sum(!keep)}/{length(ann_vals)} cells have NA {.arg annotation} and will be dropped.")
    cell_ids <- cell_ids[keep]
    ann_vals <- ann_vals[keep]
  }
  if (length(cell_ids) == 0) {
    cli::cli_abort("No cells left in {.arg annotation} after dropping NAs.")
  }

  output_dir <- fs::path_abs(fs::path_expand(output_dir))
  fs::dir_create(output_dir)

  clustering_file <- fs::path(output_dir, "clustering_results.tsv")
  nmi_file <- fs::path(output_dir, "normalized_mutual_information_scores.tsv")

  if (!overwrite && fs::file_exists(clustering_file)) {
    cli::cli_alert_info("Output already exists at {output_dir}, skipping re-run (overwrite = FALSE)")
  } else {
    annotation_file <- fs::path(output_dir, "annotation.tsv")
    writeLines(c("cell_id\tannotation", paste(cell_ids, ann_vals, sep = "\t")), annotation_file)

    cell_names_file <- fs::path(output_dir, "cell_names.txt")
    writeLines(cell_ids, cell_names_file)

    py_bin <- fs::path(reticulate::conda_python(envname = bonsai_env$env_name))
    args <- c(
      as.character(script),
      "--nwk_file", as.character(nwk_file),
      "--annotation_file", as.character(annotation_file),
      "--annotation_id", "annotation",
      "--cell_names_file", as.character(cell_names_file),
      "--results_folder", output_dir,
      "--prohibit_small_clsts", if (prohibit_small_clusters) "True" else "False",
      "--cutting_tol", as.character(cutting_tol),
      "--greedy", if (greedy) "True" else "False",
      "--verbose", "True"
    )

    cli::cli_alert_info("{py_bin} {paste(args, collapse = ' ')}")
    res <- processx::run(as.character(py_bin), args, error_on_status = FALSE, echo = TRUE)

    if (res$status != 0) {
      cli::cli_abort(c(
        "clustering_to_maximize_NMI.py exited with non-zero status ({res$status}).",
        "i" = "Check the output above for details."
      ))
    }
    if (!fs::file_exists(clustering_file)) {
      cli::cli_abort("Script exited successfully but {clustering_file} was not found.")
    }
    cli::cli_alert_success("NMI clustering written to {output_dir}")
  }

  clustering_results <- utils::read.delim(clustering_file, row.names = 1, check.names = FALSE)
  nmi_scores <- utils::read.delim(nmi_file)

  structure(
    list(
      clustering_results = clustering_results,
      nmi_scores = nmi_scores,
      output_dir = output_dir
    ),
    class = "bonsai_nmi_clustering"
  )
}

#' @export
print.bonsai_nmi_clustering <- function(x, ...) {
  cli::cli_h3("bonsai_nmi_clustering")
  cli::cli_bullets(c(
    "*" = "{nrow(x$clustering_results)} cells clustered",
    "*" = "output_dir: {x$output_dir}"
  ))
  print(x$nmi_scores)
  invisible(x)
}

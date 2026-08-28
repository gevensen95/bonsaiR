#' Average Sanity-normalized expression over groups of cells
#'
#' Shells out to \code{downstream_analyses/average_over_groups_wrapper.py}
#' to compute, for each gene and each group, a noise-aware average
#' expression estimate (accounting for each cell's individual Sanity error
#' bar, not just a plain mean) plus an uncertainty on that average, and a
#' per-gene significance score summarizing how much a gene varies across
#' the groups overall. This is the multi-group counterpart to
#' \code{\link{bonsai_marker_genes}}'s pairwise two-group comparison --
#' use this to summarize/rank genes across many groups at once (e.g. your
#' cellstates clusters, or any other cell-level grouping), not to test one
#' specific pair.
#'
#' @param sanity_output A \code{sanity_output} object, as returned by
#'   \code{\link{run_sanity}}.
#' @param group_by A named character or factor vector giving a group label
#'   per cell, named by cell ID (matching Sanity's own cell order, e.g.
#'   \code{setNames(seurat_obj$Zone, colnames(seurat_obj))}). Cells not
#'   present in \code{group_by} (or with \code{NA}) are excluded from
#'   averaging.
#' @param bonsai_env A \code{bonsai_env} object, as returned by
#'   \code{\link{bonsai_install}}.
#' @param output_dir Character. Directory to write the group-assignment
#'   file and read results from.
#' @param overwrite Logical. If \code{FALSE} (default) and output already
#'   exists in \code{output_dir}, skip re-running.
#'
#' @return An S3 object of class \code{average_over_groups_output}, a list
#'   with elements:
#'   \describe{
#'     \item{avg_activities}{Data frame, groups x genes: average expression
#'       per group.}
#'     \item{avg_deltas}{Data frame, groups x genes: uncertainty on each
#'       average (same shape as \code{avg_activities}).}
#'     \item{significance}{Data frame, one row per gene: a summary score of
#'       how much that gene's average varies across groups relative to its
#'       uncertainty (higher = more group-dependent), sorted descending.}
#'     \item{output_dir}{Path to the output directory.}
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' zone_by_cell <- setNames(seurat_obj$Zone, colnames(seurat_obj))
#' result <- run_average_over_groups(sanity_out, zone_by_cell, benv, "avg_by_zone")
#' head(result$significance)
#' }
run_average_over_groups <- function(sanity_output,
                                     group_by,
                                     bonsai_env,
                                     output_dir,
                                     overwrite = FALSE) {

  if (!inherits(sanity_output, "sanity_output")) {
    cli::cli_abort("{.arg sanity_output} must be a {.cls sanity_output} object, as returned by {.fn run_sanity}.")
  }
  if (!inherits(bonsai_env, "bonsai_env")) {
    cli::cli_abort("{.arg bonsai_env} must be a {.cls bonsai_env} object, as returned by {.fn bonsai_install}.")
  }
  if (is.null(names(group_by))) {
    cli::cli_abort("{.arg group_by} must be a named vector (names = cell IDs matching Sanity's cell order).")
  }

  script <- fs::path(bonsai_env$bonsai_repo, "downstream_analyses", "average_over_groups_wrapper.py")
  if (!fs::file_exists(script)) {
    cli::cli_abort("Could not find {script} -- check the Bonsai repo layout at {bonsai_env$bonsai_repo}.")
  }

  cellid_path <- fs::path(sanity_output$output_dir, "cellID.txt")
  delta_path <- fs::path(sanity_output$output_dir, "delta.txt")
  d_delta_path <- fs::path(sanity_output$output_dir, "d_delta.txt")
  geneid_path <- fs::path(sanity_output$output_dir, "geneID.txt")
  missing <- Filter(Negate(fs::file_exists), list(cellid_path, delta_path, d_delta_path, geneid_path))
  if (length(missing) > 0) {
    cli::cli_abort(c(
      "Expected Sanity output file(s) not found: {paste(missing, collapse = ', ')}",
      "i" = "{.arg sanity_output} should come from {.fn run_sanity}."
    ))
  }
  cell_ids <- readLines(cellid_path)

  output_dir <- fs::path_abs(fs::path_expand(output_dir))
  fs::dir_create(output_dir)

  avg_activities_file <- fs::path(output_dir, "avg_activities.csv")
  avg_deltas_file <- fs::path(output_dir, "avg_deltas.csv")
  significance_file <- fs::path(output_dir, "significance.csv")

  if (!overwrite && fs::file_exists(avg_activities_file)) {
    cli::cli_alert_info("Output already exists at {output_dir}, skipping re-run (overwrite = FALSE)")
  } else {
    # The Python script identifies excluded cells by the literal string
    # "NULL" in the group column (its own convention, not ours).
    group_vals <- as.character(group_by[cell_ids])
    n_excluded <- sum(is.na(group_vals))
    if (n_excluded > 0) {
      cli::cli_alert_info("{n_excluded}/{length(cell_ids)} cells not found in {.arg group_by} will be excluded from averaging.")
    }
    group_vals[is.na(group_vals)] <- "NULL"

    group_file <- fs::path(output_dir, "group_assignments.tsv")
    writeLines(paste(cell_ids, group_vals, sep = "\t"), group_file)

    py_bin <- fs::path(reticulate::conda_python(envname = bonsai_env$env_name))
    args <- c(
      as.character(script),
      "-x", as.character(delta_path),
      "-ue", as.character(d_delta_path),
      "-names", as.character(geneid_path),
      "-group_file", as.character(group_file),
      "-col_idx", "1",
      "-output_folder", output_dir
    )

    cli::cli_alert_info("{py_bin} {paste(args, collapse = ' ')}")
    res <- processx::run(as.character(py_bin), args, error_on_status = FALSE, echo = TRUE)

    if (res$status != 0) {
      cli::cli_abort(c(
        "average_over_groups_wrapper.py exited with non-zero status ({res$status}).",
        "i" = "Check the output above for details."
      ))
    }
    if (!fs::file_exists(avg_activities_file)) {
      cli::cli_abort("Script exited successfully but {avg_activities_file} was not found.")
    }
    cli::cli_alert_success("Group averages written to {output_dir}")
  }

  avg_activities <- utils::read.csv(avg_activities_file, row.names = 1, check.names = FALSE)
  avg_deltas <- utils::read.csv(avg_deltas_file, row.names = 1, check.names = FALSE)
  significance <- utils::read.csv(significance_file, row.names = 1, check.names = FALSE)
  names(significance) <- "significance"

  structure(
    list(
      avg_activities = avg_activities,
      avg_deltas = avg_deltas,
      significance = significance,
      output_dir = output_dir
    ),
    class = "average_over_groups_output"
  )
}

#' @export
print.average_over_groups_output <- function(x, ...) {
  cli::cli_h3("average_over_groups_output")
  cli::cli_bullets(c(
    "*" = "{nrow(x$avg_activities)} groups x {ncol(x$avg_activities)} genes",
    "*" = "output_dir: {x$output_dir}"
  ))
  invisible(x)
}

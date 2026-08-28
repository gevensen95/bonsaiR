#' Run cellstates to find statistically indistinguishable cell clusters
#'
#' Shells out to \code{scripts/run_cellstates.py} to cluster cells into
#' groups whose remaining heterogeneity is explainable by measurement
#' noise alone. At tens-of-thousands-of-cells scale, Bonsai's own
#' documentation recommends using these clusters to build a warm-start
#' tree (connecting all cells within a cellstate to a single ancestor)
#' rather than starting Bonsai's tree search from a plain star-tree, since
#' this substantially speeds up convergence and improves the final tree.
#' \code{\link{bonsai_write_config}} builds that warm-start tree
#' automatically from this function's output.
#'
#' cellstates operates on \strong{raw UMI counts} directly (not on Sanity's
#' normalized output), so this function takes the \code{sanity_input}
#' object (which still references the original count matrix) rather than
#' \code{sanity_output}.
#'
#' @param sanity_input A \code{sanity_input} object, as returned by
#'   \code{\link{bonsai_write_sanity_input}}. Used for its raw count
#'   matrix + gene/cell name files -- Sanity's normalization is not
#'   involved here.
#' @param bonsai_env A \code{bonsai_env} object, as returned by
#'   \code{\link{bonsai_install}}.
#' @param output_dir Character. Directory for cellstates' output. Created
#'   if it doesn't exist.
#' @param n_threads Integer. Number of threads. Default
#'   \code{max(1, parallel::detectCores() - 1)}.
#' @param save_intermediates Logical. If \code{TRUE}, cellstates
#'   periodically checkpoints so a long run can be resumed. Recommended at
#'   tens-of-thousands-of-cells scale, where a single run can take a
#'   while. Default \code{TRUE}.
#' @param overwrite Logical. If \code{FALSE} (default) and cellstates
#'   output already appears to exist in \code{output_dir}, skip re-running.
#'
#' @return An S3 object of class \code{cellstates_output}, a list with
#'   elements:
#'   \describe{
#'     \item{output_dir}{Path to cellstates' output directory.}
#'     \item{clusters_path}{Path to \code{optimized_clusters.txt}, giving
#'       each cell's cellstate assignment in the same order as the input
#'       cell list.}
#'     \item{n_cells}{Number of cells clustered.}
#'   }
#'   Intended to be passed to \code{bonsai_write_config()}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' cellstates_out <- run_cellstates(sanity_in, benv, output_dir = "cellstates_output")
#' }
run_cellstates <- function(sanity_input,
                            bonsai_env,
                            output_dir,
                            n_threads = max(1, parallel::detectCores() - 1),
                            save_intermediates = TRUE,
                            overwrite = FALSE) {

  if (!inherits(sanity_input, "sanity_input")) {
    cli::cli_abort("{.arg sanity_input} must be a {.cls sanity_input} object, as returned by {.fn bonsai_write_sanity_input}.")
  }
  if (!inherits(bonsai_env, "bonsai_env")) {
    cli::cli_abort("{.arg bonsai_env} must be a {.cls bonsai_env} object, as returned by {.fn bonsai_install}.")
  }
  if (!fs::file_exists(bonsai_env$cellstates_script)) {
    cli::cli_abort(c(
      "cellstates script not found at {bonsai_env$cellstates_script}.",
      "i" = "Run {.fn bonsai_install} first, and check that its C-extension build step succeeded."
    ))
  }

  # Absolute, not just tilde-expanded -- see bonsai_write_sanity_input.R for
  # why (downstream scripts that os.chdir() would silently mis-resolve a
  # relative path).
  output_dir <- fs::path_abs(fs::path_expand(output_dir))
  fs::dir_create(output_dir)

  clusters_path <- fs::path(output_dir, "optimized_clusters.txt")
  cellids_path <- fs::path(output_dir, "CellID.txt")

  if (!overwrite && fs::file_exists(clusters_path)) {
    cli::cli_alert_info("cellstates output already exists at {output_dir}, skipping re-run (overwrite = FALSE)")
    return(structure(
      list(output_dir = output_dir, clusters_path = clusters_path,
           cellids_path = cellids_path, n_cells = sanity_input$n_cells),
      class = "cellstates_output"
    ))
  }

  py_bin <- fs::path(reticulate::conda_python(envname = bonsai_env$env_name))

  args <- c(
    as.character(bonsai_env$cellstates_script),
    sanity_input$matrix_path,
    "-g", sanity_input$genes_path,
    "-c", sanity_input$cells_path,
    "-o", output_dir,
    "-t", as.character(n_threads)
  )
  if (save_intermediates) args <- c(args, "--save-intermediates")

  cli::cli_alert_info("Running cellstates on {sanity_input$n_genes} genes x {sanity_input$n_cells} cells ({n_threads} threads)")
  cli::cli_alert_info("This step can take a while at this scale -- {.code save_intermediates = TRUE} means it's safe to let it run in the background and check back.")
  cli::cli_alert_info("{py_bin} {paste(args, collapse = ' ')}")

  # bonsai_install() only runs `setup.py build_ext --inplace` for cellstates
  # (compiling its C extension in-place), not a real `pip install`/editable
  # install -- so the `cellstates` package is never placed on sys.path.
  # scripts/run_cellstates.py's own script directory (scripts/) ends up as
  # sys.path[0], which does not include the repo root that actually holds
  # the cellstates/ package, so `from cellstates.cluster import Cluster`
  # fails unless we add the repo root to PYTHONPATH ourselves.
  res <- withr::with_envvar(
    c(PYTHONPATH = paste(bonsai_env$cellstates_repo, Sys.getenv("PYTHONPATH"), sep = .Platform$path.sep)),
    processx::run(as.character(py_bin), args, error_on_status = FALSE, echo = TRUE)
  )

  if (res$status != 0) {
    cli::cli_abort(c(
      "cellstates exited with non-zero status ({res$status}).",
      "i" = "Check the output above for details.",
      "i" = "If this was interrupted partway, you can resume by passing {.code -i} and {.code --dirichlet-file} manually via a follow-up run -- this wrapper doesn't yet automate resumption, only initial checkpointing."
    ))
  }

  if (!fs::file_exists(clusters_path)) {
    cli::cli_abort(c(
      "cellstates finished but the expected output file {clusters_path} was not found.",
      "i" = "Check {output_dir} for what was actually produced -- output file naming may have changed since this package was written."
    ))
  }

  cli::cli_alert_success("cellstates output written to {output_dir}")

  structure(
    list(output_dir = output_dir, clusters_path = clusters_path,
         cellids_path = cellids_path, n_cells = sanity_input$n_cells),
    class = "cellstates_output"
  )
}

#' @export
print.cellstates_output <- function(x, ...) {
  cli::cli_h3("cellstates_output")
  cli::cli_bullets(c(
    "*" = "{x$n_cells} cells clustered",
    "*" = "clusters: {x$clusters_path}",
    "*" = "cell IDs: {x$cellids_path}",
    "*" = "output_dir: {x$output_dir}"
  ))
  invisible(x)
}

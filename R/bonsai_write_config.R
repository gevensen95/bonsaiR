#' Write a Bonsai run configuration, with a cellstates warm-start tree
#'
#' Builds the YAML configuration file that \code{run_bonsai()} needs, by
#' shelling out to Bonsai's own \code{bonsai/create_config_file.py} (rather
#' than hand-writing the YAML in R), so the schema always matches whatever
#' version of Bonsai is installed instead of a schema this package guesses
#' at and which could silently drift out of date.
#'
#' At tens-of-thousands-of-cells scale, Bonsai's documentation recommends
#' starting the tree search from a "warm-start" tree built from cellstates
#' clusters (all cells in a cellstate connected to one ancestor) rather
#' than the default star-tree, since this substantially speeds up
#' convergence. This function builds that warm-start tree by default via
#' \code{optional_preprocessing/create_cellstates_premerged_tree.py}.
#'
#' The flags used below (\code{--config_filepath}, \code{--cellstates_file},
#' \code{--premerged_folder}) were confirmed by reading the script's own
#' \code{argparse} definition directly (an earlier version of this function
#' guessed \code{--cellstates_dir}/\code{--output_dir} by analogy with
#' Bonsai's other scripts, which was wrong and fails with an
#' "unrecognized arguments" error). Note in particular that
#' \code{--cellstates_file} must be a path to cellstates' own
#' \code{optimized_clusters.txt} \strong{file} (one cluster ID per line, one
#' line per cell, in the same order as \code{cellstates_output}'s Cell IDs)
#' -- not the cellstates output directory.
#'
#' @param sanity_output A \code{sanity_output} object, as returned by
#'   \code{\link{run_sanity}}.
#' @param cellstates_output A \code{cellstates_output} object, as returned
#'   by \code{\link{run_cellstates}}, or \code{NULL} to skip the
#'   warm-start and use Bonsai's default star-tree initialization. Default
#'   is to require it (\code{NULL} is allowed but discouraged, and a
#'   warning is issued), since the warm-start is expected to matter a lot
#'   at tens-of-thousands-of-cells scale.
#' @param bonsai_env A \code{bonsai_env} object, as returned by
#'   \code{\link{bonsai_install}}.
#' @param dataset_name Character. An identifier for this dataset/run,
#'   used by Bonsai internally and in output paths.
#' @param results_folder Character. Directory where Bonsai will write its
#'   results. Created if it doesn't exist.
#' @param zscore_cutoff Numeric. Signal-to-noise cutoff for gene
#'   inclusion; Bonsai's own documentation recommends \code{1.0} (the
#'   default). Negative values keep all genes.
#' @param use_knn Integer. Candidate-pair search width for the greedy
#'   merge step; Bonsai recommends 5-20. Default \code{10}.
#' @param nnn_n_randomtrees,nnn_n_randommoves Integers controlling the
#'   thoroughness (and runtime cost) of the Nearest Neighbor Interchange
#'   search phase. Defaults follow Bonsai's own (\code{10}, \code{1000}).
#' @param verbose Logical, passed through to Bonsai. Default \code{TRUE}.
#'
#' @return An S3 object of class \code{bonsai_config}, a list with
#'   elements:
#'   \describe{
#'     \item{yaml_path}{Path to the written config YAML.}
#'     \item{results_folder}{Path Bonsai will write results into.}
#'     \item{premerged_tree_dir}{Path to the cellstates warm-start tree,
#'       or \code{NULL} if \code{cellstates_output} was \code{NULL}.}
#'   }
#'   Intended to be passed to \code{run_bonsai()}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' config <- bonsai_write_config(
#'   sanity_output = sanity_out,
#'   cellstates_output = cellstates_out,
#'   bonsai_env = benv,
#'   dataset_name = "my_dataset",
#'   results_folder = "bonsai_results"
#' )
#' }
bonsai_write_config <- function(sanity_output,
                                 cellstates_output = NULL,
                                 bonsai_env,
                                 dataset_name,
                                 results_folder,
                                 zscore_cutoff = 1.0,
                                 use_knn = 10,
                                 nnn_n_randomtrees = 10,
                                 nnn_n_randommoves = 1000,
                                 verbose = TRUE) {

  if (!inherits(sanity_output, "sanity_output")) {
    cli::cli_abort("{.arg sanity_output} must be a {.cls sanity_output} object, as returned by {.fn run_sanity}.")
  }
  if (!is.null(cellstates_output) && !inherits(cellstates_output, "cellstates_output")) {
    cli::cli_abort("{.arg cellstates_output} must be a {.cls cellstates_output} object (or NULL), as returned by {.fn run_cellstates}.")
  }
  if (is.null(cellstates_output)) {
    cli::cli_warn(c(
      "No {.arg cellstates_output} provided -- Bonsai will start from a plain star-tree.",
      "i" = "At tens-of-thousands-of-cells scale, Bonsai's own docs recommend the cellstates warm-start for reasonable runtime. Consider running {.fn run_cellstates} first."
    ))
  }

  # Absolute, not just tilde-expanded -- this is the critical one. Both
  # create_config_file.py and create_cellstates_premerged_tree.py (below)
  # os.chdir() to the Bonsai repo root before touching any path they're
  # given, so a relative results_folder gets silently written to e.g.
  # ~/bonsai_tools/Bonsai-data-representation/<relative_path>/... instead
  # of where the caller meant, while the Python script still reports
  # success -- which is exactly what produced the confusing "exited with
  # status 0" abort below (the real failure was our own file_exists()
  # check looking in the right place for the wrong reason: the file
  # existed, just not where we were looking).
  results_folder <- fs::path_abs(fs::path_expand(results_folder))
  fs::dir_create(results_folder)

  yaml_path <- fs::path(results_folder, paste0(dataset_name, "_config.yaml"))
  premerged_tree_dir <- if (!is.null(cellstates_output)) {
    fs::path(results_folder, "premerged_tree")
  } else {
    NULL
  }

  bonsai_use_conda(bonsai_env)
  py_bin <- fs::path(reticulate::conda_python(envname = bonsai_env$env_name))
  create_config_script <- fs::path(bonsai_env$bonsai_repo, "bonsai", "create_config_file.py")

  if (!fs::file_exists(create_config_script)) {
    cli::cli_abort("Could not find {create_config_script} -- check the Bonsai repo layout at {bonsai_env$bonsai_repo}.")
  }

  args <- c(
    as.character(create_config_script),
    "--new_yaml_path", yaml_path,
    "--dataset", dataset_name,
    "--data_folder", sanity_output$output_dir,
    "--verbose", if (verbose) "True" else "False",
    "--results_folder", results_folder,
    "--input_is_sanity_output", "True",
    "--zscore_cutoff", as.character(zscore_cutoff),
    "--nnn_n_randomtrees", as.character(nnn_n_randomtrees),
    "--nnn_n_randommoves", as.character(nnn_n_randommoves),
    "--use_knn", as.character(use_knn)
  )
  if (!is.null(premerged_tree_dir)) {
    args <- c(args, "--tmp_folder", premerged_tree_dir)
  }

  cli::cli_alert_info("Writing Bonsai config to {yaml_path}")
  res <- processx::run(as.character(py_bin), args, error_on_status = FALSE, echo = TRUE)
  if (res$status != 0 || !fs::file_exists(yaml_path)) {
    cli::cli_abort("Failed to write Bonsai config (create_config_file.py exited with status {res$status}).")
  }

  # ---- Build the cellstates warm-start tree, if requested ----
  if (!is.null(cellstates_output)) {
    premerge_script <- fs::path(bonsai_env$bonsai_repo, "optional_preprocessing", "create_cellstates_premerged_tree.py")
    if (!fs::file_exists(premerge_script)) {
      cli::cli_abort("Could not find {premerge_script} -- check the Bonsai repo layout at {bonsai_env$bonsai_repo}.")
    }

    fs::dir_create(premerged_tree_dir)
    premerge_args <- c(
      as.character(premerge_script),
      "--config_filepath", yaml_path,
      "--cellstates_file", cellstates_output$clusters_path,
      "--premerged_folder", premerged_tree_dir
    )

    cli::cli_alert_info("{py_bin} {paste(premerge_args, collapse = ' ')}")
    res2 <- processx::run(as.character(py_bin), premerge_args, error_on_status = FALSE, echo = TRUE)
    if (res2$status != 0) {
      cli::cli_abort(c(
        "create_cellstates_premerged_tree.py exited with non-zero status ({res2$status}).",
        "i" = "Check the output above for details."
      ))
    }
  }

  cli::cli_alert_success("Bonsai config ready at {yaml_path}")

  structure(
    list(
      yaml_path = yaml_path,
      results_folder = results_folder,
      premerged_tree_dir = premerged_tree_dir
    ),
    class = "bonsai_config"
  )
}

#' @export
print.bonsai_config <- function(x, ...) {
  cli::cli_h3("bonsai_config")
  cli::cli_bullets(c(
    "*" = "config: {x$yaml_path}",
    "*" = "results_folder: {x$results_folder}",
    "*" = if (!is.null(x$premerged_tree_dir)) "warm-start tree: {x$premerged_tree_dir}" else "warm-start: none (star-tree init)"
  ))
  invisible(x)
}

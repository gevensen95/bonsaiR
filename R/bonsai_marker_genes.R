#' Identify marker genes distinguishing two clades of a Bonsai tree
#'
#' Shells out to \code{downstream_analyses/calc_marker_genes.py} to
#' calculate, for every gene, the probability that a randomly picked cell
#' from group 1 has higher expression than a randomly picked cell from
#' group 2 (Bonsai's marker score \eqn{M}; see the Methods section of the
#' Bonsai paper). Scores near 1 or 0 indicate strong markers for group 1
#' or group 2 respectively.
#'
#' \strong{Preprocessing step:} \code{calc_marker_genes.py} does not read
#' Bonsai's raw tree output directly -- it reads a \code{bonsai_vis_data.hdf}
#' / \code{bonsai_vis_settings.json} pair that only
#' \code{bonsai_scout/bonsai_scout_preprocess.py} produces (this is the same
#' preprocessing step the Bonsai-scout Shiny app depends on). This function
#' runs that preprocessing step automatically the first time (skipped on
#' later calls once the files exist, unless \code{overwrite_preprocess =
#' TRUE}).
#'
#' \strong{JSON schema and flags (confirmed against the scripts' own source,
#' not guessed):} \code{calc_marker_genes.py} takes \code{--marker_groups_json}
#' (not \code{--groups_json}), \code{--results_folder} (not
#' \code{--tree_folder}/\code{--output_dir}), and optionally
#' \code{--marker_output_file}. The JSON file it expects contains
#' \code{cell_ids_group1}, \code{cell_ids_group2}, and \code{feature_path} --
#' the last being a path inside \code{bonsai_vis_data.hdf} to the feature
#' matrix to compare (default here: \code{"data/normalized"}, the tree's
#' posterior log-transcription-quotients with error bars, which is what
#' Bonsai's own Shiny app compares by default). An earlier version of this
#' function used a different, guessed schema
#' (\code{list(group1 = ..., group2 = ...)}); that was wrong.
#'
#' @param bonsai_tree A \code{bonsai_tree} object, as returned by
#'   \code{\link{bonsai_read_tree}}, or \code{NULL} if you already have
#'   cell ID vectors for the two groups and don't need the tree for
#'   anything else.
#' @param group1_cells,group2_cells Character vectors of cell IDs (must
#'   match the vertex/cell names Bonsai used, i.e. the Seurat object's
#'   original cell barcodes) defining the two clades to compare.
#' @param bonsai_env A \code{bonsai_env} object, as returned by
#'   \code{\link{bonsai_install}}.
#' @param bonsai_result A \code{bonsai_result} object (used both for its
#'   \code{results_folder}, where the visualization preprocessing step
#'   writes its output, and its \code{final_dir}) or a path to a
#'   \code{final_bonsai_*} directory (in which case its parent directory is
#'   assumed to be the results folder -- pass a \code{bonsai_result} object
#'   instead if that assumption doesn't hold for your layout).
#' @param output_dir Character. Directory to write the group-membership
#'   JSON and marker-score output into.
#' @param feature_path Character. Path inside \code{bonsai_vis_data.hdf} to
#'   compare cells on. Default \code{"data/normalized"} (posterior LTQs with
#'   error bars) -- only override this if you know the HDF5 layout
#'   \code{bonsai_scout_preprocess.py} produced for your run (e.g. it also
#'   writes per-annotation matrices under \code{data/<annotation_name>} when
#'   applicable).
#' @param overwrite_preprocess Logical. If \code{FALSE} (default) and
#'   \code{bonsai_vis_data.hdf}/\code{bonsai_vis_settings.json} already
#'   exist in the results folder, skip re-running
#'   \code{bonsai_scout_preprocess.py}. Set \code{TRUE} to force a re-run
#'   (e.g. after the tree changed).
#'
#' @return A data frame with columns \code{marker_genes} and
#'   \code{marker_scores}, sorted by distance of \code{marker_scores} from
#'   0.5 (strongest markers first). \code{NULL} with a warning if the
#'   expected output file isn't found.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' markers <- bonsai_marker_genes(
#'   group1_cells = clade_a_cells,
#'   group2_cells = clade_b_cells,
#'   bonsai_env = benv,
#'   bonsai_result = result,
#'   output_dir = "marker_genes"
#' )
#' }
bonsai_marker_genes <- function(bonsai_tree = NULL,
                                 group1_cells,
                                 group2_cells,
                                 bonsai_env,
                                 bonsai_result,
                                 output_dir,
                                 feature_path = "data/normalized",
                                 overwrite_preprocess = FALSE) {

  if (!inherits(bonsai_env, "bonsai_env")) {
    cli::cli_abort("{.arg bonsai_env} must be a {.cls bonsai_env} object, as returned by {.fn bonsai_install}.")
  }

  if (inherits(bonsai_result, "bonsai_result")) {
    results_folder <- bonsai_result$results_folder
  } else if (is.character(bonsai_result) && fs::dir_exists(bonsai_result)) {
    results_folder <- fs::path_dir(bonsai_result)
    cli::cli_warn("{.arg bonsai_result} was given as a raw path -- assuming its parent directory ({results_folder}) is the Bonsai results folder. Pass a {.cls bonsai_result} object instead if that's wrong.")
  } else {
    cli::cli_abort("{.arg bonsai_result} must be a {.cls bonsai_result} object or a path to a final_bonsai_* directory.")
  }
  # Absolute, not just tilde-expanded -- both bonsai_scout_preprocess.py and
  # calc_marker_genes.py (below) os.chdir() to the Bonsai repo root, which
  # silently breaks a relative path (see bonsai_write_config.R for the full
  # story). Defensive here even though bonsai_write_config() now always
  # returns an absolute results_folder, since results_folder can also come
  # from a hand-constructed bonsai_result or a raw path above.
  results_folder <- fs::path_abs(fs::path_expand(results_folder))
  if (length(group1_cells) == 0 || length(group2_cells) == 0) {
    cli::cli_abort("Both {.arg group1_cells} and {.arg group2_cells} must be non-empty.")
  }
  overlap <- intersect(group1_cells, group2_cells)
  if (length(overlap) > 0) {
    cli::cli_abort("{.arg group1_cells} and {.arg group2_cells} overlap ({length(overlap)} shared IDs) -- they must be disjoint clades.")
  }

  bonsai_use_conda(bonsai_env)
  py_bin <- fs::path(reticulate::conda_python(envname = bonsai_env$env_name))

  # ---- Step 1: build bonsai_vis_data.hdf / bonsai_vis_settings.json if needed ----
  vis_data_path <- fs::path(results_folder, "bonsai_vis_data.hdf")
  vis_settings_path <- fs::path(results_folder, "bonsai_vis_settings.json")

  if (overwrite_preprocess || !fs::file_exists(vis_data_path) || !fs::file_exists(vis_settings_path)) {
    preprocess_script <- fs::path(bonsai_env$bonsai_repo, "bonsai_scout", "bonsai_scout_preprocess.py")
    if (!fs::file_exists(preprocess_script)) {
      cli::cli_abort("Could not find {preprocess_script} -- check the Bonsai repo layout at {bonsai_env$bonsai_repo}.")
    }
    cli::cli_alert_info("Building visualization data (bonsai_scout_preprocess.py) -- needed by calc_marker_genes.py, not produced by run_bonsai() itself.")
    res_pre <- processx::run(
      as.character(py_bin),
      c(as.character(preprocess_script), "--results_folder", results_folder),
      error_on_status = FALSE, echo = TRUE
    )
    if (res_pre$status != 0 || !fs::file_exists(vis_data_path)) {
      cli::cli_abort(c(
        "bonsai_scout_preprocess.py failed to produce {vis_data_path}.",
        "i" = "Check the output above for details."
      ))
    }
  } else {
    cli::cli_alert_info("Visualization data already exists at {results_folder}, skipping bonsai_scout_preprocess.py (overwrite_preprocess = FALSE).")
  }

  # ---- Step 2: resolve requested cell IDs against the HDF5's actual node IDs ----
  # bonsai_scout_preprocess.py merges cells at (near-)zero branch-length
  # distance into single composite vertices for visualization purposes
  # (imports merge_cells_at_zero_dist) -- their node_id becomes an
  # underscore-joined concatenation of the original cell IDs (e.g.
  # "Cell1_Cell2_Cell3"), and the individual cell IDs no longer appear in
  # node_ids at all. This is expected/correct behavior (confirmed by
  # inspecting the HDF5 directly), not a bug -- it mainly shows up on very
  # clean/degenerate data (e.g. this package's own synthetic smoke-test
  # data) where many cells truly are statistically indistinguishable; real
  # datasets rarely hit exact zero distance. We resolve each requested cell
  # ID to whatever node_ids entry actually contains it (itself, or a
  # composite it was merged into) so group1_cells/group2_cells specified as
  # plain cell barcodes keep working after such a merge.
  reticulate::use_condaenv(bonsai_env$env_name, required = TRUE)
  h5py <- reticulate::import("h5py", convert = TRUE)
  hdf <- h5py$File(as.character(vis_data_path), "r")
  feature_grp <- tryCatch(
    hdf[[feature_path]],
    error = function(e) {
      hdf$close()
      cli::cli_abort("{.arg feature_path} ({.val {feature_path}}) not found in {vis_data_path}.")
    }
  )
  node_ids <- jsonlite::fromJSON(feature_grp$attrs[["node_ids"]])
  hdf$close()

  resolve_group <- function(cell_ids, group_label) {
    resolved <- character(0)
    for (cid in unique(cell_ids)) {
      if (cid %in% node_ids) {
        resolved <- c(resolved, cid)
        next
      }
      hit <- node_ids[vapply(strsplit(node_ids, "_", fixed = TRUE), function(parts) cid %in% parts, logical(1))]
      if (length(hit) == 1) {
        resolved <- c(resolved, hit)
      } else if (length(hit) == 0) {
        cli::cli_abort("Cell ID {.val {cid}} (in {group_label}) was not found among the tree's visualization node IDs at all -- check spelling/casing against bonsai_tree$phylo$tip.label.")
      } else {
        cli::cli_abort("Cell ID {.val {cid}} (in {group_label}) matched multiple composite vertices ({paste(hit, collapse = ', ')}) -- the HDF5 node_ids may be malformed.")
      }
    }
    unique(resolved)
  }

  resolved1 <- resolve_group(group1_cells, "group1_cells")
  resolved2 <- resolve_group(group2_cells, "group2_cells")

  merge_conflict <- intersect(resolved1, resolved2)
  if (length(merge_conflict) > 0) {
    cli::cli_abort(c(
      "Some requested cells from group1_cells and group2_cells were merged into the same visualization vertex by bonsai_scout_preprocess.py, so they can't be compared.",
      "x" = "Conflicting merged vertex/vertices: {paste(merge_conflict, collapse = ', ')}",
      "i" = "This happens when those cells have (near-)zero branch length distance in the reconstructed tree -- pick group boundaries that don't split a zero-distance cluster, or use a less degenerate dataset."
    ))
  }

  # ---- Step 3: write the group-membership JSON with the confirmed schema ----
  script <- fs::path(bonsai_env$bonsai_repo, "downstream_analyses", "calc_marker_genes.py")
  if (!fs::file_exists(script)) {
    cli::cli_abort("Could not find {script} -- check the Bonsai repo layout at {bonsai_env$bonsai_repo}.")
  }

  output_dir <- fs::path_abs(fs::path_expand(output_dir))
  fs::dir_create(output_dir)
  groups_json_path <- fs::path(output_dir, "marker_gene_groups.json")

  jsonlite::write_json(
    list(
      cell_ids_group1 = resolved1,
      cell_ids_group2 = resolved2,
      feature_path = jsonlite::unbox(feature_path)
    ),
    groups_json_path, auto_unbox = FALSE
  )

  marker_output_file <- fs::path(output_dir, "marker_scores.tsv")
  args <- c(
    as.character(script),
    "--marker_groups_json", groups_json_path,
    "--results_folder", results_folder,
    "--marker_output_file", marker_output_file
  )

  cli::cli_alert_info("{py_bin} {paste(args, collapse = ' ')}")
  res <- processx::run(as.character(py_bin), args, error_on_status = FALSE, echo = TRUE)

  if (res$status != 0) {
    cli::cli_abort(c(
      "calc_marker_genes.py exited with non-zero status ({res$status}).",
      "i" = "Check the output above for details."
    ))
  }

  if (!fs::file_exists(marker_output_file)) {
    cli::cli_warn(c(
      "calc_marker_genes.py exited successfully but {marker_output_file} was not found.",
      "i" = "Contents of {output_dir}: {paste(fs::path_file(fs::dir_ls(output_dir)), collapse = ', ')}"
    ))
    return(NULL)
  }

  # calc_marker_genes.py writes via pandas' to_csv() with its default
  # (unnamed) row-index column first, then marker_genes, then
  # marker_scores -- read the first column as row names rather than data.
  result <- utils::read.delim(marker_output_file, row.names = 1)
  result <- result[order(abs(result$marker_scores - 0.5), decreasing = TRUE), , drop = FALSE]
  rownames(result) <- NULL

  result
}

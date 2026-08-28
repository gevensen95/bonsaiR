#' Run an end-to-end smoke test of the full bonsaiR pipeline
#'
#' Generates a small synthetic Seurat object and runs it through every
#' stage of the pipeline (\code{bonsai_write_sanity_input} ->
#' \code{run_sanity} -> \code{run_cellstates} -> \code{bonsai_write_config}
#' -> \code{run_bonsai} -> \code{bonsai_read_tree} ->
#' \code{bonsai_marker_genes} -> \code{bonsai_to_seurat}), recording
#' pass/fail and timing for each stage.
#'
#' This deliberately does \strong{not} depend on Bonsai's own example
#' data, since this package's author could not verify the exact path,
#' format, or size of any example dataset shipped in the Bonsai repo. The
#' synthetic dataset generated here is small enough to run in a couple of
#' minutes on a laptop, and has genuine group structure (not pure noise),
#' so cellstates and Bonsai have something real to recover -- this is
#' useful for catching wiring bugs (wrong CLI flags, wrong file paths,
#' wrong output parsing) but is \strong{not} a validation of Bonsai's
#' statistical performance, and is far too small to be representative of
#' behavior at your actual tens-of-thousands-of-cells scale.
#'
#' If a stage fails, the chain stops there (each stage needs the previous
#' stage's output) but the function still returns, with remaining stages
#' marked \code{"skipped"} in the summary -- so you always get a full
#' report of exactly where things broke.
#'
#' @param bonsai_env A \code{bonsai_env} object, as returned by
#'   \code{\link{bonsai_install}}. If \code{NULL} (default), this function
#'   calls \code{bonsai_install()} itself, which can be slow the first
#'   time (conda env creation, cloning, compiling) -- pass an existing one
#'   to skip that.
#' @param work_dir Character. Directory for all intermediate files.
#'   Default a fresh \code{tempfile()} directory. Printed at the start of
#'   the run so you know where to look if something fails.
#' @param n_groups Integer. Number of synthetic cell groups. Default
#'   \code{3}.
#' @param n_cells_per_group Integer. Cells per group. Default \code{20}
#'   (total cells = \code{n_groups * n_cells_per_group}).
#' @param n_genes Integer. Total genes, split evenly into group-marker
#'   genes plus shared background genes. Default \code{150}.
#' @param use_mpi Logical. Whether \code{run_bonsai()} should use MPI.
#'   Default \code{TRUE}, to exercise the same code path you'll actually
#'   use at scale. Set \code{FALSE} if you're troubleshooting and want to
#'   rule out MPI-specific issues.
#' @param n_cores Integer. MPI processes for the smoke test's
#'   \code{run_bonsai()} call. Default \code{2} -- deliberately small
#'   since this dataset is tiny.
#' @param cleanup Logical. If \code{TRUE}, delete \code{work_dir} after
#'   the test completes (pass or fail). Default \code{FALSE} -- smoke
#'   tests exist to leave you something to inspect.
#' @param seed Integer. Random seed for synthetic data generation, for
#'   reproducibility. Default \code{1}.
#'
#' @return Invisibly, a data frame with one row per pipeline stage:
#'   \code{stage}, \code{status} (\code{"pass"}, \code{"fail"}, or
#'   \code{"skipped"}), \code{elapsed_sec}, and \code{message} (error
#'   detail, if failed). Also printed as a formatted summary at the end of
#'   the run.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Using an existing environment (recommended -- skips slow reinstall):
#' benv <- bonsai_install()
#' bonsai_smoke_test(bonsai_env = benv)
#'
#' # Or let it install fresh (slow):
#' bonsai_smoke_test()
#' }
bonsai_smoke_test <- function(bonsai_env = NULL,
                               work_dir = fs::path(tempdir(), paste0("bonsai_smoke_", format(Sys.time(), "%Y%m%d_%H%M%S"))),
                               n_groups = 3,
                               n_cells_per_group = 20,
                               n_genes = 150,
                               use_mpi = TRUE,
                               n_cores = 2,
                               cleanup = FALSE,
                               seed = 1) {

  # Absolute, not just tilde-expanded -- see bonsai_write_config.R for why
  # (downstream Bonsai scripts that os.chdir() would silently mis-resolve
  # a relative path). Harmless for the default tempfile()-based work_dir,
  # which is already absolute, but matters if a relative one is passed in.
  work_dir <- fs::path_abs(fs::path_expand(work_dir))
  fs::dir_create(work_dir)

  cli::cli_h1("bonsaiR smoke test")
  cli::cli_alert_info("Working directory: {.path {work_dir}}")
  cli::cli_alert_warning("Uses a synthetic dataset (not Bonsai's own example data -- its exact path/format was not verified). This checks pipeline wiring, not Bonsai's statistical performance.")

  results <- list()
  record <- function(stage, status, elapsed, message = NA_character_) {
    results[[stage]] <<- data.frame(stage = stage, status = status,
                                     elapsed_sec = round(elapsed, 2),
                                     message = message,
                                     stringsAsFactors = FALSE)
  }

  # A small state bag threaded through stages.
  state <- new.env()
  ok_so_far <- TRUE

  run_stage <- function(stage_name, fn) {
    if (!ok_so_far) {
      record(stage_name, "skipped", 0, "upstream stage failed")
      cli::cli_alert_info("[{stage_name}] skipped (upstream failure)")
      return(invisible())
    }
    cli::cli_h2("Stage: {stage_name}")
    t0 <- Sys.time()
    out <- tryCatch(
      { fn(); list(ok = TRUE, err = NA_character_) },
      error = function(e) list(ok = FALSE, err = conditionMessage(e))
    )
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (out$ok) {
      record(stage_name, "pass", elapsed)
      cli::cli_alert_success("[{stage_name}] passed ({round(elapsed, 1)}s)")
    } else {
      record(stage_name, "fail", elapsed, out$err)
      cli::cli_alert_danger("[{stage_name}] FAILED ({round(elapsed, 1)}s): {out$err}")
      ok_so_far <<- FALSE
    }
  }

  # ---- Stage 0: install (only if not supplied) ----
  run_stage("bonsai_install", function() {
    if (is.null(bonsai_env)) {
      cli::cli_alert_info("No bonsai_env supplied -- running bonsai_install() (this can be slow on first run).")
      state$bonsai_env <- bonsai_install()
    } else {
      if (!inherits(bonsai_env, "bonsai_env")) {
        stop("bonsai_env argument was supplied but is not a bonsai_env object.")
      }
      state$bonsai_env <- bonsai_env
    }
  })

  # ---- Stage 1: synthetic data generation ----
  run_stage("generate_synthetic_data", function() {
    set.seed(seed)
    n_cells <- n_groups * n_cells_per_group
    genes_per_group <- floor(n_genes / (n_groups + 1))  # + 1 block of shared/background genes
    stopifnot(genes_per_group > 0)

    gene_names <- paste0("Gene", seq_len(n_genes))
    cell_names <- paste0("Cell", seq_len(n_cells))
    group_id <- rep(seq_len(n_groups), each = n_cells_per_group)

    counts <- matrix(0L, nrow = n_genes, ncol = n_cells,
                      dimnames = list(gene_names, cell_names))

    background_genes <- seq_len(n_genes) > genes_per_group * n_groups

    for (g in seq_len(n_groups)) {
      cells_in_group <- which(group_id == g)
      marker_genes <- ((g - 1) * genes_per_group + 1):(g * genes_per_group)
      for (gene in seq_len(n_genes)) {
        lambda <- if (gene %in% marker_genes) 25 else if (background_genes[gene]) 5 else 1
        counts[gene, cells_in_group] <- stats::rpois(length(cells_in_group), lambda = lambda)
      }
    }

    state$true_group <- stats::setNames(paste0("group", group_id), cell_names)
    state$seurat_obj <- Seurat::CreateSeuratObject(
      counts = Matrix::Matrix(counts, sparse = TRUE),
      project = "bonsai_smoke_test",
      min.cells = 0, min.features = 0
    )
    cli::cli_alert_info("Generated {n_cells} cells x {n_genes} genes across {n_groups} synthetic groups.")
  })

  # ---- Stage 2: bonsai_write_sanity_input ----
  run_stage("bonsai_write_sanity_input", function() {
    state$sanity_input <- bonsai_write_sanity_input(
      state$seurat_obj,
      output_dir = fs::path(work_dir, "sanity_input"),
      overwrite = TRUE
    )
  })

  # ---- Stage 3: run_sanity ----
  run_stage("run_sanity", function() {
    state$sanity_output <- run_sanity(
      state$sanity_input, state$bonsai_env,
      output_dir = fs::path(work_dir, "sanity_output"),
      n_threads = 2, overwrite = TRUE
    )
  })

  # ---- Stage 4: run_cellstates ----
  run_stage("run_cellstates", function() {
    state$cellstates_output <- run_cellstates(
      state$sanity_input, state$bonsai_env,
      output_dir = fs::path(work_dir, "cellstates_output"),
      n_threads = 2, overwrite = TRUE
    )
  })

  # ---- Stage 5: bonsai_write_config ----
  run_stage("bonsai_write_config", function() {
    state$bonsai_config <- bonsai_write_config(
      state$sanity_output,
      cellstates_output = state$cellstates_output,
      bonsai_env = state$bonsai_env,
      dataset_name = "smoke_test",
      results_folder = fs::path(work_dir, "bonsai_results")
    )
  })

  # ---- Stage 6: run_bonsai ----
  run_stage("run_bonsai", function() {
    state$bonsai_result <- run_bonsai(
      state$bonsai_config, state$bonsai_env,
      n_cores = n_cores, use_mpi = use_mpi
    )
    if (is.null(state$bonsai_result$final_dir)) {
      stop("run_bonsai() completed without error but produced no final_dir -- treating as a failure for smoke-test purposes.")
    }
  })

  # ---- Stage 7: bonsai_read_tree ----
  run_stage("bonsai_read_tree", function() {
    state$bonsai_tree <- bonsai_read_tree(state$bonsai_result, state$bonsai_env)
    n_leaves <- ape::Ntip(state$bonsai_tree$phylo)
    if (n_leaves != n_groups * n_cells_per_group) {
      stop(sprintf("Tree has %d leaves, expected %d -- cell count mismatch somewhere in the pipeline.",
                    n_leaves, n_groups * n_cells_per_group))
    }
  })

  # ---- Stage 8: bonsai_marker_genes (known-shakiest function; failure here is informative, not fatal to the rest) ----
  run_stage("bonsai_marker_genes", function() {
    g1_cells <- names(state$true_group)[state$true_group == "group1"]
    g2_cells <- names(state$true_group)[state$true_group == "group2"]
    state$markers <- bonsai_marker_genes(
      bonsai_tree = state$bonsai_tree,
      group1_cells = g1_cells,
      group2_cells = g2_cells,
      bonsai_env = state$bonsai_env,
      bonsai_result = state$bonsai_result,
      output_dir = fs::path(work_dir, "marker_genes")
    )
    if (is.null(state$markers)) {
      stop("bonsai_marker_genes() returned NULL -- see warnings above; this function's I/O contract is the least-verified part of this package (see ?bonsai_marker_genes).")
    }
  })

  # ---- Stage 9: bonsai_to_seurat ----
  # NB: does not depend on stage 8 succeeding, so run it even if marker genes failed,
  # as long as we still have a tree.
  if (!is.null(state$bonsai_tree)) {
    ok_so_far_saved <- ok_so_far
    ok_so_far <- TRUE  # temporarily allow this stage even if stage 8 failed
    run_stage("bonsai_to_seurat", function() {
      state$seurat_obj <- bonsai_to_seurat(
        state$seurat_obj, state$bonsai_tree,
        cellstates_output = state$cellstates_output
      )
      if (is.null(state$seurat_obj@misc$bonsai$phylo)) {
        stop("bonsai_to_seurat() ran without error but did not attach a phylo object.")
      }
    })
    ok_so_far <- ok_so_far_saved  # restore for the final summary's own logic (unused after this point, but kept for clarity)
  } else {
    record("bonsai_to_seurat", "skipped", 0, "no tree available")
  }

  # ---- Summary ----
  summary_df <- do.call(rbind, results)
  rownames(summary_df) <- NULL

  cli::cli_h1("Smoke test summary")
  print(summary_df)

  n_pass <- sum(summary_df$status == "pass")
  n_fail <- sum(summary_df$status == "fail")
  n_skip <- sum(summary_df$status == "skipped")
  cli::cli_alert_info("{n_pass} passed, {n_fail} failed, {n_skip} skipped.")
  if (n_fail > 0) {
    cli::cli_alert_danger("First failure: {summary_df$stage[summary_df$status == 'fail'][1]}")
  } else {
    cli::cli_alert_success("All stages passed.")
  }
  cli::cli_alert_info("Artifacts left in {.path {work_dir}}{if (cleanup) ' (will be deleted now, cleanup = TRUE)' else ''}")

  if (cleanup) {
    fs::dir_delete(work_dir)
  }

  invisible(summary_df)
}

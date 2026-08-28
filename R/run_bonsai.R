#' Run the core Bonsai tree-search
#'
#' Shells out to \code{bonsai_main.py} to reconstruct the maximum
#' likelihood tree. This is the step that cannot be run through an
#' embedded \code{reticulate} interpreter (see package overview) --
#' Bonsai's own documentation describes it as designed around
#' \code{mpiexec}-launched multi-process parallelism, which is used here
#' by default given your dataset scale (tens of thousands of cells). A
#' single-core fallback is available via \code{use_mpi = FALSE}, matching
#' what Bonsai's README recommends only for smaller datasets (under
#' ~2000 cells).
#'
#' Following Bonsai's documented usage, this function invokes
#' \code{bonsai_main.py} with the working directory set to the Bonsai
#' repo root (as their README's examples do -- \code{cd
#' Bonsai-data-representation && python3 bonsai/bonsai_main.py ...}),
#' rather than assuming it is safe to run from an arbitrary directory.
#'
#' \strong{Known intermittent MPI failure:} with \code{use_mpi = TRUE},
#' Bonsai's own \code{nnnReorderRandom} step (part of its core tree search)
#' has been observed to occasionally crash with a \code{FileNotFoundError}
#' reading back one of its own \code{random_trees/random_tree_N} scratch
#' folders -- confirmed (by inspecting a failed run's leftover files) to be
#' caused by that entire scratch folder having already been deleted by
#' Bonsai's own "cleaning up intermediate datafiles" step from what appears
#' to be a separate, later round of its multi-round core calculation,
#' racing against an MPI rank still trying to read it. This is a timing-
#' dependent issue inside Bonsai's own MPI orchestration code, not
#' something this wrapper can safely patch (task-splitting across ranks
#' itself was verified working correctly). It has so far only been observed
#' on small datasets (tens of cells) -- which is also below the ~2000-cell
#' threshold Bonsai's own README recommends MPI for in the first place, so
#' affected runs are likely also the ones that gain the least from MPI. If
#' you hit this, either re-run (the race is not fully deterministic) or set
#' \code{use_mpi = FALSE}.
#'
#' @param bonsai_config A \code{bonsai_config} object, as returned by
#'   \code{\link{bonsai_write_config}}.
#' @param bonsai_env A \code{bonsai_env} object, as returned by
#'   \code{\link{bonsai_install}}.
#' @param n_cores Integer. Number of MPI processes to launch. Default
#'   \code{max(1, parallel::detectCores() - 1)}. Ignored if
#'   \code{use_mpi = FALSE}.
#' @param use_mpi Logical. Default \code{TRUE}. Set \code{FALSE} to run
#'   single-core (Bonsai's own recommendation only for small datasets --
#'   at tens of thousands of cells this will likely be very slow).
#' @param step Character. One of \code{"all"}, \code{"preprocess"},
#'   \code{"core_calc"}, \code{"metadata"}. Default \code{"all"}. Splitting
#'   into steps is mainly useful on an HPC scheduler with different
#'   resource requirements per step (see Bonsai's README); on a single
#'   workstation \code{"all"} is normally what you want.
#' @param pickup_intermediate Logical. If \code{TRUE}, resume from the
#'   latest stored intermediate result rather than starting over --
#'   useful if a previous run was interrupted. Default \code{FALSE}.
#'
#' @return An S3 object of class \code{bonsai_result}, a list with
#'   elements:
#'   \describe{
#'     \item{results_folder}{The results folder passed in via
#'       \code{bonsai_config}.}
#'     \item{final_dir}{Path to the \code{final_bonsai_*} subdirectory
#'       containing the completed tree, or \code{NULL} if \code{step} was
#'       not \code{"all"} or the run did not complete.}
#'   }
#'   Intended to be passed to \code{bonsai_read_tree()}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' result <- run_bonsai(config, benv, n_cores = 8)
#' }
run_bonsai <- function(bonsai_config,
                        bonsai_env,
                        n_cores = max(1, parallel::detectCores() - 1),
                        use_mpi = TRUE,
                        step = c("all", "preprocess", "core_calc", "metadata"),
                        pickup_intermediate = FALSE) {

  step <- match.arg(step)

  if (!inherits(bonsai_config, "bonsai_config")) {
    cli::cli_abort("{.arg bonsai_config} must be a {.cls bonsai_config} object, as returned by {.fn bonsai_write_config}.")
  }
  if (!inherits(bonsai_env, "bonsai_env")) {
    cli::cli_abort("{.arg bonsai_env} must be a {.cls bonsai_env} object, as returned by {.fn bonsai_install}.")
  }

  bonsai_main <- fs::path("bonsai", "bonsai_main.py")  # relative, per README usage
  bonsai_main_abs <- fs::path(bonsai_env$bonsai_repo, bonsai_main)
  if (!fs::file_exists(bonsai_main_abs)) {
    cli::cli_abort("Could not find {bonsai_main_abs} -- check the Bonsai repo layout at {bonsai_env$bonsai_repo}.")
  }

  py_bin <- as.character(fs::path(reticulate::conda_python(envname = bonsai_env$env_name)))

  base_args <- c(
    as.character(bonsai_main),
    "--config_filepath", bonsai_config$yaml_path,
    "--step", step
  )
  if (pickup_intermediate) base_args <- c(base_args, "--pickup_intermediate", "True")

  if (use_mpi) {
    mpiexec_bin <- Sys.which("mpiexec")
    if (mpiexec_bin == "") {
      cli::cli_abort(c(
        "{.code mpiexec} not found on PATH.",
        "i" = "This should have been installed alongside system OpenMPI in {.fn bonsai_install}. Check your PATH, or set {.code use_mpi = FALSE} to run single-core (not recommended at your dataset scale)."
      ))
    }
    cmd <- as.character(mpiexec_bin)
    args <- c("-n", as.character(n_cores), py_bin, "-m", "mpi4py", base_args)
    cli::cli_alert_info("Running Bonsai with MPI across {n_cores} processes")
  } else {
    cli::cli_warn("Running Bonsai single-core (use_mpi = FALSE). At tens-of-thousands-of-cells scale this may be very slow.")
    cmd <- py_bin
    args <- base_args
  }

  cli::cli_alert_info("(cwd: {bonsai_env$bonsai_repo}) {cmd} {paste(args, collapse = ' ')}")
  cli::cli_alert_info("This can take a long time on tens of thousands of cells -- output is streamed below.")

  # NB: when the underlying process is killed by MPI_ABORT (Open MPI prints
  # its own diagnostic message to stdout/stderr when it force-kills all
  # ranks after one crashes), that message can contain an embedded NUL
  # byte. processx's internal output-capture throws its own low-level R/C
  # error ("embedded nul in string") when it hits that byte, which
  # previously propagated up raw and completely obscured the actual Python
  # traceback (which *did* print live above via echo = TRUE, just not
  # captured into a clean R error) -- catch that specific processx failure
  # and re-point the user at the real cause instead of this R-internal
  # string-handling limitation.
  res <- tryCatch(
    processx::run(
      cmd, args,
      wd = as.character(bonsai_env$bonsai_repo),
      error_on_status = FALSE, echo = TRUE,
      # Bonsai's own core_calc step can run long; do not impose a timeout here.
      # NB: processx::run()'s "no timeout" value is Inf, not 0 -- timeout = 0
      # means "must finish within 0 seconds", which gets the process SIGKILLed
      # almost immediately (this was a real bug here: every run_bonsai() call
      # died in under a second with no captured output, which looked like a
      # crash but was actually this misconfigured timeout).
      timeout = Inf
    ),
    error = function(e) {
      if (grepl("embedded nul", conditionMessage(e), fixed = TRUE)) {
        cli::cli_abort(c(
          "The underlying Bonsai process was killed (most likely by MPI_ABORT after an internal crash on one rank), and its abort message contained a byte that R's processx can't capture into a normal error.",
          "i" = "The actual Python traceback should be visible above in the streamed output -- scroll up to find it.",
          "i" = if (use_mpi) "This specific failure (a race in Bonsai's own multi-round core-calculation step reusing shared temp folders across MPI ranks) has been observed intermittently on small datasets, which is also below the ~2000-cell threshold Bonsai's own README recommends for using MPI at all. Consider re-running (it's timing-dependent, not fully deterministic), or set use_mpi = FALSE for small/local runs where reliability matters more than the modest speedup." else NULL
        ))
      }
      cli::cli_abort(c("processx::run() failed.", "x" = conditionMessage(e)))
    }
  )

  if (res$status != 0) {
    cli::cli_abort(c(
      "Bonsai exited with non-zero status ({res$status}).",
      "i" = "Check the output above for details.",
      "i" = "If this was a resource/time-limit issue, re-run with {.code pickup_intermediate = TRUE} to resume."
    ))
  }

  final_dir <- NULL
  if (step == "all") {
    candidates <- fs::dir_ls(bonsai_config$results_folder, type = "directory", regexp = "final_bonsai_")
    if (length(candidates) == 0) {
      cli::cli_warn(c(
        "Bonsai exited successfully but no {.code final_bonsai_*} directory was found in {bonsai_config$results_folder}.",
        "i" = "Check the run output above -- the tree may not have fully converged, or the results folder naming may have changed."
      ))
    } else {
      final_dir <- as.character(candidates[[1]])
      if (length(candidates) > 1) {
        cli::cli_warn("Multiple {.code final_bonsai_*} directories found; using {final_dir}. Clean up {bonsai_config$results_folder} if this is unexpected.")
      }
      cli::cli_alert_success("Bonsai tree reconstruction complete: {final_dir}")
    }
  } else {
    cli::cli_alert_success("Bonsai step '{step}' complete.")
  }

  structure(
    list(results_folder = bonsai_config$results_folder, final_dir = final_dir),
    class = "bonsai_result"
  )
}

#' @export
print.bonsai_result <- function(x, ...) {
  cli::cli_h3("bonsai_result")
  cli::cli_bullets(c(
    "*" = "results_folder: {x$results_folder}",
    "*" = if (!is.null(x$final_dir)) "final tree: {x$final_dir}" else "final tree: not yet available"
  ))
  invisible(x)
}

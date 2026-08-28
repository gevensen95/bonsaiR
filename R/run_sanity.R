#' Run Sanity to normalize raw UMI counts with error bars
#'
#' Shells out to the compiled Sanity binary (see \code{\link{bonsai_install}})
#' to normalize a raw UMI count matrix, producing log-transcription-quotient
#' estimates with per-gene, per-cell error bars that Bonsai's likelihood
#' model requires.
#'
#' \strong{Flag choice:} Bonsai's documentation recommends running Sanity
#' with \code{-e 1 -v_m MAP} (the "Sanity 2.0+" interface) -- this function
#' hard-codes exactly that. (An earlier version of this function assumed the
#' \code{jmbreda/Sanity} build lacked \code{-v_m} and used the legacy
#' \code{-max_v true} flag instead; that was based on a stale reading of the
#' CLI help and was wrong -- the current build does have \code{-v_m}, and
#' using it produces \code{delta.txt}/\code{d_delta.txt} plus a
#' \code{sanity_command.txt} marker file, which is exactly what Bonsai's own
#' \code{bonsai_helpers.py} expects for its "2.0+" input path, as opposed to
#' its legacy \code{delta_vmax.txt}/\code{d_delta_vmax.txt} fallback for
#' pre-2.0 Sanity.)
#'
#' These flags are not exposed as user-tunable parameters, since Bonsai
#' requires them specifically -- see the Bonsai README's discussion of why
#' the MAP/single-guess gene-variance output (rather than marginalizing over
#' it, i.e. \code{MARG}) is required for Bonsai's likelihood to be
#' reconstructable.
#'
#' @param sanity_input A \code{sanity_input} object, as returned by
#'   \code{\link{bonsai_write_sanity_input}}.
#' @param bonsai_env A \code{bonsai_env} object, as returned by
#'   \code{\link{bonsai_install}}.
#' @param output_dir Character. Directory for Sanity's output. Created if
#'   it doesn't exist.
#' @param n_threads Integer. Number of threads to use. Default
#'   \code{max(1, parallel::detectCores() - 1)}.
#' @param overwrite Logical. If \code{FALSE} (default) and Sanity output
#'   already appears to exist in \code{output_dir}, skip re-running and
#'   just return the existing output's path (useful for resuming a
#'   pipeline without repeating an expensive step). Set \code{TRUE} to
#'   force a re-run.
#'
#' @return An S3 object of class \code{sanity_output}, a list with
#'   elements:
#'   \describe{
#'     \item{output_dir}{Path to Sanity's output directory.}
#'     \item{n_genes}{Number of genes in the input.}
#'     \item{n_cells}{Number of cells in the input.}
#'   }
#'   Intended to be passed to \code{run_cellstates()} and \code{run_bonsai()}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' sanity_in <- bonsai_write_sanity_input(seu, output_dir = "sanity_input")
#' sanity_out <- run_sanity(sanity_in, benv, output_dir = "sanity_output")
#' }
run_sanity <- function(sanity_input,
                        bonsai_env,
                        output_dir,
                        n_threads = max(1, parallel::detectCores() - 1),
                        overwrite = FALSE) {

  if (!inherits(sanity_input, "sanity_input")) {
    cli::cli_abort("{.arg sanity_input} must be a {.cls sanity_input} object, as returned by {.fn bonsai_write_sanity_input}.")
  }
  if (!inherits(bonsai_env, "bonsai_env")) {
    cli::cli_abort("{.arg bonsai_env} must be a {.cls bonsai_env} object, as returned by {.fn bonsai_install}.")
  }
  if (!fs::file_exists(bonsai_env$sanity_bin)) {
    cli::cli_abort(c(
      "Sanity binary not found at {bonsai_env$sanity_bin}.",
      "i" = "Run {.fn bonsai_install} first."
    ))
  }

  # Absolute, not just tilde-expanded -- see bonsai_write_sanity_input.R for
  # why (downstream scripts that os.chdir() would silently mis-resolve a
  # relative path). Sanity's own binary doesn't chdir, but this path also
  # gets passed on as create_config_file.py's --data_folder later, which does.
  output_dir <- fs::path_abs(fs::path_expand(output_dir))
  fs::dir_create(output_dir)

  # These are the files Sanity's extended output (-e 1) + -v_m MAP produce;
  # delta.txt / d_delta.txt are the ones Bonsai's "2.0+" input path (see
  # bonsai_helpers.py's SANITY_FILENAMES) actually consumes.
  key_output_file <- fs::path(output_dir, "delta.txt")

  if (!overwrite && fs::file_exists(key_output_file)) {
    cli::cli_alert_info("Sanity output already exists at {output_dir}, skipping re-run (overwrite = FALSE)")
    return(structure(
      list(output_dir = output_dir,
           n_genes = sanity_input$n_genes,
           n_cells = sanity_input$n_cells),
      class = "sanity_output"
    ))
  }

  args <- c(
    "-f", sanity_input$matrix_path,
    "-mtx_genes", sanity_input$genes_path,
    "-mtx_cells", sanity_input$cells_path,
    "-d", output_dir,
    "-n", as.character(n_threads),
    "-e", "1",
    "-v_m", "MAP"
  )

  cli::cli_alert_info("Running Sanity on {sanity_input$n_genes} genes x {sanity_input$n_cells} cells ({n_threads} threads)")
  cli::cli_alert_info("{bonsai_env$sanity_bin} {paste(args, collapse = ' ')}")

  res <- processx::run(
    as.character(bonsai_env$sanity_bin), args,
    error_on_status = FALSE, echo = TRUE
  )

  if (res$status != 0) {
    cli::cli_abort(c(
      "Sanity exited with non-zero status ({res$status}).",
      "i" = "Check the output above for details."
    ))
  }

  if (!fs::file_exists(key_output_file)) {
    cli::cli_abort(c(
      "Sanity finished but the expected output file {key_output_file} was not found.",
      "i" = "This may indicate Sanity's output file naming has changed since this package was written -- check {output_dir} for what was actually produced."
    ))
  }

  cli::cli_alert_success("Sanity output written to {output_dir}")

  structure(
    list(
      output_dir = output_dir,
      n_genes = sanity_input$n_genes,
      n_cells = sanity_input$n_cells
    ),
    class = "sanity_output"
  )
}

#' @export
print.sanity_output <- function(x, ...) {
  cli::cli_h3("sanity_output")
  cli::cli_bullets(c(
    "*" = "{x$n_genes} genes x {x$n_cells} cells",
    "*" = "output_dir: {x$output_dir}"
  ))
  invisible(x)
}

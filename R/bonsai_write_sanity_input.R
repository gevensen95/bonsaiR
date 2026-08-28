#' Write a Seurat object's raw counts as Sanity-compatible input files
#'
#' Sanity accepts two input formats: a dense tab/comma/space-separated text
#' matrix, or Matrix Market (\code{.mtx}) sparse format with optional
#' sidecar gene-name and cell-name files (the same layout cellranger writes:
#' \code{features.tsv}/\code{genes.tsv} and \code{barcodes.tsv}). At
#' tens-of-thousands-of-cells scale, writing a dense text matrix would be
#' large and slow, and Seurat already stores counts sparsely, so this
#' function always writes the \code{.mtx} format.
#'
#' Sanity expects **raw UMI counts**, not normalized or log-transformed
#' data. This function checks that the extracted matrix looks like integer
#' count data (non-negative, values close to whole numbers) and errors with
#' an explanatory message if not, since pulling the wrong Seurat layer/slot
#' (e.g. \code{"data"} instead of \code{"counts"}) is an easy mistake that
#' would otherwise fail silently deep inside Sanity.
#'
#' @param seurat_obj A \code{Seurat} object.
#' @param assay Character. Assay to pull counts from. Default \code{"RNA"}.
#' @param layer Character. Layer (Seurat v5) or slot (Seurat v4) to pull
#'   counts from. Default \code{"counts"}. Only override this if you know
#'   what you're doing -- Sanity requires raw UMI counts.
#' @param output_dir Character. Directory to write the Sanity input files
#'   into. Created if it doesn't exist.
#' @param overwrite Logical. If \code{FALSE} (default) and output files
#'   already exist, errors rather than silently overwriting them.
#'
#' @return An S3 object of class \code{sanity_input}, a list with elements:
#'   \describe{
#'     \item{matrix_path}{Path to the written \code{.mtx} file.}
#'     \item{genes_path}{Path to the gene-name sidecar file.}
#'     \item{cells_path}{Path to the cell-barcode sidecar file.}
#'     \item{n_genes}{Number of genes (rows) written.}
#'     \item{n_cells}{Number of cells (columns) written.}
#'   }
#'   This object is intended to be passed to \code{run_sanity()}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' sanity_in <- bonsai_write_sanity_input(seu, output_dir = "sanity_input")
#' }
bonsai_write_sanity_input <- function(seurat_obj,
                                       assay = "RNA",
                                       layer = "counts",
                                       output_dir,
                                       overwrite = FALSE) {

  if (!inherits(seurat_obj, "Seurat")) {
    cli::cli_abort("{.arg seurat_obj} must be a Seurat object, got class {.cls {class(seurat_obj)}}.")
  }
  if (!assay %in% Seurat::Assays(seurat_obj)) {
    cli::cli_abort(c(
      "Assay {.val {assay}} not found in this Seurat object.",
      "i" = "Available assays: {.val {Seurat::Assays(seurat_obj)}}"
    ))
  }

  # ---- Extract counts, handling both Seurat v5 (layer=) and v4 (slot=) ----
  counts <- tryCatch(
    Seurat::GetAssayData(seurat_obj, assay = assay, layer = layer),
    error = function(e) {
      tryCatch(
        Seurat::GetAssayData(seurat_obj, assay = assay, slot = layer),
        error = function(e2) {
          cli::cli_abort(c(
            "Could not extract {.val {layer}} from assay {.val {assay}}.",
            "x" = conditionMessage(e2)
          ))
        }
      )
    }
  )

  if (nrow(counts) == 0 || ncol(counts) == 0) {
    cli::cli_abort("Extracted count matrix is empty ({nrow(counts)} genes x {ncol(counts)} cells).")
  }

  # ---- Validate this actually looks like raw UMI counts ----
  # Sample rather than checking the whole matrix -- fine for tens of
  # thousands of cells x genes, this is just a sanity check (pun intended),
  # not exhaustive validation.
  sample_vals <- counts@x[seq_len(min(length(counts@x), 100000))]
  if (length(sample_vals) > 0) {
    if (any(sample_vals < 0)) {
      cli::cli_abort(c(
        "The extracted data from assay {.val {assay}}, layer/slot {.val {layer}} contains negative values.",
        "x" = "Sanity requires raw, non-negative UMI counts.",
        "i" = "You may have pulled a scaled/z-scored layer by mistake -- check {.code assay}/{.code layer}."
      ))
    }
    non_integer_frac <- mean(abs(sample_vals - round(sample_vals)) > 1e-6)
    if (non_integer_frac > 0.01) {
      cli::cli_abort(c(
        "The extracted data from assay {.val {assay}}, layer/slot {.val {layer}} does not look like integer UMI counts ({round(non_integer_frac * 100, 1)}% of sampled values are non-integer).",
        "x" = "Sanity requires raw UMI counts, not normalized or log-transformed data.",
        "i" = "Did you mean {.code layer = \"counts\"} instead of {.val {layer}}?"
      ))
    }
  } else {
    cli::cli_warn("Extracted count matrix has no non-zero entries to validate -- proceeding, but double-check {.arg assay}/{.arg layer}.")
  }

  # ---- Coerce to sparse dgCMatrix (genes x cells -- already the orientation Sanity wants) ----
  counts <- methods::as(counts, "CsparseMatrix")

  # fs::path_abs() (not just path_expand()) matters here: several Bonsai
  # scripts we shell out to downstream do os.chdir() to the Bonsai repo
  # root before touching any path they're given, which silently breaks a
  # relative path (it gets resolved against the wrong directory instead of
  # erroring). Making every output_dir/results_folder absolute at the
  # point it's first accepted means this can't happen no matter what any
  # downstream script's cwd ends up being.
  output_dir <- fs::path_abs(fs::path_expand(output_dir))
  fs::dir_create(output_dir)

  matrix_path <- fs::path(output_dir, "matrix.mtx")
  genes_path  <- fs::path(output_dir, "mtx_genes.tsv")
  cells_path  <- fs::path(output_dir, "mtx_cells.tsv")

  if (!overwrite) {
    existing <- Filter(fs::file_exists, list(matrix_path, genes_path, cells_path))
    if (length(existing) > 0) {
      cli::cli_abort(c(
        "Output file(s) already exist and {.arg overwrite} is FALSE:",
        stats::setNames(as.character(existing), rep("x", length(existing)))
      ))
    }
  }

  cli::cli_alert_info("Writing {nrow(counts)} genes x {ncol(counts)} cells to {output_dir}")

  # Matrix::writeMM() always emits triplets in column-major order (verified:
  # this holds even when starting from an RsparseMatrix, since writeMM
  # coerces to CsparseMatrix internally), i.e. row indices are NOT
  # non-decreasing across the file. The Matrix Market spec doesn't require
  # any particular order, but the Sanity binary's own .mtx reader does
  # require rows to be sorted (it errors out otherwise and points users at
  # its bundled sort_mtx_by_row.py) -- so we write the file ourselves,
  # sorted by row, rather than relying on writeMM().
  tmat <- methods::as(counts, "TsparseMatrix")
  ord <- order(tmat@i, tmat@j)
  i_sorted <- tmat@i[ord] + 1L
  j_sorted <- tmat@j[ord] + 1L
  x_sorted <- tmat@x[ord]
  is_int <- all(abs(x_sorted - round(x_sorted)) < 1e-8)

  con <- file(matrix_path, "wt")
  on.exit(close(con), add = TRUE)
  writeLines(c(
    paste("%%MatrixMarket matrix coordinate", if (is_int) "integer" else "real", "general"),
    paste(nrow(counts), ncol(counts), length(x_sorted))
  ), con)
  utils::write.table(
    data.frame(i_sorted, j_sorted, if (is_int) as.integer(round(x_sorted)) else x_sorted),
    con, row.names = FALSE, col.names = FALSE, sep = " "
  )

  writeLines(rownames(counts), genes_path)
  writeLines(colnames(counts), cells_path)

  cli::cli_alert_success("Sanity input written: {matrix_path}")

  structure(
    list(
      matrix_path = matrix_path,
      genes_path = genes_path,
      cells_path = cells_path,
      n_genes = nrow(counts),
      n_cells = ncol(counts)
    ),
    class = "sanity_input"
  )
}

#' @export
print.sanity_input <- function(x, ...) {
  cli::cli_h3("sanity_input")
  cli::cli_bullets(c(
    "*" = "{x$n_genes} genes x {x$n_cells} cells",
    "*" = "matrix: {x$matrix_path}",
    "*" = "genes: {x$genes_path}",
    "*" = "cells: {x$cells_path}"
  ))
  invisible(x)
}

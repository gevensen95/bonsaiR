#' Read a completed Bonsai tree reconstruction into R
#'
#' Reads the files Bonsai's README documents as its final output --
#' \code{.nwk} (Newick tree), \code{edgeInfo.txt}, \code{vertInfo.txt},
#' \code{metadata.json}, and the posterior log-transcription-quotient
#' arrays -- into native R objects. The \code{.npy} arrays are loaded via
#' \code{reticulate} + numpy (no dedicated R \code{.npy} reader dependency
#' is added, since \code{reticulate} is already a dependency and this is
#' exactly the single-process, function-call-style use case it's suited
#' for).
#'
#' @param bonsai_result A \code{bonsai_result} object, as returned by
#'   \code{\link{run_bonsai}}, with a non-\code{NULL} \code{final_dir}. You
#'   can also pass a path to a \code{final_bonsai_*} directory directly.
#' @param bonsai_env A \code{bonsai_env} object, as returned by
#'   \code{\link{bonsai_install}}. Needed to locate a Python/numpy
#'   interpreter for reading the \code{.npy} files.
#' @param load_posteriors Logical. If \code{TRUE} (default), also load the
#'   (potentially large) posterior LTQ arrays for internal nodes. Set
#'   \code{FALSE} to skip this and just get the tree topology + metadata,
#'   which is much faster for a first look at a large tree.
#'
#' @return An S3 object of class \code{bonsai_tree}, a list with elements:
#'   \describe{
#'     \item{phylo}{An \code{ape::phylo} object -- the tree topology and
#'       branch lengths.}
#'     \item{edge_info}{Data frame: \code{vert_ind_1}, \code{vert_ind_2},
#'       \code{edge_length}, one row per edge.}
#'     \item{vert_info}{Data frame: \code{vert_ind}, \code{vert_name}, one
#'       row per vertex (leaves are named by cell ID; internal node names
#'       are Bonsai-assigned).}
#'     \item{metadata}{Parsed \code{metadata.json} as an R list.}
#'     \item{posterior_ltqs}{Matrix of posterior LTQ point estimates
#'       (vertices x genes), or \code{NULL} if \code{load_posteriors =
#'       FALSE}.}
#'     \item{posterior_ltqs_vars}{Matrix of posterior LTQ variances
#'       (vertices x genes), or \code{NULL} if \code{load_posteriors =
#'       FALSE}.}
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' tree <- bonsai_read_tree(result, benv)
#' ape::plot.phylo(tree$phylo, show.tip.label = FALSE)
#' }
bonsai_read_tree <- function(bonsai_result, bonsai_env, load_posteriors = TRUE) {

  final_dir <- if (inherits(bonsai_result, "bonsai_result")) {
    if (is.null(bonsai_result$final_dir)) {
      cli::cli_abort("{.arg bonsai_result} has no {.field final_dir} -- the Bonsai run may not have completed with {.code step = \"all\"}.")
    }
    bonsai_result$final_dir
  } else if (is.character(bonsai_result) && fs::dir_exists(bonsai_result)) {
    bonsai_result
  } else {
    cli::cli_abort("{.arg bonsai_result} must be a {.cls bonsai_result} object or a path to a final_bonsai_* directory.")
  }

  # ---- Newick tree ----
  nwk_files <- fs::dir_ls(final_dir, glob = "*.nwk")
  if (length(nwk_files) == 0) {
    cli::cli_abort("No .nwk file found in {final_dir}.")
  }
  if (length(nwk_files) > 1) {
    cli::cli_warn("Multiple .nwk files found in {final_dir}; using {nwk_files[[1]]}.")
  }
  phylo <- ape::read.tree(as.character(nwk_files[[1]]))

  # ---- edgeInfo.txt / vertInfo.txt ----
  edge_info_path <- fs::path(final_dir, "edgeInfo.txt")
  vert_info_path <- fs::path(final_dir, "vertInfo.txt")

  edge_info <- if (fs::file_exists(edge_info_path)) {
    ei <- utils::read.delim(edge_info_path, header = FALSE)
    if (ncol(ei) >= 3) names(ei)[1:3] <- c("vert_ind_1", "vert_ind_2", "edge_length")
    ei
  } else {
    cli::cli_warn("edgeInfo.txt not found at {edge_info_path}.")
    NULL
  }

  vert_info <- if (fs::file_exists(vert_info_path)) {
    vi <- utils::read.delim(vert_info_path, header = TRUE)
    # vertInfo.txt's real header is camelCase (vertInd, nodeInd, vertName),
    # confirmed by inspecting an actual run's output -- normalize to the
    # snake_case names documented above and relied on by bonsai_to_seurat(),
    # which was otherwise silently finding no 'vert_name' column and
    # skipping posterior LTQ attachment entirely.
    rename_map <- c(vertInd = "vert_ind", nodeInd = "node_ind", vertName = "vert_name")
    matched <- intersect(names(rename_map), names(vi))
    names(vi)[match(matched, names(vi))] <- rename_map[matched]
    vi
  } else {
    cli::cli_warn("vertInfo.txt not found at {vert_info_path}.")
    NULL
  }

  # ---- metadata.json ----
  metadata_path <- fs::path(final_dir, "metadata.json")
  metadata <- if (fs::file_exists(metadata_path)) {
    jsonlite::fromJSON(metadata_path)
  } else {
    cli::cli_warn("metadata.json not found at {metadata_path}.")
    NULL
  }

  # ---- posterior LTQ arrays (.npy, via reticulate + numpy) ----
  posterior_ltqs <- NULL
  posterior_ltqs_vars <- NULL

  if (load_posteriors) {
    ltq_path <- fs::path(final_dir, "posterior_ltqs_vertByGene.npy")
    var_path <- fs::path(final_dir, "posterior_ltqsVars_vertByGene.npy")

    if (fs::file_exists(ltq_path) && fs::file_exists(var_path)) {
      reticulate::use_condaenv(bonsai_env$env_name, required = TRUE)
      np <- reticulate::import("numpy", convert = TRUE)
      posterior_ltqs <- np$load(as.character(ltq_path))
      posterior_ltqs_vars <- np$load(as.character(var_path))
      cli::cli_alert_info("Loaded posterior LTQ arrays: {nrow(posterior_ltqs)} vertices x {ncol(posterior_ltqs)} genes")
    } else {
      cli::cli_warn("Posterior LTQ .npy files not found in {final_dir} -- skipping (set load_posteriors = FALSE to silence this).")
    }
  }

  structure(
    list(
      phylo = phylo,
      edge_info = edge_info,
      vert_info = vert_info,
      metadata = metadata,
      posterior_ltqs = posterior_ltqs,
      posterior_ltqs_vars = posterior_ltqs_vars
    ),
    class = "bonsai_tree"
  )
}

#' @export
print.bonsai_tree <- function(x, ...) {
  cli::cli_h3("bonsai_tree")
  cli::cli_bullets(c(
    "*" = "{ape::Ntip(x$phylo)} leaves, {ape::Nnode(x$phylo)} internal nodes",
    "*" = if (!is.null(x$posterior_ltqs)) "posterior LTQs loaded: {nrow(x$posterior_ltqs)} x {ncol(x$posterior_ltqs)}" else "posterior LTQs: not loaded"
  ))
  invisible(x)
}

# Force reticulate to resolve conda from a known-good binary path rather
# than rediscovering it fresh (RETICULATE_CONDA env var, then a handful of
# common install locations, then PATH). That rediscovery isn't cached
# across R sessions, so a bonsai_env built in one session (or cached via
# saveRDS()/readRDS(), the pattern this package's own vignettes use) can
# silently stop resolving conda in a later session whose PATH differs --
# e.g. a cluster job launched without the interactive shell's module/PATH
# setup (verified against a real HPC failure: reticulate::conda_python()
# failing with "Unable to find conda binary" using a benv from an earlier,
# working session). bonsai_install() records the conda binary it actually
# used in bonsai_env$conda_bin for exactly this reason.
bonsai_use_conda <- function(bonsai_env) {
  if (!inherits(bonsai_env, "bonsai_env")) {
    cli::cli_abort("{.arg bonsai_env} must be a {.cls bonsai_env} object, as returned by {.fn bonsai_install}.")
  }
  if (!is.null(bonsai_env$conda_bin) && nzchar(bonsai_env$conda_bin)) {
    Sys.setenv(RETICULATE_CONDA = bonsai_env$conda_bin)
  }
  invisible(bonsai_env)
}

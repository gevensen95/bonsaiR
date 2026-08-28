#' Install the Bonsai tool stack (Sanity, cellstates, Bonsai)
#'
#' One-time setup that creates a dedicated conda environment containing
#' Python, mpi4py/OpenMPI bindings, and Bonsai's Python dependencies, then
#' clones (if needed) and builds the Sanity C++ binary and the cellstates
#' tool. This function is idempotent: steps that already appear to be done
#' are skipped unless \code{force = TRUE}.
#'
#' The core Bonsai tree-search step is designed to run as an external
#' multi-process MPI job (\code{mpiexec -n N python -m mpi4py
#' bonsai_main.py ...}). This cannot be achieved through an embedded
#' \code{reticulate} interpreter, so this package shells out to that CLI
#' rather than wrapping it as an in-process Python call. \code{reticulate}
#' is used elsewhere in this package only for single-process, function-call
#' style tasks (e.g. marker gene detection).
#'
#' @param env_name Character. Name of the conda environment to create.
#'   Default \code{"bonsai"}.
#' @param repo_path Character. Local directory into which the Bonsai,
#'   Sanity, and cellstates source repositories will be cloned. Default
#'   \code{"~/bonsai_tools"}.
#' @param python_version Character. Python version for the conda
#'   environment. Default \code{"3.11"} (Bonsai's own README currently
#'   suggests 3.14, which is unlikely to have mature package support at
#'   the time of writing; adjust if you hit dependency resolution issues).
#' @param force Logical. If \code{TRUE}, re-create the conda environment
#'   and re-clone/re-build all tools even if they already appear present.
#'   Default \code{FALSE}.
#' @param n_threads Integer. Number of threads to use when building the
#'   Sanity binary. Default \code{max(1, parallel::detectCores() - 1)}.
#'
#' @return An S3 object of class \code{bonsai_env}, a list with elements:
#'   \describe{
#'     \item{env_name}{The conda environment name.}
#'     \item{repo_path}{Path to the cloned source repositories.}
#'     \item{bonsai_repo}{Path to the Bonsai-data-representation repo.}
#'     \item{sanity_repo}{Path to the Sanity repo.}
#'     \item{sanity_bin}{Path to the built Sanity executable
#'       (\code{Sanity/bin/Sanity}).}
#'     \item{cellstates_repo}{Path to the cellstates repo.}
#'     \item{cellstates_script}{Path to cellstates' CLI entry point
#'       (\code{scripts/run_cellstates.py}).}
#'   }
#'   This object is intended to be passed as the \code{bonsai_env} argument
#'   to the other functions in this package (\code{run_sanity()},
#'   \code{run_cellstates()}, \code{run_bonsai()}, etc.) so that paths do
#'   not need to be re-specified on every call.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' benv <- bonsai_install()
#' benv <- bonsai_install(env_name = "bonsai", repo_path = "~/bonsai_tools")
#' }
bonsai_install <- function(env_name = "bonsai",
                            repo_path = "~/bonsai_tools",
                            python_version = "3.11",
                            force = FALSE,
                            n_threads = max(1, parallel::detectCores() - 1)) {

  # Absolute, not just tilde-expanded -- see bonsai_write_config.R for why.
  # Harmless for the default "~/bonsai_tools", which is already absolute
  # once expanded, but matters if a relative repo_path is passed in.
  repo_path <- fs::path_abs(fs::path_expand(repo_path))
  fs::dir_create(repo_path)

  # ---- 1. Check system OpenMPI is available before touching conda ----
  # mpi4py installation fails in confusing ways if there is no system
  # OpenMPI/mpicc to link against, so we fail fast with a clear message.
  mpicc_path <- Sys.which("mpicc")
  if (mpicc_path == "") {
    cli::cli_abort(c(
      "System OpenMPI does not appear to be installed (no {.code mpicc} on PATH).",
      "i" = "Install OpenMPI first, e.g. {.code apt-get install openmpi-bin libopenmpi-dev} on Debian/Ubuntu, or the equivalent for your system, then re-run {.fn bonsai_install}.",
      "i" = "Bonsai will technically run without mpi4py, but only on a single core, which defeats the purpose at tens-of-thousands-of-cells scale."
    ))
  }
  cli::cli_alert_success("Found system OpenMPI at {mpicc_path}")

  # ---- 2. Create (or reuse) the conda environment ----
  existing_envs <- tryCatch(reticulate::conda_list()$name, error = function(e) character(0))

  if (force || !(env_name %in% existing_envs)) {
    if (env_name %in% existing_envs) {
      cli::cli_alert_info("Removing existing conda environment {.val {env_name}} (force = TRUE)")
      reticulate::conda_remove(env_name)
    }
    cli::cli_alert_info("Creating conda environment {.val {env_name}} (python {python_version})")
    reticulate::conda_create(envname = env_name, python_version = python_version)
  } else {
    cli::cli_alert_info("Conda environment {.val {env_name}} already exists, skipping creation (force = FALSE)")
  }

  # ---- 3. Clone source repositories ----
  bonsai_repo <- fs::path(repo_path, "Bonsai-data-representation")
  sanity_repo <- fs::path(repo_path, "Sanity")
  cellstates_repo <- fs::path(repo_path, "cellstates")

  clone_if_absent <- function(url, dest) {
    if (force && fs::dir_exists(dest)) fs::dir_delete(dest)
    if (!fs::dir_exists(dest)) {
      cli::cli_alert_info("Cloning {url}")
      res <- processx::run("git", c("clone", url, dest), error_on_status = FALSE, echo = TRUE)
      if (res$status != 0) {
        cli::cli_abort("Failed to clone {url} (exit status {res$status}). See output above.")
      }
    } else {
      cli::cli_alert_info("{dest} already exists, skipping clone (force = FALSE)")
    }
  }

  clone_if_absent("https://github.com/dhdegroot/Bonsai-data-representation.git", bonsai_repo)
  clone_if_absent("https://github.com/jmbreda/Sanity.git", sanity_repo)
  clone_if_absent("https://github.com/nimwegenLab/cellstates.git", cellstates_repo)

  # ---- 4. Install Bonsai's Python requirements + mpi4py into the env ----
  requirements_file <- fs::path(bonsai_repo, "requirements.txt")
  if (!fs::file_exists(requirements_file)) {
    cli::cli_abort("Could not find {requirements_file}; the Bonsai repo clone may be incomplete or its layout has changed.")
  }

  cli::cli_alert_info("Installing Bonsai Python requirements into {.val {env_name}}")
  requirements <- readLines(requirements_file)
  requirements <- trimws(requirements)
  requirements <- requirements[nzchar(requirements) & !startsWith(requirements, "#")]

  # bonsai_marker_genes() shells out to downstream_analyses/calc_marker_genes.py,
  # which transitively imports the bonsai_scout module (for Bonvis_metadata) --
  # pulling in packages (h5py, natsort, shiny, etc.) that are not in Bonsai's
  # core requirements.txt at all. Confirmed by running it and hitting
  # ModuleNotFoundError one package at a time (h5py, then natsort). Bonsai's
  # own repo ships a separate requirements_bonsai_scout.txt for exactly this
  # -- install it alongside the core requirements rather than special-casing
  # individual missing packages as they keep surfacing.
  scout_requirements_file <- fs::path(bonsai_repo, "requirements_bonsai_scout.txt")
  if (fs::file_exists(scout_requirements_file)) {
    scout_requirements <- readLines(scout_requirements_file)
    scout_requirements <- trimws(scout_requirements)
    scout_requirements <- scout_requirements[nzchar(scout_requirements) & !startsWith(scout_requirements, "#")]
    requirements <- union(requirements, scout_requirements)
  } else {
    cli::cli_warn("Could not find {scout_requirements_file}; bonsai_marker_genes() may be missing dependencies (it needs the bonsai_scout module).")
  }

  reticulate::conda_install(
    envname = env_name,
    packages = requirements,
    pip = TRUE
  )

  cli::cli_alert_info("Installing mpi4py into {.val {env_name}}")
  # mpi4py's build resolves its C compiler by searching PATH for plain
  # "mpicc" (its bundled mpi.cfg does not honor the MPICC env var in 4.x),
  # so the compiler that wins is whichever "mpicc" is first on PATH *at
  # build time*. reticulate::conda_install() routes through a shell script
  # that sources `conda activate`, which itself unconditionally re-prepends
  # conda's own bin dirs to PATH -- so any PATH fix set before calling it is
  # clobbered by activation and never reaches the build. On some conda
  # installs (notably base Anaconda on Apple Silicon) there is a broken
  # "mpicc" wrapper in the base env's bin dir that references a
  # conda-toolchain compiler (e.g. arm64-apple-darwin20.0.0-clang) which was
  # never installed. To avoid that PATH clobbering, we bypass
  # conda_install() here and invoke the target env's own pip directly by
  # absolute path -- no activation script runs, so our PATH override holds.
  mpicc_dir <- fs::path_dir(mpicc_path)
  py_bin <- reticulate::conda_python(envname = env_name)
  withr::with_envvar(
    c(MPICC = mpicc_path,
      PATH = paste(mpicc_dir, Sys.getenv("PATH"), sep = .Platform$path.sep)),
    {
      res <- processx::run(as.character(py_bin),
                            c("-m", "pip", "install", "--upgrade", "mpi4py==4.0.0"),
                            error_on_status = FALSE, echo = TRUE)
      if (res$status != 0) {
        cli::cli_abort(c(
          "Failed to install mpi4py.",
          "i" = "Verify {mpicc_path} is a working MPI compiler wrapper and that no earlier, broken {.code mpicc} shadows it on PATH."
        ))
      }
    }
  )

  # ---- 5. Build Sanity ----
  # NB: Sanity's build must be run from the src/ subdirectory (not the repo
  # root) -- `cd Sanity/src && make` -- and it depends on OpenMP (libgomp1
  # on Linux). The compiled binary lands in Sanity/bin/Sanity regardless of
  # where `make` was invoked from, since the Makefile writes there via a
  # relative path.
  sanity_src <- fs::path(sanity_repo, "src")
  sanity_bin <- fs::path(sanity_repo, "bin", "Sanity")

  if (force || !fs::file_exists(sanity_bin)) {
    if (Sys.info()[["sysname"]] == "Linux") {
      has_libgomp <- tryCatch({
        out <- processx::run("bash", c("-c", "ldconfig -p | grep -q libgomp"),
                              error_on_status = FALSE)
        out$status == 0
      }, error = function(e) NA)
      if (isFALSE(has_libgomp)) {
        cli::cli_abort(c(
          "Sanity requires the OpenMP runtime library (libgomp1), which does not appear to be installed.",
          "i" = "Install it with {.code sudo apt-get install libgomp1} (Debian/Ubuntu) and re-run {.fn bonsai_install}."
        ))
      }
    }
    cli::cli_alert_info("Building Sanity (cd {sanity_src} && make -j{n_threads})")
    res <- processx::run("make", c(paste0("-j", n_threads)), wd = sanity_src,
                          error_on_status = FALSE, echo = TRUE)
    if (res$status != 0 || !fs::file_exists(sanity_bin)) {
      cli::cli_abort(c(
        "Failed to build Sanity.",
        "i" = "Check that a C++ compiler and OpenMP (libgomp1) are installed; see {sanity_repo}/README for macOS-specific instructions (default macOS clang lacks OpenMP support)."
      ))
    }
  } else {
    cli::cli_alert_info("Sanity binary already built at {sanity_bin}, skipping (force = FALSE)")
  }

  # ---- 6. Build cellstates C-extension and confirm the CLI entry point ----
  # cellstates is not a plain Python script -- it has a Cython/C extension
  # that needs compiling first (python setup.py build_ext --inplace), and
  # its CLI entry point is scripts/run_cellstates.py, not a top-level
  # cellstates_main.py.
  cellstates_script <- fs::path(cellstates_repo, "scripts", "run_cellstates.py")
  if (!fs::file_exists(cellstates_script)) {
    cli::cli_warn("Could not find {cellstates_script}; the cellstates repo layout may have changed since this package was written -- check before relying on run_cellstates().")
  } else {
    cli::cli_alert_info("Building cellstates C-extension in-place")
    py_bin <- fs::path(reticulate::conda_python(envname = env_name))

    # cellstates' setup.py tries to Cythonize its .pyx sources at build time,
    # but silently falls back to the .c files already checked into the repo
    # if Cython is not installed. Those checked-in .c files were generated
    # against an older CPython and reference internal headers (e.g.
    # longintrepr.h) that later CPython versions removed from the public
    # include path, so building from them fails on modern Python. Installing
    # Cython lets setup.py regenerate the extension from source against
    # whatever CPython we actually have.
    res_cython <- processx::run(as.character(py_bin), c("-m", "pip", "install", "cython"),
                                 error_on_status = FALSE, echo = TRUE)
    if (res_cython$status != 0) {
      cli::cli_abort("Failed to install Cython, needed to build cellstates' C extension against this Python version.")
    }

    # cellstates/cluster.pyx declares buffers as np.ndarray[np.int_t, ...],
    # but np.int_t was dropped entirely from numpy's Cython declarations
    # (numpy/__init__.pxd) in numpy >=2.0 -- Bonsai's requirements.txt pins
    # numpy 2.x, so Cythonizing the unmodified source fails with "Invalid
    # type". np.intp_t (still present in numpy 2.x) is a safe drop-in: both
    # are pointer-sized integers on the 64-bit platforms this package
    # targets, and np.arange()'s default int dtype matches it.
    cluster_pyx <- fs::path(cellstates_repo, "cellstates", "cluster.pyx")
    if (fs::file_exists(cluster_pyx)) {
      pyx_lines <- readLines(cluster_pyx)
      if (any(grepl("np\\.int_t", pyx_lines))) {
        cli::cli_alert_info("Patching cellstates/cluster.pyx for numpy >=2.0 (np.int_t -> np.intp_t)")
        writeLines(gsub("np\\.int_t", "np.intp_t", pyx_lines), cluster_pyx)
      }
    }

    # cellstates' setup.py build_ext resolves its C compiler via Python's
    # default sysconfig CC, which on macOS is Apple's system clang -- and
    # that clang rejects -fopenmp outright (it has no OpenMP support at
    # all, unlike Linux's gcc/clang+libgomp). Sanity's own Makefile sidesteps
    # this by invoking Homebrew's real GCC directly; do the same here via
    # CC/CXX, since setup.py has no such override built in. GCC's versioned
    # binary name (gcc-16, etc.) changes as Homebrew updates it, so detect
    # the current version rather than hardcoding it.
    build_envvars <- character(0)
    if (Sys.info()[["sysname"]] == "Darwin") {
      find_homebrew_gcc <- function() {
        for (prefix in c("/opt/homebrew/opt/gcc/bin", "/usr/local/opt/gcc/bin")) {
          gccs <- sort(fs::path_file(Sys.glob(fs::path(prefix, "gcc-*"))), decreasing = TRUE)
          gccs <- gccs[grepl("^gcc-[0-9]+$", gccs)]
          if (length(gccs) > 0) {
            ver <- sub("^gcc-", "", gccs[1])
            gxx <- fs::path(prefix, paste0("g++-", ver))
            if (fs::file_exists(gxx)) {
              return(c(CC = fs::path(prefix, gccs[1]), CXX = gxx))
            }
          }
        }
        brew_path <- Sys.which("brew")
        if (brew_path != "") {
          gcc_prefix <- tryCatch(
            trimws(processx::run(brew_path, c("--prefix", "gcc"), error_on_status = FALSE)$stdout),
            error = function(e) ""
          )
          if (nzchar(gcc_prefix)) {
            gccs <- sort(fs::path_file(Sys.glob(fs::path(gcc_prefix, "bin", "gcc-*"))), decreasing = TRUE)
            gccs <- gccs[grepl("^gcc-[0-9]+$", gccs)]
            if (length(gccs) > 0) {
              ver <- sub("^gcc-", "", gccs[1])
              gxx <- fs::path(gcc_prefix, "bin", paste0("g++-", ver))
              if (fs::file_exists(gxx)) {
                return(c(CC = fs::path(gcc_prefix, "bin", gccs[1]), CXX = gxx))
              }
            }
          }
        }
        NULL
      }
      hb_gcc <- find_homebrew_gcc()
      if (!is.null(hb_gcc)) {
        build_envvars <- hb_gcc
        cli::cli_alert_info("Using Homebrew GCC for OpenMP support: {hb_gcc[['CC']]}")
      } else {
        cli::cli_warn(c(
          "Could not find a Homebrew GCC install for building cellstates' OpenMP-enabled C extension.",
          "i" = "Apple's system clang does not support {.code -fopenmp} at all; install a real GCC with {.code brew install gcc} if the build below fails."
        ))
      }
    }

    res <- withr::with_envvar(
      build_envvars,
      processx::run(as.character(py_bin), c("setup.py", "build_ext", "--inplace"),
                    wd = cellstates_repo, error_on_status = FALSE, echo = TRUE)
    )
    if (res$status != 0) {
      cli::cli_abort(c(
        "Failed to build the cellstates C-extension.",
        "i" = "Check that a C compiler with OpenMP support is available inside the {.val {env_name}} conda environment."
      ))
    }
  }

  cli::cli_alert_success("bonsai_install() complete.")

  structure(
    list(
      env_name = env_name,
      repo_path = repo_path,
      bonsai_repo = bonsai_repo,
      sanity_repo = sanity_repo,
      sanity_bin = sanity_bin,
      cellstates_repo = cellstates_repo,
      cellstates_script = cellstates_script
    ),
    class = "bonsai_env"
  )
}

#' @export
print.bonsai_env <- function(x, ...) {
  cli::cli_h3("bonsai_env")
  cli::cli_bullets(c(
    "*" = "conda env: {.val {x$env_name}}",
    "*" = "Bonsai repo: {x$bonsai_repo}",
    "*" = "Sanity binary: {x$sanity_bin}",
    "*" = "cellstates script: {x$cellstates_script}"
  ))
  invisible(x)
}

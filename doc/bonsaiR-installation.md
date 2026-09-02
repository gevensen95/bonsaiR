## Overview

`bonsaiR` wraps a real, non-trivial systems stack – conda, git, a
compiled C++ binary (Sanity), a compiled Cython/C extension
(cellstates), and MPI – not just an R package with CRAN dependencies.
This vignette is purely about getting that stack set up correctly; see
`vignette("bonsaiR-intro")` for actually using the package once
installed.

**Platform note up front:** the macOS instructions below were worked out
by actually hitting and fixing every failure mode described, on real
Apple Silicon hardware. The Linux instructions are based on
`bonsai_install()`’s own coded logic plus one real debugging session on
an institutional HPC cluster (the shared-conda-cache issue below) –
solid, but less exhaustively battle-tested than the macOS path. If you
hit something not covered here on Linux, it’s genuinely useful
information; consider it worth reporting.

The two platforms differ in exactly one important way: **macOS’s system
compiler has no OpenMP support at all**, which Linux’s does by default.
Most of the macOS-specific steps below exist solely to work around that;
skip straight to the Linux section if that’s your platform.

## Common prerequisites (both platforms)

- **R** (4.x), with `devtools` installed
  (`install.packages("devtools")`)
- **conda or Miniconda**, with the `conda` command on your shell’s
  `PATH`
- **git**

## macOS installation

### 1. Xcode Command Line Tools

Needed for `git` and basic build tools. If you don’t already have them:

    xcode-select --install

### 2. Homebrew

If you don’t already have it: <https://brew.sh>

### 3. A real C/C++ compiler with OpenMP support

**This is the step that matters most on macOS.** Apple’s system `clang`
(the default `cc`/`gcc` you get without doing anything) has *no* OpenMP
support whatsoever – not a missing flag, not a config issue, it simply
cannot compile `-fopenmp` code. Both Sanity’s C++ build and cellstates’
Cython extension need real OpenMP. Install Homebrew’s GCC instead:

    brew install gcc

You don’t need to configure anything further yourself –
`bonsai_install()` detects Homebrew GCC automatically (by scanning for a
versioned `gcc-N`/`g++-N` pair under `/opt/homebrew/opt/gcc/bin` or
`/usr/local/opt/gcc/bin`) and uses it specifically for the steps that
need OpenMP.

### 4. System OpenMPI

Needed if you want to run the core tree search with `use_mpi = TRUE`
(recommended above ~2000 cells):

    brew install openmpi

### 5. A known gotcha: a broken `mpicc` can shadow the working one

If you also have Anaconda (not just Miniconda) installed, its base
environment can ship its own `mpicc` wrapper script that expects a
conda-toolchain compiler (e.g. `arm64-apple-darwin20.0.0-clang`) which
was never actually installed – and depending on your shell’s `PATH`
order, *that* broken `mpicc` can be found before the working Homebrew
one, causing cryptic build failures. This is a real bug that was hit and
fixed during this package’s own development, on this exact kind of
setup.

Check which `mpicc` your shell finds first, and whether it actually
works:

    command -v mpicc

    mpicc --version

If that fails with something like
`arm64-apple-darwin20.0.0-clang: command not found`, you’ve hit this
issue. `bonsai_install()` itself works around it internally (it
explicitly forces the known-good compiler to the front of `PATH` for the
specific build steps that need it), so **this doesn’t block
`bonsai_install()` from working** – but if you’re also pre-installing
packages manually (see “Pre-installing the conda environment” below),
you need to apply the same fix yourself, since a plain
`pip install mpi4py` will just use whatever `mpicc` is first on `PATH`.

### 6. Install and set up bonsaiR

    devtools::install_github("gevensen95/bonsaiR")
    library(bonsaiR)
    benv <- bonsai_install()

This clones and compiles everything (slow the first time, fast on every
subsequent call since each step is skipped if it’s already done).

### 7. Verify

    bonsai_smoke_test(bonsai_env = benv)

## Linux installation

### 1. Compiler and MPI packages

Unlike macOS, Linux’s system `gcc` supports OpenMP natively – there is
no Homebrew-GCC-equivalent workaround needed here. On Debian/Ubuntu:

    sudo apt-get install build-essential libgomp1 libopenmpi-dev openmpi-bin git

`libgomp1` is the OpenMP runtime specifically; `bonsai_install()` checks
for it before attempting to build Sanity and gives a clear error naming
this exact package if it’s missing.

On RHEL/CentOS/Rocky and other `yum`/`dnf`-based distributions:

    sudo yum install epel-release -y

    sudo yum install gcc gcc-c++ libgomp openmpi openmpi-devel git -y

(substitute `dnf` for `yum` on newer releases where `yum` is just an
alias for it – both work identically for this).

**This is not quite enough on its own.** Unlike Debian/Ubuntu’s
`openmpi-bin`, the EL family’s `openmpi`/`openmpi-devel` packages
deliberately do *not* put `mpicc`/`mpiexec` on `PATH` by default (to
allow multiple MPI implementations to coexist) – they’re normally
exposed via `module load mpi/openmpi-x86_64`. If your system doesn’t
have environment modules set up (common on a bare/minimal install, or a
shared machine where `module` isn’t configured for your account), add
the OpenMPI bin directory to `PATH` directly instead:

    find / -name "mpicc" -path "*openmpi*" 2>/dev/null

    export PATH="/usr/lib64/openmpi/bin:$PATH"  # adjust to whatever the find above returned

Verify it resolves before moving on:

    mpicc --version

Add the `export PATH=...` line to your `~/.bashrc` (or equivalent) so it
persists across sessions – otherwise `bonsai_install()` will fail its
`mpicc` check the next time you open a new shell.

On other distributions, install the equivalent of: a C/C++ compiler
(`gcc`/`g++`), the OpenMP runtime, an OpenMPI development package
(providing `mpicc`), and `git`.

### 2. conda/Miniconda

If you don’t already have conda, install Miniconda:
<https://docs.conda.io/projects/miniconda/en/latest/>

### 3. A known gotcha: EOL CentOS 7’s package repos

CentOS 7 reached end-of-life in mid-2024, and its default `yum` mirrors
were taken offline as a result – `sudo yum install ...` on an unmodified
CentOS 7 system fails with something like:

    Could not resolve host: mirrorlist.centos.org

(or a similar failure naming a `*.centos.org` mirror, sometimes for
`centos-sclo-sclo` specifically if you also have Software Collections
enabled). This isn’t a `bonsai_install()`/yum configuration issue – it
means the repo definitions still point at CentOS’s live mirror service,
which stopped serving CentOS 7 entirely. Point them at the frozen
archive instead:

    sudo sed -i \
      -e 's/mirrorlist.centos.org/vault.centos.org/g' \
      -e 's/^#.*baseurl=http/baseurl=http/g' \
      -e 's/^mirrorlist=http/#mirrorlist=http/g' \
      /etc/yum.repos.d/CentOS-*.repo

Then retry the `yum install` commands above. If you don’t have `sudo`
access to edit these repo files at all (common on a shared/managed
system), the conda-only alternative below sidesteps the system package
manager entirely.

### 4. Alternative: installing OpenMPI via conda instead of the system package manager

If you can’t get a working system OpenMPI in place – no `sudo` access, a
broken/EOL repo you can’t fix yourself, or a shared HPC login node where
you’d rather not touch system packages at all – you can install OpenMPI
from conda-forge into the `bonsai` conda environment itself, skipping
`yum`/`apt-get` for this piece entirely:

    conda create -n bonsai python=3.11 -y

    conda activate bonsai
    conda install -c conda-forge openmpi -y

**This has to happen before `bonsai_install()` runs, in the same shell/R
session.** `bonsai_install()`’s very first step checks for `mpicc` on
`PATH` and aborts immediately if it doesn’t find one – it does this
*before* it creates or looks at any conda environment, so a
conda-installed `mpicc` only satisfies that check if the `bonsai` env is
already active (and therefore already on `PATH`) when you launch R:

    conda activate bonsai
    R

    library(bonsaiR)
    benv <- bonsai_install()  # reuses the already-active/created "bonsai" env

This still needs a working C/C++ compiler with OpenMP support from the
system (`gcc`/`g++` and `libgomp`) for building Sanity and cellstates –
conda-forge’s OpenMPI only replaces the OpenMPI piece, not the whole
toolchain.

### 5. A known gotcha: shared/HPC conda installs and permission errors

If you’re on an institutional HPC cluster with a shared, multi-user
conda installation (common – look for a `conda info` output showing a
`base environment` path owned by IT/a shared group rather than you
personally), you may hit a permission error like this when
`bonsai_install()` (or a manual `conda create`) tries to write to the
shared package cache:

    PermissionError: [Errno 13] Permission denied:
      '/path/to/shared/miniconda3/pkgs/cache/....info.json'

This is not a bonsaiR bug – it means conda is trying to write a
lock/state file into a shared cache directory you don’t have write
access to. Fix it by pointing conda at directories you actually own
(this only edits your personal `~/.condarc`, not the shared system
config, so it needs no admin permissions and doesn’t affect other users
of the shared install):

    conda config --add pkgs_dirs ~/.conda/pkgs

    conda config --add envs_dirs ~/.conda/envs

Then retry. This fixes both a manual `conda create` and
`bonsai_install()` itself, since both resolve through the same
`pkgs_dirs`/`envs_dirs` search order.

### 6. Install and set up bonsaiR

    devtools::install_github("gevensen95/bonsaiR")
    library(bonsaiR)
    benv <- bonsai_install()

### 7. Verify

    bonsai_smoke_test(bonsai_env = benv)

## Pre-installing the conda environment ahead of time

On either platform, you can set up just the conda environment and Python
packages in advance (e.g. to save time later, or on a machine with
limited internet access) without running the full `bonsai_install()`.
This uses the exact same environment name, Python version, and pinned
package versions `bonsai_install()` itself expects, so a later
`bonsai_install()` call detects everything as already done and skips
straight to cloning and compiling Sanity/cellstates – the one part of
setup that has to happen fresh on every machine regardless, since it
produces compiled, architecture-specific output.

    conda create -n bonsai python=3.11 -y
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate bonsai
    pip install --upgrade matplotlib==3.11.0 numpy==2.4.6 scipy==1.17.1 pandas==3.0.3 \
      scikit-learn==1.9.0 ruamel.yaml==0.19.1 psutil==7.2.2 shiny==1.6.3 faicons==0.2.2 \
      shinyswatch==0.11.0 jinja2==3.1.6 tables==3.11.1 h5py==3.16.0 natsort==8.4.0
    export MPICC="$(command -v mpicc)"
    export PATH="$(dirname "$MPICC"):$PATH"
    pip install --upgrade mpi4py==4.0.0

(This is the same sequence walked through step-by-step earlier in this
package’s development, packaged into one block here – on macOS,
double-check `$MPICC` is the *working* `mpicc` first, per the gotcha
above, since a plain `command -v mpicc` isn’t guaranteed to find it.)

## Caching the environment object

`bonsai_install()`’s return value (`benv`) is a plain list of file paths
with nothing live or unserializable in it – safe to cache so you don’t
need to re-run the (fast, but not free) existence checks every session:

    benv_path <- "~/.bonsaiR_env.rds"
    if (file.exists(benv_path)) {
      benv <- readRDS(benv_path)
    } else {
      benv <- bonsai_install()
      saveRDS(benv, benv_path)
    }

This is only valid on the machine it was created on, and goes stale if
you later run `bonsai_install(force = TRUE)` or otherwise change the
underlying environment – if functions using a cached `benv` start
failing with path-not-found errors, regenerate it first.

## Troubleshooting checklist

If `bonsai_install()` fails partway through, it’s almost always one of:

- **No `mpicc` on `PATH` at all** – install system OpenMPI
  (`brew install openmpi` / `apt-get install openmpi-bin libopenmpi-dev`
  / `yum install openmpi openmpi-devel`), or use the conda-only
  alternative above if you can’t install system packages at all.
- **`yum install openmpi ...` succeeds, but `mpicc` still isn’t found**
  – RHEL/CentOS-specific; the EL family’s OpenMPI packages don’t put
  `mpicc`/`mpiexec` on `PATH` by default. See the gotcha above (either
  `module load mpi/openmpi-x86_64`, or add the OpenMPI bin directory to
  `PATH` manually).
- **`yum install` itself fails with a `mirrorlist.centos.org` (or
  similar) resolution error** – CentOS 7 is end-of-life and its mirrors
  are offline; see the `vault.centos.org` gotcha above.
- **An `mpicc` is found, but it’s broken** – macOS-specific; see the
  gotcha above. Check with `mpicc --version`, not just
  `command -v mpicc`.
- **`libgomp1`/OpenMP runtime missing** – Linux-specific; the error
  message from `bonsai_install()` names this package directly.
- **Permission errors writing to a shared conda install** –
  HPC-specific; see the gotcha above.
- **A stale conda environment from a previous failed attempt** – since
  every step is idempotent, first just re-run
  `bonsai_install(force = FALSE)` (the default) to retry only the step
  that failed. Only use `force = TRUE` (which rebuilds everything from
  scratch) if that doesn’t work.

For anything else, `bonsai_install()`’s own output is streamed live and
usually names the actual failing command – the fastest path to a fix is
usually reading the last ~20 lines of that output rather than guessing.

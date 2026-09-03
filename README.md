# bonsaiR

An R interface to [Bonsai](https://www.biorxiv.org/content/10.1101/2025.05.08.652944v1), a Bayesian method for reconstructing maximum-likelihood tree representations of single-cell RNA-seq data — an alternative to UMAP/t-SNE embeddings that represents cells as a tree rather than a 2D layout, with branch lengths approximating true high-dimensional distances between cells.

> **Note:** `bonsaiR` is an independent, unofficial R port of the Bonsai tool stack, created without solicitation or authorization from the original Bonsai authors. It is not affiliated with, endorsed by, or maintained by them. All credit for the underlying method and tools belongs to the original authors — see [Citation](#citation) below, and please direct questions about the algorithm itself (as opposed to this R wrapper) to the [original repository](https://github.com/dhdegroot/Bonsai-data-representation).

`bonsaiR` takes a `Seurat` object with raw UMI counts and drives it through the full Bonsai tool stack:

1. **[Sanity](https://github.com/jmbreda/Sanity)** — Bayesian normalization producing, for every gene in every cell, a point estimate of true expression *and* an error bar on that estimate.
2. **[cellstates](https://github.com/nimwegenLab/cellstates)** — pre-clustering cells into groups that are statistically indistinguishable given measurement noise, used to warm-start the tree search.
3. **[Bonsai](https://github.com/dhdegroot/Bonsai-data-representation)** — the core maximum-likelihood tree reconstruction.

Sanity and Bonsai's core search are C++/Python tools with no R bindings; `bonsaiR` shells out to them as external processes (optionally under local MPI parallelization for the tree search) rather than wrapping them in-process, and handles the environment setup for you via `bonsai_install()`.

## Requirements

This wraps a real, non-trivial systems stack — not a pure-R package. Before installing, you need:

- **conda or Miniconda** on your `PATH` (used to create an isolated Python environment for the tool stack)
- **git** (to clone the Sanity/cellstates/Bonsai source repos)
- **A C/C++ compiler with OpenMP support.** On Linux this is normally just `gcc`. **On macOS, Apple's system `clang` does not support OpenMP at all** — install a real GCC first (`brew install gcc`); `bonsai_install()` detects and uses it automatically.
- **System OpenMPI** (`brew install openmpi` on macOS, `apt-get install openmpi-bin libopenmpi-dev` on Debian/Ubuntu) if you want to run the tree search across multiple cores via MPI. A single-core fallback (`use_mpi = FALSE`) works without it and is what Bonsai's own docs recommend below ~2000 cells anyway.

## Installation

```r
# install.packages("devtools")
devtools::install_github("gevensen95/bonsaiR")
```

Then, once per machine, build the conda environment and compile the Sanity/cellstates tools:

```r
library(bonsaiR)
benv <- bonsai_install()
```

This is slow the first time (cloning repos, compiling C++/Cython) and near-instant afterward, since every step is skipped if it's already done. Pass `force = TRUE` to rebuild from scratch.

To sanity-check a fresh install end-to-end on synthetic data before pointing it at anything real:

```r
bonsai_smoke_test(bonsai_env = benv)
```

## Quick start

```r
library(bonsaiR)
library(Seurat)

benv <- bonsai_install()
work_dir <- "bonsai_output"

sanity_in    <- bonsai_write_sanity_input(my_seurat_obj, output_dir = file.path(work_dir, "sanity_input"))
sanity_out   <- run_sanity(sanity_in, benv, output_dir = file.path(work_dir, "sanity_output"))
cellstates_out <- run_cellstates(sanity_in, benv, output_dir = file.path(work_dir, "cellstates_output"))

config <- bonsai_write_config(
  sanity_output = sanity_out,
  cellstates_output = cellstates_out,
  bonsai_env = benv,
  dataset_name = "my_dataset",
  results_folder = file.path(work_dir, "bonsai_results")
)

# use_mpi = TRUE (the default) for tens of thousands of cells and up;
# use_mpi = FALSE for smaller datasets, per Bonsai's own recommendation
result <- run_bonsai(config, benv)

tree <- bonsai_read_tree(result, benv)
ape::plot.phylo(tree$phylo, type = "fan", show.tip.label = FALSE)

my_seurat_obj <- bonsai_to_seurat(my_seurat_obj, tree, cellstates_output = cellstates_out)
```

**For a full, runnable walkthrough on real data** (Seurat's built-in `pbmc_small`, including marker-gene identification between clades), see [`vignette("bonsaiR-intro")`](vignettes/bonsaiR-intro.Rmd).

## Function reference

| Function | Purpose |
|---|---|
| `bonsai_install()` | One-time setup: conda env, cloned/built Sanity + cellstates + Bonsai |
| `bonsai_smoke_test()` | Runs the full pipeline on synthetic data with per-stage pass/fail reporting |
| `bonsai_write_sanity_input()` | Extracts raw counts from a Seurat object into Sanity's input format |
| `run_sanity()` | Runs Sanity normalization |
| `run_cellstates()` | Runs cellstates pre-clustering |
| `bonsai_write_config()` | Writes Bonsai's YAML config and a cellstates-based warm-start tree |
| `run_bonsai()` | Runs the core tree search (single-core or MPI) |
| `bonsai_read_tree()` | Reads a completed tree (topology, edge/vertex info, posterior expression estimates) into R |
| `bonsai_marker_genes()` | Identifies marker genes distinguishing two clades of the tree |
| `bonsai_to_seurat()` | Attaches the tree, posterior expression estimates, and cellstate labels back onto a Seurat object |

## Known limitations

- **Scale and MPI reliability:** Bonsai's own multi-round core-calculation step has an intermittent race condition under MPI (`use_mpi = TRUE`) that has been observed on small datasets — a shared scratch folder occasionally gets cleaned up by one stage while another is still reading it. This is timing-dependent (not deterministic) and lives in Bonsai's own code, not this wrapper. It's most likely to show up exactly where Bonsai's docs already say not to use MPI anyway (below ~2000 cells); if you hit it, retry or set `use_mpi = FALSE`. See `?run_bonsai` for details.
- **Not tested at "tens of thousands of cells" scale.** This package has been validated end-to-end on synthetic data and on Seurat's small `pbmc_small` example; Bonsai's intended production scale (large real atlases, HPC/SLURM scheduling) hasn't been exercised here.
- Marker-gene comparisons will refuse groups that split a zero-branch-length cluster of cells rather than silently comparing something misleading (see `?bonsai_marker_genes`).

## License

This package is licensed under **CC-BY-NC-4.0** (non-commercial use only) — matching the license of the Bonsai/Sanity/cellstates tool stack it drives, which this package is unusable without. Check this applies to your use case before relying on it in a commercial setting. See [`LICENSE`](LICENSE) for the full text.

## Citation

`bonsaiR` is just a wrapper — if you use it, please cite the original Bonsai software and paper, not this package:

> de Groot, D.H., Morillo Leonardo, S.X., Pachkov, M., & van Nimwegen, E. (2025). *Bonsai: Tree representations for distortion-free visualization and exploratory analysis of single-cell omics data.* bioRxiv. https://doi.org/10.1101/2025.05.08.652944

This has since been peer-reviewed and published as de Groot et al., *Nature Biotechnology* (2026), https://doi.org/10.1038/s41587-026-03220-2.

> de Groot, D.H., Morillo Leonardo, S.X., Pachkov, M., & van Nimwegen, E. (2026). *Bonsai-data-representation* (v1.0.0) [Software]. Zenodo. https://doi.org/10.5281/zenodo.20370956

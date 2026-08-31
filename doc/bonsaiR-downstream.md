## Overview

This vignette covers three functions for analyzing a Bonsai tree once
you’ve reconstructed it, building on the same `pbmc_small` example used
in `vignette("bonsaiR-intro")`:

- `run_average_over_groups()` – noise-aware average expression per
  group, plus a ranking of which genes vary most across groups.
- `bonsai_cluster_by_annotation()` – cuts the tree into clusters that
  best match a known annotation, and separately reports how well the
  tree’s geometry alone (with no annotation at all) recovers those same
  labels.
- `bonsai_cluster_tree()` – cuts the tree into clusters using geometry
  alone, with no annotation involved – roughly `stats::cutree()`, but
  respecting the tree’s actual branch lengths rather than just topology.

If you’ve already worked through `vignette("bonsaiR-intro")` you’ll have
everything these functions need in your session already (`sanity_out`,
`cellstates_out`, `result`, `tree`, `pbmc_small`); this vignette
re-derives them from scratch so it also runs on its own.

## Setup

Same quick pipeline as the intro vignette – see it for a full
explanation of each step.

    library(bonsaiR)
    library(Seurat)

    benv <- bonsai_install()
    work_dir <- tempfile("bonsaiR_downstream_")
    dir.create(work_dir)

    data("pbmc_small")
    cluster_by_cell <- setNames(as.character(Idents(pbmc_small)), colnames(pbmc_small))

    sanity_in <- bonsai_write_sanity_input(pbmc_small, output_dir = file.path(work_dir, "sanity_input"))
    sanity_out <- run_sanity(sanity_in, benv, output_dir = file.path(work_dir, "sanity_output"))
    cellstates_out <- run_cellstates(sanity_in, benv, output_dir = file.path(work_dir, "cellstates_output"))

    config <- bonsai_write_config(
      sanity_output = sanity_out,
      cellstates_output = cellstates_out,
      bonsai_env = benv,
      dataset_name = "pbmc_small",
      results_folder = file.path(work_dir, "bonsai_results")
    )
    result <- run_bonsai(config, benv, use_mpi = FALSE)
    tree <- bonsai_read_tree(result, benv)

## Averaging expression over groups

`run_average_over_groups()` is the multi-group counterpart to
`bonsai_marker_genes()`’s pairwise two-group comparison – instead of
testing one pair of clades, it summarizes every gene across *all* your
groups at once, weighting each cell by its own Sanity error bar (not a
plain mean).

Here we group by `pbmc_small`’s own Seurat clusters, though any named
grouping works the same way – a Zone/cell-type annotation, cellstates
labels, or clades you’ve identified by eye in the tree plot.

    avg_result <- run_average_over_groups(
      sanity_out, cluster_by_cell, benv,
      output_dir = file.path(work_dir, "avg_by_cluster")
    )
    avg_result

    # groups x genes: average expression per cluster
    avg_result$avg_activities[, 1:5]

    # per-gene ranking of how much expression varies across clusters
    head(avg_result$significance, 10)

The top of `significance` should look familiar if you ran the
marker-genes step in the intro vignette: genes like `HLA-DPB1`,
`HLA-DRB1`, `CST3`, and `LYZ` tend to come out on top – the same
antigen-presenting/myeloid markers, just found by summarizing across all
clusters at once rather than comparing one pair.

## Clustering the tree to match a known annotation

`bonsai_cluster_by_annotation()` cuts the tree to maximize Normalized
Mutual Information (NMI) against a given annotation. It also computes a
second clustering using *only* tree geometry (minimizing within-cluster
distance, never seeing your annotation at all) for comparison.

    nmi_result <- bonsai_cluster_by_annotation(
      result, cluster_by_cell, benv,
      output_dir = file.path(work_dir, "nmi_clustering")
    )
    nmi_result

    nmi_result$nmi_scores
    head(nmi_result$clustering_results)

Read `nmi_scores` like this: the annotation-based row’s NMI just
confirms the optimizer worked – it was handed your labels and told to
match them. The **distance-based row’s NMI is the informative number**:
that clustering never saw `cluster_by_cell` at all, so a high score
there is a genuine, quantitative measure of how well the tree’s
structure alone recovers your known groups – not something you can get
from eyeballing a tree plot.

## Clustering the tree by geometry alone

`bonsai_cluster_tree()` is the same distance-based clustering, exposed
directly rather than as a side-by-side comparison – useful when you
don’t have an annotation to compare against yet, or want to explore
clusterings at several granularities at once.

    tree_clusters <- bonsai_cluster_tree(result, n_clusters = 10, bonsai_env = benv)
    str(tree_clusters)

    # cluster sizes at the finest level computed
    table(tree_clusters[[ncol(tree_clusters)]])

`tree_clusters` has one column per cut level, from a few clusters up to
`n_clusters` – a byproduct of the underlying greedy-splitting algorithm.
Pick whichever column’s granularity is useful for your analysis
(e.g. via `table()` on each column to see how cluster sizes change
across levels), the same way you’d choose a `k` when calling
`stats::cutree()`.

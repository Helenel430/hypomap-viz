# Romanov_10x

An interactive [Shiny](https://shiny.posit.co/) explorer for the **Romanov_10x**
mouse hypothalamus single-cell dataset, built with
[ShinyCell](https://github.com/SGDDNB/ShinyCell). Type a gene and see its
expression painted on the UMAP, compare cell types, view proportions, and more —
all running locally in your browser.

**Dataset:** GEO [GSE132730](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132730)
— Romanov et al. 2020, *"Design logic of hypothalamus development mapped by
single-cell RNA-seq"* (mouse, 10x Genomics, ages E15.5–P23). The label
`Romanov_10x` follows the HypoMap atlas's `[FirstAuthor]_[technique]`
dataset-naming convention.

> **Note on data:** the raw Seurat object (~1.75 GB) and the app's generated
> expression files are **not** stored in this repo (GitHub caps files at 100 MB).
> Each user regenerates them locally from their own copy of the `.rds` — see below.

## What you need

- **R** (4.4+) and ideally **RStudio**
- The dataset: `GSE132730_TractNAE_integrated.rds`
  (download from GEO, or get it from the lab's shared drive)

## One-time setup

Install the required R packages:

```r
# Core + the packages the generated app loads at runtime
install.packages(c("shiny", "Seurat", "remotes", "hdf5r",
                   "shinyhelper", "DT", "ggdendro", "ggrepel",
                   "data.table", "ggplot2", "gridExtra", "magrittr",
                   "Matrix", "RColorBrewer"))

# ShinyCell itself (installed from a direct tarball to avoid GitHub auth issues)
remotes::install_url(
  "https://github.com/SGDDNB/ShinyCell/archive/refs/heads/master.tar.gz")
```

## Run it (3 steps)

1. **Place the data** — put `GSE132730_TractNAE_integrated.rds` in the `data/`
   folder of this repo.
2. **Build the app** — from this folder, run the generator once:
   ```r
   source("prepare_app.R")
   ```
   This reads the `.rds` and writes the app into `shinyApp/` (takes a few
   minutes; needs a few GB of RAM).
3. **Launch it** — opens in your browser on `localhost`:
   ```r
   shiny::runApp("shinyApp")
   ```

## What you can explore

- **Gene expression on the UMAP / tSNE** for any gene
- **Color by** cell type (`Annotation`), broad `Class`, `cluster`, developmental
  `Age` (E15–P23), or `SampleID`
- **Two-gene co-expression**, violin plots, and a gene-vs-cell-type bubble plot
  (the interactive version of the expression-by-cell-type table)
- **Cell-type proportions** across groups

## Gene-expression lookup (lightweight, no UMAP/tSNE)

A second, much lighter app for pure expression questions — type a gene and get:

- **% of cells expressing** it (globally and per cell type) — detection rate
- **Mean expression among expressing cells** (intensity when "on")
- **Expression percentile** (genes ranked by mean-when-expressed; flags the top 5%)

**Cell cohort:** to match the HypoMap atlas, the stats are computed on **P23 cells
only** (2,415 cells, 7 cell types) — the same `RomanovDev10x` subset HypoMap used
(timepoint P23 only; embryonic + earlier postnatal cells dropped; see
`HypoMap_study_display_card.csv`). To use all 51,245 cells instead, set
`keep_ages` near the top of `compute_gene_stats.R`.

It runs off a tiny precomputed file (`gene-explorer/gene_stats.rds`, ~2 MB) — **no
1.75 GB `.rds` needed**, so it's fully self-contained in this repo:

```r
shiny::runApp("gene-explorer")
```

To regenerate the stats (only if the source data changes):

```r
source("compute_gene_stats.R")   # reads the .rds, rewrites gene-explorer/gene_stats.rds
```

## Repo layout

```
hypomap-viz/
├─ prepare_app.R         # generates the UMAP app from the .rds      (tracked)
├─ shinyApp/             # ShinyCell UMAP/tSNE explorer
│   ├─ ui.R, server.R, *.R   # app code                            (tracked)
│   └─ sc1*.rds / .h5         # generated data                     (git-ignored)
├─ compute_gene_stats.R  # precomputes per-gene stats from the .rds  (tracked)
├─ gene-explorer/        # lightweight gene-lookup app
│   ├─ app.R                  # app code                            (tracked)
│   └─ gene_stats.rds         # ~4 MB precomputed stats             (tracked)
├─ data/                 # put the .rds here                        (git-ignored)
├─ .gitignore
└─ README.md
```

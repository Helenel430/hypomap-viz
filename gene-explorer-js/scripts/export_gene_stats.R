# ============================================================
# Exports ../../gene-explorer/gene_stats.rds (R binary) into a compact
# JSON file the JS app can fetch directly — no R server involved at
# runtime. Re-run this whenever gene_stats.rds is regenerated
# (see ../../compute_gene_stats.R).
#
# Run from anywhere with:  Rscript gene-explorer-js/scripts/export_gene_stats.R
# ============================================================
library(jsonlite)

here <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))
stats_path <- file.path(here, "..", "..", "gene-explorer", "gene_stats.rds")
out_path   <- file.path(here, "..", "public", "data", "gene_stats.json")

stats <- readRDS(stats_path)
cts   <- stats$celltypes

# round to the precision the UI actually displays, to keep the JSON small
r2 <- function(x) round(x, 2)
r3 <- function(x) round(x, 3)
r1 <- function(x) round(x, 1)

# row-major: byct$<metric>[[i]] holds the per-celltype values for genes[i],
# in the same order as `celltypes` — matches how the UI looks up one gene at a time.
to_rows <- function(mat) unname(lapply(seq_len(nrow(mat)), function(i) as.numeric(mat[i, ])))

payload <- list(
  totalCells = stats$total_cells,
  celltypes  = cts,
  cellCounts = unname(as.integer(stats$cell_counts[cts])),
  genes      = stats$genes,
  global = list(
    nExpressing    = as.integer(stats$global$n_expressing),
    pctExpressing  = r2(stats$global$pct_expressing),
    meanExpressing = r3(stats$global$mean_expressing),
    percentile     = r1(stats$global$percentile)
  ),
  byct = list(
    nExpressing    = to_rows(stats$byct$n_expressing[, cts, drop = FALSE]),
    pctExpressing  = to_rows(r2(stats$byct$pct_expressing[, cts, drop = FALSE])),
    meanExpressing = to_rows(r3(stats$byct$mean_expressing[, cts, drop = FALSE])),
    percentile     = to_rows(r1(stats$byct$percentile[, cts, drop = FALSE]))
  )
)

write_json(payload, out_path, auto_unbox = TRUE, na = "null", digits = NA)
cat(sprintf("Wrote %s (%.1f MB)\n", out_path, file.info(out_path)$size / 1e6))

# ============================================================
# compute_gene_stats.R
# Precompute per-gene expression stats for the gene-lookup app:
#   - % of cells expressing the gene  (detection rate)
#   - mean expression among EXPRESSING cells (conditional mean)
#   - expression percentile, ranked by that conditional mean
# ...both GLOBALLY and within each Annotation cell type.
#
# Run once. Writes gene-explorer/gene_stats.rds (small; the app loads only this).
# ============================================================
suppressMessages({
  library(Seurat)
  library(Matrix)
})

rds_path <- "../GSE132730_TractNAE_integrated.rds"   # relative to hypomap-viz/
if (!file.exists(rds_path)) rds_path <- "GSE132730_TractNAE_integrated.rds"
cat("Loading", rds_path, "...\n")
obj <- suppressWarnings(suppressMessages(UpdateSeuratObject(readRDS(rds_path))))
DefaultAssay(obj) <- "RNA"

# --- expression matrix (genes x cells), log-normalized -------------------
expr <- as(GetAssayData(obj, assay = "RNA", layer = "data"), "CsparseMatrix")
grp  <- obj@meta.data[["Annotation"]]
keep <- !is.na(grp)
expr <- expr[, keep]
grp  <- droplevels(factor(grp[keep]))
genes <- rownames(expr)
ncell <- ncol(expr)
cat(length(genes), "genes x", ncell, "cells;", nlevels(grp), "cell types\n")

# detection matrix: every nonzero -> 1
det <- expr; det@x <- rep(1, length(det@x))

# helper: percentile rank (0-100) of meanE, computed only among "expressed"
# genes (>=1 expressing cell); never-detected genes get NA.
pctile <- function(meanE, nexp) {
  pr <- rep(NA_real_, length(meanE))
  ex <- nexp > 0
  pr[ex] <- 100 * rank(meanE[ex], ties.method = "average") / sum(ex)
  pr
}

# --- GLOBAL ----------------------------------------------------------------
g_nexp  <- as.numeric(Matrix::rowSums(det))
g_sum   <- as.numeric(Matrix::rowSums(expr))
g_meanE <- ifelse(g_nexp > 0, g_sum / g_nexp, 0)          # mean among expressing
g_pct   <- 100 * g_nexp / ncell
g_perc  <- pctile(g_meanE, g_nexp)
global <- data.frame(
  gene        = genes,
  n_expressing = g_nexp,
  pct_expressing = round(g_pct, 3),
  mean_expressing = round(g_meanE, 4),
  percentile  = round(g_perc, 2),
  top5pct     = !is.na(g_perc) & g_perc >= 95,
  row.names = NULL, stringsAsFactors = FALSE)

# --- PER CELL TYPE ---------------------------------------------------------
ind <- sparseMatrix(i = seq_along(grp), j = as.integer(grp),
                    x = 1, dims = c(length(grp), nlevels(grp)))
colnames(ind) <- levels(grp)
n_per <- as.numeric(table(grp))

nexp_ct <- as.matrix(det  %*% ind)              # genes x celltype
sum_ct  <- as.matrix(expr %*% ind)
meanE_ct <- sum_ct / nexp_ct; meanE_ct[nexp_ct == 0] <- 0
pct_ct  <- sweep(nexp_ct, 2, n_per, "/") * 100
perc_ct <- meanE_ct
for (j in seq_len(ncol(meanE_ct))) perc_ct[, j] <- pctile(meanE_ct[, j], nexp_ct[, j])
dimnames(nexp_ct) <- dimnames(meanE_ct) <- dimnames(pct_ct) <- dimnames(perc_ct) <-
  list(genes, levels(grp))

stats <- list(
  genes = genes,
  celltypes = levels(grp),
  cell_counts = setNames(as.integer(n_per), levels(grp)),
  total_cells = ncell,
  global = global,
  byct = list(n_expressing = nexp_ct, pct_expressing = round(pct_ct, 3),
              mean_expressing = round(meanE_ct, 4), percentile = round(perc_ct, 2)))

saveRDS(stats, "gene-explorer/gene_stats.rds")
cat("\nSaved gene-explorer/gene_stats.rds (",
    round(file.size("gene-explorer/gene_stats.rds")/1e6, 1), "MB )\n")
cat("Example (Pomc) global:\n"); print(global[global$gene == "Pomc", ])

# ============================================================
# prepare_app.R
# Build the "HypoMap Visualization" interactive Shiny explorer from the
# GSE132730 hypothalamus Seurat object, using ShinyCell.
#
# Run this ONCE locally to (re)generate the app into shinyApp/.
# The heavy data files it produces are git-ignored, so each user
# regenerates them from their own copy of the .rds (see README).
# ============================================================
suppressMessages({
  library(Seurat)
  library(ShinyCell)
})

# --- 1. Locate the raw Seurat object --------------------------------------
# Preferred location for sharing: data/GSE132730_TractNAE_integrated.rds
# (also falls back to the parent rstudioTest folder for the original setup).
candidates <- c(
  "data/GSE132730_TractNAE_integrated.rds",
  "../GSE132730_TractNAE_integrated.rds"
)
rds_path <- candidates[file.exists(candidates)][1]
if (is.na(rds_path)) {
  stop("Seurat .rds not found. Place it at data/GSE132730_TractNAE_integrated.rds")
}
cat("Loading", rds_path, "...\n")
obj <- suppressWarnings(suppressMessages(UpdateSeuratObject(readRDS(rds_path))))
DefaultAssay(obj) <- "RNA"   # use real expression, not the 'integrated' assay

# --- 2. Sensible default marker genes (mouse hypothalamus) ----------------
# Only genes actually present in the data are kept.
markers <- c("Pomc","Agrp","Npy","Gad1","Gad2","Slc17a6","Slc32a1",
             "Rax","Gfap","Aqp4","Mbp","Plp1","Pdgfra","Cx3cr1",
             "Sox2","Nkx2-1","Sst","Th","Hcrt","Pmch")
markers <- intersect(markers, rownames(obj))
if (length(markers) < 2) markers <- head(VariableFeatures(obj), 10)
cat("Default genes:", paste(markers, collapse = ", "), "\n")

# --- 3. Which metadata to expose ------------------------------------------
# Deliberately EXCLUDE the 100+ "prediction.score.*" columns so the app's
# dropdowns stay clean — keep the useful groupings only.
meta_keep <- intersect(
  c("Annotation", "Subclass", "Class", "cluster", "Age", "SampleID",
    "nCount_RNA", "nFeature_RNA", "percent.mito"),
  colnames(obj@meta.data))
cat("Metadata shown:", paste(meta_keep, collapse = ", "), "\n")

# --- 4. UMAP embedding column names (for the default plot) ----------------
umap_cols <- colnames(obj@reductions[["umap"]]@cell.embeddings)[1:2]

# --- 5. Generate the Shiny app --------------------------------------------
scConf <- createConfig(obj, meta.to.include = meta_keep)
makeShinyApp(
  obj, scConf,
  gex.assay         = "RNA",
  gex.slot          = "data",
  gene.mapping      = FALSE,
  shiny.title       = "HypoMap Visualization",
  shiny.dir         = "shinyApp/",
  default.gene1     = markers[1],
  default.gene2     = markers[2],
  default.multigene = head(markers, 10),
  default.dimred    = umap_cols
)

cat("\nDone. Launch the app with:\n")
cat('  shiny::runApp("shinyApp")\n')

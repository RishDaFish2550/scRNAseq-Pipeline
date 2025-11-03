suppressPackageStartupMessages({
  library(Seurat)
  library(DropletUtils)  
  library(Matrix)
  library(ggplot2)
})

.read10x_any <- function(p){
  if (dir.exists(p)) {
    Read10X(data.dir = p)
  } else {
    Read10X_h5(p)
  }
}

.tune_resolution <- function(obj, target, lo = 0.02, hi = 1.0, maxit = 12, tol = 1L){
  best <- lo
  for (i in seq_len(maxit)) {
    mid <- (lo + hi) / 2
    obj <- FindClusters(obj, resolution = mid, algorithm = 1, verbose = FALSE)
    k <- length(unique(obj$seurat_clusters))
    if (abs(k - target) <= tol) return(list(obj = obj, res = mid, k = k))
    if (k < target) { lo <- mid } else { hi <- mid }
    best <- mid
  }
  obj <- FindClusters(obj, resolution = best, algorithm = 1, verbose = FALSE)
  list(obj = obj, res = best, k = length(unique(obj$seurat_clusters)))
}

raw_path <- "data/raw_matrix"            
out_png  <- file.path("FIGURES_UMAP", "MouseBrain_S3_umap.png")
dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)

sample_id <- "MouseBrain"

targets <- c(S1=1658, S2=7495, S3=16293, S4=23762, S5=30091)

raw <- .read10x_any(raw_path)
if (is.list(raw)) {
  raw <- if ("Gene Expression" %in% names(raw)) raw[["Gene Expression"]] else raw[[1]]
}
raw <- as(raw, "dgCMatrix")

totals <- Matrix::colSums(raw)
ord    <- order(totals, decreasing = TRUE, na.last = NA)
bcs_by_rank <- colnames(raw)[ord]

n_s3 <- min(as.integer(targets[["S3"]]), length(bcs_by_rank))
bcs_s3 <- bcs_by_rank[seq_len(n_s3)]

raw_s3 <- raw[, bcs_s3, drop = FALSE]
rm(raw); gc()

obj <- CreateSeuratObject(counts = raw_s3, project = paste0(sample_id, "_S3"), min.cells = 1, min.features = 100)
rm(raw_s3); gc()

obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)
obj <- FindVariableFeatures(obj, nfeatures = 2000, selection.method = "vst", verbose = FALSE)
obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = 30, verbose = FALSE)

obj <- FindNeighbors(obj, dims = 1:30, k.param = 50, verbose = FALSE)
tuned <- .tune_resolution(obj, target = 23, lo = 0.03, hi = 1.2)
obj   <- tuned$obj
message(sprintf("Picked resolution=%.3f  -> %d clusters", tuned$res, tuned$k))

set.seed(1337)
obj <- RunUMAP(
  obj, dims = 1:30,
  seed.use = 1337,
  n.neighbors = 50,
  min.dist = 0.2,
  metric = "cosine",
  umap.method = "uwot",
  verbose = FALSE
)

k <- length(unique(obj$seurat_clusters))
p <- DimPlot(
  obj, reduction = "umap", group.by = "seurat_clusters", label = FALSE
) + ggtitle(sprintf("Graph Based (%d Clusters)", k)) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(out_png, p, width = 7.5, height = 5.0, dpi = 300)
cat("UMAP written to:", normalizePath(out_png), "\n")

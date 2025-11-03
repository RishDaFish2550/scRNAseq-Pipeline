suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
})

set.seed(1337)

raw_path  <- "data/raw_matrix"    
out_dir   <- "FIGURES_UMAP"
sample_id <- "MouseBrain"
s2_target_cells <- 7495             
target_k <- 25                     

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

.read10x_any <- function(p){
  if (dir.exists(p)) Read10X(data.dir = p) else Read10X_h5(p)
}

.pick_top_barcodes <- function(mtx, n){
  totals <- Matrix::colSums(mtx)
  ord <- order(totals, decreasing = TRUE, na.last = NA)
  colnames(mtx)[ord][seq_len(min(n, ncol(mtx)))]
}

.tune_resolution <- function(obj, target, max_iter = 10, lo = 0.05, hi = 2.0){
 
  for(i in seq_len(max_iter)){
    res <- (lo + hi) / 2
    obj <- FindClusters(obj, resolution = res, algorithm = 1, verbose = FALSE)
    k <- length(unique(Idents(obj)))
    if (k > target) hi <- res else lo <- res
  }
  list(obj = obj, res = res, k = k)
}

raw <- .read10x_any(raw_path)
if (is.list(raw)) raw <- if ("Gene Expression" %in% names(raw)) raw[["Gene Expression"]] else raw[[1]]
raw <- as(raw, "dgCMatrix")

s2_barcodes <- .pick_top_barcodes(raw, s2_target_cells)
raw_s2 <- raw[, s2_barcodes, drop = FALSE]

obj <- CreateSeuratObject(raw_s2, project = "S2", min.cells = 1, min.features = 0)
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = 30, verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:30, k.param = 30, verbose = FALSE)

tuned <- .tune_resolution(obj, target = target_k)
obj   <- tuned$obj
used_res <- tuned$res
k_now    <- length(unique(Idents(obj)))

obj <- RunUMAP(obj, dims = 1:30, seed.use = 1337, verbose = FALSE)

p <- DimPlot(obj, reduction = "umap", label = FALSE, repel = TRUE) +
  ggtitle(sprintf("Graph Based (%d Clusters)", k_now)) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "right")

out_png <- file.path(out_dir, sprintf("%s_S2_umap.png", sample_id))
ggsave(out_png, p, width = 8.5, height = 5.5, dpi = 300)
cat(sprintf("UMAP written to: %s (resolution=%.3f, clusters=%d)\n",
            normalizePath(out_png), used_res, k_now))

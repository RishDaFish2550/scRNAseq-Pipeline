suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
})

raw_path  <- "data/raw_matrix"                
sample_id <- "MouseBrain"
target_n  <- 23762                             
out_png   <- file.path("FIGURES_UMAP", sprintf("%s_S4_umap.png", sample_id))
dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)
set.seed(1337)

read10x_any <- function(p){
  x <- Read10X(p)
  if (is.list(x)) {
    if ("Gene Expression" %in% names(x)) x <- x[["Gene Expression"]] else x <- x[[1]]
  }
  as(x, "dgCMatrix")
}

raw <- read10x_any(raw_path)

totals <- Matrix::colSums(raw)
ord    <- order(totals, decreasing = TRUE, na.last = NA)
cells_S4 <- colnames(raw)[ord][seq_len(min(target_n, ncol(raw)))]

cat(sprintf("[S4] Using top %d barcodes by total counts\n", length(cells_S4)))

obj <- CreateSeuratObject(
  counts = raw[, cells_S4, drop = FALSE],
  project = paste0(sample_id, "_S4"),
  min.cells = 3,
  min.features = 200
)

obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = 30, verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:30)

target_k <- 24
res_grid <- seq(0.05, 3.0, by = 0.05)
picked_res <- NA_real_

for (r in res_grid) {
  obj <- FindClusters(obj, resolution = r, verbose = FALSE)
  k <- length(unique(obj$seurat_clusters))
  if (k == target_k) {
    picked_res <- r
    break
  }
}
if (is.na(picked_res)) {
  
  diff_vec <- sapply(res_grid, function(r){
    obj <- FindClusters(obj, resolution = r, verbose = FALSE)
    abs(length(unique(obj$seurat_clusters)) - target_k)
  })
  picked_res <- res_grid[which.min(diff_vec)]
  obj <- FindClusters(obj, resolution = picked_res, verbose = FALSE)
}

obj <- RunUMAP(obj, dims = 1:30, seed.use = 1337)

k <- length(unique(obj$seurat_clusters))
cat(sprintf("Picked resolution = %.2f -> %d clusters\n", picked_res, k))

p <- DimPlot(obj, reduction = "umap", group.by = "seurat_clusters", label = FALSE) +
  ggtitle(sprintf("Graph Based (%d Clusters)", k)) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2),
        plot.title = element_text(hjust = 0, face = "bold"))

ggsave(out_png, p, width = 7.5, height = 5.0, dpi = 300)
cat("UMAP written to:", normalizePath(out_png), "\n")

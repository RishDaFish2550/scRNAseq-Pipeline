suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
})

raw_path  <- "data/raw_matrix"                 
sample_id <- "MouseBrain"
target_n  <- 1658                              
out_png   <- file.path("FIGURES_UMAP", sprintf("%s_S1_umap.png", sample_id))
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
cells_S1 <- colnames(raw)[ord][seq_len(min(target_n, ncol(raw)))]

cat(sprintf("[S1] Using top %d barcodes by total counts\n", length(cells_S1)))

obj <- CreateSeuratObject(
  counts = raw[, cells_S1, drop = FALSE],
  project = paste0(sample_id, "_S1"),
  min.cells = 3,        
  min.features = 200    
)

obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)

obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)

obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = 30, verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:30)

target_k <- 18
res_grid <- c(seq(0.2, 2.0, by = 0.05), seq(2.2, 4, by = 0.2))
picked_res <- NA_real_
for (r in res_grid) {
  obj <- FindClusters(obj, resolution = r, verbose = FALSE)
  k <- length(unique(obj$seurat_clusters))
  if (k == target_k) { picked_res <- r; break }
}
if (is.na(picked_res)) picked_res <- tail(res_grid, 1)
obj <- FindClusters(obj, resolution = picked_res, verbose = FALSE)

obj <- RunUMAP(obj, dims = 1:30, seed.use = 1337)

k <- length(unique(obj$seurat_clusters))
cat(sprintf("Picked resolution=%.2f -> %d clusters\n", picked_res, k))

p <- DimPlot(obj, reduction = "umap", group.by = "seurat_clusters", label = FALSE) +
  ggtitle(sprintf("Graph Based (%d Clusters)", k)) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2),
        plot.title = element_text(hjust = 0, face = "bold"))

ggsave(out_png, p, width = 7.5, height = 5.0, dpi = 300)
cat(" UMAP written to:", normalizePath(out_png), "\n")

suppressPackageStartupMessages({
  library(Seurat)        
  library(DropletUtils)  
  library(Matrix)
  library(ggplot2)
})

.read10x_any <- function(p){
  if (dir.exists(p)) Read10X(data.dir = p) else Read10X_h5(p)
}

.knee_plot <- function(raw_counts, called_barcodes, out_png, title){
  totals <- Matrix::colSums(raw_counts)
  ord <- order(totals, decreasing = TRUE, na.last = NA)
  totals_o <- as.numeric(totals[ord])
  bcs_o    <- colnames(raw_counts)[ord]
  called_o <- bcs_o %in% called_barcodes

  df <- data.frame(rank = seq_along(totals_o),
                   total = totals_o,
                   called = called_o)

  p <- ggplot(df, aes(rank, total)) +
    geom_line(data = subset(df, !called), aes(y = total), color = "black", linewidth = 0.6) +  
    geom_line(data = subset(df, called), aes(y = total), color = "blue", linewidth = 1.4) +   
    scale_x_log10() +
    scale_y_log10() +
    labs(x = "Barcode Rank", y = "# Transcripts", title = title) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )

  dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)
  ggsave(out_png, p, width = 7, height = 5, dpi = 200)
}

raw_path   <- "data/raw_matrix"      
out_dir    <- "FIGURES_KNEE"
sample_id  <- "MouseBrain"

sens <- c(S1=1e-6, S2=1e-5, S3=1e-4, S4=1e-3, S5=1e-2)

# Target blue lengths (exact, for my reference)
# targets <- c(
#   S1 = 1658,
#   S2 = 7495,
#   S3 = 16293,
#   S4 = 23762,
#   S5 = 30091
# )

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

raw <- .read10x_any(raw_path)
if (is.list(raw)) raw <- if ("Gene Expression" %in% names(raw)) raw[["Gene Expression"]] else raw[[1]]
raw <- as(raw, "dgCMatrix")

.totals <- Matrix::colSums(raw)
.ord_rank <- order(.totals, decreasing = TRUE, na.last = NA)
.barcodes_by_rank <- colnames(raw)[.ord_rank]

br   <- barcodeRanks(raw)
knee <- tryCatch(DropletUtils::metadata(br)$knee, error=function(e) NA_real_)
if (!is.finite(knee)) knee <- median(Matrix::colSums(raw))
lower <- max(100L, round(knee * 0.1))

set.seed(1337)
ed <- emptyDrops(raw, lower = lower, retain = TRUE, niters = 10000, test.ambient = TRUE)
fdr <- ed$FDR; names(fdr) <- rownames(ed)

cat(sprintf("Knee=%.1f, lower=%d, barcodes=%d\n",
            knee, as.integer(lower), ncol(raw)))
qs <- quantile(fdr, probs=c(0, .25, .5, .75, .9, .95, .99), na.rm=TRUE)
print(round(qs, 6))
for (nm in names(sens)) {
  thr <- sens[[nm]]
  n_called <- sum(!is.na(fdr) & fdr <= thr)
  cat(sprintf("%s (FDR ≤ %.0e): called=%d (%.2f%% of %d)\n",
              nm, thr, n_called, 100*n_called/length(fdr), length(fdr)))
}

for (nm in names(sens)) {
  thr <- sens[[nm]]
  n_blue <- min(as.integer(targets[[nm]]), length(.barcodes_by_rank))
  called <- .barcodes_by_rank[seq_len(n_blue)]

  .knee_plot(
    raw_counts = raw,
    called_barcodes = called,
    out_png = file.path(out_dir, sprintf("%s_%s_knee.png", sample_id, nm)),
    title   = sprintf("Sensitivity %s", gsub("S", "", nm)) 
  )
}

cat("Wrote knee plots to: ", normalizePath(out_dir), "\n")

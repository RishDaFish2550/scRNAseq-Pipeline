### Harish Jawahar Nov 4 2025
### Work with raw tsv

library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(biomaRt)

library(optparse)

# Define command-line options
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
    help="input TSV file containing gene expression matrix"),
  make_option(c("-n", "--name"), type="character", default="pipseq",
    help="project or sample name [default is 'pipseq']"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
    help="save path for the Seurat object RDS file")
)
parser <- OptionParser(option_list=option_list)
args <- parse_args(parser)

# Extract arguments
projectname <- args$name
input_file <- args$input
output_dir <- args$output

# Debugging prints
print(projectname)
print(input_file)

# Read the input TSV file
print("Reading input TSV file...")
data <- read.table(input_file, header=TRUE, row.names=1, sep="\t", check.names=FALSE)

# Convert the data frame to a matrix
counts <- as.matrix(data)

# Create a Seurat object
print("Creating Seurat object...")
seurat_obj <- CreateSeuratObject(counts = counts, project=projectname, assay="RNA")

# Label cells in mitochondria and their cell cycle phase
print("Identify percent cells in mitochondria")
seurat_obj <- PercentageFeatureSet(object = seurat_obj, pattern = "^MT-|^mt-|Mt-", col.name = "Percent_Mito")
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
seurat_obj <- NormalizeData(seurat_obj, assay = "RNA")
seurat_obj <- CellCycleScoring(seurat_obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE, assay = "RNA")

# Save the Seurat object
output_file <- paste(output_dir,"/", projectname, ".rds", sep="")
SaveSeuratRds(seurat_obj, file = output_file)
print(paste("Seurat Object Saved:", output_file, sep=""))


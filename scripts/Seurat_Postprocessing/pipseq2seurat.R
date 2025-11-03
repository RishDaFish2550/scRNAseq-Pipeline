### Thuy N. Nguyen June 3, 2025
### Work with raw matricies

library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(biomaRt)

library(optparse)

option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL, # Changed action and added default=NULL
    help="input directory"),
  make_option(c("-p", "--project"), type="character", default="pipseq",
    help="project or sample name [default is 'pipseq']"),
  make_option(c("-s", "--sensitivity"), type="character", default="2",
    help="sensitivity for cell calling (default = 2)")
)
parser <- OptionParser(option_list=option_list)
args <- parse_args(parser) # Removed the testing arguments

projectname<-args$project
path2pipseq <-args$input
sensitivity<-args$sensitivity
print(projectname)
print(path2pipseq)

counts <- ReadMtx(mtx = paste(path2pipseq,"/filtered_matrix/sensitivity_",sensitivity,"/matrix.mtx.gz",sep=""), 
                  features = paste(path2pipseq,"/filtered_matrix/sensitivity_",sensitivity,"/features.tsv.gz",sep=""), 
                  cells = paste(path2pipseq,"/filtered_matrix/sensitivity_",sensitivity,"/barcodes.tsv.gz",sep=""))

seurat_obj <- CreateSeuratObject(counts = counts, project=projectname,assay="RNA")

#label cells in mitochondria and their cell cycle phase
print("Identify percent cells in mitochonrdia")
seurat_obj<-PercentageFeatureSet(object = seurat_obj, pattern = "^MT-|^mt-|Mt-", col.name = "Percent_Mito")
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
seurat_obj <- NormalizeData(seurat_obj, assay = "RNA")
seurat_obj <- CellCycleScoring(seurat_obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE, assay = "RNA")
path2outputseurat<-paste("/data/250612_pipseeker_seurat_dunnartAnnotJune2025/rds/",projectname,"_pipseq_seurat_object.rds",sep="")
SaveSeuratRds(seurat_obj, file = path2outputseurat)
print(paste("Seurat Object Saved:",path2outputseurat,sep=""))

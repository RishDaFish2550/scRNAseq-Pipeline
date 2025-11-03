# the goal here is to take each of the 8 subsetted RDS files and do SCT transform, dimension reduction, and clustering
### 
## Rscript runSCTransformAndCluster.R -i subsetted.rds -p nameOfCocktailLikeOSKMLABC
library(Seurat)
library(SeuratDisk)
library(tidyverse)
library(sceasy)
library(reticulate)
library(glmGamPoi)
library(ggplot2)
library(gridExtra)

#scObj <- seurat_obj_0<-LoadSeuratRds("250401_merged_seuratObject.rds")
options(future.globals.maxSize = 220000 * 1024^2)

library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(biomaRt)

library(optparse)

option_list <- list(
  make_option(c("-c", "--cocktail"), type="character", default="pipseq",
    help="project or sample name [default is 'pipseq']")
)
parser <- OptionParser(option_list=option_list)
args <- parse_args(parser) # Removed the testing arguments

cocktail<-args$cocktail

#cocktails = c("P2SKML", "PSKML", "PSKMNL", "PSKNL",  "OSKML",  "OSKMNL", "OSKNL")
print("SCTransform, PCA, UMAP, clustering")
figureprefix=paste("FIGURE_SUBSET_",cocktail,"_",sep="")

#seurat_obj <- LoadSeuratRds(paste("/data/250612_pipseeker_seurat_dunnartAnnotJune2025/subset_rds/seuratSubsets_",cocktails[i],".rds",sep=""))
seurat_obj <- LoadSeuratRds(paste("/data/250612_pipseeker_seurat_dunnartAnnotJune2025/subset_rds/seuratSubsets_",cocktail,".rds",sep=""))
seurat_obj <- SCTransform(seurat_obj, method = "glmGamPoi",
              assay = "RNA",
              vars.to.regress = c("nFeature_RNA", "nCount_RNA", "Percent_Mito", 'S.Score', 'G2M.Score'),
              variable.features.n = 3000)
seurat_obj <- RunPCA(seurat_obj, assay = "SCT",reduction.name = "pca",npcs = 30)
seurat_obj <- RunTSNE(seurat_obj, assay="SCT", reduction="pca", reduction.name="tsne",dims = 1:20) # Adjust dims as needed
seurat_obj<-RunUMAP(seurat_obj, assay = "SCT",reduction = "pca",reduction.name = "umap",dims = 1:20)
seurat_obj <- FindNeighbors(seurat_obj, assay = "SCT",reduction = "pca",dims = 1:20)
seurat_obj <- FindClusters(seurat_obj, resolution = seq(0.4, 1.2, 0.1),algorithm = 4)
print("Generate DimPlots of unintegrated samples")
#DimensionReductionFigures 
p_pca<- DimPlot(seurat_obj, reduction="pca", group.by="cell_type_group")
p_tsne<-DimPlot(seurat_obj, reduction="tsne", group.by="cell_type_group")
ggsave(paste(figureprefix,"pca_unintegrated.png",sep=""),plot=p_pca,dpi=220,width=12,height=8)
ggsave(paste(figureprefix,"tsne_unintegrated.png",sep=""),plot=p_tsne,dpi=220,width=12,height=8)

p1 <- DimPlot(seurat_obj, reduction = 'umap', group.by = 'cell_type_group')
p2 <- DimPlot(seurat_obj, reduction = 'umap', group.by = 'treatment')
p3 <- DimPlot(seurat_obj, reduction = 'umap', group.by = 'seurat_clusters')
p4 <- DimPlot(seurat_obj, reduction = 'umap', group.by = 'SCT_snn_res.0.4')
p5<-grid.arrange(p1,p2,p3,p4, ncol = 2, nrow = 2)
ggsave(paste(figureprefix,"umaps_unintegrated.png",sep=""),
plot=p5,dpi=220,width=12,height=8)
SaveSeuratRds(seurat_obj, file = paste("250606_seuratSubset_unintegrated_",cocktail,".rds",sep=""))

seurat_obj <- LoadSeuratRds(paste("250606_seuratSubset_unintegrated_",cocktail,".rds",sep=""))
print("Integration among species")
merged_seurat_integrated<-IntegrateLayers(seurat_obj,
                                        orig.reduction = "pca",
                                        new.reduction = "integrated_rpca",
                                        method = RPCAIntegration,
                                        normalization.method = "SCT",
                                        k.weight=5,
                                        dims=1:20)
merged_seurat_integrated<-RunTSNE(merged_seurat_integrated, 
                                assay="SCT",
                                reduction="integrated_rpca",
                                reduction.name="tsne_integrated",
                                dims=1:20)

merged_seurat_integrated<-RunUMAP(merged_seurat_integrated, 
                                assay="SCT",
                                reduction="integrated_rpca",
                                reduction.name="umap_integrated",
                                dims=1:20)

merged_seurat_integrated<-FindNeighbors(merged_seurat_integrated, 
                                  assay="SCT",
                                  reduction="integrated_rpca",
                                  dims=1:20)
merged_seurat_integrated<-FindClusters(merged_seurat_integrated,
                                      resolution = seq(0.4, 1.2, 0.1),
                                      algorithm = 4)
#SaveSeuratRds(merged_seurat_integrated, file = paste("250606_seuratSubset_",cocktails[i],".rds",sep=""))
SaveSeuratRds(merged_seurat_integrated, file = paste("250606_seuratSubset_",cocktail,".rds",sep=""))
# integrated
p1 <- DimPlot(merged_seurat_integrated, reduction = 'umap_integrated', group.by = 'cell_type_group')
p2 <- DimPlot(merged_seurat_integrated, reduction = 'umap_integrated', group.by = 'treatment')
p3 <- DimPlot(merged_seurat_integrated, reduction = 'umap_integrated', group.by = 'seurat_clusters')
p4 <- DimPlot(merged_seurat_integrated, reduction = 'umap_integrated', group.by = 'SCT_snn_res.0.4')
p5 <-grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
ggsave(paste(figureprefix,"umaps_integrated.png",sep=""),plot=p5,dpi=220,width=12,height=8)
#SaveSeuratRds(merged_seurat_integrated, file = paste("250606_seuratSubset_integrated",cocktails[i],".rds",sep=""))
SaveSeuratRds(merged_seurat_integrated, file = paste("250606_seuratSubset_integrated",cocktail,".rds",sep=""))
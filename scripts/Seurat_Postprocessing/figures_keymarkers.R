
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(gridExtra)
library(patchwork)
library(tidyverse)
library(biomaRt)
library(sceasy)
library(reticulate)
library(glmGamPoi)
py_install("leidenalg")
leidenalg<- import("leidenalg")
#options(future.globals.maxSize = 3e+09)
options(future.globals.maxSize = 4500000 * 1024^2)
library(ggplot2)
library(gridExtra)

## key markers list 
keymarkers<-c('PTPRZ1', 'EVX1', 'DPPA3', 'NCOA5', 'VGLL1', 'TRIM60', 
        'HORMAD1', 'NANOG', 'LEF1', 'OCT4', 'FGF4', 'XIST', 'HMX2', 'DNMT3L', 'CYTL1',
        'KLF4', 'DUSP6', 'MESP1', 'NLRP7', 'GREB1', 'SOX2', 'EID1', 'GATA3', 'CLDN4',
        'DPPA5', 'GATA2', 'H2AFY2', 'GATA6', 'POU5F1', 'T', 'THY1', 'ABCG2', 'SOX11',
        'KLF17', 'PDGFRA', 'SNAI1', 'ZFP729', 'SOX17', 'ESRRB', 'ERP27')

cocktails = c("P2SKML", "PSKML", "PSKMNL", "PSKNL",  "OSKML",  "OSKMNL", "OSKNL")
for (i in 1:length(cocktails)){
    cocktail=cocktails[i]
    print(cocktails[i])
    print("SCTransform, PCA, UMAP, clustering")
    figureprefix=paste("FIGURE_SUBSET_KEYMARKERS_",cocktail,"_",sep="")
    merged_seurat_integrated <- LoadSeuratRds(paste("/data/250612_pipseeker_seurat_dunnartAnnotJune2025/subset_rds/250606_seuratSubset_integrated",cocktails[i],".rds",sep=""))

    # Ensure the key markers are present in the data        
    keymarkers.present <- intersect(keymarkers, rownames(merged_seurat_integrated@assays$SCT@data)) 
    ridgeplot<-RidgePlot(merged_seurat_integrated, features = keymarkers.present, ncol = 2,log=FALSE,group.by='cell_type_group')+ggtitle(cocktail)
    ggsave(paste(figureprefix,"_ridgeplot.png",sep=""),plot=ridgeplot,dpi=220,height=12,width=20)

    vplot <- VlnPlot(merged_seurat_integrated, features = keymarkers.present,group.by='cell_type_group',same.y.lims=TRUE,ncol=8)+ggtitle(cocktail)
    ggsave(paste(figureprefix,"_violin.png",sep=""),plot=vplot,dpi=220,width=20,height=8)

    dotplotM<-DotPlot(merged_seurat_integrated, features = keymarkers.present,group.by='cell_type_group')+RotatedAxis()+ggtitle(cocktail)
    ggsave(paste(figureprefix,"_DotPlot.png",sep=""),plot=dotplotM,dpi=220)

    # Single cell heatmap of feature expression
    pheat<-DoHeatmap(merged_seurat_integrated, features = keymarkers.present, size = 3)+ggtitle(cocktail)
    ggsave(paste(figureprefix,"_heatmap.png",sep=""),plot=pheat,dpi=220)
}

### these violin plots look terrible - put them together 

keymarkers<-c('PTPRZ1', 'EVX1', 'DPPA3', 'NCOA5', 'VGLL1', 'TRIM60', 
        'HORMAD1', 'NANOG', 'LEF1', 'OCT4', 'FGF4', 'XIST', 'HMX2', 'DNMT3L', 'CYTL1',
        'KLF4', 'DUSP6', 'MESP1', 'NLRP7', 'GREB1', 'SOX2', 'EID1', 'GATA3', 'CLDN4',
        'DPPA5', 'GATA2', 'H2AFY2', 'GATA6', 'POU5F1', 'T', 'THY1', 'ABCG2', 'SOX11',
        'KLF17', 'PDGFRA', 'SNAI1', 'ZFP729', 'SOX17', 'ESRRB', 'ERP27')

sc_objects<-NA

for (i in 1:length(cocktails)){
    cocktail=cocktails[i]
    print(cocktails[i])
    print("SCTransform, PCA, UMAP, clustering")
    figureprefix=paste("FIGURE_SUBSET_KEYMARKERS_",cocktail,"_",sep="")
    sc_objects.cocktail <- LoadSeuratRds(paste("/data/250612_pipseeker_seurat_dunnartAnnotJune2025/subset_rds/250606_seuratSubset_integrated",cocktails[i],".rds",sep=""))
cocktails = c("P2SKML", "PSKML", "PSKMNL", "PSKNL",  "OSKML",  "OSKMNL", "OSKNL")

cocktail="P2SKML"
sc_object.P2SKML <- LoadSeuratRds(paste("/data/250612_pipseeker_seurat_dunnartAnnotJune2025/subset_rds/250606_seuratSubset_integrated",cocktail,".rds",sep=""))

cocktail="PSKML"
sc_object.PSKML <- LoadSeuratRds(paste("/data/250612_pipseeker_seurat_dunnartAnnotJune2025/subset_rds/250606_seuratSubset_integrated",cocktail,".rds",sep=""))

cocktail="PSKMNL"
sc_object.PSKMNL <- LoadSeuratRds(paste("/data/250612_pipseeker_seurat_dunnartAnnotJune2025/subset_rds/250606_seuratSubset_integrated",cocktail,".rds",sep=""))

cocktail="P2SKML"
sc_object.P2SKML <- LoadSeuratRds(paste("/data/250612_pipseeker_seurat_dunnartAnnotJune2025/subset_rds/250606_seuratSubset_integrated",cocktail,".rds",sep=""))

cocktail="PSKNL"
sc_object.PSKNL <- LoadSeuratRds(paste("/data/250612_pipseeker_seurat_dunnartAnnotJune2025/subset_rds/250606_seuratSubset_integrated",cocktail,".rds",sep=""))

cocktail="OSKML"
sc_object.OSKML <- LoadSeuratRds(paste("/data/250612_pipseeker_seurat_dunnartAnnotJune2025/subset_rds/250606_seuratSubset_integrated",cocktail,".rds",sep=""))

cocktail="OSKMNL"
sc_object.OSKMNL <- LoadSeuratRds(paste("/data/250612_pipseeker_seurat_dunnartAnnotJune2025/subset_rds/250606_seuratSubset_integrated",cocktail,".rds",sep=""))

cocktail="OSKNL"
sc_object.OSKNL <- LoadSeuratRds(paste("/data/250612_pipseeker_seurat_dunnartAnnotJune2025/subset_rds/250606_seuratSubset_integrated",cocktail,".rds",sep=""))

cocktail="OSKNL"

seurat_obj_list <- list(
"OSKML"=sc_object.OSKML,
"OSKMNL"=sc_object.OSKMNL,
"OSKNL"=sc_object.OSKNL,
"P2SKML"=sc_object.P2SKML,
"PSKML"=sc_object.PSKML,
"PSKMNL"=sc_object.PSKMNL,
"PSKNL"=sc_object.PSKNL)


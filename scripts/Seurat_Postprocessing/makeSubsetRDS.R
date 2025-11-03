# the goal here is to generate 8 subsetted RDS files. One for each cocktail treatment

library(Seurat)
library(SeuratDisk)
library(tidyverse)
library(sceasy)
library(reticulate)
library(glmGamPoi)
#scObj <- seurat_obj_0<-LoadSeuratRds("250401_merged_seuratObject.rds")
options(future.globals.maxSize = 320000 * 1024^2)
merged_seurat<-LoadSeuratRds("/data/250612_pipseeker_seurat_dunnartAnnotJune2025/rds/250612_merged_seuratObject.rds")

cocktails = unique(merged_seurat$treatment)
for (i in 1:length(cocktails)){
    cocktail = cocktails[i]
    print(cocktail)
    scobj_subset <- subset(merged_seurat, subset= treatment==cocktail & cell_type_group!="dEB")
    SaveSeuratRds(scobj_subset, file=paste("seuratSubsets_",cocktail,".rds",sep=''))
}

scobj_subset <- subset(merged_seurat, subset= species=="mouse" & cell_type_group!="dEB")
SaveSeuratRds(scobj_subset, file="seuratSubsets_species_mouse.rds")
scobj_subset <- subset(merged_seurat, subset= species=="dunnart" & cell_type_group!="dEB")
SaveSeuratRds(scobj_subset, file="seuratSubsets_species_dunnart.rds")


#cocktail_pairs
# OSKML P2SKML
# OSKMNL PSKMNL
# OSKNL PSKNL
cocktail = "OSKML"
cocktail2 = "P2SKML"
scobj_subset <- subset(merged_seurat, subset= treatment==cocktail &treatment==cocktail2 & cell_type_group!="dEB")
SaveSeuratRds(scobj_subset, file=paste("seuratSubsets_TWO_",cocktail,"and",cocktail2,".rds",sep=''))

cocktail = "OSKMNL"
cocktail2 = "PSKMNL"
scobj_subset <- subset(merged_seurat, subset= treatment==cocktail &treatment==cocktail2 & cell_type_group!="dEB")
SaveSeuratRds(scobj_subset, file=paste("seuratSubsets_TWO_",cocktail,"and",cocktail2,".rds",sep=''))

cocktail = "OSKNL"
cocktail2 = "PSKNL"
scobj_subset <- subset(merged_seurat, subset= treatment==cocktail &treatment==cocktail2 & cell_type_group!="dEB")
SaveSeuratRds(scobj_subset, file=paste("seuratSubsets_TWO_",cocktail,"and",cocktail2,".rds",sep=''))

scobj_subset <- subset(merged_seurat, subset= species=="mouse" & cell_type_group!="dEB")
SaveSeuratRds(scobj_subset, file="seuratSubsets_species_mouse.rds")
scobj_subset <- subset(merged_seurat, subset= species=="dunnart" & cell_type_group!="dEB")
SaveSeuratRds(scobj_subset, file="seuratSubsets_species_dunnart.rds")

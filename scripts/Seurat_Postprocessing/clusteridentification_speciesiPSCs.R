#TEST_clusteridentification.R

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


#    cocktail="P2SKML"
cocktails = c("P2SKML", "PSKML", "PSKMNL", "PSKNL",  "OSKML",  "OSKMNL")#, "OSKNL")

for (i in 1:length(cocktails)){
    cocktail=cocktails[i]
    cocktail = "OSKNL"
    sc_subset <- LoadSeuratRds(paste("/data/250612_pipseeker_seurat_dunnartAnnotJune2025/subset_rds/250606_seuratSubset_integrated",cocktail,".rds",sep=""))

    number_clusters <-unique(sc_subset$SCT_snn_res.0.4)
    # shows that there are 5 clusters
    # how much of each mipsc and dipsc are in there? 
    counts_table <- table(sc_subset$SCT_snn_res.0.4, sc_subset$species)
    counts_df <- as.data.frame(as.table(counts_table)) # as.table is important for 2-way table
    colnames(counts_df) <- c("Cluster", "CellType", "Count")
    # total number of cells is 
    total_cells<-ncol(sc_subset)
    counts_df$percentTotalCells<-100*counts_df$Count/total_cells

    write.csv(counts_df, paste("TABLE_cluster_celltype_counts",cocktail,".csv",sep=""), row.names = FALSE)
    }
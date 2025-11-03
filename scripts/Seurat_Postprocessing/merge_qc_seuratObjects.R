# first merge all of your seurat objects 

################################################################################
#mergeSeuratObjects.R
library(Seurat)
library(SeuratDisk)
library(tidyverse)
library(biomaRt)
library(sceasy)
library(reticulate)
library(glmGamPoi)
library(ggplot2)
library(gridExtra)

#scObj <- seurat_obj_0<-LoadSeuratRds("250401_merged_seuratObject.rds")
options(future.globals.maxSize = 2300000 * 1024^2)

seurat_obj_0<-LoadSeuratRds("/data/250612_pipseeker_seurat_dunnartAnnotJune2025/rds/dEB_P2SKML_0_pipseq_seurat_object.rds")
seurat_obj_1<-LoadSeuratRds("/data/250612_pipseeker_seurat_dunnartAnnotJune2025/rds/dEB_PSKML_0_pipseq_seurat_object.rds")
seurat_obj_2<-LoadSeuratRds("/data/250612_pipseeker_seurat_dunnartAnnotJune2025/rds/dEB_PSKMNL_0_pipseq_seurat_object.rds")
seurat_obj_3<-LoadSeuratRds("/data/250612_pipseeker_seurat_dunnartAnnotJune2025/rds/dEB_PSKNL_0_pipseq_seurat_object.rds")
seurat_obj_4<-LoadSeuratRds("/data/250612_pipseeker_seurat_dunnartAnnotJune2025/rds/diPSC_OSKML_1_pipseq_seurat_object.rds")
seurat_obj_5<-LoadSeuratRds("/data/250612_pipseeker_seurat_dunnartAnnotJune2025/rds/diPSC_OSKML_2_pipseq_seurat_object.rds")
seurat_obj_6<-LoadSeuratRds("/data/250612_pipseeker_seurat_dunnartAnnotJune2025/rds/diPSC_OSKMNL_1_pipseq_seurat_object.rds")
seurat_obj_7<-LoadSeuratRds("/data/250612_pipseeker_seurat_dunnartAnnotJune2025/rds/diPSC_OSKNL_1_pipseq_seurat_object.rds")
seurat_obj_8<-LoadSeuratRds("/data/250612_pipseeker_seurat_dunnartAnnotJune2025/rds/diPSC_P2SKML_1_pipseq_seurat_object.rds")
seurat_obj_9<-LoadSeuratRds("/data/250612_pipseeker_seurat_dunnartAnnotJune2025/rds/diPSC_P2SKML_2_pipseq_seurat_object.rds")
seurat_obj_10<-LoadSeuratRds("/data/250612_pipseeker_seurat_dunnartAnnotJune2025/rds/diPSC_P2SKML_3_pipseq_seurat_object.rds")
seurat_obj_11<-LoadSeuratRds("/data/250612_pipseeker_seurat_dunnartAnnotJune2025/rds/diPSC_PSKML_1_pipseq_seurat_object.rds")
seurat_obj_12<-LoadSeuratRds("/data/250612_pipseeker_seurat_dunnartAnnotJune2025/rds/diPSC_PSKML_2_pipseq_seurat_object.rds")
seurat_obj_13<-LoadSeuratRds("/data/250612_pipseeker_seurat_dunnartAnnotJune2025/rds/diPSC_PSKML_3_pipseq_seurat_object.rds")
seurat_obj_14<-LoadSeuratRds("/data/250612_pipseeker_seurat_dunnartAnnotJune2025/rds/diPSC_PSKMNL_1_pipseq_seurat_object.rds")
seurat_obj_15<-LoadSeuratRds("/data/250612_pipseeker_seurat_dunnartAnnotJune2025/rds/diPSC_PSKMNL_2_pipseq_seurat_object.rds")
seurat_obj_16<-LoadSeuratRds("/data/250612_pipseeker_seurat_dunnartAnnotJune2025/rds/diPSC_PSKMNL_3_pipseq_seurat_object.rds")
seurat_obj_17<-LoadSeuratRds("/data/250612_pipseeker_seurat_dunnartAnnotJune2025/rds/diPSC_PSKNL_1_pipseq_seurat_object.rds")
seurat_obj_18<-LoadSeuratRds("/data/250612_pipseeker_seurat_dunnartAnnotJune2025/rds/diPSC_PSKNL_2_pipseq_seurat_object.rds")
seurat_obj_19<-LoadSeuratRds("/data/250612_pipseeker_seurat_dunnartAnnotJune2025/rds/diPSC_PSKNL_3_pipseq_seurat_object.rds")
seurat_obj_20<-LoadSeuratRds("/data/250612_pipseeker_seurat_dunnartAnnotJune2025/rds/diPSC_PSKNL_4_pipseq_seurat_object.rds")
seurat_obj_21<-LoadSeuratRds("/data/250612_pipseeker_seurat_dunnartAnnotJune2025/rds/miPSC_OSKML_0_pipseq_seurat_object.rds")
seurat_obj_22<-LoadSeuratRds("/data/250612_pipseeker_seurat_dunnartAnnotJune2025/rds/miPSC_OSKMNL_0_pipseq_seurat_object.rds")
seurat_obj_23<-LoadSeuratRds("/data/250612_pipseeker_seurat_dunnartAnnotJune2025/rds/miPSC_OSKNL_1_pipseq_seurat_object.rds")
seurat_obj_24<-LoadSeuratRds("/data/250612_pipseeker_seurat_dunnartAnnotJune2025/rds/miPSC_OSKNL_2_pipseq_seurat_object.rds")
seurat_obj_25<-LoadSeuratRds("/data/250612_pipseeker_seurat_dunnartAnnotJune2025/rds/miPSC_P2SKML_0_pipseq_seurat_object.rds")
seurat_obj_26<-LoadSeuratRds("/data/250612_pipseeker_seurat_dunnartAnnotJune2025/rds/miPSC_PSKML_0_pipseq_seurat_object.rds")
seurat_obj_27<-LoadSeuratRds("/data/250612_pipseeker_seurat_dunnartAnnotJune2025/rds/miPSC_PSKMNL_0_pipseq_seurat_object.rds")
seurat_obj_28<-LoadSeuratRds("/data/250612_pipseeker_seurat_dunnartAnnotJune2025/rds/miPSC_PSKNL_0_pipseq_seurat_object.rds")

merged_seurat <- merge(seurat_obj_0, c(seurat_obj_1, seurat_obj_2, seurat_obj_3, seurat_obj_4, seurat_obj_5, seurat_obj_6, seurat_obj_7,
 seurat_obj_8, seurat_obj_9, seurat_obj_10, seurat_obj_11, seurat_obj_12, seurat_obj_13, seurat_obj_14,
 seurat_obj_15, seurat_obj_16, seurat_obj_17, seurat_obj_18, seurat_obj_19, seurat_obj_20, seurat_obj_21,
 seurat_obj_22, seurat_obj_23,seurat_obj_24,seurat_obj_25,seurat_obj_26,seurat_obj_27, seurat_obj_28))
################################################################################
path2figures = "/data/250612_pipseeker_seurat_dunnartAnnotJune2025/FIGURES/"
## Adding metadata ### 
merged_seurat$sample <- rownames(merged_seurat@meta.data)
merged_seurat$cell_type_group <- NA
merged_seurat$cell_type_group[grepl("^dEB", merged_seurat$orig.ident)] <- "dEB"
merged_seurat$cell_type_group[grepl("^diPSC", merged_seurat$orig.ident)] <- "diPSC"
merged_seurat$cell_type_group[grepl("^miPSC", merged_seurat$orig.ident)] <- "miPSC"

merged_seurat$cell_stage<- NA
merged_seurat$cell_stage[grepl("^iPSC", merged_seurat$orig.ident)] <- "iPSC"
merged_seurat$cell_stage[grepl("^dEB", merged_seurat$orig.ident)] <- "EB"

merged_seurat$species<-NA
merged_seurat$species[grepl("^diPSC",merged_seurat$orig.ident)] <- "dunnart"
merged_seurat$species[grepl("^miPSC",merged_seurat$orig.ident)] <- "mouse"
merged_seurat$species[grepl("^dEB",merged_seurat$orig.ident)] <- "dunnart"

# label each cell with the stem cell factor cocktail it was treated with
merged_seurat$treatment <- merged_seurat$orig.ident %>% sub(pattern = "^[^_]+_", replacement="", x = .) %>% sub(pattern = "_[0-9]+$", replacement = "", x = .)

#### visualizing feature counts ####
# now plot featuresxcounts and percent mitochondria in violin plots
violinplot_before <- VlnPlot(merged_seurat, features = c("nFeature_RNA", "nCount_RNA", "Percent_Mito"), ncol = 3)+ggtitle("before filter")
#ggsave(paste(path2figures,"FIGURE_pipseq_violin_beforeFilter.png"sep=""), plot=violinplot, dpi=220, width=16, height=12)
# create scatterplot between counts and features of every cell
scatterplot_before<- FeatureScatter(merged_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')+ggtitle("before filter")
 # ggsave(paste(path2figures,"FIGURE_pipseq_scatter_beforeFilter.png",sep=""),dpi=220, width=5, height=5)
# visually inspect both plots to determine cutoffs for filtering 

#merged_seurat <- subset(merged_seurat, subset = nFeature_RNA < 10000 
#                                      & nCount_RNA < 50000 
                                      #& Percent_Mito < 5)

#now plot featuresxcounts and percent mitochondria in violin plots
violinplot <- VlnPlot(merged_seurat, features = c("nFeature_RNA", "nCount_RNA", "Percent_Mito"), ncol = 3)+ggtitle("after filter")
#ggsave(paste(path2figures,"FIGURE_pipseqSensitivity5_violin_afterFilter.png",sep=""), plot=violinplot, dpi=220, width=16, height=12)
# create scatterplot between counts and features of every cell
scatterplot<- FeatureScatter(merged_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+ggtitle("after filter")
  geom_smooth(method = 'lm')
  #ggsave(paste(path2figures,"FIGURE_pipseqSensitivity5_scatter_afterFilter.png",sep=""), dpi=220, width=5, height=5)

qcplots<-grid.arrange(violinplot_before,scatterplot_before,violinplot,scatterplot, nrow = 2)
ggsave(paste(path2figures,"FIGURE_QC_grid_violin_scatter.png",sep=""),plot=qcplots, dpi=400, width=12, height=)
# visually inspect both plots to determine cutoffs for filtering 
SaveSeuratRds(merged_seurat, file="/data/250612_pipseeker_seurat_dunnartAnnotJune2025/rds/250612_merged_seuratObject.rds")

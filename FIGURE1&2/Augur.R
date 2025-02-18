rm(list=ls())
library(stringr)
library(Seurat)
library(Augur)
library(dplyr)
library(patchwork)
library(viridis)
library(qs)
library(BiocParallel)
register(MulticoreParam(workers = 4,progressbar = TRUE))
#########################################################
library(rtracklayer)
library(tibble)
library(gridExtra)
library(ggplot2)
library(Hmisc)
library(cowplot)
library(zoo)
library(scater)
library(stringr)
library(corrplot)
library(matrixStats)
library(eulerr)
library(viridis)
library(pvclust)
library(parallel)
library(dendextend)
library(future)
plan("multicore", workers = 4)
options(future.globals.maxSize= 62914560000) # 60GB (60000*1024^2)
library(SingleCellExperiment)
library(Matrix)
#####################################################
#transfer H5AD into SEURAT OBJECT
library(scater)
library(scran)
library(dplyr)
library(zellkonverter)
#BiocManager::install("scran")
#BiocManager::install("zellkonverter")
library(loomR)
library(SeuratDisk)
library(data.table)
#devtools::install_github("Bioconductor/MatrixGenerics")
#devtools::install_github("const-ae/sparseMatrixStats")
#devtools::install_github("neurorestore/Augur")
#
#install.packages("D:/R/R-4.1.3/Augur_1.0.3.tar.gz", repos = NULL,
#                 type = "source",INSTALL_opts="--no-multiarch")
############################################################
getwd()
setwd("D:/chai_WJH/Neomycin_Project/R")

adata_loom <- connect(filename = "P3_ITG_raw_scanvi_SCANPY.loom",
                      mode = "r+",skip.validate = TRUE)
matrix=adata_loom[["matrix"]][,]
matrix=t(matrix)
dim(matrix)
gene = adata_loom$row.attrs$var_names[]
barcode = adata_loom$col.attrs$obs_names[]
meta_data = read.csv('P3_ITG_raw_scanvi_SCANPY_obs.csv',row.names = 1) # as form as dataframe format
meta_feature = read.csv('P3_ITG_raw_scanvi_SCANPY_var.csv',row.names = 1)
colnames(matrix)= barcode
row.names(matrix)= gene
x_umap = adata_loom$col.attrs$X_umap[,]
x_scanvi = adata_loom$col.attrs$X_scANVI[,]
x_pca=adata_loom$col.attrs$X_pca[,]
seurat_object_p3_itg = CreateSeuratObject(counts = matrix,meta.data = meta_data,
                                          project = 'p3_itg_loom',
                                          min.cells = 0, 
                                          min.features = 0)
seurat_object_p3_itg@assays[["RNA"]]@meta.features <- meta_feature
x_umap = t(x_umap)
x_scanvi = t(x_scanvi)
x_pca = t(x_pca)
rownames(x_umap) = barcode
rownames(x_scanvi) = barcode
rownames(x_pca) = barcode
colnames(x_umap) = c('UMAP_1','UMAP_2')
colnames(x_scanvi) = c("scANVI_1","scANVI_2","scANVI_3","scANVI_4","scANVI_5","scANVI_6","scANVI_7","scANVI_8","scANVI_9","scANVI_10","scANVI_11","scANVI_12","scANVI_13","scANVI_14","scANVI_15","scANVI_16","scANVI_17","scANVI_18","scANVI_19","scANVI_20","scANVI_21","scANVI_22","scANVI_23","scANVI_24","scANVI_25","scANVI_26","scANVI_27","scANVI_28","scANVI_29","scANVI_30")
colnames(x_pca) = c('PCA_1','PCA_2','PCA_3','PCA_4','PCA_5','PCA_6','PCA_7','PCA_8',
                    'PCA_9','PCA_10','PCA_11','PCA_12','PCA_13','PCA_14','PCA_15','PCA_16',
                    'PCA_17','PCA_18','PCA_19','PCA_20','PCA_21','PCA_22','PCA_23','PCA_24',
                    'PCA_25','PCA_26','PCA_27','PCA_28','PCA_29','PCA_30','PCA_31','PCA_32',
                    'PCA_33','PCA_34','PCA_35','PCA_36','PCA_37','PCA_38','PCA_39','PCA_40',
                    'PCA_41','PCA_42','PCA_43','PCA_44','PCA_45','PCA_46','PCA_47','PCA_48',
                    'PCA_49','PCA_50')
#,'PCA_51','PCA_52','PCA_53','PCA_54','PCA_55','PCA_56','PCA_57','PCA_58',
#'PCA_59','PCA_60','PCA_61','PCA_62','PCA_63','PCA_64','PCA_65','PCA_66','PCA_67','PCA_68',
#'PCA_69','PCA_70','PCA_71','PCA_72','PCA_73','PCA_74','PCA_75','PCA_76','PCA_77','PCA_78',
#'PCA_79','PCA_80','PCA_81','PCA_82','PCA_83','PCA_84','PCA_85','PCA_86','PCA_87','PCA_88',
#'PCA_89','PCA_90','PCA_91','PCA_92','PCA_93','PCA_94','PCA_95','PCA_96','PCA_97','PCA_98',
#'PCA_99','PCA_100'
sce1 = as.SingleCellExperiment(seurat_object_p3_itg)
sce1@int_colData@listData[["reducedDims"]]@listData[["umap"]] = x_umap
sce1@int_colData@listData[["reducedDims"]]@listData[["scANVI"]] = x_scanvi
sce1@int_colData@listData[["reducedDims"]]@listData[["PCA"]] = x_pca
p3_itg_seurat <- as.Seurat(sce1)
DimPlot(p3_itg_seurat,reduction = "umap",group.by = "cell_type2")
DimPlot(p3_itg_seurat,reduction = "umap",group.by = "treatment")
p3_itg_seurat2=p3_itg_seurat
p3_itg_seurat2 <- RunUMAP(p3_itg_seurat2,reduction = "scANVI",dims = 1:30)
DimPlot(p3_itg_seurat2,reduction = "umap",group.by = "cell_type2")
DimPlot(p3_itg_seurat2,reduction = "umap",group.by = "treatment")
sce1 = as.SingleCellExperiment(p3_itg_seurat2)
sce1@int_colData@listData[["reducedDims"]]@listData[["umap"]] = x_umap
sce1@int_colData@listData[["reducedDims"]]@listData[["scANVI"]] = x_scanvi
sce1@int_colData@listData[["reducedDims"]]@listData[["PCA"]] = x_pca
#######################################augur############################
table(sce1$treatment)
table(sce1$cell_type2)
table(sce1$batch)
dat<-sce1@meta.data
#????һ????????Ϣ,ת????????
#dat$batch<-ifelse(grepl("sample2|sample4|sample6",dat$orig.ident),"1","2")
#????һ?·?????Ϣ??????????????????
#dat$group<-ifelse(grepl("sample[1-3]",dat$orig.ident),"con","treat")
#sce1@meta.data <- dat

augur <- calculate_auc(sce1,
                       cell_type_col="cell_type2",#ϸ????Ⱥ????
                       label_col="treatment",#ʵ??????????
                       n_threads=8)
qsave(augur,"augur.qs")
head(augur$AUC,5)
##A tibble:5??2
#cell_typeauc
#<fct><dbl>
# 1 NK cells 0.822
# 2 Fibroblasts 0.745
# 3 Neutrophils 0.66
# 4 B-cells 0.623
# 5 CD8+ T-cells 0.606
plot_lollipop(augur) +
  geom_segment(aes(xend = cell_type, yend = 0.5), size = 1) +
  geom_point(size = 3, aes(color = cell_type)) + # ???ӵ??Ĵ?С????ɫӳ??
  scale_color_manual(values = c(
    "HC" = "#1f77b4", # ��ɫ
    "HeC" = "#ff7f0e",# ??ɫ
    "IdC" = "#2ca02c", # ??ɫ
    "IPh_IBC" = "#d62728", # ??ɫ
    "KO" = "#9467bd", # ??ɫ
    "PC_DC" = "#8c564b", # ??ɫ
    "RMC" = "#e377c2", # ??ɫ
    "SGN" = "#7f7f7f", # ??ɫ
    "TBC" = "#EEEE00" #
  ))
# ????umap
# rankģʽ
plot_umap(augur,sce1,
          mode = "rank",
          reduction = "umap",
          palette = "cividis",
          cell_type_col = "cell_type2")
# defaultģʽ
plot_umap(augur,sce1,
          mode = "default",
          reduction = "umap",
          palette = "YlGnBu", # "viridis", "plasma", "magma", "inferno"
          cell_type_col = "cell_type2")
#########################################################
#########################################################
rm(list=ls())
getwd()
setwd("D:/chai_WJH/Neomycin_Project/R/new_p21")
adata_loom <- connect(filename = "P21_ITG_raw_scanvi_SCANPY.loom",
                      mode = "r+",skip.validate = TRUE)
matrix=adata_loom[["matrix"]][,]
matrix=t(matrix)
dim(matrix)
gene = adata_loom$row.attrs$var_names[]
barcode = adata_loom$col.attrs$obs_names[]
meta_data = read.csv('P21_ITG_raw_scanvi_SCANPY_obs.csv',row.names = 1) # as form as dataframe format
meta_feature = read.csv('P21_ITG_raw_scanvi_SCANPY_var.csv',row.names = 1)

colnames(matrix)= barcode
row.names(matrix)= gene
x_umap = adata_loom$col.attrs$X_umap[,]
x_scanvi = adata_loom$col.attrs$X_scANVI[,]
x_pca=adata_loom$col.attrs$X_pca[,]
seurat_object_p21_itg = CreateSeuratObject(counts = matrix,meta.data = meta_data,
                                          project = 'p21_itg_loom',
                                          min.cells = 0, 
                                          min.features = 0)

seurat_object_p21_itg@assays[["RNA"]]@meta.features <- meta_feature
x_umap = t(x_umap)
x_scanvi = t(x_scanvi)
x_pca = t(x_pca)
rownames(x_umap) = barcode
rownames(x_scanvi) = barcode
rownames(x_pca) = barcode
colnames(x_umap) = c('UMAP_1','UMAP_2')
colnames(x_scanvi) = c("scANVI_1","scANVI_2","scANVI_3","scANVI_4","scANVI_5","scANVI_6","scANVI_7","scANVI_8","scANVI_9","scANVI_10","scANVI_11","scANVI_12","scANVI_13","scANVI_14","scANVI_15","scANVI_16","scANVI_17","scANVI_18","scANVI_19","scANVI_20","scANVI_21","scANVI_22","scANVI_23","scANVI_24","scANVI_25","scANVI_26","scANVI_27","scANVI_28","scANVI_29","scANVI_30")
colnames(x_pca) = c('PCA_1','PCA_2','PCA_3','PCA_4','PCA_5','PCA_6','PCA_7','PCA_8',
                    'PCA_9','PCA_10','PCA_11','PCA_12','PCA_13','PCA_14','PCA_15','PCA_16',
                    'PCA_17','PCA_18','PCA_19','PCA_20','PCA_21','PCA_22','PCA_23','PCA_24',
                    'PCA_25','PCA_26','PCA_27','PCA_28','PCA_29','PCA_30','PCA_31','PCA_32',
                    'PCA_33','PCA_34','PCA_35','PCA_36','PCA_37','PCA_38','PCA_39','PCA_40',
                    'PCA_41','PCA_42','PCA_43','PCA_44','PCA_45','PCA_46','PCA_47','PCA_48',
                    'PCA_49','PCA_50')
#,'PCA_51','PCA_52','PCA_53','PCA_54','PCA_55','PCA_56','PCA_57','PCA_58',
#'PCA_59','PCA_60','PCA_61','PCA_62','PCA_63','PCA_64','PCA_65','PCA_66','PCA_67','PCA_68',
#'PCA_69','PCA_70','PCA_71','PCA_72','PCA_73','PCA_74','PCA_75','PCA_76','PCA_77','PCA_78',
#'PCA_79','PCA_80','PCA_81','PCA_82','PCA_83','PCA_84','PCA_85','PCA_86','PCA_87','PCA_88',
#'PCA_89','PCA_90','PCA_91','PCA_92','PCA_93','PCA_94','PCA_95','PCA_96','PCA_97','PCA_98',
#'PCA_99','PCA_100'
sce1 = as.SingleCellExperiment(seurat_object_p21_itg)
sce1@int_colData@listData[["reducedDims"]]@listData[["umap"]] = x_umap
sce1@int_colData@listData[["reducedDims"]]@listData[["scANVI"]] = x_scanvi
sce1@int_colData@listData[["reducedDims"]]@listData[["PCA"]] = x_pca
p21_itg_seurat <- as.Seurat(sce1)
DimPlot(p21_itg_seurat,reduction = "umap",group.by = "cell_type2")
DimPlot(p21_itg_seurat,reduction = "umap",group.by = "treatment")
p21_itg_seurat2=p21_itg_seurat
p21_itg_seurat2 <- RunUMAP(p21_itg_seurat2,reduction = "scANVI",dims = 1:30)
DimPlot(p21_itg_seurat2,reduction = "umap",group.by = "cell_type2")
DimPlot(p21_itg_seurat2,reduction = "umap",group.by = "treatment")

sce1 = as.SingleCellExperiment(p21_itg_seurat2)
sce1@int_colData@listData[["reducedDims"]]@listData[["umap"]] = x_umap
sce1@int_colData@listData[["reducedDims"]]@listData[["scANVI"]] = x_scanvi
sce1@int_colData@listData[["reducedDims"]]@listData[["PCA"]] = x_pca

#######################################augur############################
table(sce1$treatment)
table(sce1$cell_type2)
table(sce1$batch)
dat<-p21_itg_seurat2@meta.data
#????һ????????Ϣ,ת????????
#dat$batch<-ifelse(grepl("sample2|sample4|sample6",dat$orig.ident),"1","2")
#????һ?·?????Ϣ??????????????????
#dat$group<-ifelse(grepl("sample[1-3]",dat$orig.ident),"con","treat")
#sce1@meta.data <- dat

augur <- calculate_auc(sce1,
                       cell_type_col="cell_type2",#ϸ????Ⱥ????
                       label_col="treatment",#ʵ??????????
                       n_threads=8)
qsave(augur,"augur.qs")

head(augur$AUC,5)

##A tibble:5??2
#cell_typeauc
#<fct><dbl>
# 1 NK cells 0.822
# 2 Fibroblasts 0.745
# 3 Neutrophils 0.66
# 4 B-cells 0.623
# 5 CD8+ T-cells 0.606

plot_lollipop(augur) +
  geom_segment(aes(xend = cell_type, yend = 0.5), size = 1) +
  geom_point(size = 3, aes(color = cell_type)) + # ???ӵ??Ĵ?С????ɫӳ??
  scale_color_manual(values = c(
    "IHC" = "#1f77b4", # ��ɫ
    "TypeI_SGN" = "#ff7f0e",# ??ɫ
    "IDC" = "#2ca02c", # ??ɫ
    "IBC_IPh_HeC" = "#d62728", # ??ɫ
    "KO" = "#9467bd", # ??ɫ
    "DC_PC" = "#8c564b", # ??ɫ
    "RMC" = "#e377c2", # ??ɫ
    "TypeII_SGN" = "#7f7f7f", # ??ɫ
    "TBC" = "#EEEE00", #
    "Bs" = "#ADD8E6",
    "Fb" = "#008B00",
    "Is" = "#8B864E",
    "MC" = "#CD9B9B",
    "OHC" = "#E0EEEE",
    "Pericytes" = "#00FFFF",
    "Schwann.Cells" = "#000080",
    "Sp/Rt" = "#FFD700"
  ))
# ????umap
# rankģʽ
plot_umap(augur,sce1,
          mode = "rank",
          reduction = "umap",
          palette = "cividis",
          cell_type_col = "cell_type2")
# defaultģʽ
plot_umap(augur,sce1,
          mode = "default",
          reduction = "umap",
          palette = "YlGnBu", # "viridis", "plasma", "magma", "inferno"
          cell_type_col = "cell_type2")


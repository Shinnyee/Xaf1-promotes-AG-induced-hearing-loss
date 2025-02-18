rm(list = ls())
#devtools::install_github("phillipnicol/scDist",INSTALL_opts="--no-multiarch")
library(stringr)
library(Seurat)
library(scDist)
library(dplyr)
library(patchwork)
library(viridis)
library(qs)
library(BiocParallel)
library(loomR)
library(rtracklayer)
library(tibble)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(Hmisc)
library(cowplot)
library(zoo)
library(ggplot2)
#library(loomR)
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

#######################################augur############################
table(p3_itg_seurat2$treatment)
table(p3_itg_seurat2$cell_type2)
table(p3_itg_seurat2$batch)
#dat<-p3_itg_seurat2@meta.data
#增加一下批次信息,转换成数字
#dat$batch<-ifelse(grepl("sample2|sample4|sample6",dat$orig.ident),"1","2")
#增加一下分组信息，这里是随意编造的
#dat$group<-ifelse(grepl("sample[1-3]",dat$orig.ident),"con","treat")
dat <- p3_itg_seurat2@meta.data

dat <- GetAssayData(p3_itg_seurat2,slot = "counts") #其实是需要normalize data，但这样也行

out <- scDist(dat,p3_itg_seurat2@meta.data,
              fixed.effects = "treatment",
              random.effects="batch",
              clusters="cell_type2")

qsave(out,"out.qs")

out$results
DistPlot(out)

distGenes(out, cluster = "HC")
plotBetas(out,cluster = "HC")

FDRDistPlot(out)



#载入了名字空间'Matrix' 1.5-4，但需要的是>= 1.6.0
#载入了名字空间'spatstat.utils' 3.0-2，但需要的是>= 3.0.5
#载入了名字空间'ggrepel' 0.9.3，但需要的是>= 0.9.6
#matrix lme4
#########################################################
#########################################################
rm(list=ls())
getwd()
setwd("D:/chai_WJH/Neomycin_Project/R/new_p21")
set.seed(1126490984)
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

#######################################augur############################
table(p21_itg_seurat2$treatment)
table(p21_itg_seurat2$cell_type2)
table(p21_itg_seurat2$batch)
#dat<-p21_itg_seurat2@meta.data
#增加一下批次信息,转换成数字
#dat$batch<-ifelse(grepl("sample2|sample4|sample6",dat$orig.ident),"1","2")
#增加一下分组信息，这里是随意编造的
#dat$group<-ifelse(grepl("sample[1-3]",dat$orig.ident),"con","treat")
dat <- p21_itg_seurat2@meta.data

dat <- GetAssayData(p21_itg_seurat2,slot = "counts") #其实是需要normalize data，但这样也行

out <- scDist(dat,p21_itg_seurat2@meta.data,
              fixed.effects = "treatment",
              random.effects="sample",
              clusters="cell_type2")
out$results
DistPlot(out)
distGenes(out, cluster = "OHC")
distGenes(out, cluster = "IHC")
plotBetas(out,cluster = "OHC")
plotBetas(out,cluster = "IHC")
FDRDistPlot(out)
qsave(out,"out.qs")









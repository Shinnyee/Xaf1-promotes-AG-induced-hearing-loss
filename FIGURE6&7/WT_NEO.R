remove(list = ls())
#BiocManager::install("clusterProfiler")
#BiocManager::install("org.Mm.eg.db")
#install_github("husson/FactoMineR",INSTALL_opts="--no-multiarch")
library(edgeR)
library(DESeq2)
library(stringr)
library(stringi)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(pheatmap)
library(VennDiagram)
library(RColorBrewer)
library(devtools)
library(clusterProfiler)
library(FactoMineR)
library(factoextra)
library(org.Mm.eg.db)
library(ggrepel)
###设置可用内存
options(future.globals.maxSize = 10 * 1024^3)

########## 1.数据读入及表达矩阵归一化 ###########

setwd ("D:/chai_WJH/Neomycin_Project/BULK_RNA_SEQ")
rawcount <- read.csv("counts_anno_wt_neo.csv",row.names = 1)
table(!duplicated(rawcount[,1]))
rawcount$gene_id <- rownames(rawcount)
###原始矩阵归一化，获得cpm矩阵
rawcount=rawcount[!duplicated(rawcount$gene_id),] 
rawcount <- rawcount[,c(1,2,3,4,5,6,12)]
rawcount <- drop_na(rawcount) #去除含有NA值的行
rawcount <- rawcount[,!colnames(rawcount)%in%c("ENSEMBL","Ensmble","median")]#正式获得基因rawcount矩阵
keep <- rowSums(rawcount>0) >= floor(0.5*ncol(rawcount))
filter_count <- rawcount[keep,] #获得filter_count矩阵
rownames(filter_count) <- filter_count$gene_id
filter_count <- filter_count[,!colnames(filter_count)%in%"gene_id"]
express_cpm <- log2(cpm(filter_count)+ 1)
express_cpm[1:6,1:6] #获得cpm矩阵
#save(express_cpm,filter_count,file="filterGeneCountCpm.Rdata")
#########比较一：CON vs shNR2F2 ########
#filter_count1 <- filter_count[,c(1,2,3,4,5,6)]
#express_cpm1 <- express_cpm[,c(1,2,3,4,5,6)]
#colnames(filter_count1) <- c("CON_1","CON_2",'CON_3',"Oe5b_1","Oe5b_2",'Oe5b_3')#给样本命名
#colnames(express_cpm1) <- colnames(filter_count1)
group1=rep(c("WT","NEO"),each=3)
group_list1=factor(group1,levels = c("WT","NEO"))
table(group_list1)#检查一下组别数量

############## 2.查看样本数据质量，进行质检##########
## 01分别绘制两组的箱线图并且拼接起来
exprSet <- express_cpm
data <- data.frame(expression=c(exprSet),
                   sample=rep(colnames(exprSet),each=nrow(exprSet)))
p1 <- ggplot(data = data,aes(x=sample,y=expression,fill=sample))
p1 <- p1 + geom_boxplot() + theme(axis.text.x = element_text(angle = 90)) + 
  xlab(NULL) + ylab("log2(CPM+1)")+theme_bw()


## 02分别绘制两组比较的PCA图并且拼接起来

## 绘制第一组比较的PCA图
dat <- express_cpm
dat <- as.data.frame(t(dat))
dat <- na.omit(dat)
dat$group_list1 <- group_list1
dat_pca1 <- PCA(dat[,-ncol(dat)], graph = FALSE)#画图仅需数值型数据，去掉最后一列的分组信息
p2 <- fviz_pca_ind(dat_pca1,
                   geom.ind = "point", # 只显示点，不显示文字
                   col.ind = dat$group_list1, # 用不同颜色表示分组
                   palette = c("#00AFBB", "#E7B800"),
                   addEllipses = T, # 是否圈起来，少于4个样圈不起来
                   legend.title = "Groups") + theme_bw()
##绘制第二组比较的PCA图
p1
p2
############ 3.表达矩阵的差异基因分析#############

## 差异分析
exprSet1 <- filter_count
design1 <- model.matrix(~0+rev(factor(group_list1))) 
rownames(design1) <- colnames(exprSet1)
colnames(design1) <- levels(factor(group_list1))
DEG1 <- DGEList(counts=exprSet1,  #构建edgeR的DGEList对象
                group=factor(group_list1))
DEG1$samples$lib.size 
DEG1 <- calcNormFactors(DEG1)
DEG1$samples$norm.factors
## 计算线性模型的参数
DEG1 <- estimateGLMCommonDisp(DEG1,design1)
DEG1 <- estimateGLMTrendedDisp(DEG1, design1)
DEG1 <- estimateGLMTagwiseDisp(DEG1, design1)
##拟合线性模型
fit1 <- glmFit(DEG1, design1)
lrt1 <- glmLRT(fit1, contrast=c(1,-1)) 
##提取过滤差异分析结果
DEG_edgeR1 <- as.data.frame(topTags(lrt1, n=nrow(DEG1)))
head(DEG_edgeR1)
##设定阈值，筛选显著上下调差异基因
p <- 0.05
DEG_edgeR1$regulated <- ifelse(DEG_edgeR1$logFC>0.25&DEG_edgeR1$PValue<p,
                               "up",ifelse(DEG_edgeR1$logFC<(-0.25)&DEG_edgeR1$PValue<p,"down","normal"))
table(DEG_edgeR1$regulated)
write.csv(DEG_edgeR1,file="DEG_WT vs NEO.csv")

############4.差异基因可视化-绘图#############

## 01火山图的绘制及拼图
## 绘制第一组比较的火山图
data1 <- DEG_edgeR1
p3.1 <- ggplot(data=data1, aes(x=logFC, y=-log10(PValue),color=regulated)) + 
  geom_point(alpha=0.5, size=1.8) + 
  theme_set(theme_set(theme_bw(base_size=20))) + 
  xlab("log2FC") + ylab("-log10(Pvalue)") +
  scale_colour_manual(values = c('blue','black','red'))+theme_bw()
p3.1
data1$gene_name <- rownames(data1)
top_genes <- data1 %>%
  arrange(PValue) %>%
  head(30) 
p3.1 <- ggplot(data=data1, aes(x=logFC, y=-log10(PValue), color=regulated)) + 
  geom_point(alpha=0.5, size=1.8) + 
  theme_bw(base_size=20) + 
  xlab("log2FC") + ylab("-log10(Pvalue)") +
  scale_colour_manual(values = c('blue', 'black', 'red')) +
  theme_bw() +
  # 添加文本标签，避免重叠
  geom_text_repel(data=top_genes, aes(label=gene_name), box.padding=0.5, size=4, max.overlaps=30)
p3.1
top_genes <- data1 %>%
  arrange(desc(abs(logFC)), PValue) %>%
  head(30)


#####挑选感兴趣基因在火山图中可视化
## 绘制第一组比较的可视化火山图

library(EnhancedVolcano)
p3.3 <-EnhancedVolcano(data1,
                       lab = rownames(data1),
                       x = 'logFC',
                       y = 'PValue',
                       pCutoff = 0.01,
                       selectLab = c("Pcp4","Dnm3","Adamts13","Calb2","Fgf8","Nr2f2","Atoh1",#OHC/
                                     "Mab21l1","Neurod2","Ntn1","Pax8","Ret"),#sensory development
                       drawConnectors = TRUE)
p3.3
p3.4 <-EnhancedVolcano(data1,
                       lab = rownames(data1),
                       x = 'logFC',
                       y = 'PValue',
                       pCutoff = 0.01,
                       selectLab = c("Mab21l1","Neurod2","Ntn1","Pax8","Ret"),
                       drawConnectors = TRUE)
p3.4

p4.3 <-EnhancedVolcano(data1,
                       lab = rownames(data1),
                       x = 'logFC',
                       y = 'PValue',
                       selectLab = "Xaf1",
                       title = "WTvsNEO-down2006up2863",
                       drawConnectors = TRUE)
p4.3
up_genes <- rownames(DEG_edgeR1[DEG_edgeR1$regulated=="up",])
write.csv(up_genes,file = "WTvsNEO_DEG_up.csv")
down_genes <- rownames(DEG_edgeR1[DEG_edgeR1$regulated=="down",])
write.csv(down_genes,file = "WTvsNEO_DEG_down.csv")
## 02差异基因热图的绘制及拼图
## 绘制第一组比较的差异基因热图
edgeR1_sigGene1 <- DEG_edgeR1[DEG_edgeR1$regulated!="normal",]
edgeR1_sigGene1 <-rownames(edgeR1_sigGene1)
dat1 <- express_cpm[match(edgeR1_sigGene1,rownames(express_cpm)),]
dat1[1:6,1:6]
group1 <- data.frame(group=group_list1)
rownames(group1)=colnames(dat1)
library(pheatmap)
p4.1 <- pheatmap(dat1,scale = "row",show_colnames =T,show_rownames = F, 
                 cluster_cols = T, 
                 annotation_col=group1,
                 main = "edgeR's DEG")
p4.1


##########################################################
library(ComplexHeatmap)
A <- as.matrix(dat1) #将表达矩阵转化为matrix
samples <- rep(c('WT', 'NEO'), c(3, 3)) #定义样本分组信息  
for (i in 1:nrow(A)) A[i, ] <- scale(log(unlist(A[i, ] + 1), 2)) #对数据进行标准化处理
B <- Heatmap(A,#表达矩阵
             col = colorRampPalette(c("navy","white","firebrick3"))(100),#颜色定义
             show_row_names = F,#不展示行名
             top_annotation = HeatmapAnnotation(Group = samples, 
                                                simple_anno_size = unit(2, 'mm'),
                                                col = list(Group = c('WT' = '#00DAE0', 'NEO' = '#FF9289')),
                                                show_annotation_name = FALSE))#分组注释
B
genes <- c("Zbp1",
           "Xaf1", 
           "Xiap",
           "Gsdmd",
           "Ripk1",
           "Ripk3",
           "Mlkl",
           "Btk",
           "Cd14","Caspase3","Caspase7","Caspase9","Bax",
           "CCL5",
           "Tnf",
           "Tnfsf14",
           "Tab2","Ltf", "Lbp", "Nos2", "Gbp2", "Il1b", "Irgm1", "Tnfsf18")
genes <- as.data.frame(genes)
B + rowAnnotation(link = anno_mark(at = which(rownames(A) %in% genes$genes), 
                                   labels = genes$genes, labels_gp = gpar(fontsize = 10)))
######save Rdata###




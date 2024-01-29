library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(cowplot)
library(sctransform)
library(ggsci)
library(umap)
library(limma)
library(viridis)
library(patchwork)
library(SingleR)
library(xlsx)
library(harmony)
library(NMF)
library(irGSEA)
library(pheatmap)
library(RColorBrewer)
library(EnhancedVolcano)
library(CytoTRACE)
library(scCustomize)
HB.genes_total <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ") 

zl1.data <- Read_CellBender_h5_Mat("/cluster/home/yflu/UPS/cellbender/output/zl1/zl1_cellbender_filtered.h5")
zl1 <- CreateSeuratObject(counts = zl1.data, project = "zl1", min.cells = 3, min.features = 500)
d_pred <- read.table("ZL1doublets.txt",sep = ',',header = T)
zl1 <- zl1[,which(colnames(zl1) %in% d_pred$barcode[which(d_pred$predicted_doublets == 'False')])]
HB_m <- match(HB.genes_total,rownames(zl1@assays$RNA))
HB.genes <- rownames(zl1@assays$RNA)[HB_m]
HB.genes <- HB.genes[!is.na(HB.genes)]
zl1[["percent.HB"]] <- PercentageFeatureSet(zl1,features=HB.genes)
zl1[["percent.mt"]] <- PercentageFeatureSet(zl1, pattern = "^MT-")
VlnPlot(zl1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), ncol = 4)
zl1 <- subset(zl1, subset = nFeature_RNA > 250 & nFeature_RNA < 6000 & percent.mt < 20)
zl1 <- NormalizeData(zl1)
zl1 <- FindVariableFeatures(zl1, selection.method = "vst", nfeatures = 8000)
zl1 <- ScaleData(object = zl1, verbose = FALSE)
zl1 <- RunPCA(pc.genes = zl1@var.genes, npcs = 20 ,verbose = FALSE, object = zl1)
saveRDS(zl1,file = 'zl1_aPCA.rds')
gc()
zl2.data <- Read_CellBender_h5_Mat("/cluster/home/yflu/UPS/cellbender/output/zl2/zl2_cellbender_filtered.h5")
zl2 <- CreateSeuratObject(counts = zl2.data, project = "zl2", min.cells = 3, min.features = 500)
d_pred <- read.table("ZL2doublets.txt",sep = ',',header = T)
zl2 <- zl2[,which(colnames(zl2) %in% d_pred$barcode[which(d_pred$predicted_doublets == 'False')])]
HB_m <- match(HB.genes_total,rownames(zl2@assays$RNA))
HB.genes <- rownames(zl2@assays$RNA)[HB_m]
HB.genes <- HB.genes[!is.na(HB.genes)]
zl2[["percent.HB"]] <- PercentageFeatureSet(zl2,features=HB.genes)
zl2[["percent.mt"]] <- PercentageFeatureSet(zl2, pattern = "^MT-")
VlnPlot(zl2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), ncol = 4)
zl2 <- subset(zl2, subset = nFeature_RNA > 250 & nFeature_RNA < 6000 & percent.mt < 20)
zl2 <- NormalizeData(zl2)
zl2 <- FindVariableFeatures(zl2, selection.method = "vst", nfeatures = 8000)
zl2 <- ScaleData(object = zl2, verbose = FALSE)
zl2 <- RunPCA(pc.genes = zl2@var.genes, npcs = 20 ,verbose = FALSE, object = zl2)
saveRDS(zl2,file = 'zl2_aPCA.rds')
gc()
zl9.data <- Read_CellBender_h5_Mat("/cluster/home/yflu/UPS/cellbender/output/zl9/zl9_cellbender_filtered.h5")
zl9 <- CreateSeuratObject(counts = zl9.data, project = "zl9", min.cells = 3, min.features = 500)
d_pred <- read.table("ZL9doublets.txt",sep = ',',header = T)
zl9 <- zl9[,which(colnames(zl9) %in% d_pred$barcode[which(d_pred$predicted_doublets == 'False')])]
HB_m <- match(HB.genes_total,rownames(zl9@assays$RNA))
HB.genes <- rownames(zl9@assays$RNA)[HB_m]
HB.genes <- HB.genes[!is.na(HB.genes)]
zl9[["percent.HB"]] <- PercentageFeatureSet(zl9,features=HB.genes)
zl9[["percent.mt"]] <- PercentageFeatureSet(zl9, pattern = "^MT-")
VlnPlot(zl9, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), ncol = 4)
zl9 <- subset(zl9, subset = nFeature_RNA > 250 & nFeature_RNA < 6000 & percent.mt < 20)
zl9 <- NormalizeData(zl9)
zl9 <- FindVariableFeatures(zl9, selection.method = "vst", nfeatures = 8000)
zl9 <- ScaleData(object = zl9, verbose = FALSE)
zl9 <- RunPCA(pc.genes = zl9@var.genes, npcs = 20 ,verbose = FALSE, object = zl9)
saveRDS(zl9,file = 'zl9_aPCA.rds')
gc()
zl10.data <- Read_CellBender_h5_Mat("/cluster/home/yflu/UPS/cellbender/output/zl10/zl10_cellbender_filtered.h5")
zl10 <- CreateSeuratObject(counts = zl10.data, project = "zl10", min.cells = 3, min.features = 500)
d_pred <- read.table("ZL10doublets.txt",sep = ',',header = T)
zl10 <- zl10[,which(colnames(zl10) %in% d_pred$barcode[which(d_pred$predicted_doublets == 'False')])]
HB_m <- match(HB.genes_total,rownames(zl10@assays$RNA))
HB.genes <- rownames(zl10@assays$RNA)[HB_m]
HB.genes <- HB.genes[!is.na(HB.genes)]
zl10[["percent.HB"]] <- PercentageFeatureSet(zl10,features=HB.genes)
zl10[["percent.mt"]] <- PercentageFeatureSet(zl10, pattern = "^MT-")
VlnPlot(zl10, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), ncol = 4)
zl10 <- subset(zl10, subset = nFeature_RNA > 250 & nFeature_RNA < 6000 & percent.mt < 20)
zl10 <- NormalizeData(zl10)
zl10 <- FindVariableFeatures(zl10, selection.method = "vst", nfeatures = 8000)
zl10 <- ScaleData(object = zl10, verbose = FALSE)
zl10 <- RunPCA(pc.genes = zl10@var.genes, npcs = 20 ,verbose = FALSE, object = zl10)
saveRDS(zl10,file = 'zl10_aPCA.rds')
gc()
zl13.data <- Read_CellBender_h5_Mat("/cluster/home/yflu/UPS/cellbender/output/zl13/zl13_cellbender_filtered.h5")
zl13 <- CreateSeuratObject(counts = zl13.data, project = "zl13", min.cells = 3, min.features = 500)
d_pred <- read.table("ZL13doublets.txt",sep = ',',header = T)
zl13 <- zl13[,which(colnames(zl13) %in% d_pred$barcode[which(d_pred$predicted_doublets == 'False')])]
HB_m <- match(HB.genes_total,rownames(zl13@assays$RNA))
HB.genes <- rownames(zl13@assays$RNA)[HB_m]
HB.genes <- HB.genes[!is.na(HB.genes)]
zl13[["percent.HB"]] <- PercentageFeatureSet(zl13,features=HB.genes)
zl13[["percent.mt"]] <- PercentageFeatureSet(zl13, pattern = "^MT-")
VlnPlot(zl13, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), ncol = 4)
zl13 <- subset(zl13, subset = nFeature_RNA > 250 & nFeature_RNA < 6000 & percent.mt < 20)
zl13 <- NormalizeData(zl13)
zl13 <- FindVariableFeatures(zl13, selection.method = "vst", nfeatures = 8000)
zl13 <- ScaleData(object = zl13, verbose = FALSE)
zl13 <- RunPCA(pc.genes = zl13@var.genes, npcs = 20 ,verbose = FALSE, object = zl13)
saveRDS(zl13,file = 'zl13_aPCA.rds')
gc()

zl1 <- readRDS(file = "D:/R/UPS/Integrated/zl1_aPCA.rds")
zl2 <- readRDS(file = "D:/R/UPS/Integrated/zl2_aPCA.rds")
zl9 <- readRDS(file = "D:/R/UPS/Integrated/zl9_aPCA.rds")
zl10 <- readRDS(file = "D:/R/UPS/Integrated/zl10_aPCA.rds")
zl13 <- readRDS(file = "D:/R/UPS/Integrated/zl13_aPCA.rds")

zl.integrated <- merge(zl1, y = c(zl2, zl9, zl10, zl13), project = "zl.integrated")
zl.integrated <- NormalizeData(zl.integrated)
zl.integrated <- FindVariableFeatures(zl.integrated, selection.method = "vst", nfeatures = 2000)
zl.integrated <- ScaleData(object = zl.integrated)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
zl.integrated <- CellCycleScoring(zl.integrated, s.features = s.genes, g2m.features = g2m.genes)
zl.integrated <- ScaleData(zl.integrated, vars.to.regress = c("S.Score", "G2M.Score"))

zl.integrated <- RunPCA(pc.genes = zl.integrated@var.genes, npcs = 20 , object = zl.integrated)
zl.integrated <- RunHarmony(zl.integrated, "orig.ident",max.iter.harmony = 30,theta = 0.8 , plot_convergence = TRUE)
ElbowPlot(zl.integrated)
zl.integrated <- FindNeighbors(zl.integrated, reduction = "harmony", dims = 1:15)
zl.integrated <- FindClusters(zl.integrated, reduction = "harmony",resolution = 1.2)
zl.integrated <- RunUMAP(zl.integrated, reduction = "harmony", dims = 1:15)
DimPlot(zl.integrated, reduction = "umap",label = TRUE)
saveRDS(zl.integrated,"zl.integrated.cellbender.rds")
zl.integrated <- readRDS("zl.integrated.cellbender.rds")

copykat_pred_zl1 <- read.xlsx("zl1_copykat_prediction.xlsx","Sheet1")
copykat_pred_zl2 <- read.xlsx("zl2_copykat_prediction.xlsx","Sheet1")
copykat_pred_zl9 <- read.xlsx("zl9_copykat_prediction.xlsx","Sheet1")
copykat_pred_zl10 <- read.xlsx("zl10_copykat_prediction.xlsx","Sheet1")
copykat_pred_zl13 <- read.xlsx("zl13_copykat_prediction.xlsx","Sheet1")
copykat_pred <- rbind(copykat_pred_zl1,copykat_pred_zl2,copykat_pred_zl9,copykat_pred_zl10,copykat_pred_zl13)
copykat_pred_1 <- subset(copykat_pred,cell.names %in% colnames(zl.integrated))
zl.integrated <- zl.integrated[,which(colnames(zl.integrated) %in% copykat_pred$cell.names)]
copykat.pred <- factor(copykat_pred$copykat.pred)

copykat_pred_1 <- subset(copykat_pred,cell.names %in% colnames(zl.integrated))
zl.integrated@meta.data$cell.names <- rownames(zl.integrated@meta.data)
zl.integrated@meta.data <- merge(zl.integrated@meta.data,copykat_pred_1)
copykat.pred <- zl.integrated@meta.data$copykat.pred
copykat.pred <- factor(copykat.pred)
zl.integrated@meta.data$copykat.pred <- copykat.pred
rownames(zl.integrated@meta.data) <- zl.integrated@meta.data$cell.names
DimPlot(zl.integrated, reduction = "umap", group.by = "copykat.pred")
DimPlot(zl.integrated.orig, reduction = "umap", group.by = "copykat.pred")

celltype2 <- read.xlsx("celltype.integrated 1.xlsx","Sheet1")
new.cluster.ids <- celltype2$celltype2
names(new.cluster.ids) <- levels(zl.integrated.orig)
zl.integrated.orig <- RenameIdents(zl.integrated.orig, new.cluster.ids)
DimPlot(zl.integrated.orig, reduction = "umap", label = F, pt.size = 0.5) 
zl.integrated.orig@meta.data$celltype <- Idents(zl.integrated.orig)
celltype <- zl.integrated.orig@meta.data$celltype
celltype <- as.data.frame(celltype)
rownames(celltype)  <- rownames(zl.integrated.orig@meta.data)
zl.integrated <- zl.integrated[,which(colnames(zl.integrated) %in% rownames(celltype))]
celltype_1 <- celltype
celltype$cell.names <- rownames(celltype)
celltype_1 <- subset(celltype, cell.names %in% colnames(zl.integrated)) 
zl.integrated@meta.data <- merge(zl.integrated@meta.data,celltype_1)
rownames(zl.integrated@meta.data) <- zl.integrated@meta.data$cell.names
DimPlot(zl.integrated, reduction = "umap", group.by = "celltype")
zl.integrated.M.orig <- readRDS("~/UPS/cellbender/zl.integrated.M.b.score.rds")

DimPlot(zl.integrated.M.orig, reduction = "umap")
zl.integrated.M.orig@meta.data$celltype <- Idents(zl.integrated.M.orig)
zl.integrated.M <- subset(zl.integrated,cell.names %in% colnames(zl.integrated.M.orig))
DimPlot(zl.integrated.M, reduction = "umap", group.by = "celltype")

zl.integrated.M <- NormalizeData(zl.integrated.M)
zl.integrated.M <- FindVariableFeatures(zl.integrated.M, selection.method = "vst", nfeatures = 2000)
zl.integrated.M <- ScaleData(zl.integrated.M, vars.to.regress = c("S.Score", "G2M.Score"))

zl.integrated.M <- RunPCA(pc.genes = zl.integrated.M@var.genes, npcs = 20 , object = zl.integrated.M)
zl.integrated.M <- RunHarmony(zl.integrated.M, "orig.ident",max.iter.harmony = 30,theta = 0.8 , plot_convergence = TRUE)
ElbowPlot(zl.integrated.M)
zl.integrated.M <- FindNeighbors(zl.integrated.M, reduction = "harmony", dims = 1:15)
zl.integrated.M <- FindClusters(zl.integrated.M, reduction = "harmony",resolution = 1.2)
zl.integrated.M <- RunUMAP(zl.integrated.M, reduction = "harmony", dims = 1:15)
DimPlot(zl.integrated.M, reduction = "umap",label = TRUE)

celltype.M <- zl.integrated.M.orig@meta.data$celltype
celltype.M <- as.data.frame(celltype.M)
rownames(celltype.M)  <- rownames(zl.integrated.M.orig@meta.data)
celltype.M_1 <- celltype.M
celltype.M_1$cell.names <- rownames(celltype.M)
celltype.M_1 <- subset(celltype.M_1, cell.names %in% colnames(zl.integrated.M)) 
zl.integrated.M@meta.data <- merge(zl.integrated.M@meta.data,celltype.M_1)
rownames(zl.integrated.M@meta.data) <- zl.integrated.M@meta.data$cell.names
DimPlot(zl.integrated.M, reduction = "umap", group.by = "celltype.M")

copykat.pred.M <- zl.integrated.M.orig@meta.data$copykat.pred
copykat.pred.M <- as.data.frame(copykat.pred.M)
rownames(copykat.pred.M)  <- rownames(zl.integrated.M.orig@meta.data)
copykat.pred.M_1 <- copykat.pred.M
copykat.pred.M_1$cell.names <- rownames(copykat.pred.M)
copykat.pred.M_1 <- subset(copykat.pred.M_1, cell.names %in% colnames(zl.integrated.M)) 
zl.integrated.M@meta.data <- merge(zl.integrated.M@meta.data,copykat.pred.M_1)
rownames(zl.integrated.M@meta.data) <- zl.integrated.M@meta.data$cell.names
DimPlot(zl.integrated.M, reduction = "umap", group.by = "copykat.pred.M")

FeaturePlot(zl.integrated.M, features = c("AIF1","CD68","COL1A1","TAGLN","COL3A1","NME4","DLK1","DKK1"),pt.size = 0.04)
FeaturePlot(zl.integrated.M.orig, features = c("AIF1","CD68","COL1A1","TAGLN","COL3A1","NME4","DLK1","DKK1"),pt.size = 0.04)

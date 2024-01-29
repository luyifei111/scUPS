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
library(openxlsx)
library(harmony)
library(NMF)
library(irGSEA)
library(pheatmap)
library(RColorBrewer)
library(EnhancedVolcano)
library(CytoTRACE)


HB.genes_total <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ") 

zl1.data <- Read10X(data.dir = "D:/R/UPS/ZL-1/ZL-1/filtered_feature_bc_matrix")
zl1 <- CreateSeuratObject(counts = zl1.data, project = "zl1", min.cells = 3, min.features = 500)
d_pred <- read.table("zl1doublets.txt",sep = ',',header = T)
zl1 <- zl1[,which(colnames(zl1) %in% d_pred$barcode[which(d_pred$predicted_doublets == 'False')])]
HB_m <- match(HB.genes_total,rownames(zl1@assays$RNA))
HB.genes <- rownames(zl1@assays$RNA)[HB_m]
HB.genes <- HB.genes[!is.na(HB.genes)]
zl1[["percent.HB"]] <- PercentageFeatureSet(zl1,features=HB.genes)
zl1[["percent.mt"]] <- PercentageFeatureSet(zl1, pattern = "^MT-")
#VlnPlot(zl1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), ncol = 4)
zl1 <- subset(zl1, subset = nFeature_RNA > 250 & nFeature_RNA < 6000 & percent.mt < 20)
zl1 <- NormalizeData(zl1)
zl1 <- FindVariableFeatures(zl1, selection.method = "vst", nfeatures = 8000)
zl1 <- ScaleData(object = zl1, verbose = FALSE)
zl1 <- RunPCA(pc.genes = zl1@var.genes, npcs = 20 ,verbose = FALSE, object = zl1)
saveRDS(zl1,file = 'zl1_aPCA.rds')
gc()
zl2.data <- Read10X(data.dir = "D:/R/UPS/ZL-2/ZL-2/filtered_feature_bc_matrix")
zl2 <- CreateSeuratObject(counts = zl2.data, project = "zl2", min.cells = 3, min.features = 500)
d_pred <- read.table("zl2doublets.txt",sep = ',',header = T)
zl2 <- zl2[,which(colnames(zl2) %in% d_pred$barcode[which(d_pred$predicted_doublets == 'False')])]
HB_m <- match(HB.genes_total,rownames(zl2@assays$RNA))
HB.genes <- rownames(zl2@assays$RNA)[HB_m]
HB.genes <- HB.genes[!is.na(HB.genes)]
zl2[["percent.HB"]] <- PercentageFeatureSet(zl2,features=HB.genes)
zl2[["percent.mt"]] <- PercentageFeatureSet(zl2, pattern = "^MT-")
#VlnPlot(zl2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), ncol = 4)
zl2 <- subset(zl2, subset = nFeature_RNA > 250 & nFeature_RNA < 6000 & percent.mt < 20)
zl2 <- NormalizeData(zl2)
zl2 <- FindVariableFeatures(zl2, selection.method = "vst", nfeatures = 8000)
zl2 <- ScaleData(object = zl2, verbose = FALSE)
zl2 <- RunPCA(pc.genes = zl2@var.genes, npcs = 20 ,verbose = FALSE, object = zl2)
saveRDS(zl2,file = 'zl2_aPCA.rds')
gc()
zl9.data <- Read10X(data.dir = "D:/R/UPS/ZL-9/ZL-9/filtered_feature_bc_matrix")
zl9 <- CreateSeuratObject(counts = zl9.data, project = "zl9", min.cells = 3, min.features = 500)
d_pred <- read.table("zl9doublets.txt",sep = ',',header = T)
zl9 <- zl9[,which(colnames(zl9) %in% d_pred$barcode[which(d_pred$predicted_doublets == 'False')])]
HB_m <- match(HB.genes_total,rownames(zl9@assays$RNA))
HB.genes <- rownames(zl9@assays$RNA)[HB_m]
HB.genes <- HB.genes[!is.na(HB.genes)]
zl9[["percent.HB"]] <- PercentageFeatureSet(zl9,features=HB.genes)
zl9[["percent.mt"]] <- PercentageFeatureSet(zl9, pattern = "^MT-")
#VlnPlot(zl9, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), ncol = 4)
zl9 <- subset(zl9, subset = nFeature_RNA > 250 & nFeature_RNA < 6000 & percent.mt < 20)
zl9 <- NormalizeData(zl9)
zl9 <- FindVariableFeatures(zl9, selection.method = "vst", nfeatures = 8000)
zl9 <- ScaleData(object = zl9, verbose = FALSE)
zl9 <- RunPCA(pc.genes = zl9@var.genes, npcs = 20 ,verbose = FALSE, object = zl9)
saveRDS(zl9,file = 'zl9_aPCA.rds')
gc()
zl10.data <- Read10X(data.dir = "D:/R/UPS/ZL-10/ZL-10/filtered_feature_bc_matrix")
zl10 <- CreateSeuratObject(counts = zl10.data, project = "zl10", min.cells = 3, min.features = 500)
d_pred <- read.table("zl10doublets.txt",sep = ',',header = T)
zl10 <- zl10[,which(colnames(zl10) %in% d_pred$barcode[which(d_pred$predicted_doublets == 'False')])]
HB_m <- match(HB.genes_total,rownames(zl10@assays$RNA))
HB.genes <- rownames(zl10@assays$RNA)[HB_m]
HB.genes <- HB.genes[!is.na(HB.genes)]
zl10[["percent.HB"]] <- PercentageFeatureSet(zl10,features=HB.genes)
zl10[["percent.mt"]] <- PercentageFeatureSet(zl10, pattern = "^MT-")
#VlnPlot(zl10, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), ncol = 4)
zl10 <- subset(zl10, subset = nFeature_RNA > 250 & nFeature_RNA < 6000 & percent.mt < 20)
zl10 <- NormalizeData(zl10)
zl10 <- FindVariableFeatures(zl10, selection.method = "vst", nfeatures = 8000)
zl10 <- ScaleData(object = zl10, verbose = FALSE)
zl10 <- RunPCA(pc.genes = zl10@var.genes, npcs = 20 ,verbose = FALSE, object = zl10)
saveRDS(zl10,file = 'zl10_aPCA.rds')
gc()
zl13.data <- Read10X(data.dir = "D:/R/UPS/ZL-13/ZL-13/filtered_feature_bc_matrix")
zl13 <- CreateSeuratObject(counts = zl13.data, project = "zl13", min.cells = 3, min.features = 500)
d_pred <- read.table("ZL13doublets.txt",sep = ',',header = T)
zl13 <- zl13[,which(colnames(zl13) %in% d_pred$barcode[which(d_pred$predicted_doublets == 'False')])]
HB_m <- match(HB.genes_total,rownames(zl13@assays$RNA))
HB.genes <- rownames(zl13@assays$RNA)[HB_m]
HB.genes <- HB.genes[!is.na(HB.genes)]
zl13[["percent.HB"]] <- PercentageFeatureSet(zl13,features=HB.genes)
zl13[["percent.mt"]] <- PercentageFeatureSet(zl13, pattern = "^MT-")
#VlnPlot(zl13, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), ncol = 4)
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

copykat_pred_zl1 <- read.xlsx("zl1_copykat_prediction.xlsx","Sheet1")
copykat_pred_zl2 <- read.xlsx("zl2_copykat_prediction.xlsx","Sheet1")
copykat_pred_zl9 <- read.xlsx("zl9_copykat_prediction.xlsx","Sheet1")
copykat_pred_zl10 <- read.xlsx("zl10_copykat_prediction.xlsx","Sheet1")
copykat_pred_zl13 <- read.xlsx("zl13_copykat_prediction.xlsx","Sheet1")
copykat_pred <- rbind(copykat_pred_zl1,copykat_pred_zl2,copykat_pred_zl9,copykat_pred_zl10,copykat_pred_zl13)
zl.integrated <- zl.integrated[,which(colnames(zl.integrated) %in% copykat_pred$cell.names)]
copykat.pred <- factor(copykat_pred$copykat.pred)
zl.integrated@meta.data$copykat.pred <- copykat.pred

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
DimPlot(zl.integrated, reduction = "umap", split.by = "orig.ident",label = TRUE)
DimPlot(zl.integrated, reduction = "umap", split.by = "copykat.pred",label = TRUE)
DimPlot(zl.integrated,reduction = "umap", group.by = "copykat.pred")

markers <- c("RGS5","NOTCH3","ACTA2","PTPRC","NKG7","GNLY","KLRC1","TPSAB1","TPSB2","KIT",
             "PECAM1","VWF","PLVAP","CD3D","CD4","CD8A","COL3A1","COL11A1",
             "DKK1","DLK1","CRIP1","THY1","ENG","NT5E","CD68","CD163","AIF1")
indents1 <- unique(celltype2$celltype2)
DotPlot(zl.integrated,features = markers,idents = indents1, scale = F)

zl.integrated.markers <- FindAllMarkers(zl.integrated, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)
refdata <-  celldex::HumanPrimaryCellAtlasData()
refdata <- scuttle::logNormCounts(refdata)
testdata <- GetAssayData(zl.integrated,slot="data")
clusters <- zl.integrated@meta.data$seurat_clusters
cellpred <- SingleR::SingleR(testdata, 
                             refdata,
                             labels = refdata$label.main,
                             clusters = clusters,
                             assay.type.test = "logcounts",
                             assay.type.ref = "logcounts")
celltype <- data.frame(clusterID = rownames(cellpred),
                       celltype = cellpred$labels,
                       stringsAsFactors = F)
write.xlsx(celltype,"D:/R/UPS/Integrated/celltype.integrated.xlsx")
VlnPlot(zl.integrated,features = "PTPRC")
FeaturePlot(zl.integrated,features = "PTPRC")
celltype2 <- read.xlsx("celltype.integrated 1.xlsx","Sheet1")
new.cluster.ids <- celltype2$celltype2
names(new.cluster.ids) <- levels(zl.integrated)
zl.integrated <- RenameIdents(zl.integrated, new.cluster.ids)
DimPlot(zl.integrated, reduction = "umap", label = F, pt.size = 0.5) 

marker2 <- FindMarkers(zl.integrated,ident.1 = "Tumor")
marker2 <- subset(marker2, p_val_adj<0.01&abs(avg_log2FC)>1)
write.table(marker2,file="D:/R/UPS/Integrated/enrichment/symbol.xls",sep="\t",quote=F,row.names=T)

zl.integrated.tumor <-subset(zl.integrated, idents = c("Unidentified 1","Unidentified 2","Unidentified 3","Unidentified 4","UPS")) 
zl.integrated.tumor <- NormalizeData(zl.integrated.tumor)
zl.integrated.tumor <- FindVariableFeatures(zl.integrated.tumor, selection.method = "vst", nfeatures = 4000)
zl.integrated.tumor <- ScaleData(object = zl.integrated.tumor, verbose = FALSE)
zl.integrated.tumor <- RunPCA(pc.genes = zl.integrated.tumor@var.genes, npcs = 20 ,verbose = FALSE, object = zl.integrated.tumor)
zl.integrated.tumor <- RunHarmony(zl.integrated.tumor, "orig.ident",plot_convergence = TRUE)
#zl.integrated.tumor <- JackStraw(zl.integrated.tumor, dims=50, num.replicate = 100)
#zl.integrated.tumor <- ScoreJackStraw(zl.integrated.tumor, dims = 1:50)
#JackStrawPlot(zl.integrated.tumor, dims = 1:50)
ElbowPlot(zl.integrated.tumor)
zl.integrated.tumor <- FindNeighbors(zl.integrated.tumor, reduction = "harmony", dims = 1:15)
zl.integrated.tumor <- FindClusters(zl.integrated.tumor, resolution = 1.2)
zl.integrated.tumor <- RunUMAP(zl.integrated.tumor, reduction = "harmony", dims = 1:15)
DimPlot(zl.integrated.tumor, reduction = "umap",label = TRUE)
DimPlot(zl.integrated.tumor, reduction = "umap", split.by = "orig.ident",label = TRUE)
DimPlot(zl.integrated.tumor, reduction = "umap", split.by = "copykat.pred",label = TRUE)
DimPlot(zl.integrated.tumor, group.by = "copykat.pred",label = TRUE)
zl.integrated.tumor.markers <- FindAllMarkers(zl.integrated.tumor,only.pos = T ,min.pct = 0.25, logfc.threshold = 0.25)
testdata.tumor <- GetAssayData(zl.integrated.tumor,slot="data")
clusters.tumor <- zl.integrated.tumor@meta.data$seurat_clusters
cellpred.tumor <- SingleR::SingleR(testdata.tumor, 
                                   refdata, 
                                   labels = refdata$label.main,
                                   clusters = clusters.tumor,
                                   assay.type.test = "logcounts",
                                   assay.type.ref = "logcounts")
celltype.tumor <- data.frame(clusterID= rownames(cellpred.tumor),
                             celltype = cellpred.tumor$labels,
                             stringsAsFactors = FALSE)
write.xlsx(celltype.tumor,"D:/R/UPS/Integrated/tumor.figs/celltype.tumor.zl.integrated.xlsx")
celltype2.tumor <- read.xlsx("E:/UPS/Integrated/tumor.figs/celltype.tumor.zl.integrated.xlsx","Sheet1")
new.cluster.tumor.ids <- celltype2.tumor$celltype2
names(new.cluster.tumor.ids) <- levels(zl.integrated.tumor)
zl.integrated.tumor <- RenameIdents(zl.integrated.tumor, new.cluster.tumor.ids)
DimPlot(zl.integrated.tumor, reduction = "umap", label = TRUE, pt.size = 0.5)
DimPlot(zl.integrated.tumor, reduction = "umap", label = TRUE,split.by = "orig.ident" ,pt.size = 0.5)
zl.integrated.tumor <-subset(zl.integrated.tumor, idents = c("CAF","DKK1+ UPS","Rhabdomyosarcoma","MSC","NOTCH+ UPS","DKK1+ NOTCH+ UPS"))
zl.integrated.tumor.markers.D <- FindAllMarkers(zl.integrated.tumor)
markers.T <- c("DKK1","DLK1","CRIP1","COL11A1")
FeaturePlot(zl.integrated.tumor,features = markers.T,ncol = 2,keep.scale = NULL,cols = c("lightgrey", "red4") )

markers.TV <- c("COL3A1","ZFHX4","CFD","DCN","CRABP2","S100A10","CRABP1","NME4","ENG","THY1","NT5E")
VlnPlot(zl.integrated.tumor,features = markers.TV, pt.size = 0,flip = T)

zl.integrated.tumor.marker.subset <- subset(zl.integrated.tumor.markers,zl.integrated.tumor.markers$cluster %in% c(0,1))
unique <- as.data.frame(table(zl.integrated.tumor.marker.subset$gene))
unique <- subset(unique,unique$Freq == 1)
unique <- subset(zl.integrated.tumor.marker.subset, zl.integrated.tumor.marker.subset$gene %in% unique$Var1)

symbols <- paste('D:/R/UPS/Integrated/enrichment/symbol.',c("CAF","DKK1+ UPS","CRIP1+ UPS","NME4+ UPS","MSC"),'.txt',sep = '')

tumornames <- c("CAF","DKK1+ UPS","CRIP1+ UPS","DLK1+ UPS","MSC")
for (i in c(1:5)) {
  markers <- subset(zl.integrated.tumor.markers.D, p_val_adj<0.01 & abs(avg_log2FC)>1.0 & cluster==tumornames[i])
  markers <- cbind(markers[,7],markers[,2])
  colnames(markers) <- c("gene","logFC")
  write.table(markers,file=symbols[i],quote=F,row.names=F)
  gc()
}


library(monocle3)
library(tidyverse)

zl.integrated.tumor@meta.data$cell_type <- Idents(zl.integrated.tumor)

data_plot <- table(zl.integrated.tumor@meta.data$orig.ident, zl.integrated.tumor@meta.data$cell_type) %>% melt()
colnames(data_plot) <- c("Sample", "CellType","Number")
colors <- c("#F8766D","#00BF7D","#A3A500","#00B0F6","#E76BF3")
pC1 <- ggplot(data = data_plot, aes(x = Sample, y = Number, fill = CellType)) +
  geom_bar(stat = "identity", width=0.8,aes(group=CellType),position="stack")+
  scale_fill_manual(values=colors) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Average number")+
  theme(axis.text = element_text(size=12, colour = "black"))+
  theme(axis.title.y = element_text(size=12, colour = "black"))+
  theme(panel.border = element_rect(size = 1, linetype = "solid", colour = "black"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.8,  vjust = 0.6)) 


pC2 <- ggplot(data = data_plot, aes(x = Sample, y = Number, fill = CellType)) +
  geom_bar(stat = "identity", width=0.8,aes(group=CellType),position="fill")+
  scale_fill_manual(values=colors) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Cell proportion")+
  scale_y_continuous(labels = percent)+  ####用来将y轴移动位置
  theme(axis.text = element_text(size=12, colour = "black"))+
  theme(axis.title.y = element_text(size=12, colour = "black"))+
  theme(panel.border = element_rect(size = 1, linetype = "solid", colour = "black"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.8,  vjust = 0.6))      #让横轴上的标签倾斜45度

pC <- pC1 + pC2 + plot_layout(ncol = 2, widths = c(1,1),guides = 'collect')
pC2
pC



data.orig <-subset(zl.integrated.tumor, idents = tumornames)
#data.orig <-subset(data.orig, orig.ident %in% c("zl13"))
cell_metadata <- data.orig@meta.data
data <- GetAssayData(data.orig, assay = 'RNA', slot = 'counts')
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 50)
plot_pc_variance_explained(cds)
cds <- reduce_dimension(cds,preprocess_method = "PCA")
plot_cells(cds)
p1 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="cell_type") + ggtitle('cds.umap')
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(data.orig, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="cell_type") + ggtitle('int.umap')
p1|p2
cds <- cluster_cells(cds)
p1 <- plot_cells(cds, show_trajectory_graph = FALSE) + ggtitle("label by clusterID")
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE) + 
  ggtitle("label by partitionID")
wrap_plots(p1, p2)
cds <- learn_graph(cds)
plot_cells(cds, color_cells_by = "cell_type", label_groups_by_cluster = FALSE,
           label_leaves= TRUE , label_branch_points = TRUE , graph_label_size=1.5)

cds <- order_cells(cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,
           trajectory_graph_color = "green4",
           trajectory_graph_segment_size = 1.2)
Track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=1)
trace('calculateLW', edit = T, where = asNamespace("monocle3"))
Track_genes_sig <- Track_genes %>% top_n(n=10, morans_I) %>%
  pull(gene_short_name) %>% as.character()
plot_genes_in_pseudotime(cds[Track_genes_sig,], color_cells_by="cell_type", 
                         min_expr=0.5, ncol = 2)
Track_genes_sig_1 <- c("DKK1","CRIP1","NME4","COL11A1","NEAT1","THY1","APOD","DLK1","CFD","DCN","CRABP2","CRABP1")
plot_genes_in_pseudotime(cds[Track_genes_sig_1,], color_cells_by="cell_type", 
                         min_expr=0.5, ncol = 2)
plot_cells(cds, genes=Track_genes_sig, show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,  label_leaves=FALSE)
genelist <- pull(Track_genes, gene_short_name) %>% as.character()
gene_module_df <- find_gene_modules(cds[genelist,], resolution=1e-2, cores = 1)
cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                cell_group=colData(cds)$cell_type) 
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
pheatmap::pheatmap(agg_mat, scale="column", clustering_method="ward.D2")

data.orig <-zl.integrated
cell_metadata <- zl.integrated@meta.data
data <- GetAssayData(zl.integrated, assay = 'RNA', slot = 'counts')
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 50)
plot_pc_variance_explained(cds)
cds <- reduce_dimension(cds,preprocess_method = "PCA")
plot_cells(cds)
p1 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="cell_type") + ggtitle('cds.umap')
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(zl.integrated, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="cell_type") + ggtitle('int.umap')
p1|p2
cds <- cluster_cells(cds)
p1 <- plot_cells(cds, show_trajectory_graph = FALSE) + ggtitle("label by clusterID")
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE) + 
  ggtitle("label by partitionID")
wrap_plots(p1, p2)
cds <- learn_graph(cds)
plot_cells(cds, color_cells_by = "cell_type", label_groups_by_cluster = FALSE,
           label_leaves= TRUE , label_branch_points = TRUE , graph_label_size=1.5)
cds <- order_cells(cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

FeaturePlot(zl.integrated,features = c("MDM2","CDK4","MYOG","MYOD1","DES","SMN1"))

data.orig <-zl.integrated
cell_metadata <- zl.integrated@meta.data
data <- GetAssayData(zl.integrated, assay = 'RNA', slot = 'counts')
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 50)
plot_pc_variance_explained(cds)
cds <- reduce_dimension(cds,preprocess_method = "PCA")
plot_cells(cds)
p1 <- plot_cells(cds, reduction_method="UMAP") + ggtitle('cds.umap')
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(zl.integrated, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(cds, reduction_method="UMAP") + ggtitle('int.umap')
p1|p2
cds <- cluster_cells(cds)
p1 <- plot_cells(cds, show_trajectory_graph = FALSE) + ggtitle("label by clusterID")
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE) + 
  ggtitle("label by partitionID")
wrap_plots(p1, p2)
cds <- learn_graph(cds)
plot_cells(cds, label_groups_by_cluster = FALSE,
           label_leaves= TRUE , label_branch_points = TRUE , graph_label_size=1.5)
cds <- order_cells(cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
cds <- order_cells(cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
Track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=1)
Track_genes_sig <- Track_genes %>% top_n(n=10, morans_I) %>%
  pull(gene_short_name) %>% as.character()
plot_genes_in_pseudotime(cds[Track_genes_sig,], color_cells_by="cell_type", 
                         min_expr=0.5, ncol = 2)
plot_cells(cds, genes=Track_genes_sig, show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,  label_leaves=FALSE)
genelist <- pull(Track_genes, gene_short_name) %>% as.character()
gene_module_df <- find_gene_modules(cds[genelist,], resolution=1e-2, cores = 10)
cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                cell_group=colData(cds)$cell_type)
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
pheatmap::pheatmap(agg_mat, scale="column", clustering_method="ward.D2")



CRIP1_CAF_DEG <- FindMarkers(zl.integrated.tumor,ident.1 = "CRIP1+ UPS", ident.2 = "CAF")
DKK1_DLK1_DEG <- FindMarkers(zl.integrated.tumor,ident.1 = "DKK1+ UPS", ident.2 = "DLK1+ UPS")

zl.tumor.averageexp <- AverageExpression(zl.integrated.tumor)
zl.tumor.averageexp_1 <- as.data.frame(zl.tumor.averageexp$RNA)
zl.tumor.averageexp_1$feature <- rownames(zl.tumor.averageexp_1)
CRIP1_CAF_DEG$feature <- rownames(CRIP1_CAF_DEG)
DKK1_DLK1_DEG$feature <- rownames(DKK1_DLK1_DEG)
zl.tumor.averageexp_1_CAF <- zl.tumor.averageexp_1[,c(1,6)]
zl.tumor.averageexp_1_DKK1 <- zl.tumor.averageexp_1[,c(2,6)]
zl.tumor.averageexp_1_CRIP1 <- zl.tumor.averageexp_1[,c(3,6)]
zl.tumor.averageexp_1_DLK1 <- zl.tumor.averageexp_1[,c(4,6)]
CRIP1_DEG <- subset(CRIP1_CAF_DEG, subset = CRIP1_CAF_DEG$avg_log2FC > 0)
CAF_DEG <- subset(CRIP1_CAF_DEG, subset = CRIP1_CAF_DEG$avg_log2FC < 0)
DKK1_DEG <- subset(DKK1_DLK1_DEG, subset = DKK1_DLK1_DEG$avg_log2FC > 0)
DLK1_DEG <- subset(DKK1_DLK1_DEG, subset = DKK1_DLK1_DEG$avg_log2FC < 0)
zl.tumor.averageexp_1_CAF <- subset(zl.tumor.averageexp_1_CAF,subset = zl.tumor.averageexp_1_CAF$feature %in% CAF_DEG$feature)
zl.tumor.averageexp_1_DKK1 <- subset(zl.tumor.averageexp_1_DKK1,subset = zl.tumor.averageexp_1_DKK1$feature %in% DKK1_DEG$feature)
zl.tumor.averageexp_1_CRIP1 <- subset(zl.tumor.averageexp_1_CRIP1,subset = zl.tumor.averageexp_1_CRIP1$feature %in% CRIP1_DEG$feature)
zl.tumor.averageexp_1_DLK1 <- subset(zl.tumor.averageexp_1_DLK1,subset = zl.tumor.averageexp_1_DLK1$feature %in% DLK1_DEG$feature)
CRIP1_DEG <- merge(CRIP1_DEG,zl.tumor.averageexp_1_CRIP1, by = "feature")
CAF_DEG <- merge(CAF_DEG,zl.tumor.averageexp_1_CAF, by = "feature")
DKK1_DEG <- merge(DKK1_DEG,zl.tumor.averageexp_1_DKK1, by = "feature")
DLK1_DEG <- merge(DLK1_DEG,zl.tumor.averageexp_1_DLK1, by = "feature")
colnames(CRIP1_DEG)[7] <- c("average_expression")
colnames(CAF_DEG)[7] <- c("average_expression")
colnames(DKK1_DEG)[7] <- c("average_expression")
colnames(DLK1_DEG)[7] <- c("average_expression")
CRIP1_CAF_DEG_1 <- rbind(CRIP1_DEG,CAF_DEG) 
DKK1_DLK1_DEG_1 <- rbind(DKK1_DEG,DLK1_DEG)

df <- CRIP1_CAF_DEG_1
rownames(df) <- df$feature
exthreshold <- 5
fcthreshold <- 2

df$change <- as.factor(ifelse(df$average_expression > exthreshold & abs(df$avg_log2FC) > fcthreshold,
                              ifelse(df$avg_log2FC > fcthreshold,"Up","Down"),"Non"))
df$label <- ifelse(df$average_expression > exthreshold & abs(df$avg_log2FC)>fcthreshold,as.character(rownames(df)),"")

p.vol <- ggplot(data = df,
                aes(x = avg_log2FC,y = average_expression,colour = change,fill = change))+
                coord_cartesian(ylim = c(0,50))+
                scale_color_manual(values = c('green','grey','red'))+
                geom_point(alpha = 0.4,size = 3.5)+
                geom_text_repel(aes(x = avg_log2FC,y = average_expression,label = label),size = 3,
                box.padding = unit(0.6,"lines"),point.padding = unit(0.7,"lines"),
                segment.color = "black",show.legend = FALSE)+
                geom_vline(xintercept = c(-fcthreshold,fcthreshold),lty = 4,col = "black",lwd = 0.8)+
                geom_hline(yintercept = exthreshold,lty = 4,col = "black",lwd = 0.8)+
                theme_bw()+
                labs(x = "log2(FoldChange)",y = "Average Expression",title = "CRIP1 VS CAF")+
                theme(axis.text = element_text(size = 11),axis.title = element_text(size = 13), # 坐标轴标签和标题
                plot.title = element_text(hjust = 0.5,size = 15,face = "bold"), # 标题
                legend.text = element_text(size = 11),legend.title = element_text(size = 13), # 图例标签和标题
                plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm")) # 图边距
p.vol


DKK1.EXPR <- FetchData(object = zl.integrated.tumor.DKK1,vars = c("EGFR","HBEGF"))
ggplot(data = DKK1.EXPR,aes(x = EGFR,y = HBEGF))+
  geom_point()
write.xlsx(DKK1.EXPR,file = "DKK1 EGFR-HBEGF.xlsx")

phe <- zl.integrated.tumor@meta.data$cell_type
phe = as.character(phe)
names(phe) <- rownames(zl.integrated.tumor@meta.data)

mat_zl <- as.matrix(zl.integrated.tumor@assays$RNA@counts)
mat_zl[1:4,1:4]
results <- CytoTRACE(mat = mat_zl)
plotCytoGenes(results,outputDir = "E:/UPS/Integrated/CytoTRACE")
plotCytoTRACE(results,outputDir = "E:/UPS/Integrated/CytoTRACE")
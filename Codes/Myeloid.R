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
zl.integrated.M <-subset(zl.integrated, idents = c("M1 M2 Macrophage Fibroblast Doublets","M1 M2 Macrophage Tumor Doublets","M1 Macrophage Fibroblast Doublets","Monocyte","M1 M2 Macrophage","Mast","M2 Macrophage Tumor Doublets","Mast Tumor Doublets","M2 Macrophage"))
zl.integrated.M <- NormalizeData(zl.integrated.M)
zl.integrated.M <- FindVariableFeatures(zl.integrated.M, selection.method = "vst", nfeatures = 4000)
zl.integrated.M <- ScaleData(object = zl.integrated.M, verbose = FALSE,vars.to.regress = c("S.Score", "G2M.Score"))
zl.integrated.M <- RunPCA(pc.genes = zl.integrated.M@var.genes, npcs = 20 ,verbose = FALSE, object = zl.integrated.M)
zl.integrated.M <- RunHarmony(zl.integrated.M, "orig.ident",plot_convergence = TRUE,max.iter.harmony = 30)
#zl.integrated.M <- JackStraw(zl.integrated.M, dims=50, num.replicate = 100)
#zl.integrated.M <- ScoreJackStraw(zl.integrated.M, dims = 1:50)
#JackStrawPlot(zl.integrated.M, dims = 1:50)
ElbowPlot(zl.integrated.M)
zl.integrated.M <- FindNeighbors(zl.integrated.M, reduction = "harmony", dims = 1:15)
zl.integrated.M <- FindClusters(zl.integrated.M, resolution = 1.2)
zl.integrated.M <- RunUMAP(zl.integrated.M, reduction = "harmony", dims = 1:15)
DimPlot(zl.integrated.M, reduction = "umap",label = TRUE)
DimPlot(zl.integrated.M, reduction = "umap", split.by = "orig.ident",label = TRUE)
DimPlot(zl.integrated.M, reduction = "umap", split.by = "copykat.pred",label = TRUE)
DimPlot(zl.integrated.M, group.by = "copykat.pred",label = TRUE)
zl.integrated.M.markers <- FindAllMarkers(zl.integrated.M,only.pos = T)
testdata.M <- GetAssayData(zl.integrated.M,slot="data")
clusters.M <- zl.integrated.M@meta.data$seurat_clusters
cellpred.M <- SingleR::SingleR(testdata.M, 
                                   refdata, 
                                   labels = refdata$label.main,
                                   clusters = clusters.M,
                                   assay.type.test = "logcounts",
                                   assay.type.ref = "logcounts")
celltype.M <- data.frame(clusterID= rownames(cellpred.M),
                             celltype = cellpred.M$labels,
                             stringsAsFactors = FALSE)
write.xlsx(celltype.M,"D:/R/UPS/Integrated/microenvironment/M/celltype.M.zl.integrated.xlsx")

zl.integrated.M <-subset(zl.integrated.M, idents = c(8,9,11,16,20:23,25))
zl.integrated.M <- NormalizeData(zl.integrated.M)
zl.integrated.M <- FindVariableFeatures(zl.integrated.M, selection.method = "vst", nfeatures = 8000)
zl.integrated.M <- ScaleData(object = zl.integrated.M, verbose = FALSE,vars.to.regress = c("S.Score", "G2M.Score"))
zl.integrated.M <- RunPCA(pc.genes = zl.integrated.M@var.genes, npcs = 20 ,verbose = FALSE, object = zl.integrated.M)
zl.integrated.M <- RunHarmony(zl.integrated.M, "orig.ident",plot_convergence = TRUE,max.iter.harmony = 30)
#zl.integrated.M <- JackStraw(zl.integrated.M, dims=50, num.replicate = 100)
#zl.integrated.M <- ScoreJackStraw(zl.integrated.M, dims = 1:50)
#JackStrawPlot(zl.integrated.M, dims = 1:50)
ElbowPlot(zl.integrated.M)
zl.integrated.M <- FindNeighbors(zl.integrated.M, reduction = "harmony", dims = 1:15)
zl.integrated.M <- FindClusters(zl.integrated.M, resolution = 1.2)
zl.integrated.M <- RunUMAP(zl.integrated.M, reduction = "harmony", dims = 1:15)
DimPlot(zl.integrated.M, reduction = "umap",label = TRUE)
DimPlot(zl.integrated.M, reduction = "umap", split.by = "orig.ident",label = TRUE)
DimPlot(zl.integrated.M, reduction = "umap", split.by = "copykat.pred",label = TRUE)
DimPlot(zl.integrated.M, group.by = "copykat.pred",label = TRUE)
zl.integrated.M.markers <- FindAllMarkers(zl.integrated.M,only.pos = T)
testdata.M <- GetAssayData(zl.integrated.M,slot="data")
clusters.M <- zl.integrated.M@meta.data$seurat_clusters
cellpred.M <- SingleR::SingleR(testdata.M, 
                               refdata, 
                               labels = refdata$label.main,
                               clusters = clusters.M,
                               assay.type.test = "logcounts",
                               assay.type.ref = "logcounts")
celltype.M <- data.frame(clusterID= rownames(cellpred.M),
                         celltype = cellpred.M$labels,
                         stringsAsFactors = FALSE)
write.xlsx(celltype.M,"D:/R/UPS/Integrated/microenvironment/M/celltype.M.zl.integrated.xlsx")
celltype2.M <- read.xlsx("E:/UPS/Integrated/microenvironment/M/celltype.M.zl.integrated.xlsx","Sheet1")
new.cluster.M.ids <- celltype2.M$celltype2
names(new.cluster.M.ids) <- levels(zl.integrated.M)
zl.integrated.M <- RenameIdents(zl.integrated.M, new.cluster.M.ids)
DimPlot(zl.integrated.M, reduction = "umap", label = TRUE, pt.size = 0.5)
DimPlot(zl.integrated.M, reduction = "umap", label = TRUE,split.by = "orig.ident" ,pt.size = 0.5)
markers <- c("CD14","FCN1","LYZ","FCGR3A","HLA-DRA","EREG","SPP1","LYVE1","C1QC","CD1C","CLEC10A","CMTM2","MNDA","KIT","CPA3","TPSAB1")
DotPlot(zl.integrated.M,features = markers,cluster.idents = T)
VlnPlot(zl.integrated.M,features = markers)

Mnames <- c("SPP1+ LYVE1+ Macrophage","Mast","Intermediate Monocyte","CLEC10A+ cDC2","EREG+ Monocyte","Neutrophils","Nonclassical Monocyte")
zl.integrated.M.markers.D <- FindAllMarkers(zl.integrated.M,only.pos = T)
zl.integrated.M.markers.D %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(zl.integrated.M, features = top10$gene)
zl.integrated.M.data <- as.data.frame(GetAssayData(zl.integrated.M))
p=pheatmap(zl.integrated.M.data,show_rownames=T,
           scale = 'row',
           clustering_method = 'ward.D2',
           show_colnames=F,cluster_cols = TRUE,cluster_rows = FALSE,
           color=rev(colorRampPalette(brewer.pal(11,'RdBu'))(50)),
           annotation_row=top10$gene
           #clustering_distance_rows = 'binary'
)
symbols <- paste('D:/R/UPS/Integrated/microenvironment/M/enrichment/symbol.',Mnames,'.txt',sep = '')
for (i in c(1:length(Mnames))) {
  markers <- subset(zl.integrated.M.markers.D, p_val_adj<0.01 & abs(avg_log2FC)>1.0 & cluster==Mnames[i])
  markers <- cbind(markers[,7],markers[,2])
  colnames(markers) <- c("gene","logFC")
  write.table(markers,file=symbols[i],quote=F,row.names=F)
  gc()
}

library(monocle3)
library(tidyverse)

zl.integrated.M@meta.data$cell_type <- Idents(zl.integrated.M)
data.orig <-subset(zl.integrated.M, idents = Mnames)
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
           graph_label_size=1.5)
Track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=10)
trace('calculateLW', edit = T, where = asNamespace("monocle3"))
Track_genes_sig <- Track_genes %>% top_n(n=10, morans_I) %>%
  pull(gene_short_name) %>% as.character()
Track_genes_sig <- c("CD14","FCGR3A","HLA-DRA","EREG","SPP1","LYVE1","C1QC","CD1C","CLEC10A","FCN1")
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

zl.integrated.M <-subset(zl.integrated.M, idents = c(0:7,10,12:15,17:19,24))
zl.integrated.M <- NormalizeData(zl.integrated.M)
zl.integrated.M <- FindVariableFeatures(zl.integrated.M, selection.method = "vst", nfeatures = 4000)
zl.integrated.M <- ScaleData(object = zl.integrated.M, verbose = T,vars.to.regress = c("S.Score", "G2M.Score"))
zl.integrated.M <- RunPCA(pc.genes = zl.integrated.M@var.genes, npcs = 20 ,verbose = FALSE, object = zl.integrated.M)
zl.integrated.M <- RunHarmony(zl.integrated.M, "orig.ident",plot_convergence = TRUE,max.iter.harmony = 30)
#zl.integrated.M <- JackStraw(zl.integrated.M, dims=50, num.replicate = 100)
#zl.integrated.M <- ScoreJackStraw(zl.integrated.M, dims = 1:50)
#JackStrawPlot(zl.integrated.M, dims = 1:50)
ElbowPlot(zl.integrated.M)
zl.integrated.M <- FindNeighbors(zl.integrated.M, reduction = "harmony", dims = 1:15)
zl.integrated.M <- FindClusters(zl.integrated.M, resolution = 1.2)
zl.integrated.M <- RunUMAP(zl.integrated.M, reduction = "harmony", dims = 1:15)
DimPlot(zl.integrated.M, reduction = "umap",label = TRUE)
DimPlot(zl.integrated.M, reduction = "umap", split.by = "orig.ident",label = TRUE)
DimPlot(zl.integrated.M, reduction = "umap", split.by = "copykat.pred",label = TRUE)
DimPlot(zl.integrated.M, group.by = "copykat.pred",label = TRUE)
zl.integrated.M.markers <- FindAllMarkers(zl.integrated.M,only.pos = T)
testdata.M <- GetAssayData(zl.integrated.M,slot="data")
clusters.M <- zl.integrated.M@meta.data$seurat_clusters
cellpred.M <- SingleR::SingleR(testdata.M, 
                               refdata, 
                               labels = refdata$label.main,
                               clusters = clusters.M,
                               assay.type.test = "logcounts",
                               assay.type.ref = "logcounts")
celltype.M <- data.frame(clusterID= rownames(cellpred.M),
                         celltype = cellpred.M$labels,
                         stringsAsFactors = FALSE)
write.xlsx(celltype.M,"E:/UPS/Integrated/microenvironment/M/celltype.M.zl.integrated.xlsx")


celltype2.M <- read.xlsx("E:/UPS/Integrated/microenvironment/M/BACKUP1/celltype.M.zl.integrated.xlsx","Sheet1")
new.cluster.M.ids <- celltype2.M$celltype2
names(new.cluster.M.ids) <- levels(zl.integrated.M)
zl.integrated.M <- RenameIdents(zl.integrated.M, new.cluster.M.ids)
DimPlot(zl.integrated.M, reduction = "umap", label = F, pt.size = 0.5)
DimPlot(zl.integrated.M, reduction = "umap", label = TRUE,split.by = "orig.ident" ,pt.size = 0.5)
zl.integrated.M.markers.D <- FindAllMarkers(zl.integrated.M,only.pos = T)
write.xlsx(zl.integrated.M.markers.D,"E:/UPS/Integrated/microenvironment/M/celltype.M.zl.integrated.xlsx")
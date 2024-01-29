library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(cowplot)
library(sctransform)
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
zl.integrated.T <-subset(zl.integrated, idents = c("NK T + Exhausted CD8+ T","T helper","Exhausted CD8+ T")) 
zl.integrated.T <- NormalizeData(zl.integrated.T)
zl.integrated.T <- FindVariableFeatures(zl.integrated.T, selection.method = "vst", nfeatures = 4000)
zl.integrated.T <- ScaleData(object = zl.integrated.T, verbose = FALSE, vars.to.regress = c("S.Score", "G2M.Score"))
zl.integrated.T <- RunPCA(pc.genes = zl.integrated.T@var.genes, npcs = 20 ,verbose = FALSE, object = zl.integrated.T)
zl.integrated.T <- RunHarmony(zl.integrated.T, "orig.ident",plot_convergence = TRUE)
#zl.integrated.T <- JackStraw(zl.integrated.T, dims=50, num.replicate = 100)
#zl.integrated.T <- ScoreJackStraw(zl.integrated.T, dims = 1:50)
#JackStrawPlot(zl.integrated.T, dims = 1:50)
ElbowPlot(zl.integrated.T)
zl.integrated.T <- FindNeighbors(zl.integrated.T, reduction = "harmony", dims = 1:15)
zl.integrated.T <- FindClusters(zl.integrated.T, resolution = 1.2)
zl.integrated.T <- RunUMAP(zl.integrated.T, reduction = "harmony", dims = 1:15)
DimPlot(zl.integrated.T, reduction = "umap",label = TRUE)
DimPlot(zl.integrated.T, reduction = "umap", split.by = "orig.ident",label = TRUE)
DimPlot(zl.integrated.T, reduction = "umap", split.by = "copykat.pred",label = TRUE)
DimPlot(zl.integrated.T, group.by = "copykat.pred",label = TRUE)
DimPlot(zl.integrated.T, group.by = "Phase",label = TRUE)
zl.integrated.T.markers <- FindAllMarkers(zl.integrated.T,only.pos = T)
testdata.T <- GetAssayData(zl.integrated.T,slot="data")
refdata <-  celldex::HumanPrimaryCellAtlasData()
clusters.T <- zl.integrated.T@meta.data$seurat_clusters
cellpred.T <- SingleR::SingleR(testdata.T, 
                                   refdata, 
                                   labels = refdata$label.main,
                                   clusters = clusters.T,
                                   assay.type.test = "logcounts",
                                   assay.type.ref = "logcounts")
celltype.T <- data.frame(clusterID= rownames(cellpred.T),
                             celltype = cellpred.T$labels,
                             stringsAsFactors = FALSE)
write.xlsx(celltype.T,"D:/R/UPS/Integrated/microenvironment/T/celltype.T.zl.integrated.xlsx")
celltype2.T <- read.xlsx("E:/UPS/Integrated/microenvironment/T/celltype.T.zl.integrated.xlsx","Sheet1")
new.cluster.T.ids <- celltype2.T$celltype2
names(new.cluster.T.ids) <- levels(zl.integrated.T)
zl.integrated.T <- RenameIdents(zl.integrated.T, new.cluster.T.ids)
DimPlot(zl.integrated.T, reduction = "umap", label = TRUE, pt.size = 0.5)
DimPlot(zl.integrated.T, reduction = "umap", label = TRUE,split.by = "orig.ident" ,pt.size = 0.5)
markers <- c("CD3D","CD8A","PDCD1","HAVCR2","LAG3","CTLA4","TNFRSF18","TNFRSF4","GNLY","KLRG1","KLRD1","CD4","CCR7","TCF7","IL7R","CXCR3","CXCR6","IL2RA","FOXP3")
DotPlot(zl.integrated.T,features = markers,cluster.idents = T)
Tnames <- c("Exhausted CD8+ T","Pre-Exhausted CD8+ T","Naive CD4+ T","Cytotoxic T","Regulatory T","Exhausted Cytotoxic CD8+ T","DN T","Exhausted CD4+ T")
zl.integrated.T.markers.D <- FindAllMarkers(zl.integrated.T,only.pos = T)

library(monocle3)
library(tidyverse)
seuratcluster <- zl.integrated.T@meta.data$seurat_clusters
seuratcluster <- as.character(seuratcluster)
seuratcluster = c("0" = "Exhausted CD8+ T" ,"1" = "Pre-Exhausted CD8+ T","2" = "Naive CD4+ T"
                  ,"3" = "Cytotoxic T","4" = "Regulatory T","5" = "Exhausted CD8+ T"            
                  ,"6" = "Exhausted Cytotoxic CD8+ T","7" = "Regulatory T","8" =  "Exhausted CD4+ T"                         
                  ,"9" = "Pre-Exhausted CD8+ T" ,"10" =  "Exhausted Cytotoxic CD8+ T","11" = "Pre-Exhausted CD8+ T"            
                  ,"12" = "Regulatory T" ,"13" = "DN T","15" = "Exhausted Cytotoxic CD8+ T")[ as.character(seuratcluster)]
seuratcluster <- factor(seuratcluster , levels = Tnames)
zl.integrated.T@meta.data$cell_type <- seuratcluster
data.orig <-subset(zl.integrated.T, idents = Tnames)
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
Track_genes_sig <- c("LAG3","TIGIT","PDCD1","HAVCR2","CTLA4")
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

zl.integrated.T.markers.D <- FindAllMarkers(zl.integrated.T,only.pos = T)
symbols <- paste('D:/R/UPS/Integrated/microenvironment/T/enrichment/symbol.',Tnames,'.txt',sep = '')
for (i in c(1:length(Tnames))) {
  markers <- subset(zl.integrated.T.markers.D, p_val_adj<0.01 & abs(avg_log2FC)>1.0 & cluster==Tnames[i])
  markers <- cbind(markers[,7],markers[,2])
  colnames(markers) <- c("gene","logFC")
  write.table(markers,file=symbols[i],quote=F,row.names=F)
  gc()
}

exhausted.markers <- c("LAYN","LAG3","TIGIT","PDCD1","HAVCR2","CTLA4","ITGAE") 
T.score.markers <- list()
T.score.markers$Exhaustion <- exhausted.markers
zl.integrated.T <- irGSEA.score(object = zl.integrated.T, assay = "RNA", slot = "data", 
                                seeds = 123, ncores = 1,msigdb=F, 
                                custom = T, geneset = T.score.markers, method = c("AUCell", "UCell", "singscore",
                                                                                 "ssgsea"), 
                                kcdf = 'Gaussian')
scatterplot <- irGSEA.density.scatterplot(object = zl.integrated.T,
                                          method = c("AUCell", "UCell", "singscore",
                                                     "ssgsea")[4],
                                          show.geneset = c("Exhaustion"),
                                          reduction = "umap")
scatterplot
